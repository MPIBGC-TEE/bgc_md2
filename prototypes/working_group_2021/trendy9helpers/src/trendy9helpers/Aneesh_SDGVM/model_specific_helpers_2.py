import sys
import json 
from pathlib import Path
from collections import namedtuple 
import netCDF4 as nc
import numpy as np
import datetime as dt
from sympy import Symbol, symbols, var 
from typing import Callable, Tuple, Dict
from functools import reduce, partial
from importlib import import_module
from collections import OrderedDict
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.ArrayDictResult import ArrayDictResult

import bgc_md2.helper as h
from ..import general_helpers as gh
model_mod = 'bgc_md2.models.Aneesh_SDGVM'
cp_mod = import_module(f"{model_mod}.CachedParameterization")
make_func_dict = import_module(f"{model_mod}.CachedParameterization").make_func_dict
Drivers=cp_mod.Drivers
CachedParameterization=cp_mod.CachedParameterization


def spatial_mask(dataPath)->'CoorMask':
    mask=nc.Dataset(dataPath.joinpath("SDGVM_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
    sym_tr= gh.SymTransformers(
        itr=make_model_index_transforms(),
        ctr=make_model_coord_transforms()
    )
    return gh.CoordMask(
        mask,
        sym_tr
    )

def make_model_coord_transforms():
    return gh.identicalTransformers()

def make_model_index_transforms():
    return gh.transform_maker(
    lat_0 = -89.5,
    lon_0 = -179.5,
    step_lat = 1,
    step_lon = 1,
 )

Observables = namedtuple(
    'Observables',
    ["cVeg", "cRoot", "cLitter", "cSoil", "rh"]
)

# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated)
# parameters. In this case we may use only the first entry e.g. to derive startvalues.
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise

#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['SDGVM'],
        variables = Observables._fields + Drivers._fields+("gpp","ra")
    )

def get_example_site_vars(dataPath):
    # According to the netcdf metadata the datasets are not uniform
    # - npp and rh start at 360h (15 days) after 01-01-1900 and are recorded every 30 days
    #   these are interpreted as mid-monthly
    # - C_litter, C_soil, C_veg, C_root start at 4320 h = 180 days = 6 months after 01-01-1700
    #   These can be seen at midyearly values since there are 6 (midmonth) times of the npp and rh measurements after the last (midyear)
    #   measurement of C_litter, C_soil, C_veg, C_root

    # To combine these streams into a consistent array of observations we will:
    # 1. Make C_litter, C_soil, C_veg, C_root refer to hours after 01/01/1900 (as npp and rh)
    #
    # 2. cut away everything before 1900 from them (cutting of the first 200y)
    #
    # Note:
    #    We will have to adapt the costfunction and param2res later to accommodate the different
    #    resolution of the C_pool and rh observations.



    # 1.)
    # pick one of the 1700 yearly example ds to get at the times
    # convert time to refer to the same starting point (from h after 1700 to h after 1900)
    hs_from_1900=nc.Dataset(dataPath.joinpath('SDGVM_S2_cLitter.nc')).variables['time'][:]-200*12*30*24

    #2.)
    # find the index after which we are after 01/01 1900
    ind_start = 200

    # pick up 1 site   wombat state forest for the spacial selection
    s_rh  = slice(None, None, None)  # this is the same as :
    s_c  = slice(ind_start, None, None)  # this is the same as ind_start:
    #t = s, 50, 33  # [t] = [:,49,325]
    loc=(-25,16)
    t_rh = s_rh,*loc
    t_c = s_c, *loc
    print(t_c)

    # Read NetCDF data and slice out our site
    arr_dict={
        **{
            vn:nc.Dataset(str(dataPath.joinpath(fn))).variables[vn][t_c]
            for vn,fn in  {
                'cLitter': 'SDGVM_S2_cLitter.nc',
                'cSoil': 'SDGVM_S2_cSoil.nc',
                'cVeg': 'SDGVM_S2_cVeg.nc',
                'cRoot': 'SDGVM_S2_cRoot.nc',
            }.items()
        },
        **{
            vn:nc.Dataset(str(dataPath.joinpath(fn))).variables[vn][t_rh]*86400   # kg/m2/s kg/m2/day;
            for vn,fn in {
                'npp': 'SDGVM_S2_npp.nc',
                'rh': 'SDGVM_S2_rh.nc'
            }.items()
        }
    }

    return (
        Observables(*(arr_dict[k] for k in Observables._fields)),
        Drivers(*(arr_dict[k] for k in Drivers._fields))
    )

def compute_global_mean_arr_var_dict(dataPath):
    # According to the netcdf metadata the datasets are not uniform
    # - npp and rh start at 360h (15 days) after 01-01-1900 and are recorded every 30 days
    #   these are interpreted as mid-monthly
    # - C_litter, C_soil, C_veg, C_root start at 4320 h = 180 days = 6 months after 01-01-1700
    #   These can be seen at midyearly values since there are 6 (midmonth) times of the npp and rh measurements after the last (midyear)
    #   measurement of C_litter, C_soil, C_veg, C_root

    # To combine these streams into a consistent array of observations we will:
    # 1. Make C_litter, C_soil, C_veg, C_root refer to hours after 01/01/1900 (as npp and rh)
    #
    # 2. cut away everything before 1900 from them (cutting of the first 200y)

    #2.)
    # find the index after which we are after 01/01 1900
    ind_start = 200
    time_slice = slice(ind_start, None, None)

    # extract templating information from one of the files that has a mask
    ds = nc.Dataset(str(dataPath.joinpath('SDGVM_S2_cLitter.nc')))
    lats = ds.variables["latitude"].__array__()
    lons = ds.variables["longitude"].__array__()
    template = ds.variables['cLitter'][0, :, :].mask
    
    def var(vn):
        return nc.Dataset(
                str(
                    dataPath.joinpath(
                        nc_file_name(vn)
                    )
                 )
             ).variables[vn]

    
    def gm_func (var,time_slice=slice(None,None,None)):
        return gh.global_mean_var_with_resampled_mask(
            template=template,
            ctr=make_model_coord_transforms(), 
            itr=make_model_index_transforms(),
            lats=lats,
            lons=lons,
            var=var,
            time_slice=time_slice
        )

    arr_dict = {
        **{vn: gm_func(var(vn), time_slice)
            for vn in ['cLitter', 'cSoil', 'cVeg', 'cRoot']
        }, 
        **{vn: gm_func(var(vn))*86400 # kg/m2/s kg/m2/day;
            for vn in  ['npp', 'rh']
        } 
    }
    #from IPython import embed;embed()
    return arr_dict


def get_global_mean_vars(dataPath, targetPath=None, flash_cache=False):
    if targetPath is None:
        targetPath = dataPath

    arr_dict= gh.cached_var_dict(
        dataPath,
        targetPath,
        nc_global_mean_file_name,
        compute_global_mean_arr_var_dict,
        names=Observables._fields + Drivers._fields,
        flash_cache=flash_cache
    )
    obs = Observables(*(arr_dict[k] for k in Observables._fields))
    dvs = Drivers(*(arr_dict[k] for k in Drivers._fields))
    return (obs, dvs)

def get_global_mean_vars_2(conf_dict, targetPath=None):
    dataPath = Path(conf_dict["dataPath"]) 
    try:
        return get_global_mean_vars(dataPath, targetPath)    
    except FileNotFoundError:
        download_my_TRENDY_output(conf_dict)
        return get_global_mean_vars(dataPath, targetPath)    


#def make_StartVector(mvs):
#    return namedtuple(
#        "StartVector",
#        [str(v) for v in mvs.get_StateVariableTuple()]+
#        ["rh"]
#    ) 




def make_da_iterator(
        mvs,
        X_0, #: StartVector,
        par_dict,
        func_dict,
        delta_t_val=1 # defaults to 1day timestep
    ):
    # this function has to be modelspecific because some models
    # have more observables 
    mit = gh.rh_iterator(
            mvs,
            X_0,
            par_dict,
            func_dict,
            delta_t_val
    )
    # it already computes 'X', 't', 'it', 'B', 'I', 'rh', 'cVeg', 'cSoil'
    # mit.cur._fields
    # but we have to add more things we want to compare to observations 
    def numfunc(expr):
        # little helper function to compute any symbolic expression that
        # contains statevariables or time For our simple variables which are
        # just sums we could work on the results directly but this is actually
        # more easy to generalize
        return hr.numerical_func_of_t_and_Xvec(
            state_vector=mvs.get_StateVariableTuple(),
            time_symbol=mvs.get_TimeSymbol(),
            expr=expr,
            parameter_dict=par_dict,
            func_dict=func_dict,
        )

    #from IPython import embed; embed()
    var(
        """C_abvstrlit,C_abvmetlit, C_belowstrlit, C_belowmetlit, C_root,
        C_surface_microbe, C_soil_microbe, C_slowsom, C_passsom"""
    ) 
    #create the functions
    fd = OrderedDict({
        "cLitter":  numfunc(
             C_abvstrlit + C_abvmetlit + C_belowstrlit + C_belowmetlit
        ), 
        "cSoil": numfunc(C_surface_microbe + C_soil_microbe + C_slowsom + C_passsom), # redifine cSoil to keep mvs consitent with rh
        "cRoot": numfunc(C_root),
    })    
    def make_func(key):
        return lambda t,X: fd[key](t,X)

    present_step_funcs = OrderedDict(
        {
            key: make_func(key) 
            for key in fd.keys()
        }
    )
    mit.add_present_step_funcs(present_step_funcs)
    return mit

def synthetic_observables(
        mvs,
        X_0,
        par_dict,
        func_dict,
        dvs
    ):
    # - create and run the da iterator and 
    # - project the result to the size of the  osbvservables 
    #   (model specifically some variables are 
    #   yearly others monthly)
    #
    # called by the different param2res functions
    # of the  different da schemes which differ in the
    # way they produce parameters for the forward run
    # but not in how to perform it and project it to
    # the shape of the observables
    dpm = h.date.days_per_month 
    steps_per_month = 15
    delta_t_val = dpm/steps_per_month 
    bitr = ArrayDictResult(
        make_da_iterator(
            mvs,
            X_0,
            par_dict=par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
        )
    )
    number_of_months=len(dvs.npp)

    number_of_steps = number_of_months * steps_per_month
    result_dict = bitr[0: number_of_steps: steps_per_month]
    yearly_partitions = gh.partitions(0, len(result_dict['rh']), 12)
    yearly_averages = {
        key: gh.averaged_1d_array(result_dict[key],yearly_partitions)
        for key in ["cVeg", "cRoot", "cLitter", "cSoil"]
    }

    return Observables(
        cVeg=yearly_averages["cVeg"],
        cRoot=yearly_averages["cRoot"],
        cLitter=yearly_averages["cLitter"],
        cSoil=yearly_averages["cSoil"],
        rh=result_dict["rh"]#*(60*60*24)
    )


def make_weighted_cost_func(
        obs: Observables
    ) -> Callable[[Observables],np.float64]:
    def costfunction(out_simu: np.ndarray) ->np.float64:
        stream_cost = partial(gh.single_stream_cost,obs,out_simu)
        J_new = (
            stream_cost("cVeg")
            + 2*stream_cost("cRoot")
            + 2*stream_cost("cLitter")
            + 2*stream_cost("cSoil")
            + 24*stream_cost("rh")
            #+ stream_cost("ra")
        ) 
        return J_new
        
    return costfunction


def nc_file_name(nc_var_name, experiment_name="SDGVM_S2_"):
    return experiment_name+"{}.nc".format(nc_var_name)


def nc_global_mean_file_name(nc_var_name, experiment_name="SDGVM_S2_"):
    return experiment_name+"{}_gm.nc".format(nc_var_name)


def monthly_recording_times_nc_var():
    ds=nc.Dataset(str(Path(gh.confDict(Path(__file__).parent.name)['dataPath']).joinpath("SDGVM_S2_npp.nc")))
    times = ds.variables["time"]
    # times.units #claims hours since 1900-01-01 00:00:00
    # we interpret this as "trendy hours" with 
    # 1 y=360 trendy_days
    # 1 trendy_day = 24 trendy_hours
    # which is obviously different from the SI hour but
    # consistent with the assumption that 12 of the equidistant
    # values span a year which is supported by the annual
    # periodicity of the data.
    # accordingly we return the values multiplied by 24
    return times.__array__()/24

    #tm = times[0] #time of first observation in Months_since_1900-01 # print(times.units)
def start_dt():
    ## this function is important to syncronise our results
    ## because our data streams start at different times the first day of 
    ## a simulation day_ind=0 refers to different dates for different models
    ## we have to check some assumptions on which this calculation is based
    ## Here is how to get these values
    #ds=nc.Dataset(str(Path(gh.confDict(Path(__file__).parent.name)['dataPath']).joinpath("SDGVM_S2_npp.nc")))
    #times = ds.variables["time"]
    #print(times.units)
    #in days since_1900-01-01 
    return dt.datetime(
        year=1900, 
        month=1,
        day=16
    )

data_str = namedtuple( # data streams available in the model
        'data_str',
        ["cVeg", "cLitter", "cRoot", "cSoil", "gpp", "npp", "ra", "rh"]
        )

def lats_lons():
    conf_dict = gh.confDict(
        Path(__file__).parent.name
    )
    ds=nc.Dataset(Path(conf_dict["dataPath"]).joinpath("SDGVM_S2_cRoot.nc"))    
    lats=ds.variables["latitude"][:]
    lons=ds.variables["longitude"][:]
    return lats.data, lons.data

def n_months():
    mp=Path(__file__).parent
    mf=mp.name
    data_path=Path(gh.confDict(mf)["dataPath"])
    target_path = mp.joinpath("global_means")
    dvs,mvs=get_global_mean_vars(
        data_path,
        target_path
    )    
    return len(dvs.rh)
