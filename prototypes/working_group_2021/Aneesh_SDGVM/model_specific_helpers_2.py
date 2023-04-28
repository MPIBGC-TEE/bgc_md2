import sys
import json 
from pathlib import Path
from collections import namedtuple 
import netCDF4 as nc
import numpy as np
from sympy import Symbol, symbols 
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.ArrayDictResult import ArrayDictResult
from typing import Callable
import general_helpers as gh
from functools import reduce, partial
from importlib import import_module
from collections import OrderedDict

model_mod = 'bgc_md2.models.AneeshSDGVM'
cp_mod = import_module(f"{model_mod}.CachedParameterization")
make_func_dict = import_module(f"{model_mod}.CachedParameterization").make_func_dict
Drivers=cp_mod.Drivers
CachedParameterization=cp_mod.CachedParameterization

sys.path.insert(0,'..') # necessary to import general_helpers
import general_helpers as gh

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

Constants = namedtuple(
    "Constants",
    [
        "cLitter_0",
        "cSoil_0",
        "cRoot_0",
        "cVeg_0",
        "npp_0",
        "rh_0",
        "number_of_months" # necessary to prepare the output in the correct lenght
    ]
)

EstimatedParameters = namedtuple(
    "EstimatedParameters",[
         'beta_leaf',
         'beta_wood',
         #beta_root,
         'r_C_leaf2abvstrlit',
         'r_C_abvmetlit2surface_microbe',
         'r_C_abvstrlit2slowsom',
         'r_C_abvstrlit2surface_microbe',
         'r_C_belowmetlit2soil_microbe',
         'r_C_belowstrlit2slowsom',
         'r_C_belowstrlit2soil_microbe',
         #'r_C_leached',
         'r_C_leaf2abvmetlit',
         'r_C_passsom2soil_microbe',
         'r_C_root2belowmetlit',
         'r_C_root2belowstrlit',
         'r_C_slowsom2passsom',
         'r_C_slowsom2soil_microbe',
         'r_C_soil_microbe2passsom',
         'r_C_soil_microbe2slowsom',
         'r_C_surface_microbe2slowsom',
         'r_C_wood2abvmetlit',
         'r_C_wood2abvstrlit',
         'r_C_abvstrlit_rh',
         'r_C_abvmetlit_rh',
         'r_C_belowstrlit_rh',
         'r_C_belowmetlit_rh',
         'r_C_surface_microbe_rh',
         'r_C_slowsom_rh',
         'r_C_passsom_rh',
         'r_C_soil_microbe_rh',
         'C_leaf_0',
        #'C_root_0',
         'C_abvstrlit_0',
         'C_abvmetlit_0',
         'C_blwstrlit_0',
         'C_surfacemic_0',
         'C_soilmic_0',
         'C_slow_0'
    ]
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
        #flash_cash=True
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


def make_StartVector(mvs):
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()]+
        ["rh"]
    ) 



def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers
) -> Callable[[np.ndarray], np.ndarray]: 
    
    assert(cpa.number_of_months <= dvs.npp.shape[0])
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        X_0 = numeric_X_0(mvs, dvs, cpa, epa)
        dpm=30
        steps_per_month = 2
        delta_t_val = dpm/steps_per_month 

        par_dict = gh.make_param_dict(mvs, cpa, epa)
        func_dict = make_func_dict(dvs , cpa=cpa, epa=epa)
        bitr = ArrayDictResult(
            make_da_iterator(
                mvs,
                X_0,
                par_dict=par_dict,
                func_dict=func_dict,
                delta_t_val=delta_t_val
            )
        )
        number_of_steps = int(cpa.number_of_months/delta_t_val)
        steps_per_month = int(dpm / delta_t_val)
        result_dict = bitr[0: number_of_steps: steps_per_month]
        steps_per_year = steps_per_month*12
        yearly_partitions = gh.partitions(0, number_of_steps, steps_per_year)
        yearly_averages = {
            key: gh.averaged_1d_array(result_dict[key],yearly_partitions)
            for key in ["cVeg", "cRoot", "cLitter", "cSoil"]
        }

        return Observables(
            cVeg=yearly_averages["cVeg"],
            cRoot=yearly_averages["cRoot"],
            cLitter=yearly_averages["cLitter"],
            cSoil=yearly_averages["cSoil"],
            rh=result_dict["rh"]#/(60*60*24)
        )
    return param2res


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
    C_abvstrlit, C_abvmetlit, C_belowstrlit, C_belowmetlit, C_root = symbols(
        "C_abvstrlit,C_abvmetlit, C_belowstrlit, C_belowmetlit, C_root"
    ) 
    #create the functions
    fd = {
        "cLitter":  numfunc(
             C_abvstrlit + C_abvmetlit + C_belowstrlit + C_belowmetlit
        ), 
        "cRoot": numfunc(C_root)
    }    
    present_step_funcs = OrderedDict(
        {
            key: lambda t,X: fd[key](t,X)
            for key in fd.keys()
        }
    )
    mit.add_present_step_funcs(present_step_funcs)
    return mit

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
            + 100.0*stream_cost("rh")
            #+ stream_cost("ra")
        ) 
        return J_new
        
        #number_of_ys=out_simu.cVeg.shape[0]
        #number_of_ms=out_simu.rh.shape[0]

        #J_obj1 = np.mean (( out_simu.cVeg - obs.cVeg )**2)/(2*np.var(obs.cVeg))
        #J_obj2 = np.mean (( out_simu.cRoot - obs.cRoot )**2)/(2*np.var(obs.cRoot))
        #J_obj3 = np.mean (( out_simu.cLitter - obs.cLitter )**2)/(2*np.var(obs.cLitter))
        #J_obj4 = np.mean (( out_simu.cSoil -  obs.cSoil )**2)/(2*np.var(obs.cSoil))

        #J_obj5 = np.mean (( out_simu.rh - obs.rh )**2)/(2*np.var(obs.rh))

        #J_new = (J_obj1 + J_obj2 + J_obj3 + J_obj4)/200 + J_obj5/4
        ## to make this special costfunction comparable (in its effect on the
        ## acceptance rate) to the general costfunction proposed by Feng we
        ## rescale it by a factor
        #return J_new*400
    return costfunction


def make_param_filter_func(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters, 
        cpa: Constants,
        ) -> Callable[[np.ndarray], bool]:

    # find position of beta_leaf and beta_wood
        
    
    def isQualified(c):
        def value(field_name):
            try:
                return c[EstimatedParameters._fields.index(field_name)]
            except Exception as e:
                print("###########################")
                print(e)
                print(field_name)
                raise e
        conds=[
            (c >= c_min).all(), 
            (c <= c_max).all(), 
            sum(map(value, ["beta_leaf", "beta_wood"])) <= 0.99,
            value("C_leaf_0") <= cpa.cVeg_0-cpa.cRoot_0, 
            sum(
                map(
                    value,
                    ["C_abvstrlit_0","C_abvmetlit_0","C_blwstrlit_0"]
                )
            ) <= cpa.cVeg_0,
            sum(
                map(
                    value,
                    ["C_surfacemic_0", "C_soilmic_0", "C_slow_0"]
                )
            ) <= cpa.cSoil_0  
        ]    
        res=all(conds)
        if not res:
            print(conds)
        return res
        
    return isQualified

def numeric_X_0(mvs,dvs,cpa,epa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "C_leaf": apa['C_leaf_0'],     
        "C_root": apa['cRoot_0'],     
        "C_wood": apa['cVeg_0'] - (apa['C_leaf_0'] +  apa['cRoot_0']),  
        "C_abvstrlit": apa['C_abvstrlit_0'],
        "C_abvmetlit": apa['C_abvmetlit_0'],
        "C_belowstrlit": apa["C_blwstrlit_0"],
        "C_belowmetlit":apa["cLitter_0"]- apa["C_abvstrlit_0"] - apa["C_abvmetlit_0"] - apa["C_blwstrlit_0"],
        "C_surface_microbe":apa["C_surfacemic_0"],
        "C_soil_microbe": apa["C_soilmic_0"],
        "C_slowsom":apa["C_slow_0"],
        "C_passsom":apa["cSoil_0"] - apa["C_surfacemic_0"] - apa["C_soilmic_0"] - apa["C_slow_0"]
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    return X_0


def nc_file_name(nc_var_name, experiment_name="SDGVM_S2_"):
    return experiment_name+"{}.nc".format(nc_var_name)


def nc_global_mean_file_name(nc_var_name, experiment_name="SDGVM_S2_"):
    return experiment_name+"{}_gm.nc".format(nc_var_name)


def start_date():
    ## this function is important to syncronise our results
    ## because our data streams start at different times the first day of 
    ## a simulation day_ind=0 refers to different dates for different models
    ## we have to check some assumptions on which this calculation is based
    ## Here is how to get these values
    #ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("SDGVM_S2_npp.nc")))
    #times = ds.variables["time"]
    #tm = times[0] #time of first observation in Months_since_1900-01 # print(times.units)
    #td = ts/24  #in days since_1900-01-01 
    #import datetime as dt
    #ad = dt.date(1, 1, 1) # first of January of year 1 
    #sd = dt.date(1900, 1, 1)
    #td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    ## from td_aD (days since 1-1-1) we can compute the year month and day
    return gh.date(
        year=1900, 
        month=1,
        day=16
    )

data_str = namedtuple( # data streams available in the model
        'data_str',
        ["cVeg", "cLitter", "cRoot", "cSoil", "gpp", "npp", "ra", "rh"]
        )


# fixme mm 7-4 2023:
def get_global_mean_vars_all(experiment_name):
    print("""
    deprecation warning:
    model_folder="Aneesh_SDGVM" is self referential and would
    break if we rename the model folder. This is a design flaw 
    If there is any model specific information it should be computed
    by model specific functions. 
    """
    )
    return(
        gh.get_global_mean_vars_all(model_folder="Aneesh_SDGVM", 
                            experiment_name=experiment_name,
                            lat_var="latitude",
                            lon_var="longitude",
                            ) 
        )

