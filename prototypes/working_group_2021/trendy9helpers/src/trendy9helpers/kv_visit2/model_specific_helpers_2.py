import sys
import json
import inspect
import os
import netCDF4 as nc
import numpy as np
import datetime as dt
import pandas as pd
from pathlib import Path
from collections import namedtuple, OrderedDict
from sympy import Symbol, symbols
from importlib import import_module
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.ArrayDictResult import ArrayDictResult
from typing import Dict,Tuple, Callable
from functools import lru_cache, reduce, partial
from copy import copy
from typing import Callable
from scipy.interpolate import interp1d

from ComputabilityGraphs.CMTVS import  CMTVS
import CompartmentalSystems.helpers_reservoir as hr
from CompartmentalSystems.ArrayDict import ArrayDict
import bgc_md2.helper as h
from .. import general_helpers as gh

model_mod = 'bgc_md2.models.kv_visit2'
cp_mod=import_module(f"{model_mod}.CachedParameterization")

mvs=import_module(f"{model_mod}.source").mvs

# some classes (not instances) as arguments for 
# functions that can produce output of different type
Drivers=cp_mod.Drivers 
OrgDrivers=cp_mod.OrgDrivers 
CachedParameterization=cp_mod.CachedParameterization

#from bgc_md2.models.kv_visit2.source import  mvs
#from bgc_md2.models.kv_visit2.CachedParameterization import  Drivers, OrgDrivers

def spatial_mask(dataPath) -> "CoorMask":
    mask = (
        nc.Dataset(dataPath.joinpath("VISIT_S2_cSoil.nc"))
        .variables["cSoil"][0, :, :]
        .mask
    )
    sym_tr = gh.SymTransformers(
        itr=make_model_index_transforms(), ctr=make_model_coord_transforms()
    )
    return gh.CoordMask(mask, sym_tr)


def make_model_coord_transforms():
    """This function can is used to achieve a target grid LAT,LON with
    - LAT ==   0 at the equator with
    - LAT == -90 at the south pole,
    - LAT== +90 at the north pole,
    - LON ==   0 at Greenich and
    - LON is counted positive eastwards from -180 to 180
    """
    return gh.CoordTransformers(
        lat2LAT=lambda lat: lat,
        LAT2lat=lambda LAT: LAT,
        lon2LON=lambda lon: lon,
        LON2lon=lambda LON: LON,
    )


def make_model_index_transforms():
    return gh.transform_maker(
        lat_0=89.75,
        lon_0=-179.75,
        step_lat=-0.5,
        step_lon=0.5,
    )


# we will use the trendy output names directly in other parts of the output
Observables = namedtuple("Observables", ["cVeg", "cLitter", "cSoil", "rh"])  # ,"ra"]

# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated)
# parameters. In this case we may use only the first entry e.g. to derive startvalues.
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise
# to guard agains accidentally changed order we use a namedtuple again. Since B_func and u_func rely
# on the correct ordering of the statevariables we build V dependent on this order

# create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(
            conf_dict["dataPath"]
        ),  # platform independent path desc. (Windows vs. linux)
        models=["VISIT"],
        variables=Observables._fields + OrgDrivers._fields,
    )


def get_example_site_vars(dataPath):
    # pick up 1 site
    s = slice(None, None, None)  # this is the same as :
    t = s, 50, 33  # [t] = [:,49,325]

    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them
        ds = nc.Dataset(str(path))
        return ds.variables[vn][t]

    o_names = [(f, "VISIT_S2_{}.nc".format(f)) for f in Observables._fields]
    d_names = [(f, "VISIT_S2_{}.nc".format(f)) for f in OrgDrivers._fields]

    # we want to drive with npp and can create it from gpp and ra
    # observables
    odvs = OrgDrivers(*map(f, d_names))
    obss = Observables(*map(f, o_names))

    dvs = Drivers(
        npp=odvs.gpp - odvs.ra,
        mrso=odvs.mrso,
        tas=odvs.tas  # ,
        # xi=odvs.xi_t*odvs.xi_w
    )
    return (obss, dvs)


def nc_file_name(nc_var_name, experiment_name="VISIT_S2_"):
    return experiment_name + "{}.nc".format(nc_var_name)


def nc_global_mean_file_name(nc_var_name, experiment_name="VISIT_S2_"):
    return experiment_name + "{}_gm.nc".format(nc_var_name)


def compute_global_mean_arr_var_dict(dataPath):
    ds = nc.Dataset(dataPath.joinpath("VISIT_S2_gpp.nc"))
    vs = ds.variables
    lats = vs["lat"].__array__()
    lons = vs["lon"].__array__()
    template = vs['gpp'][0,:,:].mask

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

    
    scaled=["gpp","rh","ra"]
    arr_dict = {
        **{vn: gm_func(var(vn) )
            for vn in set(Observables._fields + OrgDrivers._fields).difference(scaled)
        }, 
        **{vn: gm_func(var(vn))*86400 # kg/m2/s kg/m2/day;
            for vn in scaled 
        } 
    }
    return arr_dict

def get_global_mean_vars(dataPath,targetPath=None,flash_cache=False):
    if targetPath is None:
        targetPath = dataPath

    arr_dict = gh.cached_var_dict(
        dataPath,
        targetPath,
        nc_global_mean_file_name,
        compute_global_mean_arr_var_dict,
        names=Observables._fields + OrgDrivers._fields,
        flash_cache=flash_cache
    )
    obs = Observables(*(arr_dict[k] for k in Observables._fields))
    odvs = OrgDrivers(*(arr_dict[k] for k in OrgDrivers._fields))
    dvs = Drivers(npp=odvs.gpp - odvs.ra, mrso=odvs.mrso, tas=odvs.tas)
    
    return (obs,dvs)


def get_global_mean_vars_2(conf_dict,targetPath=None):
    o_names = Observables._fields
    d_names = OrgDrivers._fields
    names = o_names + d_names
    dataPath = Path(conf_dict["dataPath"])
    if not(all([dataPath.joinpath(nc_file_name(vn)).exists() for vn in names])):
        download_my_TRENDY_output(conf_dict)
    return get_global_mean_vars(dataPath,targetPath)    


def ivp_from_cache_dir(path):
    par_dict = gh.load_dict_from_json_path(dirPath.joinpath("parameter_dict.json"))
    func_dict = load_funcdict(path)
    #x0 =load(
    return 

def key2path(path,k):
    return path.joinpath(f"{k}.nc")

def make_da_iterator(
        mvs,
        X_0, #: StartVector,
        par_dict,
        func_dict,
        delta_t_val=1 # defaults to 1day timestep
    ):
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

    C_leaf_litter, C_wood_litter, C_root_litter = symbols(
        "C_leaf_litter, C_wood_litter, C_root_litter"
    )
    #create the functions
    fd = {
        "cLitter":  numfunc(C_leaf_litter + C_wood_litter + C_root_litter) 
    }    
    present_step_funcs = OrderedDict(
        {
            key: lambda t,X: fd[key](t,X)
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
    steps_per_month = 2
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

    return Observables(
        cVeg=result_dict["cVeg"],
        cLitter=result_dict["cLitter"],
        cSoil=result_dict["cSoil"],
        rh=result_dict["rh"]
    )

def make_weighted_cost_func(
        obs: Observables
    ) -> Callable[[Observables],np.float64]:
    
    # Define cost function 
    def costfunction(out_simu: Observables) -> np.float64:
        stream_cost = partial(gh.single_stream_cost,obs,out_simu)
        # sum costs with possibly manually ajusted weights
        J_new = (
            stream_cost("cVeg")
            + stream_cost("cLitter")
            + stream_cost("cSoil")
            + stream_cost("rh")
            #+ stream_cost("ra")
        ) 
        return J_new

    return costfunction



def monthly_recording_times_nc_var()->nc.Variable:
    ds=nc.Dataset(str(Path(gh.confDict(Path(__file__).parent.name)['dataPath']).joinpath("VISIT_S2_gpp.nc")))
    times = ds.variables["time"]
    # times.units #claims months since 1860-01-01
    # we interpret this as months of equal length
    # consistent with the assumption that 12 of the equidistant
    # values span a year which is supported by the annual
    # periodicity of the data.
    return times

def start_dt()->dt.datetime:
    # this function is important to syncronise our results because the trendy9
    # data streams start at different times the first day of a simulation
    # day_ind=0 refers to different dates for different models we have to check
    # some assumptions on which this calculation is based Here is how to get
    # these values
    #ds=nc.Dataset(str(Path(gh.confDict(Path(__file__).parent.name)['dataPath']).joinpath("VISIT_S2_gpp.nc")))
    #times = ds.variables["time"]
    # tm = times[0] #time of first observation in Months_since_1860-01 # print(times.units)
    # td = int(tm *30)  #in days since_1860-01-01
    return dt.datetime(year=1860, month=1, day=1)


data_str = namedtuple("data_str", ["cVeg", "cLitter", "cSoil", "gpp", "ra", "rh"])

def lats_lons():
    conf_dict = gh.confDict(
        Path(__file__).parent.name
    ) 
    ds=nc.Dataset(Path(conf_dict["dataPath"]).joinpath("VISIT_S2_cLitter.nc"))    
    lats=ds.variables["lat"][:]
    lons=ds.variables["lon"][:]
    return lats.data, lons.data


def n_months():
    mp = Path(__file__).parent
    mf = mp.name
    data_path=Path(gh.confDict(mf)["dataPath"])
    target_path = mp.joinpath("global_means")
    svs, dvs =get_global_mean_vars(
        data_path,
        target_path
    )
    return len(dvs.npp)
