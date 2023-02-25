import sys
import json
from pathlib import Path
from collections import namedtuple
import netCDF4 as nc
import numpy as np
from sympy import Symbol
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
    TimeStepIterator2,
)
from typing import Dict
from functools import lru_cache
from copy import copy
from typing import Callable
from functools import reduce
from collections import OrderedDict
from scipy.interpolate import interp1d

from ComputabilityGraphs.CMTVS import  CMTVS
import CompartmentalSystems.helpers_reservoir as hr
from CompartmentalSystems.ArrayDict import ArrayDict

sys.path.insert(0, "..")  # necessary to import general_helpers
import general_helpers as gh


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
OrgDrivers = namedtuple("OrgDrivers", ["gpp", "ra", "mrso", "tas"])
Drivers = namedtuple("Drivers", ("npp",) + OrgDrivers._fields[2:])
# As a safety measure we specify those parameters again as 'namedtuples', which are like a mixture of dictionaries and tuples
# They preserve order as numpy arrays which is great (and fast) for the numeric computations
# and they allow to access values by keys (like dictionaries) which makes it difficult to accidentally mix up values.

Constants = namedtuple(
    "Constants",
    [
        "cLitter_0",
        "cSoil_0",
        "cVeg_0",
        "npp_0",
        "rh_0",
        # "ra_0",
        "number_of_months",  # necessary to prepare the output in the correct lenght
    ],
)
# Note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated)
# parameters. In this case we may use only the first entry e.g. to derive startvalues and their length (number of months)
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise
# It is better to start with only a few

EstimatedParameters = namedtuple(
    "EstimatedParameters",
    [
        "beta_leaf",
        "beta_wood",
        "T_0",
        "E",
        "KM",
        "r_C_leaf_litter_rh",
        "r_C_wood_litter_rh",
        "r_C_root_litter_rh",
        "r_C_soil_fast_rh",
        "r_C_soil_slow_rh",
        "r_C_soil_passive_rh",
        "r_C_leaf_2_C_leaf_litter",
        "r_C_wood_2_C_wood_litter",
        "r_C_root_2_C_root_litter",
        "r_C_leaf_litter_2_C_soil_fast",
        "r_C_leaf_litter_2_C_soil_slow",
        "r_C_leaf_litter_2_C_soil_passive",
        "r_C_wood_litter_2_C_soil_fast",
        "r_C_wood_litter_2_C_soil_slow",
        "r_C_wood_litter_2_C_soil_passive",
        "r_C_root_litter_2_C_soil_fast",
        "r_C_root_litter_2_C_soil_slow",
        "r_C_root_litter_2_C_soil_passive",
        "C_leaf_0",  # for the trendy data also the startvalues have to be estimated but
        "C_wood_0",
        # C_root_0 can be inferred as cVeg_0-(C_leaf_0+C_wood_0)
        "C_leaf_litter_0",
        "C_wood_litter_0",
        # C_root_litter_0 can be inferred
        "C_soil_fast_0",
        "C_soil_slow_0",
        # C_soil_passive_0 can be inferred
    ],
)
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




def compute_and_write_global_mean_var(dataPath, targetPath, index_mask, vn):
    path = dataPath.joinpath(nc_file_name(vn))
    ds = nc.Dataset(str(path))
    vs = ds.variables
    lats = vs["lat"].__array__()
    lons = vs["lon"].__array__()
    print(vn)
    var = ds.variables[vn]
    # check if we have a cached version (which is much faster)
    gm_path = targetPath.joinpath(nc_global_mean_file_name(vn))

    gm = gh.global_mean_var(
        lats,
        lons,
        index_mask,
        # combined_mask,
        var,
    )
    gh.write_timeline(gm_path, gm, vn)
    return gm * 86400 if vn in ["gpp", "rh", "ra"] else gm

def compute_and_write_global_mean_vars(dataPath,targetPath):
    o_names = Observables._fields
    d_names = OrgDrivers._fields
    
    # load an example file with mask
    gm = gh.globalMask()
    template = (
        nc.Dataset(dataPath.joinpath("VISIT_S2_cSoil.nc"))
        .variables["cSoil"][0, :, :]
        .mask
    )
    gcm = gh.project_2(
        source=gm,
        target=gh.CoordMask(
            index_mask=np.zeros_like(template),
            tr=gh.SymTransformers(
                ctr=make_model_coord_transforms(), itr=make_model_index_transforms()
            ),
        ),
    )

    print("computing means, this may take some minutes...")
    def f(name):
        return compute_and_write_global_mean_var(
            dataPath, 
            targetPath, 
            gcm.index_mask, 
            name
        )

    odvs = OrgDrivers(*map(f, d_names))
    obss = Observables(*map(f, o_names))
    dvs = Drivers(npp=odvs.gpp - odvs.ra, mrso=odvs.mrso, tas=odvs.tas)

    return (obss, dvs)


def get_cached_global_mean_drivers(targetPath):
    o_names = Observables._fields
    d_names = OrgDrivers._fields
    def f(vn):
        gm = gh.get_nc_array(
            targetPath.joinpath(nc_global_mean_file_name(vn)), vn
        )
        return gm * 86400 if vn in ["gpp", "rh", "ra"] else gm

    # map variables to data
    odvs = OrgDrivers(*map(f, d_names))
    dvs = Drivers(npp=odvs.gpp - odvs.ra, mrso=odvs.mrso, tas=odvs.tas)
    return dvs

def get_cached_global_mean_observables(targetPath):
    o_names = Observables._fields
    def f(vn):
        gm = gh.get_nc_array(
            targetPath.joinpath(nc_global_mean_file_name(vn)), vn
        )
        return gm * 86400 if vn in ["gpp", "rh", "ra"] else gm

    # map variables to data
    obss = Observables(*map(f, o_names))
    return obss

def get_cached_global_mean_vars(targetPath):
    return (
        get_cached_global_mean_observables(targetPath), 
        get_cached_global_mean_drivers(targetPath)
    )    

def get_global_mean_vars(dataPath,targetPath=None):
    if targetPath is None:
        targetPath = dataPath 
    
    o_names = Observables._fields
    d_names = OrgDrivers._fields
    names = o_names + d_names
        

    if all([targetPath.joinpath(nc_global_mean_file_name(vn)).exists() for vn in names]):
        print(
            """ Found cached global mean files. If you want to recompute the global means
            remove the following files: """
        )
        for vn in names:
            print(targetPath.joinpath(nc_global_mean_file_name(vn)))

        return get_cached_global_mean_vars(targetPath)

    else:
        dataPath = Path(conf_dict["dataPath"])
        return msh.compute_and_write_global_mean_vars(dataPath, targetPath)
        

def get_global_mean_vars_2(conf_dict,targetPath=None):
    o_names = Observables._fields
    d_names = OrgDrivers._fields
    names = o_names + d_names
    dataPath = Path(conf_dict["dataPath"])
    if not(all([dataPath.joinpath(nc_file_name(vn)).exists() for vn in names])):
        download_my_TRENDY_output(conf_dict)
    return get_global_mean_vars(dataPath,targetPath)    

def make_func_dict(dvs,**kwargs):
    f_d = {
        "TAS": gh.make_interpol_of_t_in_days(dvs.tas),
        "mrso": gh.make_interpol_of_t_in_days(dvs.mrso),
        "NPP": gh.make_interpol_of_t_in_days(dvs.npp),
    }
    return f_d

# deprecated
def make_func_dict_old(mvs, dvs, cpa, epa):
    def npp_func(day):
        month = gh.day_2_month_index(day)
        # kg/m2/s kg/m2/day;
        return dvs.npp[month]  # * 86400

    def tas_num(day):
        month = gh.day_2_month_index(day)
        return dvs.tas[month]

    def mrso_num(day):
        month = gh.day_2_month_index(day)
        return dvs.mrso[month]

    f_d = {
        "TAS": tas_num,
        "mrso": mrso_num,
        "NPP": npp_func,
        # "xi": make_xi_func(dvs, epa)
    }
    return f_d


# def make_traceability_iterator(mvs,dvs,cpa,epa):
#    par_dict={
#    Symbol(k): v for k,v in {
#        "beta_leaf": epa.beta_leaf,
#        "beta_wood": epa.beta_wood,
#        "T_0": epa.T_0,
#        "E": epa.E,
#        "KM": epa.KM,
#        "r_C_leaf_litter_rh": epa.r_C_leaf_litter_rh,
#        "r_C_wood_litter_rh": epa.r_C_wood_litter_rh,
#        "r_C_root_litter_rh": epa.r_C_root_litter_rh,
#        "r_C_soil_fast_rh": epa.r_C_soil_fast_rh,
#        "r_C_soil_slow_rh": epa.r_C_soil_slow_rh,
#        "r_C_soil_passive_rh": epa.r_C_soil_passive_rh,
#        "r_C_leaf_2_C_leaf_litter": epa.r_C_leaf_2_C_leaf_litter,
#        "r_C_wood_2_C_wood_litter": epa.r_C_wood_2_C_wood_litter,
#        "r_C_root_2_C_root_litter": epa.r_C_root_2_C_root_litter,
#        "r_C_leaf_litter_2_C_soil_fast": epa.r_C_leaf_litter_2_C_soil_fast,
#        "r_C_leaf_litter_2_C_soil_slow": epa.r_C_leaf_litter_2_C_soil_slow,
#        "r_C_leaf_litter_2_C_soil_passive": epa.r_C_leaf_litter_2_C_soil_passive,
#        "r_C_wood_litter_2_C_soil_fast": epa.r_C_wood_litter_2_C_soil_fast,
#        "r_C_wood_litter_2_C_soil_slow": epa.r_C_wood_litter_2_C_soil_slow,
#        "r_C_wood_litter_2_C_soil_passive": epa.r_C_wood_litter_2_C_soil_passive,
#        "r_C_root_litter_2_C_soil_fast": epa.r_C_root_litter_2_C_soil_fast,
#        "r_C_root_litter_2_C_soil_slow": epa.r_C_root_litter_2_C_soil_slow,
#        "r_C_root_litter_2_C_soil_passive": epa.r_C_root_litter_2_C_soil_passive
#    }.items()
# }
#    X_0_dict={
#        "C_leaf": epa.C_leaf_0,
#        "C_wood": epa.C_wood_0,
#        "C_root": cpa.cVeg_0-(epa.C_leaf_0 + epa.C_wood_0),
#        "C_leaf_litter": epa.C_leaf_litter_0,
#        "C_wood_litter": epa.C_wood_litter_0,
#        "C_root_litter": cpa.cLitter_0-(epa.C_leaf_litter_0 + epa.C_wood_litter_0),
#        "C_soil_fast": epa.C_soil_fast_0,
#        "C_soil_slow": epa.C_soil_slow_0,
#        "C_soil_passive": cpa.cSoil_0-(epa.C_soil_fast_0 + epa.C_soil_slow_0)
#    }
#    X_0= np.array(
#        [
#            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
#        ]
#    ).reshape(9,1)
#    fd = make_func_dict(dvs, cpa, epa)
#    V_init = gh.make_InitialStartVectorTrace(
#            X_0,mvs,
#            par_dict=par_dict,
#            func_dict=fd
#    )
#    it_sym_trace = gh.make_daily_iterator_sym_trace(
#        mvs,
#        V_init=V_init,
#        par_dict=par_dict,
#        func_dict=fd
#    )
#    return it_sym_trace

# We now build the essential object to run the model forward. We have a
# - startvector $V_0$ and
# - a function $f$ to compute the next value: $V_{it+1} =f(it,V_{it})$
#   the dependence on the iteration $it$ allows us to represent drivers that
#   are functions of time
#
# So we can compute all values:
#
# $V_1=f(0,V_0)$
#
# $V_2=f(1,V_1)$
#
# ...
#
# $V_n+1=f(n,V_n)$
#
# Technically this can be implemented as an `iterator` object with a `__next__()` method to move our system one step forward in time.
#
# What we want to build is an object `it_sym` that can use as follows.
# ```python
# for i in range(10):
#     print(it_sym.__next__())
# ```
# to print out the first 10 values.
#
# If iterators had not been invented yet we would invent them now, because they capture exactly the mathematical concept of an initial value system,
# without all the nonessential technical details of e.g. how many steps we want to make or which of the results we want to store.
# This is essential if we want to test and use the iterator later in different scenarios but avoid reimplemtation of the basic timestep.
#
# Remark:
#
# If we were only interested in the timeseries of the pool contents `bgc_md2` could compute the solution automatically without the need to build an iterator ourselves.
# In our case we are also interested in tracking the autotrophic and heterotrophic respiration and therefore have to recreate and extend the code `bgc_md2` would use in the background.
# We will let `bgc_md2` do part of the work and derive numeric functions for the Compartmental matrix $B$ and the input $u$ and the Outfluxes - from which to compute $ra$ $rh$ - from our symbolic description but we will build our own iterator by combining these functions.
# We will proceed as follows:
# - create $V_0$
# - build the function $f$

# +
def make_iterator_sym(
    mvs,
    V_init,  #: StartVector,
    par_dict,
    func_dict,
    delta_t_val=1,  # defaults to 1day timestep
):
    B_func, u_func = gh.make_B_u_funcs_2(mvs, par_dict, func_dict, delta_t_val)
    sv = mvs.get_StateVariableTuple()
    # mass production of output functions

    n = len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr = np.array(
        [V_init.__getattribute__(str(v)) for v in sv]
        + [V_init.rh]
        # [V_init.ra,V_init.rh]
    ).reshape(
        n + 1, 1
    )  # reshaping is neccessary for matmul (the @ in B @ X)

    # To compute the ra and rh we have to some up the values for autotrophic and heterotrophic respiration we have to sum up outfluxes.
    # We first create numerical functions for all the outfluxes from our symbolic description.
    numOutFluxesBySymbol = {
        k: gh.numfunc(
            expr_cont,
            mvs,
            delta_t_val=delta_t_val,
            par_dict=par_dict,
            func_dict=func_dict,
        )
        for k, expr_cont in mvs.get_OutFluxesBySymbol().items()
    }

    def f(it, V):
        X = V[0:n]  # .reshape(-1)
        b = u_func(it, X)
        B = B_func(it, X)
        X_new = X + b + B @ X
        # from IPython import embed; embed()
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep

        #         ra = np.sum(
        #             [
        #               numOutFluxesBySymbol[Symbol(k)](it,*X)
        #                 for k in ["C_leaf","C_wood","C_root"]
        #                 if Symbol(k) in numOutFluxesBySymbol.keys()
        #             ]
        #         )
        rh = np.sum(
            [
                numOutFluxesBySymbol[Symbol(k)](it, *X)
                for k in [
                    "C_leaf_litter",
                    "C_wood_litter",
                    "C_root_litter",
                    "C_soil_fast",
                    "C_soil_slow",
                    "C_soil_passive",
                ]
                if Symbol(k) in numOutFluxesBySymbol.keys()
            ]
        )
        # V_new = np.concatenate((X_new.reshape(n,1),np.array([ra,rh]).reshape(2,1)), axis=0)
        V_new = np.concatenate(
            (X_new.reshape(n, 1), np.array(rh).reshape(1, 1)), axis=0
        )

        return V_new

    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


# -


def make_StartVector(mvs):
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()] + ["rh"]
        # ["ra","rh"]
    )


def make_param2res_sym(
    mvs, cpa: Constants, dvs: Drivers
) -> Callable[[np.ndarray], np.ndarray]:
    # To compute numeric solutions we will need to build and iterator
    # as we did before. As it will need numeric values for all the parameters
    # we will have to create a complete dictionary for all the symbols
    # exept those for the statevariables and time.
    # This set of symbols does not change during the mcmc procedure, since it only
    # depends on the symbolic model.
    # Therefore we create it outside the mcmc loop and bake the result into the
    # param2res function.
    # The iterator does not care if they are estimated or constant, it only
    # wants a dictionary containing key: value pairs for all
    # parameters that are not state variables or the symbol for time
    StartVector = make_StartVector(mvs)
    srm = mvs.get_SmoothReservoirModel()
    model_par_dict_keys = srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))] + list(mvs.get_StateVariableTuple())
    )

    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test
    seconds_per_day = 86400
    # def npp_func(day):
    #    month=gh.day_2_month_index(day)
    #    return dvs.npp[month] #* seconds_per_day   # kg/m2/s kg/m2/day;

    def param2res(pa):
        epa = EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
        V_init = StartVector(
            C_leaf=epa.C_leaf_0,
            C_wood=epa.C_wood_0,
            C_root=cpa.cVeg_0 - (epa.C_leaf_0 + epa.C_wood_0),
            C_leaf_litter=epa.C_leaf_litter_0,
            C_wood_litter=epa.C_wood_litter_0,
            C_root_litter=cpa.cLitter_0 - (epa.C_leaf_litter_0 + epa.C_wood_litter_0),
            C_soil_fast=epa.C_soil_fast_0,
            C_soil_slow=epa.C_soil_slow_0,
            C_soil_passive=cpa.cSoil_0 - (epa.C_soil_fast_0 + epa.C_soil_slow_0),
            # ra=cpa.ra_0,
            rh=cpa.rh_0,
        )
        # next we create the parameter dict for the iterator
        # The iterator does not care if they are estimated or not so we look for them
        # in the combination
        apa = {**cpa._asdict(), **epa._asdict()}
        model_par_dict = {
            Symbol(k): v for k, v in apa.items() if Symbol(k) in model_par_dict_keys
        }

        # print(model_par_dict)
        # from IPython import embed;embed()

        # Beside the par_dict the iterator also needs the python functions to replace the symbolic ones with
        # our fake xi_func could of course be defined outside of param2res but in general it
        # could be build from estimated parameters and would have to live here...
        # def xi_func(day):
        #    month=gh.day_2_month_index(day)
        #    TS = (dvs.tas[month]-273.15) # # convert from Kelvin to Celcius
        #    TS = 0.748*TS + 6.181 # approximate soil T at 20cm from air T (from https://doi.org/10.1155/2019/6927045)
        #    if TS > epa.T_0:
        #        xi_out = np.exp(epa.E*(1/(10-epa.T_0)-1/(TS-epa.T_0))) * dvs.mrso[month]/(epa.KM+dvs.mrso[month])
        #    else:
        #        xi_out=0
        #    return(xi_out)
        #    # return 1.0 # preliminary fake for lack of better data...
        #
        # func_dict={
        #    'NPP':npp_func,
        #     'xi':xi_func
        # }
        func_dict = make_func_dict(dvs)

        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2)
        delta_t_val = 15
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val,
        )

        # Now that we have the iterator we can start to compute.
        # the first thing to notice is that we don't need to store all values (daili yi delta_t_val=1)
        # since the observations are recorded monthly while our iterator possibly has a smaller timestep.
        # - for the pool contents we only need to record the last day of the month
        # - for the respiration(s) ra and rh we want an accumulated value (unit mass)
        #   have to sum up the products of the values*delta_t over a month
        #
        # Note: check if TRENDY months are like this...
        # days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        # sols=[]

        # empty arrays for saving data
        cVeg_arr = np.zeros(cpa.number_of_months)
        cLitter_arr = np.zeros(cpa.number_of_months)
        cSoil_arr = np.zeros(cpa.number_of_months)
        rh_arr = np.zeros(cpa.number_of_months)
        # ra_arr = np.zeros(cpa.number_of_months)
        im = 0
        dpm = 30
        steps_per_month = int(dpm / delta_t_val)
        # n=len(V_init)
        # forward simulation by year
        for y in range(int(cpa.number_of_months / 12)):
            for m in range(12):
                cVeg_avg = 0
                cLitter_avg = 0
                cSoil_avg = 0
                rh_avg = 0
                # ra_avg = 0
                for d in range(steps_per_month):
                    V = StartVector(*it_sym.__next__())
                    rh_avg += V.rh
                    # ra_avg += V.ra
                    cVeg_avg += float(V.C_leaf + V.C_wood + V.C_root)
                    cLitter_avg += float(
                        V.C_leaf_litter + V.C_wood_litter + V.C_root_litter
                    )
                    cSoil_avg += float(V.C_soil_fast + V.C_soil_slow + V.C_soil_passive)
                rh_arr[im] = rh_avg / steps_per_month
                # ra_arr[im] = ra_avg / steps_per_month
                cVeg_arr[im] = cVeg_avg / steps_per_month
                cLitter_arr[im] = cLitter_avg / steps_per_month
                cSoil_arr[im] = cSoil_avg / steps_per_month
                im += 1
            # if y == 100:
            #    print(V)
        return Observables(
            cVeg=cVeg_arr, cLitter=cLitter_arr, cSoil=cSoil_arr, rh=rh_arr
        )  # ,
        # ra=ra_arr)

    return param2res


# +
Full_output = namedtuple(
    "Full_output",
    [
        "C_leaf",
        "C_wood",
        "C_root",
        "C_leaf_litter",
        "C_wood_litter",
        "C_root_litter",
        "C_soil_fast",
        "C_soil_slow",
        "C_soil_passive",
        "rh",
    ],  # ,"ra"]
)


def make_param2res_full_output(
    mvs, cpa: Constants, dvs: Drivers
) -> Callable[[np.ndarray], np.ndarray]:
    # To compute numeric solutions we will need to build and iterator
    # as we did before. As it will need numeric values for all the parameters
    # we will have to create a complete dictionary for all the symbols
    # exept those for the statevariables and time.
    # This set of symbols does not change during the mcmc procedure, since it only
    # depends on the symbolic model.
    # Therefore we create it outside the mcmc loop and bake the result into the
    # param2res function.
    # The iterator does not care if they are estimated or constant, it only
    # wants a dictionary containing key: value pairs for all
    # parameters that are not state variables or the symbol for time
    StartVector = make_StartVector(mvs)
    srm = mvs.get_SmoothReservoirModel()
    model_par_dict_keys = srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))] + list(mvs.get_StateVariableTuple())
    )

    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test
    #     seconds_per_day = 86400
    #     def npp_func(day):
    #         month=gh.day_2_month_index_vm(day)
    #         return dvs.npp[month] #* seconds_per_day   # kg/m2/s kg/m2/day;

    def param2res_full_output(pa):
        epa = EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
        V_init = StartVector(
            C_leaf=epa.C_leaf_0,
            C_wood=epa.C_wood_0,
            C_root=cpa.cVeg_0 - (epa.C_leaf_0 + epa.C_wood_0),
            C_leaf_litter=epa.C_leaf_litter_0,
            C_wood_litter=epa.C_wood_litter_0,
            C_root_litter=cpa.cLitter_0 - (epa.C_leaf_litter_0 + epa.C_wood_litter_0),
            C_soil_fast=epa.C_soil_fast_0,
            C_soil_slow=epa.C_soil_slow_0,
            C_soil_passive=cpa.cSoil_0 - (epa.C_soil_fast_0 + epa.C_soil_slow_0),
            # ra=cpa.ra_0,
            rh=cpa.rh_0,
        )
        # next we create the parameter dict for the iterator
        # The iterator does not care if they are estimated or not so we look for them
        # in the combination
        apa = {**cpa._asdict(), **epa._asdict()}
        model_par_dict = {
            Symbol(k): v for k, v in apa.items() if Symbol(k) in model_par_dict_keys
        }

        # print(model_par_dict)
        # from IPython import embed;embed()

        # Beside the par_dict the iterator also needs the python functions to replace the symbolic ones with
        # our fake xi_func could of course be defined outside of param2res but in general it
        # could be build from estimated parameters and would have to live here...
        #         def xi_func(day):
        #             month=gh.day_2_month_index_vm(day)
        #             TS = (dvs.tas[month]-273.15) # # convert from Kelvin to Celcius
        #             TS = 0.748*TS + 6.181 # approximate soil T at 20cm from air T (from https://doi.org/10.1155/2019/6927045)
        #             if TS > epa.T_0:
        #                 xi_out = np.exp(epa.E*(1/(10-epa.T_0)-1/(TS-epa.T_0))) * dvs.mrso[month]/(epa.KM+dvs.mrso[month])
        #             else:
        #                 xi_out=0
        #             return(xi_out)
        #             # return 1.0 # preliminary fake for lack of better data...

        #         func_dict={
        #             'NPP':npp_func,
        #              'xi':xi_func
        #         }
        func_dict = make_func_dict(dvs)
        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2)
        delta_t_val = 15
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val,
        )

        # Now that we have the iterator we can start to compute.
        # the first thing to notice is that we don't need to store all values (daili yi delta_t_val=1)
        # since the observations are recorded monthly while our iterator possibly has a smaller timestep.
        # - for the pool contents we only need to record the last day of the month
        # - for the respiration(s) ra and rh we want an accumulated value (unit mass)
        #   have to sum up the products of the values*delta_t over a month
        #
        # Note: check if TRENDY months are like this...
        # days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        # sols=[]

        # empty arrays for saving data
        C_leaf_arr = np.zeros(cpa.number_of_months)
        C_wood_arr = np.zeros(cpa.number_of_months)
        C_root_arr = np.zeros(cpa.number_of_months)
        C_leaf_litter_arr = np.zeros(cpa.number_of_months)
        C_wood_litter_arr = np.zeros(cpa.number_of_months)
        C_root_litter_arr = np.zeros(cpa.number_of_months)
        C_soil_fast_arr = np.zeros(cpa.number_of_months)
        C_soil_slow_arr = np.zeros(cpa.number_of_months)
        C_soil_passive_arr = np.zeros(cpa.number_of_months)
        rh_arr = np.zeros(cpa.number_of_months)
        # ra_arr = np.zeros(cpa.number_of_months)
        im = 0
        dpm = 30
        steps_per_month = int(dpm / delta_t_val)
        # n=len(V_init)
        # forward simulation by year
        for y in range(int(cpa.number_of_months / 12)):
            for m in range(12):
                C_leaf_avg = 0
                C_wood_avg = 0
                C_root_avg = 0
                C_leaf_litter_avg = 0
                C_wood_litter_avg = 0
                C_root_litter_avg = 0
                C_soil_fast_avg = 0
                C_soil_slow_avg = 0
                C_soil_passive_avg = 0
                rh_avg = 0
                # ra_avg = 0
                for d in range(steps_per_month):
                    V = StartVector(*it_sym.__next__())
                    rh_avg += V.rh
                    # ra_avg += V.ra
                    C_leaf_avg += float(V.C_leaf)
                    C_wood_avg += float(V.C_wood)
                    C_root_avg += float(V.C_root)
                    C_leaf_litter_avg += float(V.C_leaf_litter)
                    C_wood_litter_avg += float(V.C_wood_litter)
                    C_root_litter_avg += float(V.C_root_litter)
                    C_soil_fast_avg += float(V.C_soil_fast)
                    C_soil_slow_avg += float(V.C_soil_slow)
                    C_soil_passive_avg += float(V.C_soil_passive)
                rh_arr[im] = rh_avg / steps_per_month
                # ra_arr[im] = ra_avg / steps_per_month
                C_leaf_arr[im] = C_leaf_avg / steps_per_month
                C_wood_arr[im] = C_wood_avg / steps_per_month
                C_root_arr[im] = C_root_avg / steps_per_month
                C_leaf_litter_arr[im] = C_leaf_litter_avg / steps_per_month
                C_wood_litter_arr[im] = C_wood_litter_avg / steps_per_month
                C_root_litter_arr[im] = C_root_litter_avg / steps_per_month
                C_soil_fast_arr[im] = C_soil_fast_avg / steps_per_month
                C_soil_slow_arr[im] = C_soil_slow_avg / steps_per_month
                C_soil_passive_arr[im] = C_soil_passive_avg / steps_per_month
                im += 1
            # if y == 100:
            #    print(V)
        return Full_output(
            C_leaf=C_leaf_arr,
            C_wood=C_wood_arr,
            C_root=C_root_arr,
            C_leaf_litter=C_leaf_litter_arr,
            C_wood_litter=C_wood_litter_arr,
            C_root_litter=C_root_litter_arr,
            C_soil_fast=C_soil_fast_arr,
            C_soil_slow=C_soil_slow_arr,
            C_soil_passive=C_soil_passive_arr,
            rh=rh_arr,
        )  # ,
        # ra=ra_arr)

    return param2res_full_output


# -


def make_param_filter_func(
    c_max: EstimatedParameters, c_min: EstimatedParameters
) -> Callable[[np.ndarray], bool]:

    # find position of beta_leaf and beta_wood
    beta_leaf_ind = EstimatedParameters._fields.index("beta_leaf")
    beta_wood_ind = EstimatedParameters._fields.index("beta_wood")

    def isQualified(c):
        beta_leaf_ind
        cond1 = (c >= c_min).all()
        cond2 = (c <= c_max).all()
        cond3 = c[beta_leaf_ind] + c[beta_wood_ind] < 1
        return cond1 and cond2 and cond3

    return isQualified


def numeric_X_0(mvs, dvs, cpa, epa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict = gh.make_param_dict(mvs, cpa, epa)
    X_0_dict = {
        "C_leaf": apa["C_leaf_0"],
        "C_wood": apa["C_wood_0"],
        "C_root": apa["cVeg_0"] - (apa["C_leaf_0"] + apa["C_wood_0"]),
        "C_leaf_litter": apa["C_leaf_litter_0"],
        "C_wood_litter": apa["C_wood_litter_0"],
        "C_root_litter": apa["cLitter_0"]
        - (apa["C_leaf_litter_0"] + apa["C_wood_litter_0"]),
        "C_soil_fast": apa["C_soil_fast_0"],
        "C_soil_slow": apa["C_soil_slow_0"],
        "C_soil_passive": apa["cSoil_0"]
        - (apa["C_soil_fast_0"] + apa["C_soil_slow_0"]),
    }
    X_0 = np.array([X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()]).reshape(
        len(X_0_dict), 1
    )
    return X_0


# def days_per_month():
#    #dpm= [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
#    dpm= [30 for i in range(12)]
#    return dpm


def start_date():
    ## this function is important to syncronise our results
    ## because our data streams start at different times the first day of
    ## a simulation day_ind=0 refers to different dates for different models
    ## we have to check some assumptions on which this calculation is based
    ## Here is how to get these values
    # ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("VISIT_S2_gpp.nc")))
    # times = ds.variables["time"]
    # tm = times[0] #time of first observation in Months_since_1860-01 # print(times.units)
    # td = int(tm *30)  #in days since_1860-01-01
    # import datetime as dt
    # ad = dt.date(1, 1, 1) # first of January of year 1
    # sd = dt.date(1860, 1, 1)
    # td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    ## from td_aD (days since 1-1-1) we can compute the year month and day
    return gh.date(year=1860, month=1, day=1)


data_str = namedtuple("data_str", ["cVeg", "cLitter", "cSoil", "gpp", "ra", "rh"])


def get_global_mean_vars_all(experiment_name):
    return gh.get_global_mean_vars_all(
        model_folder="kv_visit2",
        experiment_name=experiment_name,
        lat_var="lat",
        lon_var="lon",
    )


# ################ function for computing global mean for custom data streams ###################

# def get_global_mean_vars_all(experiment_name="VISIT_S2_"):
# def nc_file_name(nc_var_name, experiment_name="VISIT_S2_"):
# return experiment_name+"{}.nc".format(nc_var_name)

# def nc_global_mean_file_name(experiment_name="VISIT_S2_"):
# return experiment_name+"mean_vars.nc"

# data_str = namedtuple(
# 'data_str',
# ["cVeg", "cLitter", "cSoil", "gpp", "ra", "rh"]
# )
# names = data_str._fields
# conf_dict = gh.confDict("kv_visit2")
# dataPath=Path(conf_dict["dataPath"])

# if dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name)).exists():
# print(""" Found cached global mean files. If you want to recompute the global means
# remove the following files: """
# )
# print( dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name)))

# def get_cached_global_mean(vn):
# gm = gh.get_nc_array(dataPath.joinpath(
# nc_global_mean_file_name(experiment_name=experiment_name)),vn)
# return gm

# #map variables to data
# output=gh.data_streams(*map(get_cached_global_mean, gh.data_streams._fields))
# return (
# output
# )

# else:
# gm=gh.globalMask()
# # load an example file with mask
# template = nc.Dataset(dataPath.joinpath("VISIT_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
# gcm=gh.project_2(
# source=gm,
# target=gh.CoordMask(
# index_mask=np.zeros_like(template),
# tr=gh.SymTransformers(
# ctr=make_model_coord_transforms(),
# itr=make_model_index_transforms()
# )
# )
# )

# print("computing means, this may take some minutes...")

# def compute_and_cache_global_mean(vn):
# path = dataPath.joinpath(nc_file_name(vn, experiment_name=experiment_name))
# ds = nc.Dataset(str(path))
# vs=ds.variables
# lats= vs["lat"].__array__()
# lons= vs["lon"].__array__()
# print(vn)
# var=ds.variables[vn]
# # check if we have a cached version (which is much faster)
# gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))

# gm=gh.global_mean_var(
# lats,
# lons,
# gcm.index_mask,
# var
# )
# return gm * 86400 if vn in ["gpp", "npp", "rh", "ra"] else gm

# #map variables to data
# output=data_str(*map(compute_and_cache_global_mean, data_str._fields))
# output_final=gh.data_streams(
# cVeg=output.cVeg,
# cSoil=output.cLitter+output.cSoil,
# gpp=output.gpp,
# npp=output.gpp-output.ra,
# ra=output.ra,
# rh=output.rh,
# )
# gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))
# gh.write_data_streams_cache(gm_path, output_final)
# return (
# output_final
# )

class ConsistentTestArgs():
    def __init__(
            self,
            confDict: Dict[str,str],
            mvs: CMTVS
        ):
        """A model (and data) dependent class that expresses
        dependencies between certain parts of the parameterization as
        dictated by a particular choice of estimated and constant
        paramters for parameter estimation """
        self.confDict = confDict
        self.mvs = mvs
    
    @property
    def dataPath(self):
        return Path(self.confDict["dataPath"])

    @property
    @lru_cache(maxsize=None)
    def timelines(self):
        svs, dvs = get_cached_global_mean_vars(self.dataPath)
        return svs, dvs

    @property
    def cpa(self):
        svs, dvs = self.timelines
        svs_0 = Observables(*map(lambda v: v[0],svs))
        cpa = Constants(
            cVeg_0=svs_0.cVeg,
            cLitter_0=svs_0.cLitter,
            cSoil_0=svs_0.cSoil,
            rh_0=svs_0.rh,   # kg/m2/s kg/m2/day
            npp_0=dvs.npp[0],   # kg/m2/s kg/m2/day
            #ra_0=svs_0.ra,   # kg/m2/s kg/m2/day
            #r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
            #r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
            number_of_months=len(svs.rh)
            #number_of_months=24 # for testing and tuning mcmc
        )
        return cpa


    @property
    def epa_opt(self):
        # this function should compute epa_opt by parameter estimation
        # this is a fake until this is implemented
        return EstimatedParameters(
            beta_leaf=0.6102169482865195, 
            beta_wood=0.26331553815787545, 
            T_0=1.9560345980471245, 
            E=6.951145421284498, 
            KM=12.73895376887386, 
            r_C_leaf_litter_rh=0.0012830039575098323, 
            r_C_wood_litter_rh=0.0010536416454437036, 
            r_C_root_litter_rh=0.00022271755326847413, 
            r_C_soil_fast_rh=0.00013707839288781872, 
            r_C_soil_slow_rh=3.228645064482276e-05, 
            r_C_soil_passive_rh=3.8079656062059605e-06, 
            r_C_leaf_2_C_leaf_litter=0.011755309034589333, 
            r_C_wood_2_C_wood_litter=0.00012990716959685548, 
            r_C_root_2_C_root_litter=8.243205281709114e-05, 
            r_C_leaf_litter_2_C_soil_fast=0.0014521759026031634, 
            r_C_leaf_litter_2_C_soil_slow=0.000200225210999593, 
            r_C_leaf_litter_2_C_soil_passive=8.380707345301035e-05, 
            r_C_wood_litter_2_C_soil_fast=3.19128931685429e-05, 
            r_C_wood_litter_2_C_soil_slow=7.278721749448471e-05, 
            r_C_wood_litter_2_C_soil_passive=3.275165103336979e-06, 
            r_C_root_litter_2_C_soil_fast=0.00044055481426693227, 
            r_C_root_litter_2_C_soil_slow=3.1019188662910444e-05, 
            r_C_root_litter_2_C_soil_passive=0.00012243099600679218, 
            C_leaf_0=0.042596582017273114, 
            C_wood_0=2.4493874052342517, 
            C_leaf_litter_0=0.16251047110808622, 
            C_wood_litter_0=0.17601405444541945, 
            C_soil_fast_0=1.8746682268250323, 
            C_soil_slow_0=1.9807766341505468
        )

    @property
    def global_mean_drivers(self):
        return get_cached_global_mean_drivers(self.dataPath)

    @property
    def par_dict(self):
        return  gh.make_param_dict(
            self.mvs, 
            self.cpa,
            self.epa_opt
        )

    def write_data(self, targetPath):
        # write global mean var *.nc files in 
        # the local directory
        #compute_and_write_global_mean_vars(self.dataPath, targetPath)
        gh.dump_dict_to_json_path(
            self.par_dict,
            targetPath.joinpath('parameter_dict.json')
        )
        dvs_dict=self.global_mean_drivers._asdict()
        for k,v in dvs_dict.items():
            gh.write_timeline_to_nc_file(targetPath.joinpath(f"{k}.nc"),v,k)

        
        #gh.dump_named_tuple_to_json_path(
        #    self.cpa,
        #    targetPath.joinpath('cpa.json')
        #)
        #gh.dump_named_tuple_to_json_path(
        #    self.epa_opt,
        #    targetPath.joinpath('epa_opt.json')
        #)
        
