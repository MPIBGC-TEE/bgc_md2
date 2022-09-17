from cmath import pi, sin
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
from copy import copy
from typing import Callable
from functools import reduce

sys.path.insert(0,'..') # necessary to import general_helpers
import general_helpers as gh

def spatial_mask(dataPath)->'CoorMask':
    f_mask=nc.Dataset(dataPath.joinpath("IBIS_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
    
    print("computing masks to exclude pixels with nan entries, this may take some minutes...")
    
    def f(vn):
        path = dataPath.joinpath(nc_file_name(vn))
        ds = nc.Dataset(str(path))
        var =ds.variables[vn]
        ##return after assessing NaN data values
        return gh.get_nan_pixel_mask(var)

    o_names=Observables._fields
    d_names=Drivers._fields
    names = o_names + d_names 

    masks=[ f(name)    for name in names ]
    # We compute the common mask so that it yields valid pixels for ALL variables 
    combined_mask= reduce(lambda acc,m: np.logical_or(acc,m),masks,f_mask)
    print("found additional {} NaN pixels".format(combined_mask.sum()-f_mask.sum()))
    #from IPython import embed;embed() 
    sym_tr= gh.SymTransformers(
        itr=make_model_index_transforms(),
        ctr=make_model_coord_transforms()
    )
    return gh.CoordMask(
        combined_mask,
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

# we will use the trendy output names directly in other parts of the output
Observables_annual = namedtuple(
    'Observables_annual',
    ["cVeg","cLitter","cSoil"]
)
Observables_monthly = namedtuple(
    'Observables_monthly',
    #["rh","ra"]
    ["rh"]
)
Observables = namedtuple(
    "Observables",
    Observables_annual._fields+Observables_monthly._fields
)
# OrgDrivers=namedtuple(
#     "OrgDrivers",
#     ["gpp", "mrso", "tas"]
# )    
Drivers=namedtuple(
    "Drivers",
    ["npp","mrso","tas"] #
)      
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
        #"ra_0",
        
        # "k_C_mll",
        # "k_C_mwl",
        # "k_C_mrl",
        # "k_C_sll",
        # "k_C_swl",
        # "k_C_srl"test_get_global_mean_vars,
        # "k_C_lll",
        # "k_C_lwl",
        # "k_C_lrl",
        
        # "r_C_mll_2_C_mic", # f13_4
        # "r_C_mwl_2_C_mic", # f13_5
        # "r_C_mrl_2_C_mic", # f13_6
        # "r_C_sll_2_C_mic", # f13_7
        # "r_C_swl_2_C_mic", # f13_8
        # "r_C_srl_2_C_mic", # f13_9
        # "r_C_pass_2_C_mic", # f13_16

        # "r_C_lll_2_C_prot", # f14_10
        # "r_C_lwl_2_C_prot", # f14_11
        # "r_C_lrl_2_C_prot", # f14_12
        # "r_C_lll_2_C_nonprot", # f15_10
        # "r_C_lwl_2_C_nonprot", # f15_11
        # "r_C_lrl_2_C_nonprot", # f15_12
        
        "number_of_months" # necessary to prepare the output in the correct lenght 
    ]
)
# Note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues and their length (number of months)
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 
# It is better to start with only a few

EstimatedParameters = namedtuple(
    "EstimatedParameters",[ 

        "beta_leaf",
        "beta_wood",
        
        "r_C_leaf_2_C_mll",     # 2,  f4_1
        "r_C_wood_2_C_mwl",     # 3,  f5_2
        "r_C_root_2_C_mrl",     # 4,  f6_3
        "r_C_leaf_2_C_sll",     # 5,  f7_1
        "r_C_wood_2_C_swl",     # 6,  f8_2
        "r_C_root_2_C_srl",     # 7,  f9_3
        
        "r_C_leaf_2_C_lll",
        "r_C_wood_2_C_lwl",
        "r_C_root_2_C_lrl",
        
        "r_C_prot_2_C_mic",     # 8， f13_14
        "r_C_nonprot_2_C_mic",  # 9， f13_15
        "r_C_mic_2_C_prot",     # 10, f14_13
        "r_C_mic_2_C_nonprot",  # 11, f15_13
        "r_C_prot_2_C_pass",    # 12, f16_14
        "r_C_nonprot_2_C_pass", # 13, f16_15

        # "k_C_mic",              # 17
        # "k_C_protsom",          # 18
        # "k_C_nonprotsom",       # 19
        # "k_C_passsom",          # 20

        "C_wood_0",           # 23
        "C_root_0",           # 24

        "C_mll_0",            # 25
        "C_mwl_0",            # 26
        "C_sll_0",            # 27
        "C_swl_0",            # 28
        "C_lll_0",            # 29

        "C_mrl_0",            # 30
        "C_srl_0",            # 31
        "C_lrl_0",            # 32 
        "C_mic_0",            # 33 
        "C_prot_0",           # 34
        "C_nonprot_0",        # 35
        
         "r_C_mll_rh",
         "r_C_mwl_rh",
         "r_C_mrl_rh",
         "r_C_sll_rh",
         "r_C_swl_rh",
         "r_C_srl_rh",          
         "r_C_lll_rh",
         "r_C_lwl_rh",
         "r_C_lrl_rh",            
         "r_C_mic_rh",           
         "r_C_prot_rh",
         "r_C_nonprot_rh",        
         "r_C_pass_rh",
         
         
         
        "r_C_mll_2_C_mic", # f13_4
        "r_C_mwl_2_C_mic", # f13_5
        "r_C_mrl_2_C_mic", # f13_6
        "r_C_sll_2_C_mic", # f13_7
        "r_C_swl_2_C_mic", # f13_8
        "r_C_srl_2_C_mic", # f13_9
        "r_C_pass_2_C_mic", # f13_16

        "r_C_lll_2_C_prot", # f14_10
        "r_C_lwl_2_C_prot", # f14_11
        "r_C_lrl_2_C_prot", # f14_12
        "r_C_lll_2_C_nonprot", # f15_10
        "r_C_lwl_2_C_nonprot", # f15_11
        "r_C_lrl_2_C_nonprot", # f15_12

    ]
)
# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues. 
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise         
# to guard agains accidentally changed order we use a namedtuple again. Since B_func and u_func rely 
# on the correct ordering of the statevariables we build V dependent on this order 

#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['IBIS'],
        variables = Observables._fields + Drivers._fields
    )

experiment_name="IBIS_S2_"
def nc_file_name(nc_var_name):
    return experiment_name+"{}.nc".format(nc_var_name)


def nc_global_mean_file_name(nc_var_name):
    return experiment_name+"{}_gm.nc".format(nc_var_name)


# this function is deprecated: please use get_global_mean_vars

# def get_global_vars(dataPath):
#     #define function to average variables
#     def f(tup):
#         #define parts of function from nc file
#         vn, fn = tup
#         path = dataPath.joinpath(fn)
#         ds = nc.Dataset(str(path))
#         lats = ds.variables["latitude"]
#         lons = ds.variables["longitude"]

#         #check for npp/gpp/rh/ra to convert from kg/m2/s to kg/m2/day
#         if vn in ["npp","gpp","rh","ra"]:
#             #for name, variable in ds.variables.items():            
#             #    for attrname in variable.ncattrs():
#             #        print("{} -- {}".format(attrname, getattr(variable, attrname)))
#             return (gh.global_mean(lats, lons,ds.variables[vn].__array__())*24*60*60)
#         else:
#             #for name, variable in ds.variables.items():            
#             #    for attrname in variable.ncattrs():
#             #        print("{} -- {}".format(attrname, getattr(variable, attrname)))
#             return (gh.global_mean(lats, lons, ds.variables[vn].__array__()))

#     # Link symbols and data:
#     # IBIS has annual vs monthly file names so they are linked separately
#     # If all your data is similarly named you can do this in one step

#     # Create annual file names (single step if files similarly named)
#     o_names=[(f,"IBIS_S2_{}.nc".format(f)) for f in Observables_annual._fields]

#     # Create monthly file names (can remove if done in one step above)
#     monthly_names=[(f,"IBIS_S2_{}.nc".format(f)) for f in Observables_monthly._fields]
#     # Extend name list with monthly names
#     o_names.extend(monthly_names)

#     # create file names for Drivers
#     d_names=[(f,"IBIS_S2_{}.nc".format(f)) for f in Drivers._fields]

#     # we want to drive with npp and can create it from gpp and ra 
#     # observables
#     odvs=OrgDrivers(*map(f,d_names))
#     obss=Observables(*map(f, o_names))

#     dvs=Drivers(
#         npp=odvs.gpp-obss.ra,
#         mrso=odvs.mrso,
#         tas=odvs.tas
#     )

#     # Link symbols and data for Observables/Drivers
#     # return (Observables(*map(f, o_names)),Drivers(*map(f,d_names)))
#     return (obss, dvs)
# -

def get_example_site_vars(dataPath):
    # pick up 1 site
    s = slice(None, None, None)  # this is the same as :
    t = s, 72, 117  # [t] = [:,49,325]

    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them
        ds = nc.Dataset(str(path))
        return ds.variables[vn][t]

    # note that there use S2 for IBIS model
    o_names = [(f, "IBIS_S2_{}.nc".format(f)) for f in Observables._fields]
    d_names = [(f, "IBIS_S2_{}.nc".format(f)) for f in Drivers._fields]

    # we want to drive with npp and can create it from gpp and ra 
    # observables
    dvs = Drivers(*map(f, d_names))
    obss = Observables(*map(f, o_names))

    return (obss, dvs)


def get_global_mean_vars(dataPath):
    o_names=Observables._fields
    d_names=Drivers._fields
    names = o_names + d_names
    print("names")
    print(names)
    print("Observables")
    print(Observables._fields)
    if all([dataPath.joinpath(nc_global_mean_file_name(vn)).exists() for vn in names]):
        print(""" Found cached global mean files. If you want to recompute the global means
            remove the following files: """
        )
        for vn in names:
            print( dataPath.joinpath(nc_global_mean_file_name(vn)))

        def get_cached_global_mean(vn):
            gm = gh.get_cached_global_mean(dataPath.joinpath(nc_global_mean_file_name(vn)),vn)
            return gm * 86400 if vn in ["gpp", "rh", "ra", "npp"] else gm
    
        return (
            Observables(*map(get_cached_global_mean, o_names)),
            Drivers(*map(get_cached_global_mean,d_names))
        )

    else:
        # we now check if any of the arrays has a time lime containing nan values 
        # APART FROM values that are already masked by the fillvalue
        
        # print("computing masks to exclude pixels with nan entries, this may take some minutes...")
        # def f(vn):
            # path = dataPath.joinpath(nc_file_name(vn))
            # ds = nc.Dataset(str(path))
            # #scale fluxes vs pools
            # var =ds.variables[vn]
            # return gh.get_nan_pixel_mask(var)

        # masks=[ f(name)    for name in names ]
        # # We compute the common mask so that it yields valid pixels for ALL variables 
        # combined_mask= reduce(lambda acc,m: np.logical_or(acc,m),masks)
        
        # now we use common mask for all models
        gm=gh.globalMask()
        # load an example file with mask
        template = nc.Dataset(dataPath.joinpath("IBIS_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
        gcm=gh.project_2(
                source=gm,
                target=gh.CoordMask(
                    index_mask=np.zeros_like(template),
                    tr=gh.SymTransformers(
                        ctr=make_model_coord_transforms(),
                        itr=make_model_index_transforms()
                    )
                )
        )
       
        print("computing means, this may also take some minutes...")

        def compute_and_cache_global_mean(vn):
            path = dataPath.joinpath(nc_file_name(vn))
            ds = nc.Dataset(str(path))
            vs=ds.variables
            #print(vs)
            lats= vs["latitude"].__array__()
            lons= vs["longitude"].__array__()
            var=ds.variables[vn]
            # check if we have a cached version (which is much faster)
            gm_path = dataPath.joinpath(nc_global_mean_file_name(vn))
            
            # ## THIS IS TEMPORARY. GLOBAL MASK DOES NOT WORK FOR RH AND RA, THEREFORE COMPUTING MODEL-SPECIFIC MASK HERE
            # print("computing masks to exclude pixels with nan entries, this may take some minutes...")
            # def f(vn):
                # path = dataPath.joinpath(nc_file_name(vn))
                # ds = nc.Dataset(str(path))
                # #scale fluxes vs pools
                # var =ds.variables[vn]
                # return gh.get_nan_pixel_mask(var)

            # masks=[ f(name)    for name in names ]
            # # We compute the common mask so that it yields valid pixels for ALL variables 
            # combined_mask= reduce(lambda acc,m: np.logical_or(acc,m),masks)
            

            gm=gh.global_mean_var(
                    lats,
                    lons,
                    gcm.index_mask,
                    #combined_mask,
                    var
            )                                
           
            gh.write_global_mean_cache(
                    gm_path,
                    gm,
                    vn
            )
            return gm * 86400 if vn in ["gpp", "rh", "ra", "npp"] else gm
    
        #map variables to data
        return (
            Observables(*map(compute_and_cache_global_mean, o_names)),
            Drivers(*map(compute_and_cache_global_mean, d_names))
        )

def make_npp_func(dvs):
    def func(day):
        month=gh.day_2_month_index(day)
        # kg/m2/s kg/m2/day;
        return (dvs.npp[month]) # * 86400

    return func


import math
def make_xi_func(dvs):
    def func(day):
        
        month = gh.day_2_month_index(day)
        
        tconst = 344.0 # constant for Lloyd and Taylor (1994) function
        bconst = 10.0  # base temperrature used for carbon decompositon
        btemp = 288.16 # maximum value of decomposition factor
        
        T = dvs.tas[month] # do not have soil temp so we use air temp to replace
        
        # temp regulates factor
        if (T > 237.13):
            factor = min(math.exp(tconst * ((1.0 /(btemp-227.13)) - (1.0 /(T-227.13)) )), bconst)
        else:
            factor = math.exp(tconst * ((1.0 /(btemp-227.13)) - (1.0 /(237.13-227.13)) ))
        
        wfps = 55.0 #
        moist = math.exp((wfps - 60.0)**2 /-(800.0))   # moisture regulates factor
        
        factor = max(0.001, min(bconst, factor * moist))
        
        if (factor > 1.0):
            factor = 1
                        
        #print(factor)
        
        return factor # preliminary fake for lack of better data... 
    return func


def make_func_dict(mvs,dvs,cpa,epa):
    return {
        "NPP": make_npp_func(dvs),
        "xi": make_xi_func(dvs)
    }

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

def make_iterator_sym(
        mvs,
        V_init, #: StartVector,
        par_dict,
        func_dict,
        delta_t_val=1 # defaults to 1day timestep
    ):
    B_func, u_func = gh.make_B_u_funcs_2(mvs,par_dict,func_dict,delta_t_val)

    

    sv=mvs.get_StateVariableTuple()
    #mass production of output functions

    n=len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct 
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr=np.array(
        [V_init.__getattribute__(str(v)) for v in sv]+
        #[V_init.ra,V_init.rh]
        [V_init.rh]
    ).reshape(n+1,1) #reshaping is neccessary for matmul (the @ in B @ X)
    
    
    # To compute the ra and rh we have to some up the values for autotrophic and heterotrophic respiration we have to sum up outfluxes.
    # We first create numerical functions for all the outfluxes from our symbolic description.
    numOutFluxesBySymbol={
        k:gh.numfunc(expr_cont, mvs, delta_t_val=delta_t_val, par_dict=par_dict, func_dict=func_dict) 
        for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
    } 
    def f(it,V):
        X = V[0:n]
        b = u_func(it,X)
        B = B_func(it,X)
        X_new = X + b + B @ X
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
        
        # ra = np.sum(
            # [
              # numOutFluxesBySymbol[Symbol(k)](it,*X)
                # for k in ["C_leaf","C_wood","C_root"] 
                # if Symbol(k) in numOutFluxesBySymbol.keys()
            # ]
        # )
        rh = np.sum(
            [
                numOutFluxesBySymbol[Symbol(k)](it,*X)
                for k in ["C_mll","C_mwl","C_mrl","C_sll","C_swl","C_srl","C_lll","C_lwl","C_lrl","C_mic","C_prot","C_nonprot","C_pass"] 
                if Symbol(k) in numOutFluxesBySymbol.keys()
            ]
        )
        V_new = np.concatenate((X_new.reshape(n,1),np.array([rh]).reshape(1,1)), axis=0)
        
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )

def make_StartVector(mvs):
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()]+
        #["ra","rh"]
        ["rh"]
    ) 


# +
def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers
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
    # StartVector=make_StartVector(mvs)
    
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test
    
#     seconds_per_day = 86400
#     def npp_func(day):
#         month=gh.day_2_month_index(day)
#         return dvs.npp[month] * seconds_per_day   # kg/m2/s kg/m2/day;
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
        StartVector=make_StartVector(mvs)
        V_init = StartVector(
            C_leaf=cpa.cVeg_0-(epa.C_wood_0 + epa.C_root_0),
            C_wood=epa.C_wood_0,            
            C_root=epa.C_root_0,
            
            C_mll=epa.C_mll_0,
            C_mwl=epa.C_mwl_0,
            C_mrl=epa.C_mrl_0,
            C_sll=epa.C_sll_0,
            C_swl=epa.C_swl_0,
            C_srl=epa.C_srl_0,
            C_lll=epa.C_lll_0,
            C_lwl=cpa.cLitter_0-(epa.C_mll_0 + epa.C_mwl_0 + epa.C_sll_0 + epa.C_swl_0 + epa.C_lll_0),
            C_lrl=epa.C_lrl_0,
            
            C_mic=epa.C_mic_0,
            C_prot=epa.C_prot_0,
            C_nonprot=epa.C_nonprot_0,
            C_pass=cpa.cSoil_0-(epa.C_mrl_0 + epa.C_srl_0 + epa.C_lrl_0 + epa.C_mic_0 + epa.C_prot_0 + epa.C_nonprot_0),

            #ra=cpa.ra_0,
            rh=cpa.rh_0
        )
        # next we create the parameter dict for the iterator
        # The iterator does not care if they are estimated or not so we look for them
        # in the combination
        apa = {**cpa._asdict(),**epa._asdict()}
        model_par_dict = {
            Symbol(k):v for k,v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }
        
        # commented out by Kostia: par_dict should be computed from epa and cpa. 
        # I changed epa and cpa to be consistent and use only flux rates (r), no turnover rates (k) and transfer coeffitients (f) 
        # added by cybian for fix some errors
        # model_par_dict = {
            # "r_C_mll_rh": cpa.k_C_mll*(1 - cpa.r_C_mll_2_C_mic),
            # "r_C_mwl_rh": cpa.k_C_mwl*(1 - cpa.r_C_mwl_2_C_mic),
            # "r_C_mrl_rh": cpa.k_C_mrl*(1 - cpa.r_C_mrl_2_C_mic),
            # "r_C_sll_rh": cpa.k_C_sll*(1 - cpa.r_C_sll_2_C_mic),
            # "r_C_swl_rh": cpa.k_C_swl*(1 - cpa.r_C_swl_2_C_mic),
            # "r_C_srl_rh": cpa.k_C_srl*(1 - cpa.r_C_srl_2_C_mic),          
            # "r_C_lll_rh": cpa.k_C_lll*(-cpa.r_C_lll_2_C_nonprot - cpa.r_C_lll_2_C_prot + 1),
            # "r_C_lwl_rh": cpa.k_C_lwl*(-cpa.r_C_lwl_2_C_nonprot - cpa.r_C_lwl_2_C_prot + 1),
            # "r_C_lrl_rh": cpa.k_C_lrl*(-cpa.r_C_lrl_2_C_nonprot - cpa.r_C_lrl_2_C_prot + 1),            
            # "r_C_mic_rh": epa.k_C_mic*(-epa.r_C_mic_2_C_nonprot - epa.r_C_mic_2_C_prot + 1),           
            # "r_C_prot_rh": epa.k_C_protsom*(-epa.r_C_prot_2_C_mic - epa.r_C_prot_2_C_pass + 1),
            # "r_C_nonprot_rh": epa.k_C_nonprotsom*(-epa.r_C_nonprot_2_C_mic - epa.r_C_nonprot_2_C_pass + 1),            
            # "r_C_pass_rh": epa.k_C_passsom*(1 - cpa.r_C_pass_2_C_mic),  
            
            # "r_C_mll_rh": cpa.k_C_mll*(1 - cpa.r_C_mll_2_C_mic),
            # "r_C_mwl_rh": cpa.k_C_mwl*(1 - cpa.r_C_mwl_2_C_mic),
            # "r_C_mrl_rh": cpa.k_C_mrl*(1 - cpa.r_C_mrl_2_C_mic),
            # "r_C_sll_rh": cpa.k_C_sll*(1 - cpa.r_C_sll_2_C_mic),
            # "r_C_swl_rh": cpa.k_C_swl*(1 - cpa.r_C_swl_2_C_mic),
            # "r_C_srl_rh": cpa.k_C_srl*(1 - cpa.r_C_srl_2_C_mic),          
            # "r_C_lll_rh": cpa.k_C_lll*(-cpa.r_C_lll_2_C_nonprot - cpa.r_C_lll_2_C_prot + 1),
            # "r_C_lwl_rh": cpa.k_C_lwl*(-cpa.r_C_lwl_2_C_nonprot - cpa.r_C_lwl_2_C_prot + 1),
            # "r_C_lrl_rh": cpa.k_C_lrl*(-cpa.r_C_lrl_2_C_nonprot - cpa.r_C_lrl_2_C_prot + 1),            
            # "r_C_mic_rh": epa.k_C_mic*(-epa.r_C_mic_2_C_nonprot - epa.r_C_mic_2_C_prot + 1),           
            # "r_C_prot_rh": epa.k_C_protsom*(-epa.r_C_prot_2_C_mic - epa.r_C_prot_2_C_pass + 1),
            # "r_C_nonprot_rh": epa.k_C_nonprotsom*(-epa.r_C_nonprot_2_C_mic - epa.r_C_nonprot_2_C_pass + 1),            
            # "r_C_pass_rh": epa.k_C_passsom*(1 - cpa.r_C_pass_2_C_mic),  

            # #"r_C_leaf_2_C_lll": 1.0-(epa.r_C_leaf_2_C_mll + epa.r_C_leaf_2_C_sll),
            # #"r_C_wood_2_C_lwl": 1.0-(epa.r_C_wood_2_C_mwl + epa.r_C_wood_2_C_swl),
            # #"r_C_root_2_C_lrl": 1.0-(epa.r_C_root_2_C_mrl + epa.r_C_root_2_C_srl),
            # "r_C_leaf_2_C_lll":epa.r_C_leaf_2_C_lll,
            # "r_C_wood_2_C_lwl":epa.r_C_wood_2_C_lwl,
            # "r_C_root_2_C_lrl":epa.r_C_root_2_C_lrl,

            # "r_C_mll_2_C_mic": cpa.r_C_mll_2_C_mic,
            # "r_C_mwl_2_C_mic": cpa.r_C_mwl_2_C_mic,
            # "r_C_mrl_2_C_mic": cpa.r_C_mrl_2_C_mic,
            # "r_C_sll_2_C_mic": cpa.r_C_sll_2_C_mic,
            # "r_C_swl_2_C_mic": cpa.r_C_swl_2_C_mic,
            # "r_C_srl_2_C_mic": cpa.r_C_srl_2_C_mic,
            # "r_C_pass_2_C_mic": cpa.r_C_pass_2_C_mic,
            # "r_C_lll_2_C_prot": cpa.r_C_lll_2_C_prot,
            # "r_C_lwl_2_C_prot": cpa.r_C_lwl_2_C_prot,
            # "r_C_lrl_2_C_prot": cpa.r_C_lrl_2_C_prot,
            # "r_C_lll_2_C_nonprot": cpa.r_C_lll_2_C_nonprot,
            # "r_C_lwl_2_C_nonprot": cpa.r_C_lwl_2_C_nonprot,
            # "r_C_lrl_2_C_nonprot": cpa.r_C_lrl_2_C_nonprot,
            
            # "beta_leaf": epa.beta_leaf,
            # "beta_wood": epa.beta_wood,

            # "r_C_leaf_2_C_mll": epa.r_C_leaf_2_C_mll,
            # "r_C_wood_2_C_mwl": epa.r_C_wood_2_C_mwl,
            # "r_C_root_2_C_mrl": epa.r_C_root_2_C_mrl,
            # "r_C_leaf_2_C_sll": epa.r_C_leaf_2_C_sll,
            # "r_C_wood_2_C_swl": epa.r_C_wood_2_C_swl,
            # "r_C_root_2_C_srl": epa.r_C_root_2_C_srl,
            # "r_C_prot_2_C_mic": epa.r_C_prot_2_C_mic,
            # "r_C_nonprot_2_C_mic": epa.r_C_nonprot_2_C_mic,
            # "r_C_mic_2_C_prot": epa.r_C_mic_2_C_prot,
            # "r_C_mic_2_C_nonprot": epa.r_C_mic_2_C_nonprot,
            # "r_C_prot_2_C_pass": epa.r_C_prot_2_C_pass,
            # "r_C_nonprot_2_C_pass": epa.r_C_nonprot_2_C_pass
        # }


        #print('model_par_dict:',model_par_dict)
        #from IPython import embed;embed()
        
        # Beside the par_dict the iterator also needs the python functions to replace the symbolic ones with
        # our fake xi_func could of course be defined outside of param2res but in general it
        # could be build from estimated parameters and would have to live here...
#         def xi_func(day):
#             return 1.0 # preliminary fake for lack of better data... 
    
#         func_dict={
#             'NPP':npp_func,
#              'xi':xi_func
#         }

        func_dict=make_func_dict(mvs,dvs,cpa,epa)
        
        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2) 
        delta_t_val=30
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
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
        sols=[]
        dpm=30 # 
        n=len(V_init)
        
        # added by cybian 
        rhs=np.zeros(cpa.number_of_months)
        #ras=np.zeros(cpa.number_of_months)
        number_of_years=int(cpa.number_of_months/12)
        #print('number_of_years:',number_of_years)
        
        cVegs=np.zeros(number_of_years)
        cLitters=np.zeros(number_of_years)
        cSoils=np.zeros(number_of_years)
        dpm=30
        dpy=dpm*12
        m_id=0
        steps_per_month=int(dpm/delta_t_val)
        steps_per_year=int(dpm/delta_t_val)*12
        
        for y in range(number_of_years):
            #print('y:',y)
            cVeg_ave=0
            cLitter_ave=0
            cSoil_ave=0
            for m in range(12):
                #print('y:',y,'m:',m,'m_id:',m_id)
                #dpm = days_per_month[ m % 12]  
                #mra_ave=0
                mrh_ave=0
                for d in range(steps_per_month): #int(dpm/delta_t_val)
                    v = it_sym.__next__()#.reshape(n,)
                    # actually the original datastream seems to be a flux per area (kg m^-2 ^-s )
                    # at the moment the iterator also computes a flux but in kg^-2 ^day                    
                    V=StartVector(*v)
                    
                    #cVeg_ave=np.array(cVeg_ave, dtype=object)+float(V.C_leaf+V.C_wood+V.C_root),
                    #cLitter_ave=np.array(cLitter_ave, dtype=object)+float(V.C_mll + V.C_mwl + V.C_sll + V.C_swl + V.C_lll + V.C_lwl),
                    #cSoil_ave=np.array(cSoil_ave, dtype=object)+float(V.C_mrl + V.C_srl + V.C_lrl + V.C_mic + V.C_prot + V.C_nonprot + V.C_pass),
                    #mrh_ave=np.array(mrh_ave)+V.rh
                    #mra_ave=np.array(mra_ave)+V.ra
                    cVeg_ave+=float(V.C_leaf+V.C_wood+V.C_root)
                    cLitter_ave+=float(V.C_mll + V.C_mwl + V.C_sll + V.C_swl + V.C_lll + V.C_lwl)
                    cSoil_ave+=float(V.C_mrl + V.C_srl + V.C_lrl + V.C_mic + V.C_prot + V.C_nonprot + V.C_pass)
                    mrh_ave+=V.rh
                    #mra_ave+=V.ra                    
                #print('here:m_id:',m_id)
                rhs[m_id]=mrh_ave/steps_per_month
                #ras[m_id]=mra_ave/steps_per_month
                m_id=m_id+1
            
            #print('Here:y:',y)
            cVegs[y]=cVeg_ave/steps_per_year
            cLitters[y]=cLitter_ave/steps_per_year
            cSoils[y]=cSoil_ave/steps_per_year
            
        return Observables(
            cVeg=cVegs,
            cSoil=cSoils,
            cLitter=cLitters,
            rh=rhs,
            #ra=ras
        )
    
            #comment by cybian
#                 #from IPython import embed;embed()
#                 o=Observables(
#                     cVeg=float(V.C_leaf+V.C_wood+V.C_root),
#                     cLitter=float(V.C_mll + V.C_mwl + V.C_sll + V.C_swl + V.C_lll + V.C_lwl),
#                     cSoil=float(V.C_mrl + V.C_srl + V.C_lrl + V.C_mic + V.C_prot + V.C_nonprot + V.C_pass),
#                     ra=V.ra/seconds_per_day,
#                     rh=V.rh/seconds_per_day # the data is per second while the time units for the iterator refer to days
#                 )
#                 sols.append(o)
            
        #sol=np.stack(sols)       
        #convert to yearly output if necessary (the monthly pool data looks very "yearly")
        #sol_yr=np.zeros(int(cpa.number_of_months/12)*sol.shape[1]).reshape([int(cpa.number_of_months/12),sol.shape[1]])  
        #for i in range(sol.shape[1]):
        #   sol_yr[:,i]=monthly_to_yearly(sol[:,i])
        #sol=sol_yr
        #return sol 
        
    return param2res


# -

def make_param_filter_func(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters 
        ) -> Callable[[np.ndarray], bool]:

    # find position of beta_leaf and beta_wood
    beta_leaf_ind=EstimatedParameters._fields.index("beta_leaf")
    beta_wood_ind=EstimatedParameters._fields.index("beta_wood")

    def isQualified(c):
        beta_leaf_ind
        cond1 =  (c >= c_min).all() 
        cond2 =  (c <= c_max).all() 
        cond3 =  c[beta_leaf_ind]+c[beta_wood_ind] <=1  
        return (cond1 and cond2 and cond3)
        
    
    return isQualified


def make_weighted_cost_func(
        obs: Observables
    ) -> Callable[[Observables],np.float64]:
    
    # first unpack the observation array into its parts
    #cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
    def costfunction(out_simu: np.ndarray) ->np.float64:
        # fixme
        #   as indicated by the fact that the function lives in this
        #   model-specific module it is not apropriate for (all) other models.
        #   There are model specific properties:
        #   1.) The weight for the different observation streams
        #
        number_of_ys=out_simu.cVeg.shape[0]
        number_of_ms=out_simu.rh.shape[0]

        J_obj1 = np.mean (( out_simu.cVeg - obs.cVeg )**2)/(2*np.var(obs.cVeg))
        J_obj2 = np.mean (( out_simu.cLitter - obs.cLitter )**2)/(2*np.var(obs.cLitter))
        J_obj3 = np.mean (( out_simu.cSoil -  obs.cSoil )**2)/(2*np.var(obs.cSoil))

        J_obj4 = np.mean (( out_simu.rh - obs.rh )**2)/(2*np.var(obs.rh))

        J_new = (J_obj1 + J_obj2 + J_obj3)/200 + J_obj4/4
        # to make this special costfunction comparable (in its effect on the
        # acceptance rate) to the general costfunction proposed by Feng we
        # rescale it by a factor
        return J_new*400
    return costfunction


def pesudo_yearly_to_monthly(yearly):
    # Added by cybian just for extending the data from yearly to monthly for matching all variables' length
    # not used in InspectModel.py
    yearly_len = len(yearly)
    months_per_year = 12
    monthly_data = np.zeros(yearly_len*months_per_year)+1

    for iyear in range(0, yearly_len):
        for imonth in range(0, months_per_year):
            monthly_data[iyear*12:(iyear+1)*12] = monthly_data[iyear*12:(iyear+1)*12] * yearly[iyear]#(yearly[iyear]-yearly[iyear]*sin(2*pi/12*imonth - pi/2))

    return monthly_data

# this function is deprecated - see general helpers traceability_iterator
# def make_traceability_iterator(mvs,dvs,cpa,epa):

def numeric_X_0(mvs,dvs,cpa,epa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    #apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "C_leaf":cpa.cVeg_0-(epa.C_wood_0 + epa.C_root_0),
        "C_wood":epa.C_wood_0,            
        "C_root":epa.C_root_0,

        "C_mll":epa.C_mll_0,
        "C_mwl":epa.C_mwl_0,
        "C_mrl":epa.C_mrl_0,
        "C_sll":epa.C_sll_0,
        "C_swl":epa.C_swl_0,
        "C_srl":epa.C_srl_0,
        "C_lll":epa.C_lll_0,
        "C_lwl":cpa.cLitter_0-(epa.C_mll_0 + epa.C_mwl_0 + epa.C_sll_0 + epa.C_swl_0 + epa.C_lll_0),
        "C_lrl":epa.C_lrl_0,

        "C_mic":epa.C_mic_0,
        "C_prot":epa.C_prot_0,
        "C_nonprot":epa.C_nonprot_0,
        "C_pass":cpa.cSoil_0-(epa.C_mrl_0 + epa.C_srl_0 + epa.C_lrl_0 + epa.C_mic_0 + epa.C_prot_0 + epa.C_nonprot_0),
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    return X_0

def start_date():
    ## this function is important to syncronise our results
    ## because our data streams start at different times the first day of 
    ## a simulation day_ind=0 refers to different dates for different models
    ## we have to check some assumptions on which this calculation is based
    ## Here is how to get these values
    #ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("VISIT_S2_gpp.nc")))
    #times = ds.variables["time"]
    #tm = times[0] #time of first observation in Months_since_1860-01 # print(times.units)
    #td = int(tm *30)  #in days since_1860-01-01 
    #import datetime as dt
    #ad = dt.date(1, 1, 1) # first of January of year 1 
    #sd = dt.date(1860, 1, 1)
    #td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    ## from td_aD (days since 1-1-1) we can compute the year month and day
    return gh.date(
        year=1700, 
        month=1,
        day=1
    )

################ function for computing global mean for custom data streams ###################
    
def get_global_mean_vars_all(experiment_name="IBIS_S2_"):
    
    def nc_file_name(nc_var_name, experiment_name="IBIS_S2_"):
        return experiment_name+"{}.nc".format(nc_var_name)

    def nc_global_mean_file_name(nc_var_name, experiment_name="IBIS_S2_"):
        return experiment_name+"{}_gm_all.nc".format(nc_var_name)

    data_str = namedtuple( # data streams available in the model
        'data_str',
        ["cVeg", "cLitter", "cSoil", "gpp", "npp", "ra", "rh"]
        )
        
    names = data_str._fields
    conf_dict = gh.confDict("bian_ibis2")
    # with Path('config.json').open(mode='r') as f:
        # conf_dict = frozendict(json.load(f))
    dataPath=Path(conf_dict["dataPath"])    
    
    if all([dataPath.joinpath(nc_global_mean_file_name(vn, experiment_name=experiment_name)).exists() for vn in names]):
        print(""" Found cached global mean files. If you want to recompute the global means
            remove the following files: """
        )
        for vn in names:
            print( dataPath.joinpath(nc_global_mean_file_name(vn,experiment_name=experiment_name)))

        def get_cached_global_mean(vn):
            gm = gh.get_cached_global_mean(dataPath.joinpath(nc_global_mean_file_name(vn,experiment_name=experiment_name)),vn)
            return gm * 86400 if vn in ["gpp", "npp", "rh", "ra"] else gm

        #map variables to data
        output=gh.data_streams(*map(get_cached_global_mean, data_str._fields))
        return (
            output
        )

    else:
        gm=gh.globalMask()
        # load an example file with mask
        template = nc.Dataset(dataPath.joinpath("IBIS_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
        gcm=gh.project_2(
                source=gm,
                target=gh.CoordMask(
                    index_mask=np.zeros_like(template),
                    tr=gh.SymTransformers(
                        ctr=make_model_coord_transforms(),
                        itr=make_model_index_transforms()
                    )
                )
        )

        print("computing means, this may take some minutes...")

        def compute_and_cache_global_mean(vn):
            path = dataPath.joinpath(nc_file_name(vn, experiment_name=experiment_name))
            ds = nc.Dataset(str(path))
            vs=ds.variables
            lats= vs["latitude"].__array__()
            lons= vs["longitude"].__array__()
            print(vn)
            var=ds.variables[vn]
            # check if we have a cached version (which is much faster)
            gm_path = dataPath.joinpath(nc_global_mean_file_name(vn, experiment_name=experiment_name))

            gm=gh.global_mean_var(
                    lats,
                    lons,
                    gcm.index_mask,
                    var
            )
            gh.write_global_mean_cache(
                    gm_path,
                    gm,
                    vn
            )
            return gm * 86400 if vn in ["gpp", "npp", "rh", "ra"] else gm
        
        #map variables to data
        output=data_str(*map(compute_and_cache_global_mean, data_str._fields))
        return (
            gh.data_streams( # required data streams
                cVeg=output.cVeg,
                cSoil=output.cLitter+output.cSoil,
                gpp=output.gpp,
                npp=output.npp,
                ra=output.ra,
                rh=output.rh,
            )
        )