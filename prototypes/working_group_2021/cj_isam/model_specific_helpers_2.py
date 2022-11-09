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
    mask=nc.Dataset(dataPath.joinpath("ISAM_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
    sym_tr= gh.SymTransformers(
        itr=make_model_index_transforms(),
        ctr=make_model_coord_transforms()
    )
    return gh.CoordMask(
        mask,
        sym_tr
    )

def make_model_coord_transforms():
    """ This function can is used to achieve a target grid LAT,LON with
    - LAT ==   0 at the equator with 
    - LAT == -90 at the south pole,
    - LAT== +90 at the north pole,
    - LON ==   0 at Greenich and 
    - LON is counted positive eastwards from -180 to 180
    """

    return gh.CoordTransformers(
            lat2LAT=lambda lat: lat,
            LAT2lat=lambda LAT: LAT,
            lon2LON=lambda lon: -180+ lon-180 if lon > 180 else lon,
            LON2lon=lambda LON: 360+LON if LON < 0 else LON
    )


def make_model_index_transforms():
    return gh.transform_maker(
    lat_0 = -89.75,
    lon_0 = 0.25,
    step_lat = 0.5,
    step_lon = 0.5,
 )
# we will use the trendy output names directly in other parts of the output
Observables = namedtuple(
    'Observables',
    ["cVeg","cLitter","cSoil","rh","ra"]
)
#OrgDrivers=namedtuple(
#    "OrgDrivers",
#    ["gpp", "mrso", "tas"]
#)    
#Drivers=namedtuple(
#    "Drivers",
#    ("npp",) + OrgDrivers._fields[1:]
#)    
Drivers=namedtuple(
    "Drivers",
    ["npp", "mrso", "tas"] #[3840,36,720]
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
        "ra_0",
        #"r_C_NWT_rh",
        #"r_C_AGWT_rh",
        #"r_C_TR_rh",
        #"r_C_GVF_rh",
        #"r_C_GVR_rh",
        "r_C_AGML_rh",
        "r_C_AGSL_rh",
        "r_C_AGMS_rh",
        "r_C_YHMS_rh",
        "k_C_BGDL",
        "k_C_BGRL",
        "k_C_BGMS",
        "k_C_SHMS",
        "r_C_AGML_2_C_AGMS",
        "r_C_AGMS_2_C_YHMS",
        "r_C_YHMS_2_C_AGMS",
        "r_C_YHMS_2_C_SHMS",
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
        "fwt",
        "fgv",
        "fco",
        "fml",
        "fd",
        "k_C_NWT",
        "k_C_AGWT",
        "k_C_TR",
        "k_C_GVF",
        "k_C_GVR",
        "f_C_AGSL_2_C_AGMS",
        "f_C_BGRL_2_C_SHMS",
        'C_NWT_0',#for the trendy data also the startvalues have to be estimated but 
        'C_AGWT_0',
        'C_GVF_0',
        'C_GVR_0',
        'C_AGML_0',
        'C_AGSL_0',
        'C_BGDL_0',
        'C_AGMS_0',
        'C_YHMS_0',
        'C_SHMS_0',
    ]
)

# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues. 
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise         
# to guard agains accidentally changed order we use a namedtuple again. Since B_func and u_func rely 
# on the correct ordering of the statevariables we build V dependent on this order 

#create a small model specific function that will later be stored in the file model_specific_helpers.py
#def download_my_TRENDY_output(conf_dict):
#    gh.download_TRENDY_output(
#        username=conf_dict["username"],
#        password=conf_dict["password"],
#        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
#        models=['ISAM'],
#        variables = Observables._fields + OrgDrivers._fields
#    )
def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['ISAM'],
        variables = Observables._fields + Drivers._fields
    )

def get_example_site_vars(dataPath):
    # pick up 1 site
    s = slice(None, None, None)  # this is the same as :
    lat=180
    lon=200
    t = s, lat, lon  # [t] = [:,lat,lon]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them
        ds = nc.Dataset(str(path))
        return ds.variables[vn][t]

    o_names=[(f,"ISAM_S2_{}.nc".format(f)) for f in Observables._fields]
    #d_names=[(f,"ISAM_S2_{}.nc".format(f)) for f in OrgDrivers._fields]
    d_names=[(f,"ISAM_S2_{}.nc".format(f)) for f in Drivers._fields]

    # we want to drive with npp and can create it from gpp and ra 
    # observables
    #odvs=OrgDrivers(*map(f,d_names))
    dvs=Drivers(*map(f,d_names))
    obss=Observables(*map(f, o_names))

    #dvs=Drivers(
    #    npp=odvs.gpp-obss.ra,
    #    mrso=odvs.mrso,
    #    tas=odvs.tas
    #)
    return (obss, dvs)

def nc_file_name(nc_var_name,experiment_name="ISAM_S2_"):
    return experiment_name+"{}.nc".format(nc_var_name)

def nc_global_mean_file_name(nc_var_name,experiment_name="ISAM_S2_"):
    return experiment_name+"{}_gm.nc".format(nc_var_name)

# +
def get_global_mean_vars(dataPath):
    o_names=Observables._fields
    d_names=Drivers._fields
    names = o_names + d_names 
    
    if all([dataPath.joinpath(nc_global_mean_file_name(vn)).exists() for vn in names]):
        print(""" Found cached global mean files. If you want to recompute the global means
            remove the following files: """
        )
        for vn in names:
            print( dataPath.joinpath(nc_global_mean_file_name(vn)))

        def get_cached_global_mean(vn):
            return gh.get_cached_global_mean(dataPath.joinpath(nc_global_mean_file_name(vn)),vn)
    
        return (
            Observables(*map(get_cached_global_mean, o_names)),
            Drivers(*map(get_cached_global_mean,d_names))
        )

    else:
        # we now check if any of the arrays has a time lime containing nan values 
        # APART FROM values that are already masked by the fillvalue
#         print("computing masks to exclude pixels with nan entries, this may take some minutes...")
#         def f(vn):
#             path = dataPath.joinpath(nc_file_name(vn))
#             ds = nc.Dataset(str(path))
#             #scale fluxes vs pools
#             var =ds.variables[vn]
#             return gh.get_nan_pixel_mask(var)

#         masks=[ f(name)    for name in names ]
#         # We compute the common mask so that it yields valid pixels for ALL variables 
#         combined_mask= reduce(lambda acc,m: np.logical_or(acc,m),masks)

        gm=gh.globalMask()
        # load an example file with mask
        template = nc.Dataset(dataPath.joinpath("ISAM_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
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
            lats= vs["lat"].__array__()
            lons= vs["lon"].__array__()
            print(vn)
            var=ds.variables[vn]
            # check if we have a cached version (which is much faster)
            gm_path = dataPath.joinpath(nc_global_mean_file_name(vn))

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
            return gm # * 86400 if vn in ["npp", "rh"] else gm
    
        #map variables to data
        return (
            Observables(*map(compute_and_cache_global_mean, o_names)),
            Drivers(*map(compute_and_cache_global_mean, d_names))
        )


# -

def make_npp_func(dvs):
    def func(day):
        month=gh.day_2_month_index(day)
        #print(day,month)
        # kg/m2/s kg/m2/day;
        return (dvs.npp[month]) * 86400

    return func


# +
#def make_xi_func(dvs):
#    def func(day):
#        return 1.0 # preliminary fake for lack of better data... 
#    return func
# -

import math
def make_xi_func(dvs):
    def xi_func(day):
        mi = gh.day_2_month_index(day)
        # alternative FT
        FT = 0.08 * math.exp(0.095 * (dvs.tas[mi]-273.15)) # temperature rate modifier
        FW = 1 #/ (1 + 30 * math.exp(-8.5 * dvs.mrso[mi])) # water rate modifier
        #print("FT,FW", FT, FW)
        rh_factor = FT * FW
        return rh_factor # 1.0     # Set to 1 if no scaler implemented
        # return 1.0

    return xi_func


def make_func_dict(mvs,dvs,cpa,epa):
    return {
        "NPP": make_npp_func(dvs),
        "xi": make_xi_func(dvs)
    }

# this function is deprecated - see general helpers traceability_iterator
# def make_traceability_iterator(mvs,dvs,cpa,epa):

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
        [V_init.ra,V_init.rh]
    ).reshape(n+2,1) #reshaping is neccessary for matmul (the @ in B @ X)
    

    
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
        
        ra = np.sum(
            [
              numOutFluxesBySymbol[Symbol(k)](it,*X)
                for k in ['C_NWT','C_AGWT','C_TR','C_GVF','C_GVR'] 
                if Symbol(k) in numOutFluxesBySymbol.keys()
            ]
        )
        rh = np.sum(
            [
                numOutFluxesBySymbol[Symbol(k)](it,*X)
                for k in ['C_AGML','C_AGSL','C_BGDL','C_BGRL','C_AGMS','C_YHMS','C_SHMS','C_BGMS'] 
                if Symbol(k) in numOutFluxesBySymbol.keys()
            ]
        )
        V_new = np.concatenate((X_new.reshape(n,1),np.array([ra,rh]).reshape(2,1)), axis=0)
        
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


def make_StartVector(mvs):
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()]+
        ["ra","rh"]
    ) 


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
    StartVector=make_StartVector(mvs)
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test
    seconds_per_day = 86400
#     def npp_func(day):
#         month=gh.day_2_month_index(day)
#         return dvs.npp[month] * seconds_per_day   # kg/m2/s kg/m2/day;
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
        V_init = StartVector(
            C_NWT=epa.C_NWT_0,
            C_AGWT=epa.C_AGWT_0,
            C_GVF=epa.C_GVF_0,
            C_GVR=epa.C_GVR_0,
            C_TR=cpa.cVeg_0-(epa.C_NWT_0 + epa.C_AGWT_0 + epa.C_GVF_0 + epa.C_GVR_0),
            C_AGML=epa.C_AGML_0,
            C_AGSL=epa.C_AGSL_0,
            C_BGDL=epa.C_BGDL_0,
            C_BGRL=cpa.cLitter_0-(epa.C_AGML_0 + epa.C_AGSL_0 + epa.C_BGDL_0),
            C_AGMS=epa.C_AGMS_0,
            C_YHMS=epa.C_YHMS_0,
            C_SHMS=epa.C_SHMS_0,
            C_BGMS=cpa.cSoil_0-(epa.C_AGMS_0 + epa.C_YHMS_0 + epa.C_SHMS_0),
            ra=cpa.ra_0,
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
        #for k,v in model_par_dict.items():
        #    print('"{}": {},'.format(k,v))

#         model_par_dict = {
#             'r_C_AGMS_rh':cpa.r_C_AGMS_rh,
#             'r_C_AGML_2_C_AGMS':cpa.r_C_AGML_2_C_AGMS,
#             'beta_TR':1-epa.fgv-epa.fwt,
#             'r_C_GVF_2_C_AGML':epa.k_C_GVF*(1-epa.fml),
#             'r_C_AGMS_2_C_YHMS':cpa.r_C_AGMS_2_C_YHMS,
#             'r_C_GVR_2_C_BGDL':epa.k_C_GVR*epa.fd,
#             'r_C_GVR_2_C_BGRL':epa.k_C_GVR*(1-epa.fd),
#             'r_C_YHMS_2_C_AGMS':cpa.r_C_YHMS_2_C_AGMS,
#             'r_C_BGMS_2_C_SHMS':cpa.r_C_YHMS_2_C_SHMS,
#             'r_C_AGSL_2_C_AGMS':cpa.r_C_AGSL_rh/0.7*epa.f_C_AGSL_2_C_AGMS,
#             'r_C_YHMS_rh':cpa.r_C_YHMS_rh,
#             'beta_GVF':epa.fgv*0.5,
#             'r_C_AGML_rh':cpa.r_C_AGML_rh,
#             'r_C_BGRL_rh':cpa.k_C_BGRL*epa.fco,
#             'r_C_AGSL_2_C_YHMS':cpa.r_C_AGSL_rh/0.7*(1-epa.f_C_AGSL_2_C_AGMS),
#             'r_C_AGWT_2_C_AGSL':epa.k_C_AGWT*1,
#             'r_C_TR_2_C_BGRL':epa.k_C_TR*(1-epa.fd),
#             'beta_NWT':epa.fwt*0.5,
#             'r_C_BGMS_rh':cpa.k_C_BGMS*epa.fco,
#             'r_C_GVF_2_C_AGSL':epa.k_C_GVF*epa.fml,
#             'r_C_BGRL_2_C_SHMS':cpa.k_C_BGRL*(1-epa.fco)*epa.f_C_BGRL_2_C_SHMS,
#             'r_C_BGDL_rh':cpa.k_C_BGDL*epa.fco,
#             'r_C_BGRL_2_C_BGMS':cpa.k_C_BGRL*(1-epa.fco)*(1-epa.f_C_BGRL_2_C_SHMS),
#             'r_C_SHMS_rh':cpa.k_C_SHMS*epa.fco,
#             'r_C_NWT_2_C_AGSL':epa.k_C_NWT*(1-epa.fml),
#             'r_C_TR_2_C_BGDL':epa.k_C_TR*epa.fd,
#             'r_C_AGSL_rh':cpa.r_C_AGSL_rh,
#             'beta_AGWT':epa.fwt*0.5,
#             'r_C_SHMS_2_C_BGMS':cpa.k_C_SHMS*(1-epa.fco),
#             'r_C_NWT_2_C_AGML':epa.k_C_NWT*epa.fml,
#             'r_C_BGDL_2_C_SHMS':cpa.k_C_BGDL*(1-epa.fco)
#         }
        
        func_dict=make_func_dict(mvs,dvs,cpa,epa)

        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2) 
        delta_t_val=15 
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
        for m in range(cpa.number_of_months):
            #dpm = days_per_month[ m % 12]  
            mra=0
            mrh=0
            for d in range(int(dpm/delta_t_val)):
                v = it_sym.__next__().reshape(n,)
                # actually the original datastream seems to be a flux per area (kg m^-2 ^-s )
                # at the moment the iterator also computes a flux but in kg^-2 ^day
            V=StartVector(*v)
            #from IPython import embed;embed()
            o=Observables(
                cVeg=float(V.C_NWT+V.C_AGWT+V.C_GVF+V.C_TR+V.C_GVR),
                cLitter=float(V.C_AGML+V.C_AGSL+V.C_BGDL+V.C_BGRL),
                cSoil=float(V.C_AGMS+V.C_YHMS+V.C_BGMS+V.C_SHMS),
                ra=V.ra/seconds_per_day,
                rh=V.rh/seconds_per_day # the data is per second while the time units for the iterator refer to days
            )
            sols.append(o)
            
        sol=np.stack(sols)       
        #convert to yearly output if necessary (the monthly pool data looks very "yearly")
        #sol_yr=np.zeros(int(cpa.number_of_months/12)*sol.shape[1]).reshape([int(cpa.number_of_months/12),sol.shape[1]])  
        #for i in range(sol.shape[1]):
        #   sol_yr[:,i]=monthly_to_yearly(sol[:,i])
        #sol=sol_yr
        return sol 
        
    return param2res

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


def numeric_X_0(mvs,dvs,cpa,epa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "C_NWT": epa.C_NWT_0,
        "C_AGWT": epa.C_AGWT_0,
        "C_GVF": epa.C_GVF_0,
        "C_GVR": epa.C_GVR_0,
        "C_TR": cpa.cVeg_0-(epa.C_NWT_0 + epa.C_AGWT_0 + epa.C_GVF_0 + epa.C_GVR_0),
        "C_AGML": epa.C_AGML_0,
        "C_AGSL": epa.C_AGSL_0,
        "C_BGDL": epa.C_BGDL_0,
        "C_BGRL": cpa.cLitter_0-(epa.C_AGML_0 + epa.C_AGSL_0 + epa.C_BGDL_0),
        "C_AGMS": epa.C_AGMS_0,
        "C_YHMS": epa.C_YHMS_0,
        "C_SHMS": epa.C_SHMS_0,
        "C_BGMS": cpa.cSoil_0-(epa.C_AGMS_0 + epa.C_YHMS_0 + epa.C_SHMS_0),
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)

    return X_0

# +
# # Define start and end dates of the simulation
# import datetime as dt
# start_date=dt.date(1700, 1, 16)
# end_date = dt.date(2019, 12, 16)

# +
# def make_sim_day_2_day_since_a_D(conf_dict):
#     # this function is extremely important to syncronise our results
#     # because our data streams start at different times the first day of 
#     # a simulation day_ind=0 refers to different dates for different models
#     # we have to check some assumptions on which this calculation is based
#     # for jules the data points are actually spaced monthly with different numbers of days
#     ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("ISAM_S2_cVeg.nc")))
#     times = ds.variables["time"]

#     # we have to check some assumptions on which this calculation is based

#     tm = times[0] #time of first observation in Months_since_1860-01 # print(times.units)
#     td = int(tm *31)  #in days since_1700-01-01 
#     #NOT assuming a 30 day month...
#     import datetime as dt
#     ad = dt.date(1, 1, 1) # first of January of year 1 
#     sd = dt.date(1700, 1, 16)
#     td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    

#     def f(day_ind: int)->int:
#         return day_ind+td_aD

#     return f
# -

def start_date():
    ## this function is important to syncronise our results
    ## because our data streams start at different times the first day of 
    ## a simulation day_ind=0 refers to different dates for different models
    ## we have to check some assumptions on which this calculation is based
    ## Here is how to get these values
    #ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("ISAM_S2_npp.nc")))
    #times = ds.variables["time"]
    #tm = times[0] #time of first observation in Months_since_1700-01 # print(times.units)
    #td = int(tm *30)  #in days since_1700-01-01 
    #import datetime as dt
    #ad = dt.date(1, 1, 1) # first of January of year 1 
    #sd = dt.date(1700, 1, 1)
    #td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    ## from td_aD (days since 1-1-1) we can compute the year month and day
    return gh.date(
        year=1700, 
        month=1,
        day=16
    )

################ function for computing global mean for custom data streams ###################

data_str = namedtuple( # data streams available in the model
    'data_str',
    ["cVeg", "cLitter", "cSoil", "gpp", "npp", "ra", "rh"]
    )
    
def get_global_mean_vars_all(experiment_name):
        return(
            gh.get_global_mean_vars_all(model_folder="cj_isam", 
                            experiment_name=experiment_name,
                            lat_var="lat",
                            lon_var="lon",
                            ) 
        )       
    
# def get_global_mean_vars_all(experiment_name="ISAM_S2_"):
    
    # def nc_file_name(nc_var_name, experiment_name="ISAM_S2_"):
        # return experiment_name+"{}.nc".format(nc_var_name)

    # def nc_global_mean_file_name(nc_var_name, experiment_name="ISAM_S2_"):
        # return experiment_name+"{}_gm_all.nc".format(nc_var_name)

    # data_str = namedtuple( # data streams available in the model
        # 'data_str',
        # ["cVeg", "cLitter", "cSoil", "gpp", "npp", "ra", "rh"]
        # )
        
    # names = data_str._fields
    # conf_dict = gh.confDict("cj_isam")
    # # with Path('config.json').open(mode='r') as f:
        # # conf_dict = frozendict(json.load(f))
    # dataPath=Path(conf_dict["dataPath"])    
    
    # if all([dataPath.joinpath(nc_global_mean_file_name(vn, experiment_name=experiment_name)).exists() for vn in names]):
        # print(""" Found cached global mean files. If you want to recompute the global means
            # remove the following files: """
        # )
        # for vn in names:
            # print( dataPath.joinpath(nc_global_mean_file_name(vn,experiment_name=experiment_name)))

        # def get_cached_global_mean(vn):
            # gm = gh.get_cached_global_mean(dataPath.joinpath(nc_global_mean_file_name(vn,experiment_name=experiment_name)),vn)
            # return gm * 86400 if vn in ["gpp", "npp", "rh", "ra"] else gm

        # #map variables to data
        # output=gh.data_streams(*map(get_cached_global_mean, data_str._fields))
        # return (
            # output
        # )

    # else:
        # gm=gh.globalMask()
        # # load an example file with mask
        # template = nc.Dataset(dataPath.joinpath("ISAM_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
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
            # gm_path = dataPath.joinpath(nc_global_mean_file_name(vn, experiment_name=experiment_name))

            # gm=gh.global_mean_var(
                    # lats,
                    # lons,
                    # gcm.index_mask,
                    # var
            # )
            # gh.write_global_mean_cache(
                    # gm_path,
                    # gm,
                    # vn
            # )
            # return gm * 86400 if vn in ["gpp", "npp", "rh", "ra"] else gm
        
        # #map variables to data
        # output=data_str(*map(compute_and_cache_global_mean, data_str._fields))
        # return (
            # gh.data_streams( # required data streams
                # cVeg=output.cVeg,
                # cSoil=output.cLitter+output.cSoil,
                # gpp=output.gpp,
                # npp=output.npp,
                # ra=output.ra,
                # rh=output.rh,
            # )
        # )