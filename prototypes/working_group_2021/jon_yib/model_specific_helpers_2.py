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

Observables_annual = namedtuple(
    'Observables_annual',
    ["cVeg", "cSoil"]
)
Observables_monthly = namedtuple(
    'Observables_monthly',
    ["rh" ]
)
Observables = namedtuple(
    "Observables",
    Observables_annual._fields+Observables_monthly._fields
)
Drivers=namedtuple(
    "Drivers",
    ["npp"]
)

Constants = namedtuple(
    "Constants",
    [
        'npp_0',       #Initial input/pools
        'rh_0',
        'c_veg_0',
        'c_soil_0',
        'clay',        #Constants like clay
        'silt',
        'nyears'       #Run time (years for my model)
    ]
)
EstimatedParameters = namedtuple(
    'EstimatedParameters', 
    [
        'c_leaf_0',               #Names: c_poolname_0
        'c_root_0',               #Only initial pools that are estimated
        'c_lit_cwd_0',
        'c_lit_met_0',
        'c_lit_str_0',
        'c_lit_mic_0',
        'c_soil_met_0',
        'c_soil_str_0',
        'c_soil_mic_0',
        'c_soil_slow_0',
        'beta_leaf',
        'beta_root',
        'r_c_leaf_rh',
        'r_c_root_rh',
        'r_c_wood_rh',
        'r_c_lit_cwd_rh',
        'r_c_lit_met_rh',
        'r_c_lit_str_rh',
        'r_c_lit_mic_rh',
        'r_c_soil_met_rh',
        'r_c_soil_str_rh',
        'r_c_soil_mic_rh',
        'r_c_soil_slow_rh',
        'r_c_soil_passive_rh',
        'r_c_leaf_2_c_lit_met',
        'r_c_leaf_2_c_lit_str',
        'r_c_root_2_c_soil_met',
        'r_c_root_2_c_soil_str',
        'r_c_wood_2_c_lit_cwd',
        'r_c_lit_cwd_2_c_lit_mic',
        'r_c_lit_cwd_2_c_soil_slow',
        'r_c_lit_met_2_c_lit_mic',
        'r_c_lit_str_2_c_lit_mic',
        'r_c_lit_str_2_c_soil_slow',
        'r_c_lit_mic_2_c_soil_slow',
        'r_c_soil_met_2_c_soil_mic',
        'r_c_soil_str_2_c_soil_mic',
        'r_c_soil_str_2_c_soil_slow',
        'r_c_soil_mic_2_c_soil_slow',
        'r_c_soil_mic_2_c_soil_passive',
        'r_c_soil_slow_2_c_soil_mic',
        'r_c_soil_slow_2_c_soil_passive',
        'r_c_soil_passive_2_c_soil_mic'
    ]
)


#when downloading data make sure model names match TRENDY server names:
#"CABLE-POP","CLASSIC","CLM5","DLEM","IBIS","ISAM","ISBA_CTRIP",
#"JSBACH","JULES-ES-1.0","LPJ-GUESS","LPJwsl","LPX-Bern","OCN",
#"ORCHIDEE","ORCHIDEE-CNP","ORCHIDEEv3","ORCHIDEEv3_0.5deg"
#"SDGVM","VISIT","YIBs"

#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['YIBs'],
        variables =Observables._fields + Drivers._fields
    )
#call it to test that the download works the data
#download_my_TRENDY_output()

def get_example_site_vars(dataPath):
    # pick up 1 site   
    s = slice(None, None, None)  # this is the same as :
    t = s, 74, 118  # [t] = [:,49,325]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them 
        ds = nc.Dataset(str(path))
        #check for npp/gpp/rh/ra to convert from kg/m2/s to kg/m2/day
        if vn in ["npp","gpp","rh","ra"]:
            for name, variable in ds.variables.items():            
                for attrname in variable.ncattrs():
                    print("{} -- {}".format(attrname, getattr(variable, attrname)))
            return ds.variables[vn][t]
        else:
            for name, variable in ds.variables.items():            
                for attrname in variable.ncattrs():
                    print("{} -- {}".format(attrname, getattr(variable, attrname)))
            return ds.variables[vn][t]

    # Link symbols and data:
    # YIBS has annual vs monthly file names so they are linked separately
    # If all your data is similarly named you can do this in one step

    # Create annual file names (single step if files similarly named)
    o_names=[(f,"YIBs_S2_Annual_{}.nc".format(f)) for f in Observables_annual._fields]

    # Create monthly file names (can remove if done in one step above)
    monthly_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Observables_monthly._fields]
    # Extend name list with monthly names
    o_names.extend(monthly_names)

    # create file names for Drivers
    d_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Drivers._fields]

    # Link symbols and data for Observables/Drivers
    return (Observables(*map(f, o_names)),Drivers(*map(f,d_names)))


def make_iterator_sym(
        mvs,
        V_init, 
        par_dict,
        func_dict,
        delta_t_val=1 # defaults to 1day timestep
    ):
    B_func, u_func = gh.make_B_u_funcs_2(mvs,par_dict,func_dict,delta_t_val)  

    sv=mvs.get_StateVariableTuple()
    n=len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct 
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr=np.array(
        [V_init.__getattribute__(str(v)) for v in sv]+
        [V_init.rh]
    ).reshape(n+1,1) #reshaping is neccessary for matmul (the @ in B @ X)

    numOutFluxesBySymbol={
        k: gh.numfunc(expr_cont, mvs, delta_t_val=delta_t_val, par_dict=par_dict, func_dict=func_dict) 
    
        for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
    } 
    def f(it,V):
        X = V[0:n]
        b = u_func(it,X)
        B = B_func(it,X)
        X_new = X + b + B @ X
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
        
        l=[
                numOutFluxesBySymbol[Symbol(k)](it,*X.reshape(n,))
                for k in ["c_lit_cwd","c_lit_met","c_lit_str","c_lit_mic","c_soil_met","c_soil_str","c_soil_mic","c_soil_slow","c_soil_passive"] 
                if Symbol(k) in numOutFluxesBySymbol.keys()
        ]
        rh = np.array(l).sum()
        V_new = np.concatenate(
            (
                X_new.reshape(n,1),
                np.array([rh]).reshape(1,1)
            )
            , axis=0
        )
        return V_new

    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )

def make_StartVector(mvs):
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()]+
        ["rh"]
    ) 


def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers,
    ) -> Callable[[np.ndarray], np.ndarray]:
    
    # Build iterator 
    # Need dictionary of numeric values for all parameters that are not state variables/time
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    
    # Time dependent driver function does not change with the estimated parameters
    # Defined once outside param2res function
    #seconds_per_day = 86400
    
    def npp_func(day):
        month=gh.day_2_month_index(day)
        return dvs.npp[month]*86400
    
    # Build environmental scaler function
    def xi_func(day):
        return 1.0     # Set to 1 if no scaler implemented 

    # Define function dictionary
    func_dict={
        'NPP':npp_func,
        'xi':xi_func
    }
    
    # Create namedtuple for initial values
    StartVector=make_StartVector(mvs)
    
    # Define actual forward simulation function
    def param2res(pa):
        
        # Parameter vector
        epa=EstimatedParameters(*pa)
         
        # Create a startvector for the iterator 
        V_init = StartVector(
            c_leaf=epa.c_leaf_0,
            c_root=epa.c_root_0,
            c_wood=cpa.c_veg_0-(epa.c_leaf_0 + epa.c_root_0),
            c_lit_cwd=epa.c_lit_cwd_0,
            c_lit_met=epa.c_lit_met_0,
            c_lit_str=epa.c_lit_str_0,
            c_lit_mic=epa.c_lit_mic_0,
            c_soil_met=epa.c_soil_met_0,
            c_soil_str=epa.c_soil_str_0,
            c_soil_mic=epa.c_soil_mic_0,
            c_soil_slow=epa.c_soil_slow_0,
            c_soil_passive=cpa.c_soil_0-(epa.c_lit_cwd_0+epa.c_lit_met_0+epa.c_lit_str_0+epa.c_lit_mic_0+epa.c_soil_met_0+epa.c_soil_str_0+epa.c_soil_mic_0+epa.c_soil_slow_0),
            rh=cpa.rh_0
        )
        
        # Parameter dictionary for the iterator
        apa = {**cpa._asdict(),**epa._asdict()}
        model_par_dict = {
            Symbol(k):v for k,v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }
        
        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2) 
        delta_t_val=10 #time step length in days
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
        )
        ######################################################################
        # calling the iterator and projecting the results
        def cVegF(V):
            return float(V.c_leaf+V.c_wood+V.c_root)
        
        def cSoilF(V): 
            return float(
                V.c_lit_cwd+
                V.c_lit_met+
                V.c_lit_str+
                V.c_lit_mic+ 
                V.c_soil_met+
                V.c_soil_str+
                V.c_soil_mic+
                V.c_soil_slow+
                V.c_soil_passive
            )
        
        #empty arrays for saving data
        cVeg_arr=np.zeros(cpa.nyears)
        cSoil_arr=np.zeros(cpa.nyears)
        rh_arr=np.zeros(cpa.nyears*12)
        
        dpm = 30                             # Set days per month
        im = 0
        steps_per_month=int(dpm/delta_t_val)
        steps_per_year=int(dpm/delta_t_val)*12
        # forward simulation by year
        for y in range(cpa.nyears):
            # reset the yearly average values to 0
            cVeg_avg= 0  
            cSoil_avg = 0
            for m in range(12):
                rh_avg=0
                for d in range(steps_per_month):    
                    V = StartVector(*it_sym.__next__())                  
                    rh_avg=V.rh/steps_per_month
                    cVeg_avg += cVegF(V)/steps_per_year
                    cSoil_avg += cSoilF(V)/steps_per_year
                rh_arr[im] = rh_avg/(24*60*60) #convert to kg/s from kg/day
                im += 1
                
            cVeg_arr[y] = cVeg_avg
            cSoil_arr[y] = cSoil_avg
        return Observables(cVeg=cVeg_arr,cSoil=cSoil_arr,rh=rh_arr)
    return param2res

def make_weighted_cost_func(
        obs: Observables
    ) -> Callable[[Observables],np.float64]:
    # first unpack the observation array into its parts
    #cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
    def costfunction(out_simu: Observables) -> np.float64:
        # fixme
        #   as indicated by the fact that the function lives in this
        #   model-specific module it is NOT apropriate for (all) other models.
        #   There are model specific properties:
        #   1.) The weight for the different observation streams
        #
        number_of_ys=out_simu.cVeg.shape[0]
        number_of_ms=out_simu.rh.shape[0]

        J_obj1 = np.mean (( out_simu.cVeg - obs.cVeg )**2)/(2*np.var(obs.cVeg))
        J_obj2 = np.mean (( out_simu.cSoil -  obs.cSoil )**2)/(2*np.var(obs.cSoil))

        J_obj3 = np.mean (( out_simu.rh - obs.rh )**2)/(2*np.var(obs.rh))

        J_new = (J_obj1 + J_obj2)+ J_obj3/12 #the 12 is purely conjectural
        return J_new
    return costfunction

def make_param_filter_func(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters 
        ) -> Callable[[np.ndarray], bool]:

    # find position of beta_leaf and beta_wood
    beta_leaf_ind=EstimatedParameters._fields.index("beta_leaf")
    beta_wood_ind=EstimatedParameters._fields.index("beta_root")

    def isQualified(c):
        beta_leaf_ind
        cond1 =  (c >= c_min).all() 
        cond2 =  (c <= c_max).all() 
        cond3 =  c[beta_leaf_ind]+c[beta_wood_ind] <=1  
        return (cond1 and cond2 and cond3)
        
    
    return isQualified


def make_traceability_iterator(mvs,dvs,cpa,epa):
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "c_leaf": epa.c_leaf_0,     
        "c_root": epa.c_root_0,     
        "c_wood": cpa.c_veg_0 - (epa.c_leaf_0 +  epa.c_root_0),  
        "c_lit_cwd": epa.c_lit_cwd_0,
        "c_lit_met": epa.c_lit_met_0,
        "c_lit_str": epa.c_lit_str_0,
        "c_lit_mic": epa.c_lit_mic_0,
        "c_soil_met": epa.c_soil_met_0,
        "c_soil_str": epa.c_soil_str_0,
        "c_soil_mic": epa.c_soil_mic_0,
        "c_soil_slow": epa.c_soil_slow_0,
        "c_soil_passive": cpa.c_soil_0 - (
                              epa.c_lit_cwd_0 
                            + epa.c_lit_met_0 
                            + epa.c_lit_str_0 
                            + epa.c_lit_mic_0 
                            + epa.c_soil_met_0 
                            + epa.c_soil_str_0 
                            + epa.c_soil_mic_0 
                            + epa.c_soil_slow_0
                        )
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    fd=make_func_dict(mvs,dvs)
    V_init = gh.make_InitialStartVectorTrace(
            X_0,mvs,
            par_dict=par_dict,
            func_dict=fd
    )
    it_sym_trace = gh.make_daily_iterator_sym_trace(
        mvs,
        V_init=V_init,
        par_dict=par_dict,
        func_dict=fd
    )
    return it_sym_trace

def make_npp_func(dvs):
    def func(day):
        month=gh.day_2_month_index(day)
        # kg/m2/s kg/m2/day;
        return (dvs.npp[month]) * 86400

    return func


def make_xi_func(dvs):
    def func(day):
        return 1.0 # preliminary fake for lack of better data... 
    return func


def make_func_dict(mvs,dvs):
    return {
        "NPP": make_npp_func(dvs),
        "xi": make_xi_func(dvs)
    }


