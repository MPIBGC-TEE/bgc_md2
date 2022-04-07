# +
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
# -

# we will use the trendy output names directly in other parts of the output
Observables = namedtuple(
    'Observables',
    ["cVeg", "cLitter", "cSoil", "rh" ]
)
Drivers=namedtuple(
    "Drivers",
    ["npp","mrso"]
)    
  

Constants = namedtuple(
    "Constants",
    [
        "cVeg_0",
        "cLitter_0",
        "cSoil_0",
        "npp_0",
        "rh_0",
       # "mrso",
       # "tsl",
        "number_of_years" # necessary to prepare the output in the correct length 
    ]
)

EstimatedParameters = namedtuple(
    "EstimatedParameters",
    [
        "beta_leaf",
        "beta_wood",
        "Theta_sat",
        "Theta_fc",
        "r_C_leaf_rh",
        "r_C_wood_rh",
        "r_C_root_rh",
        "r_C_aom1_rh",
        "r_C_aom2_rh",
        "r_C_smb1_rh",
        "r_C_smb2_rh",
        "r_C_smr_rh",
        "r_C_nom_rh",
        "r_C_dom_rh",
        "r_C_psom_rh",
        "r_C_leaf_2_C_aom1",
        "r_C_leaf_2_C_aom2",
        "r_C_wood_2_C_aom1",
        "r_C_wood_2_C_aom2",
        "r_C_root_2_C_aom1",
        "r_C_root_2_C_aom2",
        "r_C_aom1_2_C_smb1",
        "r_C_aom1_2_C_smb2",
        "r_C_aom1_2_C_nom",
        "r_C_aom1_2_C_dom",
        "r_C_aom2_2_C_smb1",
        "r_C_aom2_2_C_smb2",
        "r_C_aom2_2_C_dom",
        "r_C_smb1_2_C_nom",
        "r_C_smb1_2_C_psom",
        "r_C_smb2_2_C_smr",
        "r_C_smr_2_C_smb1",
        "r_C_nom_2_C_smb1",
        "r_C_nom_2_C_dom",
        "r_C_nom_2_C_psom",
        "r_C_dom_2_C_smb1",
        "r_C_dom_2_C_nom",
        "r_C_psom_2_C_smb1",
        "C_leaf_0",
        "C_wood_0",
        "C_aom1_0",
        "C_smb1_0",
        "C_smb2_0",
        "C_smr_0",
        "C_nom_0",
        "C_dom_0"
    ]
)


def download_my_TRENDY_output():
    download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['DLEM'],
        variables = Observables._fields + Drivers._fields
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
        if vn in ["npp", "rh"]:
            return ds.variables[vn][t]*86400
        else:
            return ds.variables[vn][t]

    o_names=[(f,"DLEM_S2_{}.nc".format(f)) for f in Observables._fields]
    d_names=[(f,"DLEM_S2_{}.nc".format(f)) for f in Drivers._fields]
    return (Observables(*map(f, o_names)),Drivers(*map(f,d_names)))


def get_globalmean_vars(dataPath):
    # pick up 1 site   
    #s = slice(None, None, None)  # this is the same as :
    #t = s, 50, 33  # [t] = [:,49,325]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them 
        ds = nc.Dataset(str(path))
        #print(ds.variables)
        # access lat/long of netCDF file
        lats= ds.variables["lat"]
        lons= ds.variables["lon"]
        #scale fluxes vs pools
        if vn in ["npp", "rh"]:
            return gh.global_mean_JULES(lats,lons,ds.variables[vn])*86400
        else:
            return gh.global_mean_JULES(lats,lons,ds.variables[vn])
    #map variables to data
    o_names=[(f,"DLEM_S2_{}.nc".format(f)) for f in Observables._fields]
    d_names=[(f,"DLEM_S2_{}.nc".format(f)) for f in Drivers._fields]
    return (Observables(*map(f, o_names)),Drivers(*map(f,d_names)))


def make_StartVector(mvs):
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()]+
        ["rh"]
    ) 


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
        k:gh.numfunc(expr_cont, mvs, delta_t_val=delta_t_val, par_dict=par_dict, func_dict=func_dict) 
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
                for k in ["C_aom1","C_aom2","C_smb1","C_smb2","C_smr","C_nom","C_dom","C_psom"] 
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


def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers,
    ) -> Callable[[np.ndarray], np.ndarray]: 
    # To compute numeric solutions we will need to build and iterator 
    # as we did before. As it will need numeric values for all the parameters 
    # we will have to create a complete dictionary for all the symbols
    # exept those for the statevariables and time.
    # This set of symbols does not change during the mcmc procedure, since it only
    # depends on the symbolic model.
    # Therefore we create it outside the mcmc loop and bake the result into the 
    # param2res function.
    # The iterator does not care if they are estimated or not it only 
    # wants a dictionary containing key: value pairs for all
    # parameters that are not state variables or the symbol for time
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    # Create namedtuple for initial values
    StartVector=make_StartVector(mvs)
    
    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test
    def npp_func(day):
        month=gh.day_2_month_index(day)
        return dvs.npp[month]   # kg/m2/s kg/m2/day
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
         # Parameter dictionary for the iterator
        apa = {**cpa._asdict(),**epa._asdict()}
        model_par_dict = {
            Symbol(k):v for k,v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }
        
        # Create a startvector for the iterator 
        V_init = StartVector(
            C_leaf = apa['C_leaf_0'],
            C_wood = apa['C_wood_0'],
            C_root = apa['cVeg_0'] - (
                apa['C_leaf_0'] + 
                apa['C_wood_0']
            ),
            C_aom1 = apa['C_aom1_0'],
            C_aom2 = apa['cLitter_0'] - apa['C_aom1_0'],
            C_smb1 = apa['C_smb1_0'],
            C_smb2 = apa['C_smb2_0'],
            C_smr = apa['C_smr_0'],
            C_nom = apa['C_nom_0'],
            C_dom = apa['C_dom_0'],
            C_psom = apa['cSoil_0'] - (
                apa['C_smb1_0'] +
                apa['C_smb2_0'] +
                apa['C_smr_0'] +
                apa['C_nom_0'] +
                apa['C_dom_0'] 
            ),
            rh = apa['rh_0']
        )

        # Beside the par_dict the iterator also needs the python functions to replace the symbolic ones with
        # our fake xi_func could of course be defined outside of param2res but in general it
        # could be build from estimated parameters and would have to live here...
        def xi_func(day):
            return 1.0 # preliminary fake for lack of better data... 
    
        func_dict={
            'NPP':npp_func,
             'xi':xi_func
        }
        
        delta_t_val=15
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
        )
        
        def cVegF(V):
            return float(V.C_leaf+V.C_wood+V.C_root)
        
        def cLitF(V):
            return float(V.C_aom1+V.C_aom2)
        
        def cSoilF(V): 
            return float(V.C_smb1+V.C_smb2+V.C_smr+V.C_nom+V.C_dom+V.C_psom)
        
        # Now that we have the iterator we can start to compute.
        # the first thing to notice is that we don't need to store all daily values,
        # since the observations are recorded monthly while our iterator has a
        # daily timestep.
        # - for the pool contents we only need to record the last day of the month
        # - for the respiration(s) ra and rh we have to sum up the daily values 
        #   over a month
        # 
        nyears = int(cpa.number_of_years)
        #empty arrays for saving data
        cVeg_arr=np.zeros(nyears)
        cLit_arr=np.zeros(nyears)
        cSoil_arr=np.zeros(nyears)
        rh_arr=np.zeros(nyears*12)

        #constants for forward simulation                   
        im = 0
        dpm = 30
        steps_per_month=int(dpm/delta_t_val)
        steps_per_year=int((dpm/delta_t_val)*12)
        
        #forward sim
        for y in range(nyears):
            cVeg_avg= 0  
            cSoil_avg = 0
            cLit_avg = 0
            for m in range(12): 
                rh_avg = 0
                for d in range(steps_per_month):    
                    V = StartVector(*it_sym.__next__())                  
                    rh_avg += V.rh
                    cVeg_avg += cVegF(V)
                    cLit_avg += cLitF(V)
                    cSoil_avg += cSoilF(V)
                rh_arr[im] = rh_avg/steps_per_month
                im += 1
            cVeg_arr[y] = cVeg_avg/steps_per_year
            cLit_arr[y] = cLit_avg/steps_per_year
            cSoil_arr[y] = cSoil_avg/steps_per_year
            #if y == 100:
            #    print(V)
        return Observables(
            cVeg = cVeg_arr,
            cLitter = cLit_arr,
            cSoil = cSoil_arr,
            rh = rh_arr)
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


def make_weighted_cost_func(
        obs: Observables
    ) -> Callable[[Observables],np.float64]:
    
    # Define cost function 
    def costfunction(out_simu: Observables) -> np.float64:
   
        #calculate costs for each data stream
        J_obj1 = (100/obs.cVeg.shape[0]) * np.sum((out_simu.cVeg - obs.cVeg)**2, axis=0) / (obs.cVeg.mean(axis=0)**2)
        J_obj2 = (100/obs.cLitter.shape[0]) * np.sum((out_simu.cLitter -  obs.cLitter)**2, axis=0) / (obs.cLitter.mean(axis=0)**2)
        J_obj3 = (100/obs.cSoil.shape[0]) * np.sum((out_simu.cSoil - obs.cSoil)**2, axis=0) / (obs.cSoil.mean(axis=0)**2)
        J_obj4 = (100/obs.rh.shape[0]) * np.sum((out_simu.rh - obs.rh)**2, axis=0) / (obs.rh.mean(axis=0)**2)
        
        # sum costs
        J_new = (J_obj1 + J_obj2 + J_obj3 + J_obj4)
        return J_new
    return costfunction
