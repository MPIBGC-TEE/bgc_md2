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
from general_helpers import month_2_day_index, monthly_to_yearly
from functools import reduce

sys.path.insert(0,'..') # necessary to import general_helpers
from general_helpers import download_TRENDY_output, day_2_month_index, make_B_u_funcs_2

# we will use the trendy output names directly in other parts of the output
Observables = namedtuple(
    'Observables',
    ["cVeg","cLitter","cSoil","rh","ra"]
)
Drivers=namedtuple(
    "Drivers",
    ["gpp", "mrso", "tas"]
)    
from source import mvs
StartVector = namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()]+
        ["ra","rh"]
    ) 
# As a safety measure we specify those parameters again as 'namedtuples', which are like a mixture of dictionaries and tuples
# They preserve order as numpy arrays which is great (and fast) for the numeric computations
# and they allow to access values by keys (like dictionaries) which makes it difficult to accidentally mix up values.

UnEstimatedParameters = namedtuple(
    "UnEstimatedParameters",
    [
        "cLitter_0",
        "cSoil_0",
        "cVeg_0",
        "gpp_0",
        "rh_0",
        "ra_0",
        "r_C_root_litter_2_C_soil_passive",# here  we pretend to know these two rates 
        "r_C_root_litter_2_C_soil_slow",# it would be much better to know more  
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
        "T_0",
        "E",
        "KM",
        #"r_C_leaf_rh",
        #"r_C_wood_rh",
        #"r_C_root_rh",
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
        'C_leaf_0',#for the trendy data also the startvalues have to be estimated but 
        'C_wood_0',
        #C_root_0 can be inferred as cVeg_0-(C_leaf_0+C_wood_0)
        'C_leaf_litter_0',
        'C_wood_litter_0',
        #C_root_litter_0 can be inferred
        'C_soil_fast_0',
        'C_soil_slow_0',
        #C_soil_passive_0 can be inferred 
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
    download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['VISIT'],
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
        return ds.variables[vn][t]

    o_names=[(f,"VISIT_S2_{}.nc".format(f)) for f in Observables._fields]
    d_names=[(f,"VISIT_S2_{}.nc".format(f)) for f in Drivers._fields]
    return (Observables(*map(f, o_names)),Drivers(*map(f,d_names)))

def make_npp_func(dvs,svs):
    def func(day):
        month=day_2_month_index(day)
        # kg/m2/s kg/m2/day;
        return (dvs.gpp[month]-svs.ra[month]) * 86400

    return func


def make_xi_func(dvs,svs):
    def func(day):
        return 1.0 # preliminary fake for lack of better data... 
    return func


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
    B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict,delta_t_val)  
    sv=mvs.get_StateVariableTuple()
    #mass production of output functions

    def numfunc(expr_cont,delta_t_val):
        # build the discrete expression (which depends on it,delta_t instead of
        # the continius one that depends on t (TimeSymbol))
        it=Symbol("it")           #arbitrary symbol for the step index )
        t=mvs.get_TimeSymbol()
        delta_t=Symbol('delta_t')
        expr_disc = expr_cont.subs({t:delta_t*it})
        return hr.numerical_function_from_expression(
            expr=expr_disc.subs({delta_t:delta_t_val}),
            tup=(it, *mvs.get_StateVariableTuple()),
            parameter_dict=par_dict,
            func_set=func_dict
        )

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
        k:numfunc(expr_cont,delta_t_val=delta_t_val) 
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
                for k in ["C_leaf","C_wood","C_root"] 
                if Symbol(k) in numOutFluxesBySymbol.keys()
            ]
        )
        rh = np.sum(
            [
                numOutFluxesBySymbol[Symbol(k)](it,*X)
                for k in ["C_leaf_litter","C_wood_litter","C_root_litter","C_soil_fast","C_soil_slow","C_soil_passive"] 
                if Symbol(k) in numOutFluxesBySymbol.keys()
            ]
        )
        V_new = np.concatenate((X_new.reshape(n,1),np.array([ra,rh]).reshape(2,1)), axis=0)
        
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )

def make_param2res_sym(
        cpa: UnEstimatedParameters,
        dvs: Drivers,
        svs: Observables
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
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test
    seconds_per_day = 86400
    def npp_func(day):
        month=day_2_month_index(day)
        return (dvs.gpp[month]-svs.ra[month]) * seconds_per_day   # kg/m2/s kg/m2/day;
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
        V_init = StartVector(
            C_leaf=epa.C_leaf_0,
            C_wood=epa.C_wood_0,
            C_root=cpa.cVeg_0-(epa.C_leaf_0 + epa.C_wood_0),
            C_leaf_litter=epa.C_leaf_litter_0,
            C_wood_litter=epa.C_wood_litter_0,
            C_root_litter=cpa.cLitter_0-(epa.C_leaf_litter_0 + epa.C_wood_litter_0),
            C_soil_fast=epa.C_soil_fast_0,
            C_soil_slow=epa.C_soil_slow_0,
            C_soil_passive=cpa.cSoil_0-(epa.C_soil_fast_0 + epa.C_soil_slow_0),
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

        #print(model_par_dict)
        #from IPython import embed;embed()
        
        # Beside the par_dict the iterator also needs the python functions to replace the symbolic ones with
        # our fake xi_func could of course be defined outside of param2res but in general it
        # could be build from estimated parameters and would have to live here...
        def xi_func(day):
            return 1.0 # preliminary fake for lack of better data... 
    
        func_dict={
            'NPP':npp_func,
             'xi':xi_func
        }
        
        ##mass production of output functions
        #def numfunc(expr_cont,delta_t_val):
        #    # build the discrete expression (which depends on it,delta_t instead of
        #    # the continius one that depends on t (TimeSymbol))
        #    it=Symbol("it")           #arbitrary symbol for the step index )
        #    t=mvs.get_TimeSymbol()
        #    delta_t=Symbol('delta_t')
        #    expr_disc = expr_cont.subs({t:delta_t*it})
        #    return hr.numerical_function_from_expression(
        #        expr=expr_disc.subs({delta_t:delta_t_val}),
        #        tup=(it, *mvs.get_StateVariableTuple()),
        #        parameter_dict=par_dict,
        #        func_set=func_dict
        #    )
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
        for m in range(cpa.number_of_months):
            #dpm = days_per_month[ m % 12]  
            dpm=30 # 
            mra=0
            mrh=0
            for d in range(int(dpm/delta_t_val)):
                v = it_sym.__next__()
                #mra +=v[9,0] 
                #mrh +=v[10,0]
                # actually the original datastream seems to be a flux per area (kg m^-2 ^-s )
                # at the moment the iterator also computes a flux but in kg^-2 ^day
            V=StartVector(*v)
            o=Observables(
                cVeg=float(V.C_leaf+V.C_wood+V.C_root),
                cLitter=float(V.C_leaf_litter+V.C_wood_litter+V.C_root_litter),
                cSoil=float(V.C_soil_fast+V.C_soil_slow+V.C_soil_passive),
                ra=v[9,0]/seconds_per_day,
                rh=v[10,0]/seconds_per_day # the data is per second while the time units for the iterator refer to days
                #ra=mra/dpm, # monthly respiration back to kg/m2/day units
                #rh=mrh/dpm, # monthly respiration back to kg/m2/day units
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
