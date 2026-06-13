from importlib import import_module
from collections import namedtuple
from ComputabilityGraphs import CMTVS
from typing import Callable, Dict, Tuple
from pathlib import Path
import numpy as np

from CompartmentalSystems.ArrayDictResult import ArrayDictResult
import CompartmentalSystems.helpers_reservoir as hr
import bgc_md2.helper as h

"""
Here we implement a parameter estimation method based on the assumption that we
start in "equilibrium" to estimate parameters from  the trendy9 S2 experiments.
While the definition of equilibrium is not further specified in the trendy9
description we take it literally here in the mathematical sense and formalize it in the following way:
We (have to) make the  assumption that all drivers (time dependent parameters
in the model description) had been constant for an infinately long past before
t_0 when we start the reconstructed model.  While this is wildly unrealistic
for meausured data, it is of course a possible scenario for the creation of the
trendy S2 output.(although the trendy description is not specific enough to be
sure)

The nonautonomous system becomes thus autonomous and the mathematical concept
of equilibrium becomes applicable.  While in general (for non linear models)
mathematics does not quarantee the existence of equilibria even for this
autonomous case, we are in a better position  for linear output connected
compartmental systems (to which this method should be restricted) A unique
equilibrium can be shown to exist and computed by the inversion of the
compartmental matrix multiplied by the (also constant) Input B(t_0)^{-1}*I(t_0).
(Rasmussen 2016)

Note that:
- This assumption allows us to remove all startvalues from the estimated
  parameters (They become deterministically dependent on those parameters that
  control input and compartmental matrix at t=t0).  
- This also implies that the startvalues might actually contradict the
  "observables" of aggregated values cLitter,cSoil,cVeg at the starting point
  t_0 
- The question remains as to WHICH VALUES we set the drivers for the infinite
  past preceeding t_0 Heuristically the constant driver value should be as
  representative as possible for the infinite past before we start the model".
  To capture at least the yearly cycle we could compute the yearly averages of
  the drivers and assume a model start t_0 when these values were actually
  reached for the first time in the first year of output probably at some time
  in spring and hopefully not at too different times for different drivers)

"""
from  ... import general_helpers as gh
from .. import model_specific_helpers_2 as msh  

model_mod = 'bgc_md2.models.kv_visit2'
cp_mod=import_module(f"{model_mod}.CachedParameterization")
mvs=import_module(f"{model_mod}.source").mvs
make_func_dict = cp_mod.make_func_dict
#Drivers=cp_mod.Drivers

FreeConstants = namedtuple(
    "FreeConstants",
    [] 
    #funny you actually don't have any free constants, but you 
    # could transfer some of the FreeEstimatedParameters here if 
    # you have a good guess for their value
)    

# this is what is read from cpa.json 
FreeConstants = namedtuple(
    "FreeConstants",
    []
    # funny! 
    # You actually don't have any free constants, but you 
    # could transfer some of the FreeEstimatedParameters here if 
    # you have a good guess for their value
)    
DerivedConstants = namedtuple(
    "DerivedConstants",
    [
        "npp_0",
        "rh_0"
    ],
)
Constants = namedtuple(
    "Constants",
    [*FreeConstants._fields,*DerivedConstants._fields] 
)

def cpa(
        fcpa: FreeConstants,
        dvs: msh.Drivers,
        svs: msh.Observables
    ) -> Constants:
    dcpa = DerivedConstants(
        npp_0=dvs.npp[0],
        rh_0=svs.rh[0],
    )
    return Constants(*fcpa, *dcpa)

# Note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated)
# parameters. 
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise
# It is better to start with only a few
FreeEstimatedParameters = namedtuple(
    "FreeEstimatedParameters",
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
    ],
)
DerivedEstimatedParameters = namedtuple(
    "DerivedEstimatedParameters",
    []
)
EstimatedParameters = namedtuple(
    "EstimatedParameters",
    [*FreeEstimatedParameters._fields,*DerivedEstimatedParameters._fields] 
)

def epa_min(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # this function can obviously be changed
    # it avoids contradictions between the startvalues and 
    # data 
    depa = DerivedEstimatedParameters()
    return EstimatedParameters(*fepa,*depa)

def epa_max(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # this function can obviously be changed
    # it avoids contradictions between the startvalues and 
    # data 
    depa = DerivedEstimatedParameters()
    return EstimatedParameters(*fepa,*depa)

def epa_0(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # this function can obviously be changed
    # it avoids contradictions between the startvalues and 
    # data 
    depa = DerivedEstimatedParameters()
    return EstimatedParameters(*fepa,*depa)

def make_param_filter_func(
    c_max: EstimatedParameters,
    c_min: EstimatedParameters,
    fcpa: FreeConstants,
    dvs: msh.Drivers,
    svs: msh.Observables
) -> Callable[[np.ndarray], bool]:

    def isQualified(
            c,
            print_conds=False
        ):
        epa=EstimatedParameters(*c)
        apa = {**fcpa._asdict(), **epa._asdict()}
        
        def value(field_name):
            return apa[field_name]
        
        conds=[
            (c >= c_min).all(), 
            (c <= c_max).all(), 
            sum(map(value, ["beta_leaf", "beta_wood"])) <= 0.99,
        ]
        res = all(conds)

        if print_conds:
            if not res:
                print(conds)

        return res
        
    return isQualified

def make_proposer(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters, 
        fcpa: FreeConstants,
        dvs: msh.Drivers,
        svs: msh.Observables,
    ) -> Callable[[np.ndarray,float,bool], np.ndarray]:
    
    cpa_val = cpa(fcpa,dvs,svs)
    """Returns a function that will be used by the mcmc algorithm to propose
    a new parameter value tuple based on a given one.
    The two arrays c_max and c_min define the boundaries
    of the n-dimensional rectangular domain for the parameters and must be of
    the same shape.  After a possible parameter value has been sampled the
    filter_func will be applied to it to either accept or discard it.  So
    filter func must accept parameter array and return either True or False
    :param c_max: array of maximum parameter values
    :param c_min: array of minimum parameter values
    :param D: a damping parameter to regulate the proposer step 
     larger D means smaller step size. If the maximum distance for a new value is
     max_dist then the proposer will make a max_dist/ 
    :param filter_func: model-specific function to filter out impossible 
    parameter combinations
    """


    g = np.random.default_rng()
    # filter out the startvalues and use a dirichlet distribution 
    # for them  
    dirichlet_tups= [
        (["beta_leaf", "beta_wood"],1)
    ]
    #return GenerateParamValues
    return gh.make_dirichlet_uniform_proposer(
        dirichlet_tups=dirichlet_tups,
        EstimatedParameters=EstimatedParameters,
        c_min=c_min,
        c_max=c_max
    )

def make_param2res_sym(
    mvs: CMTVS,
    fcpa: FreeConstants,
    dvs: msh.Drivers,
    svs: msh.Observables,
) -> Callable[[np.ndarray], np.ndarray]:

    # select the first value for the drivers (to represent the past)
    # (alternatively compute average of the drivers for the first year or chose a midyear value)
    func_dict_const={ 
        "TAS": lambda t: dvs.tas[0],
        "mrso": lambda t: dvs.mrso[0],
        "NPP": lambda t: dvs.npp[0],
    }
    cpa_v=cpa(fcpa, dvs, svs)
    def param2res(pa):

        epa = EstimatedParameters(*pa)

        par_dict = gh.make_param_dict(mvs, cpa_v, epa)
        func_dict = make_func_dict(dvs , cpa=cpa_v, epa=epa)
        X_0 = numeric_X_0_internal(mvs, par_dict, func_dict_const)
        return msh.synthetic_observables(
            mvs,
            X_0,
            par_dict=par_dict,
            func_dict=func_dict,
            dvs=dvs
        )

    return param2res

# fixme mm 5-15 2023
# hopefully disposable in the future
def numeric_X_0(mvs, dvs, apa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    par_dict = gh.make_param_dict2(mvs, apa)
    func_dict=make_func_dict(dvs)
    return numeric_X_0_internal(mvs,par_dict,func_dict)


def numeric_X_0_internal(mvs, par_dict, func_dict,t_0=0):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    state_vector = mvs.get_StateVariableTuple()
    t = mvs.get_TimeSymbol()

    B_func = hr.numerical_array_func(
        state_vector=state_vector,
        time_symbol=t,
        expr=mvs.get_CompartmentalMatrix(),
        parameter_dict=par_dict,
        func_dict=func_dict,
    )
    I_func = hr.numerical_1d_vector_func(
        state_vector=state_vector,
        time_symbol=t,
        expr=mvs.get_InputTuple(),
        parameter_dict=par_dict,
        func_dict=func_dict,
    )
    # pseudo x0 =(0,0,...,0) is ok since for linear models
    # B_func and I_func do not really depend on x
    px0=np.zeros(shape=(len(state_vector),))
    #t_0=90 # for now should actually be 30 march (should be a value when the actual driver value should be close to the yearly average of the drivers 
    B_0=B_func(t_0, px0)
    I_0=I_func(t_0, px0)
    
    X_0 = (np.linalg.inv(B_0) @ (-I_0)).reshape(
        len(state_vector), 1
    )
    return X_0
