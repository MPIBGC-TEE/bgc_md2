from importlib import import_module
from collections import namedtuple
from ComputabilityGraphs import CMTVS
from typing import Callable, Dict, Tuple
from pathlib import Path
import numpy as np

from CompartmentalSystems.ArrayDictResult import ArrayDictResult
import bgc_md2.helper as h

from  ... import general_helpers as gh
from .. import model_specific_helpers_2 as msh  

model_mod = 'bgc_md2.models.kv_visit2'
cp_mod=import_module(f"{model_mod}.CachedParameterization")
make_func_dict = cp_mod.make_func_dict
mvs=import_module(f"{model_mod}.source").mvs
Drivers=cp_mod.Drivers

# this is what is read from cpa.json 
FreeConstants = namedtuple(
    "FreeConstants",
    []
)    
# this is what is inferred but stays constant during parameter estimation 
DerivedConstants = namedtuple(
    "Derived",
    [
        "cLitter_0",
        "cSoil_0",
        "cVeg_0",
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
        cVeg_0=svs.cVeg[0],
        cSoil_0=svs.cSoil[0],
        cLitter_0=svs.cLitter[0],
    )
    return Constants(*fcpa, *dcpa)

# Note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated)
# parameters. In this case we may use only the first entry e.g. to derive startvalues and their length (number of months)
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
DerivedEstimatedParameters = namedtuple(
     # These should actually not be set by hand but by a function because there
     # are conditions which could be easily violated (resulting in the 
     # wrong kind of proposals (the betas could go here) 
    "DerivedEstimatedParameters",
    [ ]
)
EstimatedParameters = namedtuple(
    "EstimatedParameters",
    [*FreeEstimatedParameters._fields,*DerivedEstimatedParameters._fields] 
)

def epa_min(
    fepa: FreeEstimatedParameters,
    dvs: msh.Drivers,
    svs: msh.Observables
) -> EstimatedParameters:
    # this function could set the betas 
    # but in this scheme it is only used to pass on the
    # values from the epa_0.json file
    depa = DerivedEstimatedParameters()
    return EstimatedParameters(*fepa,*depa)

def epa_max(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # this function could set the betas 
    # but in this scheme it is only used to pass on the
    # values from the epa_0.json file
    depa = DerivedEstimatedParameters()
    return EstimatedParameters(*fepa,*depa)

def epa_0(
    fepa: FreeEstimatedParameters,
    dvs: msh.Drivers,
    svs: msh.Observables
) -> EstimatedParameters:
    # this function could set the betas 
    # but in this scheme it is only used to pass on the
    # values from the epa_0.json file
    depa = DerivedEstimatedParameters()
    return EstimatedParameters(*fepa, *depa)


def make_proposer(
    c_max: EstimatedParameters,
    c_min: EstimatedParameters, 
    fcpa: FreeConstants,
    dvs: msh.Drivers,
    svs: msh.Observables
) -> Callable[[np.ndarray, float, bool], np.ndarray]:
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
        (["beta_leaf", "beta_wood"],1),
    ]
    #return GenerateParamValues
    return gh.make_dirichlet_uniform_proposer(
        dirichlet_tups=dirichlet_tups,
        EstimatedParameters=EstimatedParameters,
        c_min=c_min,
        c_max=c_max
    )


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
            sum(map(value, ["C_leaf_0", "C_wood_0"])) <= svs.cVeg[0],
            sum(map(value, ["C_leaf_litter_0", "C_wood_litter_0"])) <= svs.cLitter[0],
            sum(
                map(
                    value, ['C_soil_fast_0', 'C_soil_slow_0']
                )
            ) <= svs.cSoil[0],
        ]
        res = all(conds)
        if print_conds:
            if not res:
                print(conds)
        return res
        
    return isQualified


def make_param2res_sym(
    mvs: CMTVS,
    fcpa: FreeConstants,
    dvs: Drivers,
    svs: msh.Observables
) -> Callable[[np.ndarray], np.ndarray]:

    cpa_v = cpa(fcpa, dvs, svs)
    def param2res(pa):

        epa = EstimatedParameters(*pa)
        apa = {**cpa_v._asdict(), **epa._asdict()}
        X_0 = numeric_X_0(mvs, dvs, apa)
        #dpm = h.date.days_per_month
        #steps_per_month = 2
        #delta_t_val = dpm/steps_per_month 

        par_dict = gh.make_param_dict(mvs, cpa_v, epa)
        func_dict = make_func_dict(dvs , cpa=cpa_v, epa=epa)
        return msh.synthetic_observables(
            mvs,
            X_0,
            par_dict=par_dict,
            func_dict=func_dict,
            dvs=dvs
        )

    return param2res

def numeric_X_0(mvs, dvs, apa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    par_dict = gh.make_param_dict2(mvs, apa)
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

