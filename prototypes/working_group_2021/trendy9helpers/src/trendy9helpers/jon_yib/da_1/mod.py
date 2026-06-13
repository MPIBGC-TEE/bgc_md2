from importlib import import_module
from collections import namedtuple
from ComputabilityGraphs import CMTVS
from typing import Callable, Dict, Tuple
from pathlib import Path
import numpy as np
from copy import copy

from CompartmentalSystems.ArrayDictResult import ArrayDictResult
import bgc_md2.helper as h

from  ... import general_helpers as gh
from .. import model_specific_helpers_2 as msh  

model_mod = 'bgc_md2.models.jon_yib'
cp_mod=import_module(f"{model_mod}.CachedParameterization")
mvs=import_module(f"{model_mod}.source").mvs
make_func_dict = cp_mod.make_func_dict
Drivers=cp_mod.Drivers
FreeConstants = namedtuple(
    "FreeConstants",
    [
        'clay',        #Constants set in the cpa.json files
        'silt',
    ]
)
DerivedConstants = namedtuple(
    "DerivedConstants",
    [
        #'ra_0',
        'c_veg_0',
        'c_soil_0',
    ]
)

Constants = namedtuple(
    "Constants",
    [*FreeConstants._fields,*DerivedConstants._fields] 
)

FreeEstimatedParameters = namedtuple(
    'FreeEstimatedParameters', 
    [
        'beta_leaf',
        'beta_root',
        #'r_c_leaf_rh',
        #'r_c_root_rh',
        #'r_c_wood_rh',
        'r_c_leaf_2_c_lit_met',
        'r_c_leaf_2_c_lit_str',
        'r_c_root_2_c_soil_met',
        'r_c_root_2_c_soil_str',
        'r_c_wood_2_c_lit_cwd',
        'r_c_lit_cwd_rh',
        'r_c_lit_met_rh',
        'r_c_lit_str_rh',
        'r_c_lit_mic_rh',
        'r_c_soil_met_rh',
        'r_c_soil_str_rh',
        'r_c_soil_mic_rh',
        'r_c_soil_slow_rh',
        'r_c_soil_passive_rh',
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
        'r_c_soil_passive_2_c_soil_mic',
    ]
)
DerivedEstimatedParameters = namedtuple(
    'DerivedEstimatedParameters', 
    [
        'c_leaf_0',               
        'c_root_0',
        'c_lit_cwd_0',
        'c_lit_met_0',
        'c_lit_str_0',
        'c_lit_mic_0',
        'c_soil_met_0',
        'c_soil_str_0',
        'c_soil_mic_0',
        'c_soil_slow_0',
    ]
)
EstimatedParameters = namedtuple(
    "EstimatedParameters",
    [*FreeEstimatedParameters._fields, *DerivedEstimatedParameters._fields] 
)
def cpa(
        fcpa: FreeConstants,
        dvs: msh.Drivers,
        svs: msh.Observables
    ) -> Constants:
    dcpa = DerivedConstants(
        c_veg_0=svs.cVeg[0],
        c_soil_0=svs.cSoil[0],
    )
    return Constants(*fcpa,*dcpa)

def epa_0(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # this function can obviously be changed
    # it avoids contradictions between the startvalues and 
    # data 
    cv0=svs.cVeg[0]/2
    cs0=svs.cSoil[0]/9
    depa = DerivedEstimatedParameters(
        c_leaf_0=cv0,               
        c_root_0=cv0,
        c_lit_cwd_0=cs0,
        c_lit_met_0=cs0,
        c_lit_str_0=cs0,
        c_lit_mic_0=cs0,
        c_soil_met_0=cs0,
        c_soil_str_0=cs0,
        c_soil_mic_0=cs0,
        c_soil_slow_0=cs0,
    )
    return EstimatedParameters(*fepa,*depa)

def epa_min(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # this function can obviously be changed
    # it avoids contradictions between the startvalues and 
    # data 
    depa = DerivedEstimatedParameters(
        c_leaf_0=0,               
        c_root_0=0,
        c_lit_cwd_0=0,
        c_lit_met_0=0,
        c_lit_str_0=0,
        c_lit_mic_0=0,
        c_soil_met_0=0,
        c_soil_str_0=0,
        c_soil_mic_0=0,
        c_soil_slow_0=0
    )
    return EstimatedParameters(*fepa,*depa)

def epa_max(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # this function can obviously be changed
    # it avoids contradictions between the startvalues and 
    # data 
    cv0=svs.cVeg[0]
    cs0=svs.cSoil[0]
    depa = DerivedEstimatedParameters(
        c_leaf_0=cv0,               
        c_root_0=cv0,
        c_lit_cwd_0=cs0,
        c_lit_met_0=cs0,
        c_lit_str_0=cs0,
        c_lit_mic_0=cs0,
        c_soil_met_0=cs0,
        c_soil_str_0=cs0,
        c_soil_mic_0=cs0,
        c_soil_slow_0=cs0,
    )
    return EstimatedParameters(*fepa,*depa)


def make_param2res_sym(
        mvs,
        fcpa: Constants,
        dvs: Drivers,
        svs: msh.Observables
    ) -> Callable[[np.ndarray], np.ndarray]:
    
    cpa_v = cpa(fcpa, dvs, svs)
    # Define actual forward simulation function
    def param2res(pa):
        
        # Parameter vector
        epa = EstimatedParameters(*pa)
        
        apa = {**cpa_v._asdict(), **epa._asdict()}
        X_0 = numeric_X_0(mvs, dvs, apa)
        #
        par_dict = gh.make_param_dict2(mvs,apa)
        func_dict = make_func_dict(dvs)
        return msh.synthetic_observables(
            mvs,
            X_0,
            par_dict=par_dict,
            func_dict=func_dict,
            dvs=dvs
        )

    return param2res


def make_proposer(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters, 
        fcpa: FreeConstants,
        dvs: msh.Drivers,
        svs: msh.Observables,
    ) -> Callable[[np.ndarray,float,bool], np.ndarray]:
    cpa_v = cpa(fcpa,dvs,svs)
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
        (["beta_leaf", "beta_root"], 1), 
        (["c_leaf_0", "c_root_0"], cpa_v.c_veg_0),
        (
            [
                'c_lit_cwd_0',
                'c_lit_met_0',
                'c_lit_str_0',
                'c_lit_mic_0',
                'c_soil_met_0',
                'c_soil_str_0',
                'c_soil_mic_0',
                'c_soil_slow_0',
            ]
            ,cpa_v.c_soil_0
        )
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
        dvs: Drivers,
        svs: msh.Observables
    ) -> Callable[[np.ndarray], bool]:

    def isQualified(
            c,
            print_conds=False
        ):
        epa=EstimatedParameters(*c)
        apa = {**fcpa._asdict(), **epa._asdict()}
        def value(field_name):
            try:
                return apa[field_name]
            except Exception as e:
                print("###########################")
                print(e)
                print(field_name)
                raise e
        conds=[
            (c >= c_min).all(), 
            (c <= c_max).all(), 
            sum(map(value, ["beta_leaf", "beta_root"])) <= 0.99,
            sum(map(value, ["c_leaf_0", "c_root_0"])) <= svs.cVeg[0],
            sum(
                map(
                    value,
                    [   
                        'c_lit_cwd_0',
                        'c_lit_met_0',
                        'c_lit_str_0',
                        'c_lit_mic_0',
                        'c_soil_met_0',
                        'c_soil_str_0',
                        'c_soil_mic_0',
                        'c_soil_slow_0',
                    ]
                )
            ) <= svs.cSoil[0] # so that c_soil_passive >0
        ]
        res = all(conds)
        if print_conds:
            if not res:
                print(conds)
            return res
        
    return isQualified

def numeric_X_0(mvs,dvs,apa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    par_dict=gh.make_param_dict2(mvs, apa)
    X_0_dict={
        "c_leaf": apa['c_leaf_0'],     
        "c_root": apa['c_root_0'],     
        "c_wood": apa['c_veg_0'] - (apa['c_leaf_0'] +  apa['c_root_0']),  
        "c_lit_cwd": apa['c_lit_cwd_0'],
        "c_lit_met": apa['c_lit_met_0'],
        "c_lit_str": apa['c_lit_str_0'],
        "c_lit_mic": apa['c_lit_mic_0'],
        "c_soil_met": apa['c_soil_met_0'],
        "c_soil_str": apa['c_soil_str_0'],
        "c_soil_mic": apa['c_soil_mic_0'],
        "c_soil_slow": apa['c_soil_slow_0'],
        "c_soil_passive": apa['c_soil_0'] - (
                              apa['c_lit_cwd_0'] 
                            + apa['c_lit_met_0'] 
                            + apa['c_lit_str_0'] 
                            + apa['c_lit_mic_0'] 
                            + apa['c_soil_met_0'] 
                            + apa['c_soil_str_0'] 
                            + apa['c_soil_mic_0'] 
                            + apa['c_soil_slow_0']
                        )
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    return X_0
