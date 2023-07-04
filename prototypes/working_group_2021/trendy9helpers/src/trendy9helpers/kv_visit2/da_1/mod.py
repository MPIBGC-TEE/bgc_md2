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
Constants = namedtuple(
    "Constants",
    [
        "cLitter_0",
        "cSoil_0",
        "cVeg_0",
        "npp_0",
        "rh_0",
        # "ra_0",
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

def make_param_filter_func(
    c_max: EstimatedParameters,
    c_min: EstimatedParameters,
    cpa: Constants
) -> Callable[[np.ndarray], bool]:

    def isQualified(
            c,
            print_conds=False
        ):
        epa=EstimatedParameters(*c)
        apa = {**cpa._asdict(), **epa._asdict()}
        
        def value(field_name):
            return apa[field_name]
        
        conds=[
            (c >= c_min).all(), 
            (c <= c_max).all(), 
            sum(map(value, ["beta_leaf", "beta_wood"])) <= 0.99,
            sum(map(value, ["C_leaf_0", "C_wood_0"])) <= value("cVeg_0"),
            sum(map(value, ["C_leaf_litter_0", "C_wood_litter_0"])) <= value("cLitter_0"),
            sum(
                map(
                    value, ['C_soil_fast_0', 'C_soil_slow_0']
                )
            ) <= value('cSoil_0'),
        ]
        res = all(conds)
        if print_conds:
            if not res:
                print(conds)
        return res
        
    return isQualified

def make_param2res_sym(
    mvs: CMTVS,
    cpa: Constants,
    dvs: Drivers,
) -> Callable[[np.ndarray], np.ndarray]:

    def param2res(pa):

        epa = EstimatedParameters(*pa)
        X_0 = numeric_X_0(mvs, dvs, cpa, epa)
        dpm = h.date.days_per_month
        steps_per_month = 2
        delta_t_val = dpm/steps_per_month 

        par_dict = gh.make_param_dict(mvs, cpa, epa)
        func_dict = make_func_dict(dvs , cpa=cpa, epa=epa)
        bitr = ArrayDictResult(
            msh.make_da_iterator(
                mvs,
                X_0,
                par_dict=par_dict,
                func_dict=func_dict,
                delta_t_val=delta_t_val
            )
        )
        #from IPython import embed; embed()
        number_of_months=dvs.npp.shape[0]

        number_of_steps = number_of_months* steps_per_month
        result_dict = bitr[0: number_of_steps: steps_per_month]

        return msh.Observables(
            cVeg=result_dict["cVeg"],
            cLitter=result_dict["cLitter"],
            cSoil=result_dict["cSoil"],
            rh=result_dict["rh"]
        )
    return param2res

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

