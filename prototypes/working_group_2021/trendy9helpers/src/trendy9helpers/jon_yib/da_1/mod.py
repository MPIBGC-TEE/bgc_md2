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

model_mod = 'bgc_md2.models.jon_yib'
cp_mod=import_module(f"{model_mod}.CachedParameterization")
mvs=import_module(f"{model_mod}.source").mvs
make_func_dict = cp_mod.make_func_dict
Drivers=cp_mod.Drivers

Constants = namedtuple(
    "Constants",
    [
        'npp_0',       #Initial input/pools
        'rh_0',
        #'ra_0',
        'c_veg_0',
        'c_soil_0',
        'clay',        #Constants like clay
        'silt',
    ]
)
EstimatedParameters = namedtuple(
    'EstimatedParameters', 
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
        'c_leaf_0',               
        'c_root_0',
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

def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers,
    ) -> Callable[[np.ndarray], np.ndarray]:
    
    
    # Define actual forward simulation function
    def param2res(pa):
        
        # Parameter vector
        epa = EstimatedParameters(*pa)
        X_0 = numeric_X_0(mvs, dvs, cpa, epa)
        dpm = h.date.days_per_month
        steps_per_month = 2
        delta_t_val = dpm/steps_per_month 
        
        par_dict = gh.make_param_dict(mvs, cpa, epa)
        func_dict = make_func_dict(dvs)
        
        #from IPython import embed; embed()
        number_of_months=dvs.npp.shape[0]
        steps_per_month = 2
        number_of_steps = number_of_months*steps_per_month
        bitr = ArrayDictResult(
            msh.make_da_iterator(
                mvs,
                X_0,
                par_dict=par_dict,
                func_dict=func_dict,
                delta_t_val=delta_t_val
            )
        )
        result_dict = bitr[0: number_of_steps: steps_per_month]
        steps_per_year = steps_per_month*12
        yearly_partitions = gh.partitions(0, number_of_steps, steps_per_year)
        yearly_averages = {
            key: gh.averaged_1d_array(result_dict[key],yearly_partitions)
            for key in ["cVeg", "cSoil"]
        }

        return msh.Observables(
            cVeg=yearly_averages["cVeg"],
            cSoil=yearly_averages["cSoil"],
            rh=result_dict["rh"]#/(60*60*24)
        )

    return param2res

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
            sum(map(value, ["c_leaf_0", "c_root_0"])) <= value("c_veg_0"),
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
            ) <= value('c_soil_0')
        ]
        res = all(conds)
        if print_conds:
            if not res:
                print(conds)
            return res
        
    return isQualified

def numeric_X_0(mvs,dvs,cpa,epa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
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
