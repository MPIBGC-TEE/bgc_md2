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

model_mod = 'bgc_md2.models.AneeshSDGVM'
cp_mod=import_module(f"{model_mod}.CachedParameterization")
mvs=import_module(f"{model_mod}.source").mvs
make_func_dict = cp_mod.make_func_dict
Drivers=cp_mod.Drivers

Constants = namedtuple(
    "Constants",
    [
        "cLitter_0",
        "cSoil_0",
        "cRoot_0",
        "cVeg_0",
        "npp_0",
        "rh_0"
    ]
)

EstimatedParameters = namedtuple(
    "EstimatedParameters",[
         'beta_leaf',
         'beta_wood',
         #beta_root,
         'r_C_leaf2abvstrlit',
         'r_C_abvmetlit2surface_microbe',
         'r_C_abvstrlit2slowsom',
         'r_C_abvstrlit2surface_microbe',
         'r_C_belowmetlit2soil_microbe',
         'r_C_belowstrlit2slowsom',
         'r_C_belowstrlit2soil_microbe',
         #'r_C_leached',
         'r_C_leaf2abvmetlit',
         'r_C_passsom2soil_microbe',
         'r_C_root2belowmetlit',
         'r_C_root2belowstrlit',
         'r_C_slowsom2passsom',
         'r_C_slowsom2soil_microbe',
         'r_C_soil_microbe2passsom',
         'r_C_soil_microbe2slowsom',
         'r_C_surface_microbe2slowsom',
         'r_C_wood2abvmetlit',
         'r_C_wood2abvstrlit',
         'r_C_abvstrlit_rh',
         'r_C_abvmetlit_rh',
         'r_C_belowstrlit_rh',
         'r_C_belowmetlit_rh',
         'r_C_surface_microbe_rh',
         'r_C_slowsom_rh',
         'r_C_passsom_rh',
         'r_C_soil_microbe_rh',
         'C_leaf_0',
        #'C_root_0',
         'C_abvstrlit_0',
         'C_abvmetlit_0',
         'C_blwstrlit_0',
         'C_surfacemic_0',
         'C_soilmic_0',
         'C_slow_0'
    ]
)

def make_param_filter_func(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters, 
        cpa: Constants,
        ) -> Callable[[np.ndarray], bool]:

    # find position of beta_leaf and beta_wood
        
    
    def isQualified(c):
        def value(field_name):
            try:
                return c[EstimatedParameters._fields.index(field_name)]
            except Exception as e:
                print("###########################")
                print(e)
                print(field_name)
                raise e
        conds=[
            (c >= c_min).all(), 
            (c <= c_max).all(), 
            sum(map(value, ["beta_leaf", "beta_wood"])) <= 0.99,
            value("C_leaf_0") <= cpa.cVeg_0-cpa.cRoot_0, 
            sum(
                map(
                    value,
                    ["C_abvstrlit_0","C_abvmetlit_0","C_blwstrlit_0"]
                )
            ) <= cpa.cVeg_0,
            sum(
                map(
                    value,
                    ["C_surfacemic_0", "C_soilmic_0", "C_slow_0"]
                )
            ) <= cpa.cSoil_0  
        ]    
        res=all(conds)
        if not res:
            print(conds)
        return res
        
    return isQualified


def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers
) -> Callable[[np.ndarray], np.ndarray]: 
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        X_0 = numeric_X_0(mvs, dvs, cpa, epa)
        dpm=30
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
        number_of_steps = int(number_of_months/delta_t_val)
        result_dict = bitr[0: number_of_steps: steps_per_month]
        steps_per_year = steps_per_month*12
        yearly_partitions = gh.partitions(0, number_of_steps, steps_per_year)
        yearly_averages = {
            key: gh.averaged_1d_array(result_dict[key],yearly_partitions)
            for key in ["cVeg", "cRoot", "cLitter", "cSoil"]
        }

        return msh.Observables(
            cVeg=yearly_averages["cVeg"],
            cRoot=yearly_averages["cRoot"],
            cLitter=yearly_averages["cLitter"],
            cSoil=yearly_averages["cSoil"],
            rh=result_dict["rh"]#/(60*60*24)
        )
    return param2res


def numeric_X_0(mvs,dvs,cpa,epa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "C_leaf": apa['C_leaf_0'],     
        "C_root": apa['cRoot_0'],     
        "C_wood": apa['cVeg_0'] - (apa['C_leaf_0'] +  apa['cRoot_0']),  
        "C_abvstrlit": apa['C_abvstrlit_0'],
        "C_abvmetlit": apa['C_abvmetlit_0'],
        "C_belowstrlit": apa["C_blwstrlit_0"],
        "C_belowmetlit":apa["cLitter_0"]- apa["C_abvstrlit_0"] - apa["C_abvmetlit_0"] - apa["C_blwstrlit_0"],
        "C_surface_microbe":apa["C_surfacemic_0"],
        "C_soil_microbe": apa["C_soilmic_0"],
        "C_slowsom":apa["C_slow_0"],
        "C_passsom":apa["cSoil_0"] - apa["C_surfacemic_0"] - apa["C_soilmic_0"] - apa["C_slow_0"]
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    return X_0


