
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

model_mod = 'bgc_md2.models.yz_jules'
cp_mod=import_module(f"{model_mod}.CachedParameterization")
mvs=import_module(f"{model_mod}.source").mvs
make_func_dict = cp_mod.make_func_dict
Drivers=cp_mod.Drivers

Constants = namedtuple(
    "Constants",
    [
        'npp_0',  # Initial input/pools
        'rh_0',
        'c_veg_0',
        'c_soil_0',
        'fVegSoil_0',  # Total carbon mass from vegetation directly into the soil
        'number_of_months'  # Run time 
    ]
)

EstimatedParameters = namedtuple(
    'EstimatedParameters',
    [
        'c_leaf_0',  # Names: c_poolname_0
        'c_wood_0',  # Only initial pools that are estimated
        'c_DPM_0',
        'c_RPM_0',
        'c_BIO_0',

        'beta_leaf',
        'beta_wood',
        'Mw',
        'Ms',
        'Topt',
        'Tcons',

        'r_c_DPM_rh',
        'r_c_RPM_rh',
        'r_c_BIO_rh',
        'r_c_HUM_rh',

        'r_c_leaf_2_c_DPM',
        'r_c_leaf_2_c_RPM',
        'r_c_wood_2_c_DPM',
        'r_c_wood_2_c_RPM',
        'r_c_root_2_c_DPM',
        'r_c_root_2_c_RPM',
        'r_c_DPM_2_c_BIO',
        'r_c_DPM_2_c_HUM',
        'r_c_RPM_2_c_BIO',
        'r_c_RPM_2_c_HUM',
        'r_c_BIO_2_c_HUM',
        'r_c_HUM_2_c_BIO',
    ]
)

def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers,
) -> Callable[[np.ndarray], np.ndarray]:
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        X_0 = numeric_X_0(mvs, dvs, cpa, epa)
        dpm=30
        steps_per_month = 2
        delta_t_val = dpm/steps_per_month 

        par_dict = gh.make_param_dict(mvs, cpa, epa)
        func_dict = make_func_dict(
            dvs,
            epa.Mw,
            epa.Ms,
            epa.Topt,
            epa.Tcons,
        )
        bitr = ArrayDictResult(
            msh.make_da_iterator(
                mvs,
                X_0,
                par_dict=par_dict,
                func_dict=func_dict,
                delta_t_val=delta_t_val
            )
        )
        number_of_steps= int(cpa.number_of_months * 30 / delta_t_val)
        steps_per_month = int(dpm / delta_t_val)
        result_dict = bitr[0: number_of_steps :steps_per_month]
        #steps_per_year = steps_per_month*12
        #yearly_partitions = gh.partitions(0, number_of_steps, steps_per_year)
        #yearly_averages = {
        #    key: gh.averaged_1d_array(result_dict[key],yearly_partitions)
        #    for key in ["cVeg", "cRoot", "cLitter", "cSoil"]
        #}

        return msh.Observables(
            cVeg=result_dict["cVeg"],
            cSoil=result_dict["cSoil"],
            fVegSoil=result_dict["fVegSoil"],
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
        "c_leaf": apa['c_leaf_0'],     
        "c_wood": apa['c_wood_0'],     
        "c_root": apa['c_veg_0'] - (apa['c_leaf_0'] +  apa['c_wood_0']),  
        "c_DPM": apa['c_DPM_0'],
        "c_RPM": apa['c_RPM_0'],
        "c_BIO": apa['c_BIO_0'],
        "c_HUM": apa['c_soil_0'] - (apa['c_DPM_0'] + apa['c_RPM_0'] + apa['c_BIO_0'])
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    return X_0

# fixme mm 5-18
# there are some conditions missing to 
# prevent negative startvalues
def make_param_filter_func(
    c_max: EstimatedParameters,
    c_min: EstimatedParameters,
    cpa: Constants,
) -> Callable[[np.ndarray], bool]:

    # find position of beta_leaf and beta_wood
    beta_leaf_ind = EstimatedParameters._fields.index("beta_leaf")
    beta_wood_ind = EstimatedParameters._fields.index("beta_wood")

    def isQualified(c):
        beta_leaf_ind
        cond1 = (c >= c_min).all()
        cond2 = (c <= c_max).all()
        cond3 = c[beta_leaf_ind] + c[beta_wood_ind] < 1
        print(
            "cond1",cond1,
            "cond2",cond2,
            "cond3",cond3,
        )
        return cond1 and cond2 and cond3

    return isQualified
