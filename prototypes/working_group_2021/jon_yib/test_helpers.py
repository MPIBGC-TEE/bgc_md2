# this file provides some variables for testing.
from collections import namedtuple
from sympy import Symbol, Function
from pathlib import Path
import json 
import numpy as np
from functools import lru_cache

@lru_cache
def make_test_args(conf_dict,msh,mvs):
    TestArgs=namedtuple(
        "TestArgs",
        [
            "V_init",
            "par_dict",
            "func_dict",
            "mvs",
            "dvs",
            "svs",
            "epa_0",
            "epa_min",
            "epa_max",
            "epa_opt",
            "cpa"
        ]
    )
    svs,dvs=msh.get_example_site_vars(
        dataPath=Path(conf_dict["dataPath"])
    )
    svs_0=msh.Observables(*map(lambda v: v[0],svs))
    #Define the constant values of parameters 
    # that are NOT affected by data assimilation
    cpa = msh.Constants(             
        npp_0 = dvs.npp[0],
        rh_0 = svs.rh[0],
        ra_0 = svs.ra[0],
        c_veg_0 = svs_0.cVeg,
        c_soil_0 = svs_0.cSoil,
        clay = 0.2028,
        silt = 0.2808,
        #nyears = 320
        nyears = 5
    )
    par_dict = {
        'beta_leaf': 0.3,
        'beta_root': 0.3,
        'r_c_leaf_rh': 0,
        'r_c_root_rh': 0,
        'r_c_wood_rh': 0,
        'r_c_lit_cwd_rh': 0.00150000000000000,
        'r_c_lit_met_rh': 0.0400000000000000,
        'r_c_lit_str_rh': 0.0300000000000000,
        'r_c_lit_mic_rh': 0.0400000000000000,
        'r_c_soil_met_rh': 0.0320000000000000,
        'r_c_soil_str_rh': 0.0234000000000000,
        'r_c_soil_mic_rh': 0.0300000000000000,
        'r_c_soil_slow_rh': 7.99034960000000e-5,
        'r_c_soil_passive_rh': 0.000261600000000000,
        'r_c_leaf_2_c_lit_met': 0.0100000000000000,
        'r_c_leaf_2_c_lit_str': 0.0100000000000000,
        'r_c_root_2_c_soil_met': 0.00500000000000000,
        'r_c_root_2_c_soil_str': 0.00500000000000000,
        'r_c_wood_2_c_lit_cwd': 0.00700000000000000,
        'r_c_lit_cwd_2_c_lit_mic': 0.00350000000000000,
        'r_c_lit_cwd_2_c_soil_slow': 0.00500000000000000,
        'r_c_lit_met_2_c_lit_mic': 0.0100000000000000,
        'r_c_lit_str_2_c_lit_mic': 0.0100000000000000,
        'r_c_lit_str_2_c_soil_slow': 0.0100000000000000,
        'r_c_lit_mic_2_c_soil_slow': 0.0100000000000000,
        'r_c_soil_met_2_c_soil_mic': 0.00800000000000000,
        'r_c_soil_str_2_c_soil_mic': 0.00780000000000000,
        'r_c_soil_str_2_c_soil_slow': 0.00780000000000000,
        'r_c_soil_mic_2_c_soil_slow': 0.0100000000000000,
        'r_c_soil_mic_2_c_soil_passive': 0.0100000000000000,
        'r_c_soil_slow_2_c_soil_mic': 2.00000000000000e-5,
        'r_c_soil_slow_2_c_soil_passive': 9.65040000000000e-8,
        'r_c_soil_passive_2_c_soil_mic': 6.54000000000000e-5
    }
    epa_0 = msh.EstimatedParameters(
    **{
        "c_leaf_0": svs_0.cVeg/3,          #set inital pool values to svs values
        "c_root_0": svs_0.cVeg/3,          #you can set numerical values here directly as well
        "c_lit_cwd_0": svs_0.cSoil/35,
        "c_lit_met_0": svs_0.cSoil/35,
        "c_lit_str_0": svs_0.cSoil/35,
        "c_lit_mic_0": svs_0.cSoil/35,
        "c_soil_met_0": svs_0.cSoil/20,
        "c_soil_str_0": svs_0.cSoil/15,
        "c_soil_mic_0": svs_0.cSoil/10,
        "c_soil_slow_0": svs_0.cSoil/3
        
    },
    **par_dict
    )

    # Jon please change this to the optimized parameters if you changed them  
    epa_opt=epa_0
    
    apa = {**cpa._asdict(),**epa_0._asdict()}
    StartVector = msh.make_StartVector(mvs)
    V_init = StartVector(
            c_leaf = apa['c_leaf_0'],
            c_root = apa['c_root_0'],
            c_wood = apa['c_veg_0'] - (
                apa['c_leaf_0'] + 
                apa['c_root_0']
            ),
            c_lit_cwd = apa['c_lit_cwd_0'],
            c_lit_met = apa['c_lit_met_0'],
            c_lit_str = apa['c_lit_str_0'],
            c_lit_mic = apa['c_lit_mic_0'],
            c_soil_met = apa['c_soil_met_0'],
            c_soil_str = apa['c_soil_str_0'],
            c_soil_mic = apa['c_soil_mic_0'],
            c_soil_slow = apa['c_soil_slow_0'],
            c_soil_passive = apa['c_soil_0'] - (
                apa['c_lit_cwd_0'] +
                apa['c_lit_met_0'] +
                apa['c_lit_str_0'] +
                apa['c_lit_mic_0'] +
                apa['c_soil_met_0'] +
                apa['c_soil_str_0'] +
                apa['c_soil_mic_0'] +
                apa['c_soil_slow_0']
            ),
            rh = apa['rh_0'],
            ra = apa['ra_0']
    )

    # set min/max parameters to +- 100 times initial values
    epa_min=msh.EstimatedParameters._make(tuple(np.array(epa_0)*0.01))
    epa_max=msh.EstimatedParameters._make(tuple(np.array(epa_0)*100))

    # fix values that are problematic from calculation
    #epa_max = epa_max._replace(beta_leaf = 0.9)
    #epa_max = epa_max._replace(beta_root = 0.9)
    #epa_max = epa_max._replace(c_leaf_0 = svs_0.cVeg)
    #epa_max = epa_max._replace(c_root_0 = svs_0.cVeg)
    epa_max = epa_max._replace(c_lit_cwd_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_lit_met_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_lit_str_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_lit_mic_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_soil_met_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_soil_str_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_soil_mic_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_soil_slow_0 = svs_0.cSoil)
    
    return TestArgs(
        V_init=V_init,
        par_dict=par_dict,
        func_dict=msh.make_func_dict(mvs,dvs,cpa,epa_0),
        dvs=dvs,
        svs=svs,
        mvs=mvs,
        epa_0=epa_0,
        epa_min=epa_min,
        epa_max=epa_max,
        epa_opt=epa_opt,
        cpa=cpa
    )
