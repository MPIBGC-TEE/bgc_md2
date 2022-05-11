# this file provides some variables for testing.
from collections import namedtuple
from sympy import Symbol, Function
from pathlib import Path
import numpy as np
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
    svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
    
    svs_0=msh.Observables(*map(lambda v: v[0],svs))
    #Define the constant values of parameters 
    # that are NOT affected by data assimilation
    cpa = msh.Constants(             #use Constants namedtuple to define constant values
        npp_0 = dvs.npp[0],
        rh_0 = svs.rh[0],
        #ra_0 = svs.ra[0],
        c_veg_0 = svs.cVeg[0],
        c_soil_0 = svs.cSoil[0],
        clay = 0.2028,
        silt = 0.2808,
        nyears = 320,
    )

    par_dict = {
        'beta_leaf': 0.37152535661667285,
        'beta_root': 0.2118738332472721,
        #'r_c_leaf_rh': 0.0022972292016441116,
        #'r_c_root_rh': 0.0015470633697005037,
        #'r_c_wood_rh':0.0003981642399033648,
        'r_c_leaf_2_c_lit_met': 0.0008419144443122888, 
        'r_c_leaf_2_c_lit_str': 7.253712507163508e-05,
        'r_c_root_2_c_soil_met': 0.0007599224861792184,
        'r_c_root_2_c_soil_str': 0.0007161706404910827,
        'r_c_wood_2_c_lit_cwd': 0.0009217945194693122,
        'c_leaf_0': 0.11328379866881665,
        'c_root_0': 0.14464613373390392,
        'r_c_lit_cwd_rh': 0.02026318476587012, 
        'r_c_lit_met_rh': 0.00340079410753037, 
        'r_c_lit_str_rh': 0.008989119944533677,
        'r_c_lit_mic_rh': 0.011276949417831122,
        'r_c_soil_met_rh': 0.0006741622348146495,
        'r_c_soil_str_rh': 0.00017592886085999286,
        'r_c_soil_mic_rh': 0.000519741477608671,
        'r_c_soil_slow_rh': 1.0255263440555624e-06,
        'r_c_soil_passive_rh': 3.881935738016802e-07,
        'r_c_lit_cwd_2_c_lit_mic': 1.3188464625334016e-05,
        'r_c_lit_cwd_2_c_soil_slow': 1.6316549662914743e-05,
        'r_c_lit_met_2_c_lit_mic': 2.9433144645429653e-06,
        'r_c_lit_str_2_c_lit_mic': 0.00010298015064924245,
        'r_c_lit_str_2_c_soil_slow': 0.0016579805745133146,
        'r_c_lit_mic_2_c_soil_slow': 0.0011840494205249575,
        'r_c_soil_met_2_c_soil_mic': 7.861811338124696e-05,
        'r_c_soil_str_2_c_soil_mic': 2.578967926776423e-05,
        'r_c_soil_str_2_c_soil_slow': 1.7394627034766953e-06,
        'r_c_soil_mic_2_c_soil_slow': 0.00021605360652605818,
        'r_c_soil_mic_2_c_soil_passive': 4.569266267503945e-05,
        'r_c_soil_slow_2_c_soil_mic': 4.1146075824754925e-07,
        'r_c_soil_slow_2_c_soil_passive': 2.9993396188473066e-08,
        'r_c_soil_passive_2_c_soil_mic': 2.751360714464457e-06,
        'c_lit_cwd_0': 0.011122590276073926, 
        'c_lit_met_0': 0.04563448012195457,
        'c_lit_str_0': 0.022083588329899793,
        'c_lit_mic_0': 0.011910319433275054,
        'c_soil_met_0': .048208986458370635,
        'c_soil_str_0': 0.6643525311241724,
        'c_soil_mic_0': 0.05837121211447685,
        'c_soil_slow_0': 0.3228602860446373
    }

    
    epa_0 = msh.EstimatedParameters(**par_dict)
    epa_opt = epa_0
    
    # set min/max parameters to +- 100 times initial values
    epa_min=msh.EstimatedParameters._make(tuple(np.array(epa_0)*0.01))
    epa_max=msh.EstimatedParameters._make(tuple(np.array(epa_0)*100))

    # fix values that are problematic from calculation
    epa_max = epa_max._replace(beta_leaf = 0.9)
    epa_max = epa_max._replace(beta_root = 0.9)
    epa_max = epa_max._replace(c_leaf_0 = svs_0.cVeg)
    epa_max = epa_max._replace(c_root_0 = svs_0.cVeg)
    epa_max = epa_max._replace(c_lit_cwd_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_lit_met_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_lit_str_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_lit_mic_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_soil_met_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_soil_str_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_soil_mic_0 = svs_0.cSoil)
    epa_max = epa_max._replace(c_soil_slow_0 = svs_0.cSoil)
    
    # Create namedtuple for initial values
    StartVector = msh.make_StartVector(mvs)
    
    # Parameter dictionary for the iterator
    apa = {**cpa._asdict(),**epa_0._asdict()}
    
    # Build dictionary of model parameters
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    model_par_dict = {
        Symbol(k):v for k,v in apa.items()
        if Symbol(k) in model_par_dict_keys
    }

    # Create a startvector for the iterator 
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
        rh = apa['rh_0'] #,
        #ra = apa['ra_0']
    )
        
    return TestArgs(
        V_init=V_init,
        par_dict=model_par_dict,
        func_dict=msh.make_func_dict(mvs,dvs,cpa,epa_opt),
        dvs=dvs,
        svs=svs,
        mvs=mvs,
        epa_0=epa_0,
        epa_min=epa_min,
        epa_max=epa_max,
        epa_opt=epa_opt,
        cpa=cpa
    )
