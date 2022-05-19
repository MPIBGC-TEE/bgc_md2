# this file provides some variables for testing.
from collections import namedtuple
from sympy import Symbol, Function
from pathlib import Path
import numpy as np
import json 
import numpy as np
import netCDF4 as nc
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
            "cpa",
            "lats",
            "lons"
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
        'beta_leaf': 0.7682384553229892, 
        'beta_root': 0.14204586101167424, 
        #'r_c_leaf_rh': 0.0022972292016441116,
        #'r_c_root_rh': 0.0015470633697005037,
        #'r_c_wood_rh':0.0003981642399033648,
        'r_c_leaf_2_c_lit_met': 0.0007120305898251202, 
        'r_c_leaf_2_c_lit_str': 0.0007883254343687162, 
        'r_c_root_2_c_soil_met': 0.0012363546536425828,
        'r_c_root_2_c_soil_str': 0.0004067071571929328, 
        'r_c_wood_2_c_lit_cwd': 0.0002446124661809959, 
        'c_leaf_0': 0.5615104265472759, 
        'c_root_0': 0.6495807787373601, 
        'r_c_lit_cwd_rh': 0.012210598784500313, 
        'r_c_lit_met_rh': 0.032329371473925055, 
        'r_c_lit_str_rh': 0.0068962597828840545, 
        'r_c_lit_mic_rh': 0.012288969300236126,
        'r_c_soil_met_rh': 0.0006758167922701092, 
        'r_c_soil_str_rh': 0.00015605646641441882, 
        'r_c_soil_mic_rh': 0.0019289112487673949, 
        'r_c_soil_slow_rh': 4.189139209275766e-05, 
        'r_c_soil_passive_rh': 1.35207047030553e-06,
        'r_c_lit_cwd_2_c_lit_mic': 1.0045960071963091e-05,
        'r_c_lit_cwd_2_c_soil_slow': 2.3397158076915994e-05,
        'r_c_lit_met_2_c_lit_mic': 4.002638739947475e-06,
        'r_c_lit_str_2_c_lit_mic': 0.00015293789003356162,
        'r_c_lit_str_2_c_soil_slow': 0.0006227680438540832,
        'r_c_lit_mic_2_c_soil_slow': 0.0009027280237572068,
        'r_c_soil_met_2_c_soil_mic': 6.415659555529971e-05, 
        'r_c_soil_str_2_c_soil_mic': 4.210824237095305e-05, 
        'r_c_soil_str_2_c_soil_slow': 2.286834247048494e-06,
        'r_c_soil_mic_2_c_soil_slow': 0.0002479841662654259,
        'r_c_soil_mic_2_c_soil_passive': 4.462823543110515e-05,
        'r_c_soil_slow_2_c_soil_mic': 2.1991855605003557e-06,
        'r_c_soil_slow_2_c_soil_passive': 1.2525484949464374e-06,
        'r_c_soil_passive_2_c_soil_mic': 2.7413297461931206e-06,
        'c_lit_cwd_0': 0.17874982599134848,
        'c_lit_met_0': 0.0780909757802213,
        'c_lit_str_0': 0.11870057582418095,
        'c_lit_mic_0': 0.07406262306023444,
        'c_soil_met_0': 0.3098341636788049,
        'c_soil_str_0': 0.5370201126137074,
        'c_soil_mic_0': 0.16367672062228367, 
        'c_soil_slow_0': 1.6344026460653003
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
    # actually we want a dataset with mask but we don't have one at 
    # the moment
    ds=nc.Dataset(Path(conf_dict["dataPath"]).joinpath("YIBs_S0_Monthly_npp.nc"))    
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
        cpa=cpa,
        lats=ds.variables["latitude"][:],
        lons=ds.variables["longitude"][:],
    )
