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
            "lons",
            "start_date"
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
    'beta_leaf':0.7930800797551596, 
    'beta_root':0.09594854231998003, 
    'r_c_leaf_2_c_lit_met':0.000787334849363968, 
    'r_c_leaf_2_c_lit_str':0.0005841023155216069, 
    'r_c_root_2_c_soil_met':0.0012742764626961106, 
    'r_c_root_2_c_soil_str':0.0004462365679457925, 
    'r_c_wood_2_c_lit_cwd':0.00018001438095217938, 
    'c_leaf_0':0.44009643074602856, 
    'c_root_0':0.6678484442904881, 
    'r_c_lit_cwd_rh':0.012953231396982123, 
    'r_c_lit_met_rh':0.05393965198634569, 
    'r_c_lit_str_rh':0.007017095898872865, 
    'r_c_lit_mic_rh':0.00983547116108986, 
    'r_c_soil_met_rh':0.0005261300243452956, 
    'r_c_soil_str_rh':0.000167029792112257, 
    'r_c_soil_mic_rh':0.001256426932787047, 
    'r_c_soil_slow_rh':2.9300467335491405e-05, 
    'r_c_soil_passive_rh':9.62437139608528e-07, 
    'r_c_lit_cwd_2_c_lit_mic':1.281457434214037e-05,
    'r_c_lit_cwd_2_c_soil_slow':2.2140304421250735e-05, 
    'r_c_lit_met_2_c_lit_mic':4.765056198078948e-06,
    'r_c_lit_str_2_c_lit_mic':0.0001484429124541755, 
    'r_c_lit_str_2_c_soil_slow':0.0008938597443093927, 
    'r_c_lit_mic_2_c_soil_slow':0.0007298297320239218, 
    'r_c_soil_met_2_c_soil_mic':5.3768071140100535e-05, 
    'r_c_soil_str_2_c_soil_mic':4.069550359283175e-05,
    'r_c_soil_str_2_c_soil_slow':2.2987957568951882e-06, 
    'r_c_soil_mic_2_c_soil_slow':0.00019623354207997005,
    'r_c_soil_mic_2_c_soil_passive':3.5528012678573316e-05,
    'r_c_soil_slow_2_c_soil_mic':1.6856472203388676e-06, 
    'r_c_soil_slow_2_c_soil_passive':1.0091955670820996e-06,
    'r_c_soil_passive_2_c_soil_mic':2.096068650869897e-06, 
    'c_lit_cwd_0':0.1069874406303686, 
    'c_lit_met_0':0.06272469710387509, 
    'c_lit_str_0':0.102891111018907, 
    'c_lit_mic_0':0.10278739781806061, 
    'c_soil_met_0':0.3201141255719188, 
    'c_soil_str_0':0.6667514376729721, 
    'c_soil_mic_0':0.1745864850357272, 
    'c_soil_slow_0':1.1416674862662826
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
    ds=nc.Dataset(Path(conf_dict["dataPath"]).joinpath("YIBs_S2_Monthly_npp.nc"))    
    return TestArgs(
        V_init=V_init,
        par_dict=model_par_dict,
        func_dict=msh.make_func_dict(dvs),
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
        start_date=msh.start_date()
    )
