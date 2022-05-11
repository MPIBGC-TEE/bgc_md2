# this file provides some variables for testing.
from collections import namedtuple
from sympy import Symbol, Function
from pathlib import Path
import numpy as np
import json 

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
    
    cpa=msh.Constants(
     cVeg_0=svs_0.cVeg,
     cLitter_0=svs_0.cLitter,
     cSoil_0=svs_0.cSoil,
     npp_0=dvs.npp[0],   # kg/m2/s kg/m2/day
     rh_0=svs_0.rh,   # kg/m2/s kg/m2/day
     #mrso_0=dvs.mrso[0],
     #tsl_0=dvs.tsl[0],
     number_of_years=int(320) # for testing and tuning mcmc
     #number_of_months=len(svs.rh)
    )
    
    par_dict= {
        'beta_leaf': 0.3028851224272647,
        'beta_wood': 0.4564030561461778,
        'Theta_sat': 0.09469315267122164,
        'Theta_fc': 0.0795237425542056,
        'r_C_leaf_rh': 0,
        'r_C_wood_rh': 0,
        'r_C_root_rh': 0,
        'r_C_aom1_rh': 0.0003037126928450523,
        'r_C_aom2_rh': 0.00029780017612134036,
        'r_C_smb1_rh': 1.0689106224513612e-05,
        'r_C_smb2_rh': 0.00017382741384128404,
        'r_C_smr_rh': 7.29448370877164e-06,
        'r_C_nom_rh': 1.0164098602070474e-05,
        'r_C_dom_rh': 1.1500691016736846e-05,
        'r_C_psom_rh': 1.7644513592821266e-05,
        'r_C_leaf_2_C_aom1': 0.01970013741245839,
        'r_C_leaf_2_C_aom2': 0.015317429471496852,
        'r_C_wood_2_C_aom1': 0.0004097261277630772,
        'r_C_wood_2_C_aom2': 0.0003689152669226263,
        'r_C_root_2_C_aom1': 0.0006603079638088481,
        'r_C_root_2_C_aom2': 0.00043562236684542056,
        'r_C_aom1_2_C_smb1': 0.00023145662814734428,
        'r_C_aom1_2_C_smb2': 0.0002977928058255672,
        'r_C_aom1_2_C_nom': 0.00015287835816890072,
        'r_C_aom1_2_C_dom': 0.0001781924058628594,
        'r_C_aom2_2_C_smb1': 2.0810976139531716e-05,
        'r_C_aom2_2_C_smb2': 2.8349756784116157e-05,
        'r_C_aom2_2_C_dom': 2.14750752962314e-05,
        'r_C_smb1_2_C_nom': 2.3404172674296687e-05,
        'r_C_smb1_2_C_psom': 3.3408243666981594e-05,
        'r_C_smb2_2_C_smr': 1.2338527499527314e-05,
        'r_C_smr_2_C_smb1': 4.801573811093888e-05,
        'r_C_nom_2_C_smb1': 1.0587791011076786e-05,
        'r_C_nom_2_C_dom': 2.551756974855787e-05,
        'r_C_nom_2_C_psom': 1.307830769309431e-05,
        'r_C_dom_2_C_smb1': 1.0950958379980324e-05,
        'r_C_dom_2_C_nom': 1.4109489107131771e-05,
        'r_C_psom_2_C_smb1': 2.222172994302753e-06,
        'C_leaf_0': 0.0034735199394882624,
        'C_wood_0': 0.23539047059980747,
        'C_aom1_0': 0.2312532697753018,
        'C_smb1_0': 2.267901320473235,
        'C_smb2_0': 0.5198209480949317,
        'C_smr_0': 0.11596093546336099,
        'C_nom_0': 2.3560315471577784,
        'C_dom_0': 3.309888037557665
    }
    
    epa_0 = msh.EstimatedParameters(**par_dict)
    epa_opt = epa_0
    
    # set min/max parameters to +- 100 times initial values
    epa_min=msh.EstimatedParameters._make(tuple(np.array(epa_0)*0.01))
    epa_max=msh.EstimatedParameters._make(tuple(np.array(epa_0)*100))

    # fix values that are problematic from calculation
    epa_max = epa_max._replace(beta_leaf = 0.9)
    epa_max = epa_max._replace(beta_wood = 0.9)
    epa_max = epa_max._replace(C_leaf_0 = svs_0.cVeg)
    epa_max = epa_max._replace(C_wood_0 = svs_0.cVeg)
    epa_max = epa_max._replace(C_aom1_0 = svs_0.cLitter)
    epa_max = epa_max._replace(C_smb1_0 = svs_0.cSoil)
    epa_max = epa_max._replace(C_smb2_0 = svs_0.cSoil)
    epa_max = epa_max._replace(C_smr_0 = svs_0.cSoil)
    epa_max = epa_max._replace(C_nom_0 = svs_0.cSoil)
    epa_max = epa_max._replace(C_dom_0 = svs_0.cSoil)
    
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
        C_leaf = apa['C_leaf_0'],
        C_wood = apa['C_wood_0'],
        C_root = apa['cVeg_0'] - (
            apa['C_leaf_0'] + 
            apa['C_wood_0']
        ),
        C_aom1 = apa['C_aom1_0'],
        C_aom2 = apa['cLitter_0'] - apa['C_aom1_0'],
        C_smb1 = apa['C_smb1_0'],
        C_smb2 = apa['C_smb2_0'],
        C_smr = apa['C_smr_0'],
        C_nom = apa['C_nom_0'],
        C_dom = apa['C_dom_0'],
        C_psom = apa['cSoil_0'] - (
            apa['C_smb1_0'] +
            apa['C_smb2_0'] +
            apa['C_smr_0'] +
            apa['C_nom_0'] +
            apa['C_dom_0'] 
        ),
        rh = apa['rh_0']
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
