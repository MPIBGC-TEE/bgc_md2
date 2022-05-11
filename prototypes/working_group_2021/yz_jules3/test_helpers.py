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
    par_dict={
        Symbol(k): v for k,v in {           
            'beta_leaf': 0.35,
            'beta_wood': 0.30,
            # 'r_c_leaf_rh': 0,
            # 'r_c_wood_rh': 0,
            # 'r_c_root_rh': 0,
            'r_c_DPM_rh':0.0218855218855219,
            'r_c_RPM_rh':0.000866666666666667,
            'r_c_BIO_rh':0.00174841269841270,
            'r_c_HUM_rh':5.87450980392157e-5,
            'r_c_leaf_2_c_DPM':0.000152777777777778,
            'r_c_leaf_2_c_RPM':0.000541666666666667,
            'r_c_wood_2_c_DPM':2.00364298724954e-5,
            'r_c_wood_2_c_RPM':7.10382513661202e-5,
            'r_c_root_2_c_DPM':0.000152777777777778,
            'r_c_root_2_c_RPM':0.000541666666666667,
            'r_c_DPM_2_c_BIO':0.00283950617283951,
            'r_c_DPM_2_c_HUM':0.00333333333333333,
            'r_c_RPM_2_c_BIO':0.000112444444444444,
            'r_c_RPM_2_c_HUM':0.000132000000000000,
            'r_c_BIO_2_c_HUM':0.000235714285714286,
            'r_c_HUM_2_c_BIO':6.61437908496732e-6
        }.items()
    }
    svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
    # func_dict={
    #     Function(k):v
    #     for k,v in {
    #         'NPP':msh.make_npp_func(dvs),
    #         'xi':msh.make_xi_func(dvs)
    #     }.items()
    # }
    svs_0=msh.Observables(*map(lambda v: v[0],svs))
    # # create a start parameter tuple for the mcmc.
    epa_0=msh.EstimatedParameters(
    **{
            'c_leaf_0': svs_0.cVeg * 0.12,  # set inital pool values to svs values
            'c_wood_0': svs_0.cVeg * 0.76,  # you can set numerical values here directly as well
            'c_DPM_0': svs_0.cSoil * 0.0025,  # set according to QY's single site results: 0.0025 DPM, 0.22 RPM, 0.02 BIO, 0.7575 HUM
            'c_RPM_0': svs_0.cSoil * 0.248,
            'c_BIO_0': svs_0.cSoil * 0.022
         },
         **{
            'beta_leaf': 0.35,
            'beta_wood': 0.3,
            'Mw': 0.1,
            'Ms': np.max(dvs.mrsos) + 500, #, may need add a condition here ## ASK MARKUS
            'Topt': 18.32,
            'Tcons': 47.91,
            # 'r_c_leaf_rh': 0,
            # 'r_c_wood_rh': 0,
            # 'r_c_root_rh': 0,
            'r_c_DPM_rh':0.0218855218855219,
            'r_c_RPM_rh':0.000866666666666667,
            'r_c_BIO_rh':0.00174841269841270,
            'r_c_HUM_rh':5.87450980392157e-5,
            'r_c_leaf_2_c_DPM':0.000152777777777778,
            'r_c_leaf_2_c_RPM':0.000541666666666667,
            'r_c_wood_2_c_DPM':2.00364298724954e-5,
            'r_c_wood_2_c_RPM':7.10382513661202e-5,
            'r_c_root_2_c_DPM':0.000152777777777778,
            'r_c_root_2_c_RPM':0.000541666666666667,
            'r_c_DPM_2_c_BIO':0.00283950617283951,
            'r_c_DPM_2_c_HUM':0.00333333333333333,
            'r_c_RPM_2_c_BIO':0.000112444444444444,
            'r_c_RPM_2_c_HUM':0.000132000000000000,
            'r_c_BIO_2_c_HUM':0.000235714285714286,
            'r_c_HUM_2_c_BIO':6.61437908496732e-6
         }
    )
    # Yu please change this to the optimized parameters if you changed them  
    epa_opt=epa_0
    # set min/max parameters to +- 100 times initial values
    epa_min = msh.EstimatedParameters(* tuple(np.array(epa_0) * 0.001)) # for fluxes 
    epa_max = msh.EstimatedParameters(* tuple(np.array(epa_0) * 200))

    # fix values that are problematic from calculation
    epa_max = epa_max._replace(beta_leaf = epa_0.beta_leaf * 1.5)
    epa_max = epa_max._replace(beta_wood = epa_0.beta_wood * 1.5)
    epa_max = epa_max._replace(c_leaf_0 = svs_0.cVeg * epa_0.c_leaf_0 * 1.5)
    epa_max = epa_max._replace(c_wood_0 = svs_0.cVeg * epa_0.c_wood_0 * 1.5)
    epa_max = epa_max._replace(c_DPM_0 = svs_0.cSoil * epa_0.c_DPM_0 * 10)
    epa_max = epa_max._replace(c_RPM_0 = svs_0.cSoil * epa_0.c_RPM_0 * 10)
    epa_max = epa_max._replace(c_BIO_0 = svs_0.cSoil * epa_0.c_BIO_0 * 10)
    epa_max = epa_max._replace(Mw = 0.5)
    epa_max = epa_max._replace(Ms = (max(dvs.mrsos) + 500) * 2) # max(dvs.mrso) * 2
    epa_min = epa_min._replace(Topt = epa_0.Topt - 5)
    epa_max = epa_max._replace(Topt = epa_0.Topt + 5)
    epa_min = epa_min._replace(Tcons = epa_0.Tcons / 1.2)
    epa_max = epa_max._replace(Tcons = epa_0.Tcons * 1.2)
    
    cpa = msh.Constants(
        # use Constants namedtuple to define constant values #Define the constant values of parameters NOT affected by data assimilation
        npp_0=dvs.npp[0],
        rh_0=svs_0.rh,
        c_veg_0=svs_0.cVeg,
        c_soil_0=svs_0.cSoil,
        fVegSoil_0=svs_0.fVegSoil,  # add the fraction
        nyears=320
    )
    # Assign values to initial pools using InitialPools named tuple
    StartVector = msh.make_StartVector(mvs)
    V_init = StartVector(
        c_leaf=svs_0.cVeg * 0.12,  # set inital pool values to svs values
        c_root=svs_0.cVeg * 0.12,  # you can set numerical values here directly as well
        c_wood=svs_0.cVeg * 0.76,
        c_DPM=svs_0.cSoil * 0.0025,
        c_RPM=svs_0.cSoil * 0.248,
        c_BIO=svs_0.cSoil * 0.022,
        c_HUM=svs_0.cSoil * 0.7275,
        rh=svs_0.rh,
        fVegSoil=svs_0.fVegSoil
        # f_veg2soil=svs_0.f_veg2soil# add the fraction
    )
    
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
