# this file provides some variables for testing.
from collections import namedtuple
from sympy import Symbol, Function
from pathlib import Path
import json
import numpy as np
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
            "cpa"
        ]
    )
    par_dict={
        Symbol(k): v for k,v in {
            'beta_leaf': 0.3333333333333333,
            'beta_root': 0.3333333333333333,
            'r_c_DPM_rh': 1e-4,
            'r_c_RPM_rh': 1e-4,
            'r_c_BIO_rh': 9.04109589041096e-5,
            'r_c_HUM_rh': 1.36849315068493e-5,
            'r_c_leaf_2_c_DPM': 0.000136986301369863,
            'r_c_leaf_2_c_RPM': 0.00123287671232877,
            'r_c_wood_2_c_DPM': 4.56621004566210e-7,
            'r_c_wood_2_c_RPM': 4.52054794520548e-5,
            'r_c_root_2_c_DPM': 9.13242009132420e-8,
            'r_c_root_2_c_RPM': 9.12328767123288e-5,
            'r_c_DPM_2_c_BIO': 4.56621004566210e-6,
            'r_c_DPM_2_c_HUM': 4.10958904109589e-5,
            'r_c_RPM_2_c_BIO': 9.13242009132420e-8,
            'r_c_RPM_2_c_HUM': 9.12328767123288e-5,
            'r_c_BIO_2_c_HUM': 9.13242009132420e-7,
            'r_c_HUM_2_c_BIO': 1.36986301369863e-8
        }.items()
    }
    svs,dvs=msh.get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))
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
             'c_leaf_0': svs_0.cVeg / 3,  # set inital pool values to svs values
             'c_wood_0': svs_0.cVeg / 3,  # you can set numerical values here directly as well
             'c_DPM_0': svs_0.cSoil / 4,
             'c_RPM_0': svs_0.cSoil / 4,
             'c_BIO_0': svs_0.cSoil / 4
         },
         **{
             'beta_leaf': 0.3333333333333333,
             'beta_wood': 0.3333333333333333,
             'Mw': 0.1,
             'Ms': np.max(dvs.mrso)+500, #, may need add a condition here ## ASK MARKUS
             # 'r_c_leaf_rh': 0,
             # 'r_c_wood_rh': 0,
             # 'r_c_root_rh': 0,
             'r_c_DPM_rh': 1e-4,
             'r_c_RPM_rh': 1e-4,
             'r_c_BIO_rh': 9.04109589041096e-5,
             'r_c_HUM_rh': 1.36849315068493e-5,
             'r_c_leaf_2_c_DPM': 0.000136986301369863,
             'r_c_leaf_2_c_RPM': 0.00123287671232877,
             'r_c_wood_2_c_DPM': 4.56621004566210e-7,
             'r_c_wood_2_c_RPM': 4.52054794520548e-5,
             'r_c_root_2_c_DPM': 9.13242009132420e-8,
             'r_c_root_2_c_RPM': 9.12328767123288e-5,
             'r_c_DPM_2_c_BIO': 4.56621004566210e-6,
             'r_c_DPM_2_c_HUM': 4.10958904109589e-5,
             'r_c_RPM_2_c_BIO': 9.13242009132420e-8,
             'r_c_RPM_2_c_HUM': 9.12328767123288e-5,
             'r_c_BIO_2_c_HUM': 9.13242009132420e-7,
             'r_c_HUM_2_c_BIO': 1.36986301369863e-8
         }
    )
    # set min/max parameters to +- 100 times initial values
    epa_min = msh.EstimatedParameters(* tuple(np.array(epa_0) * 0.01))
    epa_max = msh.EstimatedParameters(* tuple(np.array(epa_0) * 100))
    
    # fix values that are problematic from calculation
    epa_max = epa_max._replace(beta_leaf = 0.9)
    epa_max = epa_max._replace(beta_wood = 0.9)
    epa_max = epa_max._replace(c_leaf_0 = svs_0.cVeg * 0.9)
    epa_max = epa_max._replace(c_wood_0 = svs_0.cVeg * 0.9)
    epa_max = epa_max._replace(c_DPM_0=svs_0.cSoil)
    epa_max = epa_max._replace(c_RPM_0=svs_0.cSoil)
    epa_max = epa_max._replace(c_BIO_0=svs_0.cSoil)
    epa_max = epa_max._replace(Mw = 0.8)
    epa_max = epa_max._replace(Ms = (max(dvs.mrso)+500) * 2) # max(dvs.mrso) * 2

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
        c_leaf=svs_0.cVeg / 3,  # set inital pool values to svs values
        c_root=svs_0.cVeg / 3,  # you can set numerical values here directly as well
        c_wood=svs_0.cVeg / 3,
        c_DPM=svs_0.cSoil / 4,
        c_RPM=svs_0.cSoil / 4,
        c_BIO=svs_0.cSoil / 4,
        c_HUM=svs_0.cSoil / 4,
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
        cpa=cpa
    )
