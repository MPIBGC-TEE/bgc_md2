# this file provides some variables for testing.
from collections import namedtuple
from sympy import Symbol, Function
from pathlib import Path
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
            "cpa"
        ]
    )
    par_dict={
        Symbol(k): v for k,v in {
            'beta_leaf': 0.6,
            'beta_wood': 0.25,
            'T_0': 2,
            'E': 4,
            'KM': 10,
            'r_C_leaf_litter_rh': 0.000415110004151100,
            'r_C_wood_litter_rh': 0.000124533001245330,
            'r_C_root_litter_rh': 0.000122042341220423,
            'r_C_soil_fast_rh': 0.000152207001522070,
            'r_C_soil_slow_rh': 2.73972602739726e-5,
            'r_C_soil_passive_rh': 7.82778864970646e-6,
            'r_C_leaf_2_C_leaf_litter': 0.00833333333333333,
            'r_C_wood_2_C_wood_litter': 9.13242009132420e-5,
            'r_C_root_2_C_root_litter': 0.000124533001245330,
            'r_C_leaf_litter_2_C_soil_fast': 0.000340390203403902,
            'r_C_leaf_litter_2_C_soil_slow': 5.81154005811540e-5,
            'r_C_leaf_litter_2_C_soil_passive': 1.66044001660440e-5,
            'r_C_wood_litter_2_C_soil_fast': 7.47198007471980e-5,
            'r_C_wood_litter_2_C_soil_slow': 2.98879202988792e-5,
            'r_C_wood_litter_2_C_soil_passive': 1.99252801992528e-5,
            'r_C_root_litter_2_C_soil_fast': 7.47198007471980e-5,
            'r_C_root_litter_2_C_soil_slow': 3.48692403486924e-5,
            'r_C_root_litter_2_C_soil_passive': 1.74346201743462e-5
        }.items()
    }
    svs,dvs=msh.get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))
    func_dict={
        Function(k):v
        for k,v in {
            'NPP':msh.make_npp_func(dvs),
            'xi':msh.make_xi_func(dvs)
        }.items()
    }
    svs_0=msh.Observables(*map(lambda v: v[0],svs))
    # create a start parameter tuple for the mcmc.
    epa_old=msh.EstimatedParameters(
            beta_leaf=0.15,
	    beta_root=0.2,
	    lig_leaf=0.15,
	    f_leaf2metlit=0.28,
	    f_root2metlit=0.6,
	    k_leaf=1/365,
	    k_root=1/(365*5),
	    k_wood=1/(365*40),
	    k_metlit=0.5/(365*0.1),
	    k_mic=0.3/(365*0.137),
	    k_slowsom=0.3/(365*5),
	    k_passsom=0.3/(222.22*365),
	    C_metlit_0=0.05,
	    C_CWD_0=0.1,
	    C_mic_0=1,
	    C_passom_0=5,
        )
    )
    #epa_min=msh.EstimatedParameters(
    #    beta_leaf=0,
    #    beta_wood=0,
    #    T_0=-20,
    #    E=.1,
    #    KM=1,
    #    r_C_leaf_litter_rh=epa_0.r_C_leaf_litter_rh/100,
    #    r_C_wood_litter_rh=epa_0.r_C_wood_litter_rh/100,
    #    r_C_root_litter_rh=epa_0.r_C_root_litter_rh/100,
    #    r_C_soil_fast_rh=epa_0.r_C_soil_fast_rh/100,
    #    r_C_soil_slow_rh=epa_0.r_C_soil_slow_rh/100,
    #    r_C_soil_passive_rh=epa_0.r_C_soil_passive_rh/100,
    #    r_C_leaf_2_C_leaf_litter=epa_0.r_C_leaf_2_C_leaf_litter/100,       
    #    r_C_wood_2_C_wood_litter=epa_0.r_C_wood_2_C_wood_litter/100,
    #    r_C_root_2_C_root_litter=epa_0.r_C_root_2_C_root_litter/100,
    #    r_C_leaf_litter_2_C_soil_fast=epa_0.r_C_leaf_litter_2_C_soil_fast/100,
    #    r_C_leaf_litter_2_C_soil_slow=epa_0.r_C_leaf_litter_2_C_soil_slow/100,
    #    r_C_leaf_litter_2_C_soil_passive=epa_0.r_C_leaf_litter_2_C_soil_passive/100,
    #    r_C_wood_litter_2_C_soil_fast=epa_0.r_C_wood_litter_2_C_soil_fast/100,
    #    r_C_wood_litter_2_C_soil_slow=epa_0.r_C_wood_litter_2_C_soil_slow/100,
    #    r_C_wood_litter_2_C_soil_passive=epa_0.r_C_wood_litter_2_C_soil_passive/100,
    #    r_C_root_litter_2_C_soil_fast=epa_0.r_C_root_litter_2_C_soil_fast/100,
    #    r_C_root_litter_2_C_soil_slow=epa_0.r_C_root_litter_2_C_soil_slow/100,
    #    r_C_root_litter_2_C_soil_passive=epa_0.r_C_root_litter_2_C_soil_passive/100,
    #    C_leaf_0=0,
    #    C_wood_0=0,
    #    C_leaf_litter_0=0,
    #    C_wood_litter_0=0,
    #    C_soil_fast_0=0,
    #    C_soil_slow_0=0,
    #)


    #epa_max=msh.EstimatedParameters(
    #    beta_leaf=0.99,
    #    beta_wood=0.99,
    #    T_0=10,
    #    E=100,
    #    KM=100,
    #    r_C_leaf_litter_rh=epa_0.r_C_leaf_litter_rh*100,
    #    r_C_wood_litter_rh=epa_0.r_C_wood_litter_rh*100,
    #    r_C_root_litter_rh=epa_0.r_C_root_litter_rh*100,
    #    r_C_soil_fast_rh=epa_0.r_C_soil_fast_rh*100,
    #    r_C_soil_slow_rh=epa_0.r_C_soil_slow_rh*100,
    #    r_C_soil_passive_rh=epa_0.r_C_soil_passive_rh*100,
    #    r_C_leaf_2_C_leaf_litter=epa_0.r_C_leaf_2_C_leaf_litter*100,       
    #    r_C_wood_2_C_wood_litter=epa_0.r_C_wood_2_C_wood_litter*100,
    #    r_C_root_2_C_root_litter=epa_0.r_C_root_2_C_root_litter*100,
    #    r_C_leaf_litter_2_C_soil_fast=epa_0.r_C_leaf_litter_2_C_soil_fast*100,
    #    r_C_leaf_litter_2_C_soil_slow=epa_0.r_C_leaf_litter_2_C_soil_slow*100,
    #    r_C_leaf_litter_2_C_soil_passive=epa_0.r_C_leaf_litter_2_C_soil_passive*100,
    #    r_C_wood_litter_2_C_soil_fast=epa_0.r_C_wood_litter_2_C_soil_fast*100,
    #    r_C_wood_litter_2_C_soil_slow=epa_0.r_C_wood_litter_2_C_soil_slow*100,
    #    r_C_wood_litter_2_C_soil_passive=epa_0.r_C_wood_litter_2_C_soil_passive*100,
    #    r_C_root_litter_2_C_soil_fast=epa_0.r_C_root_litter_2_C_soil_fast*100,
    #    r_C_root_litter_2_C_soil_slow=epa_0.r_C_root_litter_2_C_soil_slow*100,
    #    r_C_root_litter_2_C_soil_passive=epa_0.r_C_root_litter_2_C_soil_passive*100,
    #    C_leaf_0=svs_0.cVeg,
    #    C_wood_0=svs_0.cVeg,
    #    C_leaf_litter_0=svs_0.cLitter,
    #    C_wood_litter_0=svs_0.cLitter,
    #    C_soil_fast_0=svs_0.cSoil,
    #    C_soil_slow_0=svs_0.cSoil,
    #)    
    cpa = msh.Constants(
        C_leaf_0 = svs.cleaf[0],
        C_root_0 = svs.croot[0],
        C_wood_0 = svs.cwood[0],
        cLitter_0 = svs.clitter[0],
        csoil_0 = svs.csoil[0],
        rh_0 = svs.rh[0],
        clay=0.2028,
        silt=0.2808,
        lig_wood=0.4,
        f_wood2CWD=1,
        f_metlit2mic=0.45,
        npp=npp, #array!!!
        number_of_months=10
    )

    #StartVector = msh.make_StartVector(mvs) 
    V_init= StartVector(
        C_leaf=svs_0.cVeg/3,
        C_wood=svs_0.cVeg/3,
        C_root=svs_0.cVeg/3,
        C_leaf_litter=svs_0.cLitter/3,
        C_wood_litter=svs_0.cLitter/3,
        C_root_litter=svs_0.cLitter/3,
        C_soil_fast=svs_0.cSoil/3,
        C_soil_slow=svs_0.cSoil/3,
        C_soil_passive=svs_0.cSoil/3,
        ra=svs_0.ra*86400,   # kg/m2/s kg/m2/day;,
        rh=svs_0.rh*86400   # kg/m2/s kg/m2/day;
    )
    return TestArgs(
        V_init=V_init,
        par_dict=par_dict,
        func_dict=func_dict,
        dvs=dvs,
        svs=svs,
        mvs=mvs,
        epa_0=epa_0,
        epa_min=epa_min,
        epa_max=epa_max,
        cpa=cpa
    )


