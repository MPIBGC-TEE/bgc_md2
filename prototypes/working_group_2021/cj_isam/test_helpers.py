# this file provides some variables for testing.
import json 
import numpy as np
from collections import namedtuple
from sympy import Symbol, Function
from pathlib import Path
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
            "cpa",
            "obs_arr"
        ]
    )
    par_dict={
        Symbol(k): v for k,v in {
            "beta_NWT": 0.3,
            "beta_AGWT": 0.3,
            "beta_TR": 0.10000000000000003,
            "beta_GVF": 0.15,
            "beta_GVR": 0.15,
            "r_C_NWT_rh": 0,
            "r_C_AGWT_rh": 0,
            "r_C_TR_rh": 0,
            "r_C_GVF_rh": 0,
            "r_C_GVR_rh": 0,
            "r_C_AGML_rh": 0.00678082191780822,
            "r_C_AGSL_rh": 0.0354794520547945,
            "r_C_AGMS_rh": 0.00800000000000000,
            "r_C_YHMS_rh": 0.00246575342465753,
            "r_C_BGDL_rh": 0.0200000000000000,
            "r_C_BGRL_rh": 0.000600000000000000,
            "r_C_BGMS_rh": 0.00132000000000000,
            "r_C_SHMS_rh": 4.00000000000000e-5,
            "r_C_NWT_2_C_AGML": 0.00116438356164384,
            "r_C_NWT_2_C_AGSL": 0.000205479452054795,
            "r_C_AGWT_2_C_AGSL": 9.13242009132420e-5,
            "r_C_TR_2_C_BGDL": 9.21544209215442e-5,
            "r_C_TR_2_C_BGRL": 3.23785803237858e-5,
            "r_C_GVF_2_C_AGML": 7.76255707762557e-5,
            "r_C_GVF_2_C_AGSL": 1.36986301369863e-5,
            "r_C_GVR_2_C_BGDL": 9.21544209215442e-5,
            "r_C_GVR_2_C_BGRL": 3.23785803237858e-5,
            "r_C_AGML_2_C_AGMS": 0.00554794520547945,
            "r_C_AGSL_2_C_AGMS": 0.00760273972602740,
            "r_C_AGSL_2_C_YHMS": 0.00760273972602740,
            "r_C_AGMS_2_C_YHMS": 0.0120000000000000,
            "r_C_YHMS_2_C_AGMS": 0.00271232876712329,
            "r_C_YHMS_2_C_SHMS": 0.000301369863013699,
            "r_C_BGDL_2_C_SHMS": 0.00739726027397260,
            "r_C_BGRL_2_C_BGMS": 0.000110958904109589,
            "r_C_BGRL_2_C_SHMS": 0.000110958904109589,
            "r_C_BGMS_2_C_SHMS": 0.000488219178082192,
            "r_C_SHMS_2_C_BGMS": 1.47945205479452e-5
        }.items()
    }
    svs,dvs=msh.get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))
    obs_arr= np.column_stack((np.repeat(svs.cVeg, 12),np.repeat(svs.cLitter, 12),np.repeat(svs.cSoil, 12),svs.rh,svs.ra))
    func_dict={
        Function(k):v
        for k,v in {
            'NPP':msh.make_npp_func(dvs),
            'xi':msh.make_xi_func(dvs)
        }.items()
    }
    svs_0=msh.Observables(*map(lambda v: v[0],svs))
    # create a start parameter tuple for the mcmc.
    epa_0=msh.EstimatedParameters(
        fwt=0.62,
	fgv=0.3,
	fco=0.95,
	fml=0.8,
	fd=0.75,
	k_C_NWT=0.0013698630136986301,
	k_C_AGWT=9.132420091324201e-05,
	k_C_TR=0.000136986301369863,
	k_C_GVF=9.132420091324201e-05,
	k_C_GVR=0.000136986301369863,
	f_C_AGSL_2_C_AGMS=0.15,
	f_C_BGRL_2_C_SHMS=0.5,
	C_NWT_0=0.6298059149319745,
	C_AGWT_0=9.447088723979736,
	C_GVF_0=4.723544361989868,
	C_GVR_0=3.4639325321259533,
	C_AGML_0=0.08922250461536606,
	C_AGSL_0=0.02085168231869287,
	C_BGDL_0=0.01941901571040178,
	C_AGMS_0=0.0446000716845726,
	C_YHMS_0=0.08792074042111502,
	C_SHMS_0=8.152881730528613
    )
    epa_min=msh.EstimatedParameters(
        fwt=0.5,
        fgv=0.1,
        fco=0.6,
        fml=0.6,
        fd=0.6,
        k_C_NWT=1/(365*10),
        k_C_AGWT=1/(365*40),
        k_C_TR=1/(365*40),
        k_C_GVF=1/(365*40),
        k_C_GVR=1/(365*40),
        f_C_AGSL_2_C_AGMS=0.1*0.3,
        f_C_BGRL_2_C_SHMS=0.1,
        C_NWT_0=0,
        C_AGWT_0=0,
        C_GVF_0=0,
        C_GVR_0=0,
        C_AGML_0=0,
        C_AGSL_0=0,
        C_BGDL_0=0,
        C_AGMS_0=0,
        C_YHMS_0=0,
        C_SHMS_0=0,
    )


    epa_max=msh.EstimatedParameters(
        fwt=0.8,
        fgv=0.3,
        fco=0.99,
        fml=0.9,
        fd=0.9,
        k_C_NWT=1/(365*1),
        k_C_AGWT=1/(365*10),
        k_C_TR=1/(365*10),
        k_C_GVF=1/(365*10),
        k_C_GVR=1/(365*10),
        f_C_AGSL_2_C_AGMS=0.9*0.3,
        f_C_BGRL_2_C_SHMS=0.9,
        C_NWT_0=svs_0.cVeg,
        C_AGWT_0=svs_0.cVeg,
        C_GVF_0=svs_0.cVeg,
        C_GVR_0=svs_0.cVeg,
        C_AGML_0=svs_0.cLitter,
        C_AGSL_0=svs_0.cLitter,
        C_BGDL_0=svs_0.cLitter,
        C_AGMS_0=svs_0.cSoil,
        C_YHMS_0=svs_0.cSoil,
        C_SHMS_0=svs_0.cSoil,
    )
    
    cpa = msh.Constants(
        cVeg_0=svs_0.cVeg,
        cLitter_0=svs_0.cLitter,
        cSoil_0=svs_0.cSoil,
        npp_0=dvs.npp[0] * 86400,   # kg/m2/s kg/m2/day
        rh_0=svs_0.rh * 86400,   # kg/m2/s kg/m2/day
        ra_0=svs_0.ra * 86400,   # kg/m2/s kg/m2/day
        r_C_NWT_rh=0,
        r_C_AGWT_rh=0,
        r_C_TR_rh=0,
        r_C_GVF_rh=0,
        r_C_GVR_rh=0,
        r_C_AGML_rh=0.55*4.5/365,
        r_C_AGSL_rh=0.7*18.5/365,
        r_C_AGMS_rh=0.4*7.3/365,
        r_C_YHMS_rh=0.45*2.0/365,
        k_C_BGDL=10/365,
        k_C_BGRL=0.3/365,
        k_C_BGMS=0.66/365,
        k_C_SHMS=0.02/365,
        r_C_AGML_2_C_AGMS=0.45*4.5/365,
        r_C_AGMS_2_C_YHMS=0.6*7.3/365,
        r_C_YHMS_2_C_AGMS=0.9*0.55*2.0/365,
        r_C_YHMS_2_C_SHMS=0.1*0.55*2.0/365,
        #number_of_months=len(svs.rh)
        number_of_months=120 # for testing and tuning mcmc
    )
    StartVector = msh.make_StartVector(mvs) 
    V_init= StartVector(
        C_NWT=svs_0.cVeg/5,
        C_AGWT=svs_0.cVeg/5,
        C_TR=svs_0.cVeg/5,
        C_GVF=svs_0.cVeg/5,
        C_GVR=svs_0.cVeg/5,
        C_AGML=svs_0.cLitter/4,
        C_AGSL=svs_0.cLitter/4,
        C_BGDL=svs_0.cLitter/4,
        C_BGRL=svs_0.cLitter/4,
        C_AGMS=svs_0.cSoil/4,
        C_YHMS=svs_0.cSoil/4,
        C_SHMS=svs_0.cSoil/4,
        C_BGMS=svs_0.cSoil/4,
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
        cpa=cpa,
        obs_arr=obs_arr
    )


