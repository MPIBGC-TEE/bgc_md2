# this file provides some variables for testing.
import json 
import numpy as np
import netCDF4 as nc
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
            "epa_opt",
            "cpa",
            "lats",
            "lons",
            "obs_arr",
            "start_date"
        ]
    )
    par_dict={
        Symbol(k): v for k,v in {
            "r_C_AGML_rh": 0.006780821917808219,
            "r_C_AGSL_rh": 0.03547945205479452,
            "r_C_AGMS_rh": 0.008,
            "r_C_YHMS_rh": 0.002465753424657534,
            "k_C_BGDL": 0.0273972602739726,
            "k_C_BGRL": 0.000821917808219178,
            "k_C_BGMS": 0.0018082191780821918,
            "k_C_SHMS": 5.479452054794521e-05,
            "r_C_AGML_2_C_AGMS": 0.0055479452054794515,
            "r_C_AGMS_2_C_YHMS": 0.012,
            "r_C_YHMS_2_C_AGMS": 0.002712328767123288,
            "r_C_YHMS_2_C_SHMS": 0.00030136986301369865,
            "fwt": 0.62,
            "fgv": 0.3,
            "fco": 0.95,
            "fml": 0.8,
            "fd": 0.75,
            "k_C_NWT": 0.0013698630136986301,
            "k_C_AGWT": 9.132420091324201e-05,
            "k_C_TR": 0.000136986301369863,
            "k_C_GVF": 9.132420091324201e-05,
            "k_C_GVR": 0.000136986301369863,
            "f_C_AGSL_2_C_AGMS": 0.15,
            "f_C_BGRL_2_C_SHMS": 0.5,
        }.items()
    }
    svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
    obs_arr= np.column_stack((np.repeat(svs.cVeg, 12),np.repeat(svs.cLitter, 12),np.repeat(svs.cSoil, 12),svs.rh,svs.ra))
    func_dict = msh.make_func_dict(dvs)
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
    epa_opt=msh.EstimatedParameters(
        fwt=0.5907770914828289, 
        fgv=0.10708374044873868, 
        fco=0.9502719613629499, 
        fml=0.6985590765466911, 
        fd=0.8108017779961694, 
        k_C_NWT=0.0018600810916478165, 
        k_C_AGWT=0.00017354142452106252, 
        k_C_TR=0.00016065843641210772, 
        k_C_GVF=0.00022102017216433633, 
        k_C_GVR=0.00017926856125131916, 
        f_C_AGSL_2_C_AGMS=0.20853426509202325, 
        f_C_BGRL_2_C_SHMS=0.24638112975102788, 
        C_NWT_0=svs_0.cVeg * 0.04,#0.39641121927763323, 
        C_AGWT_0=svs_0.cVeg * 0.4,#1.0098899271611432, 
        C_GVF_0=svs_0.cVeg * 0.06,#0.1784893310039542, 
        C_GVR_0=svs_0.cVeg * 0.07,#2.1680315400436174, 
        C_AGML_0=svs_0.cLitter * 0.17,#0.1251689278629053, 
        C_AGSL_0=svs_0.cLitter * 0.06,#0.005800531050824444, 
        C_BGDL_0=svs_0.cLitter * 0.11,#0.0484130929152639, 
        C_AGMS_0=svs_0.cSoil * 0.006,#0.10074331291791151, 
        C_YHMS_0=svs_0.cSoil * 0.018,#0.5036084965444287, 
        C_SHMS_0=svs_0.cSoil * 0.975,#8.080067914918454,
)
    
    cpa = msh.Constants(
        cVeg_0=svs_0.cVeg,
        cLitter_0=svs_0.cLitter,
        cSoil_0=svs_0.cSoil,
        npp_0=dvs.npp[0] * 86400,   # kg/m2/s kg/m2/day
        rh_0=svs_0.rh * 86400,   # kg/m2/s kg/m2/day
        ra_0=svs_0.ra * 86400,   # kg/m2/s kg/m2/day
        #r_C_NWT_rh=0,
        #r_C_AGWT_rh=0,
        #r_C_TR_rh=0,
        #r_C_GVF_rh=0,
        #r_C_GVR_rh=0,
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
        number_of_months=320*12 # for testing and tuning mcmc
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
    ds=nc.Dataset(Path(conf_dict["dataPath"]).joinpath("ISAM_S2_cLitter.nc"))    
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
        epa_opt=epa_opt,
        cpa=cpa,
        lats=ds.variables["lat"][:],
        lons=ds.variables["lon"][:],
        obs_arr=obs_arr,
        start_date=msh.start_date()        
    )


