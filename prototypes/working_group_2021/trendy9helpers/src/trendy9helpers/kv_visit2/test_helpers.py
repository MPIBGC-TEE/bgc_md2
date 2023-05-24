# this file provides some variables for testing.
from collections import namedtuple
from sympy import Symbol, Function
from pathlib import Path
import netCDF4 as nc
import json 
from functools import lru_cache
from frozendict import frozendict

# deprecated
@lru_cache
def make_test_args(conf_dict,msh,mvs):
    TestArgs=namedtuple(
        "TestArgs",
        [
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
    par_dict={
        Symbol(k): v for k,v in {
            'beta_leaf': 0.6,
            'beta_wood': 0.25,
            'T_0': 2,
            'E': 6.5,
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
    svs,dvs=msh.get_global_mean_vars_2(conf_dict)
    svs_0=msh.Observables(*map(lambda v: v[0],svs))
    # create a start parameter tuple for the mcmc.
    epa_0=msh.EstimatedParameters(
        beta_leaf=0.6,
        beta_wood=0.25,
        T_0=2,
        E=6.5,
        KM=10,
        r_C_leaf_litter_rh=0.0004151100041511,
        r_C_wood_litter_rh=0.00012453300124533,
        r_C_root_litter_rh=0.000122042341220423,
        r_C_soil_fast_rh=0.00015220700152207,
        r_C_soil_slow_rh=3.73972602739726e-05,
        r_C_soil_passive_rh=8.82778864970646e-06,
        r_C_leaf_2_C_leaf_litter=0.00833333333333333,
        r_C_wood_2_C_wood_litter=9.1324200913242e-05,
        r_C_root_2_C_root_litter=0.00012453300124533,
        r_C_leaf_litter_2_C_soil_fast=0.000340390203403902,
        r_C_leaf_litter_2_C_soil_slow=5.8115400581154e-05,
        r_C_leaf_litter_2_C_soil_passive=1.6604400166044e-05,
        r_C_wood_litter_2_C_soil_fast=7.4719800747198e-05,
        r_C_wood_litter_2_C_soil_slow=2.98879202988792e-05,
        r_C_wood_litter_2_C_soil_passive=1.99252801992528e-05,
        r_C_root_litter_2_C_soil_fast=7.4719800747198e-05,
        r_C_root_litter_2_C_soil_slow=3.48692403486924e-05,
        r_C_root_litter_2_C_soil_passive=1.74346201743462e-05,
        C_leaf_0=0.051828761170322826,
        C_wood_0=1.970572690329994,
        C_leaf_litter_0=0.1202311902470766,
        C_wood_litter_0=0.2225433197876749,
        C_soil_fast_0=1.7309510511856925,
        C_soil_slow_0=2.4435101360092473        
    )
    #func_dict={
    #    Function(k):v
    #    for k,v in {
    #        'NPP':msh.make_npp_func(dvs),
    #        'xi':msh.make_xi_func(dvs,epa_0)
    #    }.items()
    #}
    epa_min=msh.EstimatedParameters(
        beta_leaf=0,
        beta_wood=0,
        T_0=-10,
        E=1,
        KM=1,
        r_C_leaf_litter_rh=epa_0.r_C_leaf_litter_rh/100,
        r_C_wood_litter_rh=epa_0.r_C_wood_litter_rh/100,
        r_C_root_litter_rh=epa_0.r_C_root_litter_rh/100,
        r_C_soil_fast_rh=epa_0.r_C_soil_fast_rh/100,
        r_C_soil_slow_rh=epa_0.r_C_soil_slow_rh/100,
        r_C_soil_passive_rh=epa_0.r_C_soil_passive_rh/100,
        r_C_leaf_2_C_leaf_litter=epa_0.r_C_leaf_2_C_leaf_litter/100,       
        r_C_wood_2_C_wood_litter=epa_0.r_C_wood_2_C_wood_litter/100,
        r_C_root_2_C_root_litter=epa_0.r_C_root_2_C_root_litter/100,
        r_C_leaf_litter_2_C_soil_fast=epa_0.r_C_leaf_litter_2_C_soil_fast/100,
        r_C_leaf_litter_2_C_soil_slow=epa_0.r_C_leaf_litter_2_C_soil_slow/100,
        r_C_leaf_litter_2_C_soil_passive=epa_0.r_C_leaf_litter_2_C_soil_passive/100,
        r_C_wood_litter_2_C_soil_fast=epa_0.r_C_wood_litter_2_C_soil_fast/100,
        r_C_wood_litter_2_C_soil_slow=epa_0.r_C_wood_litter_2_C_soil_slow/100,
        r_C_wood_litter_2_C_soil_passive=epa_0.r_C_wood_litter_2_C_soil_passive/100,
        r_C_root_litter_2_C_soil_fast=epa_0.r_C_root_litter_2_C_soil_fast/100,
        r_C_root_litter_2_C_soil_slow=epa_0.r_C_root_litter_2_C_soil_slow/100,
        r_C_root_litter_2_C_soil_passive=epa_0.r_C_root_litter_2_C_soil_passive/100,
        C_leaf_0=0,
        C_wood_0=0,
        C_leaf_litter_0=0,
        C_wood_litter_0=0,
        C_soil_fast_0=0,
        C_soil_slow_0=0,
    )


    epa_max=msh.EstimatedParameters(
        beta_leaf=0.99,
        beta_wood=0.99,
        T_0=5,
        E=15,
        KM=100,
        r_C_leaf_litter_rh=epa_0.r_C_leaf_litter_rh*100,
        r_C_wood_litter_rh=epa_0.r_C_wood_litter_rh*100,
        r_C_root_litter_rh=epa_0.r_C_root_litter_rh*100,
        r_C_soil_fast_rh=epa_0.r_C_soil_fast_rh*100,
        r_C_soil_slow_rh=epa_0.r_C_soil_slow_rh*100,
        r_C_soil_passive_rh=epa_0.r_C_soil_passive_rh*100,
        r_C_leaf_2_C_leaf_litter=epa_0.r_C_leaf_2_C_leaf_litter*100,       
        r_C_wood_2_C_wood_litter=epa_0.r_C_wood_2_C_wood_litter*100,
        r_C_root_2_C_root_litter=epa_0.r_C_root_2_C_root_litter*100,
        r_C_leaf_litter_2_C_soil_fast=epa_0.r_C_leaf_litter_2_C_soil_fast*100,
        r_C_leaf_litter_2_C_soil_slow=epa_0.r_C_leaf_litter_2_C_soil_slow*100,
        r_C_leaf_litter_2_C_soil_passive=epa_0.r_C_leaf_litter_2_C_soil_passive*100,
        r_C_wood_litter_2_C_soil_fast=epa_0.r_C_wood_litter_2_C_soil_fast*100,
        r_C_wood_litter_2_C_soil_slow=epa_0.r_C_wood_litter_2_C_soil_slow*100,
        r_C_wood_litter_2_C_soil_passive=epa_0.r_C_wood_litter_2_C_soil_passive*100,
        r_C_root_litter_2_C_soil_fast=epa_0.r_C_root_litter_2_C_soil_fast*100,
        r_C_root_litter_2_C_soil_slow=epa_0.r_C_root_litter_2_C_soil_slow*100,
        r_C_root_litter_2_C_soil_passive=epa_0.r_C_root_litter_2_C_soil_passive*100,
        C_leaf_0=svs_0.cVeg,
        C_wood_0=svs_0.cVeg,
        C_leaf_litter_0=svs_0.cLitter,
        C_wood_litter_0=svs_0.cLitter,
        C_soil_fast_0=svs_0.cSoil,
        C_soil_slow_0=svs_0.cSoil,
    ) 
    epa_opt=msh.EstimatedParameters(
        beta_leaf=0.6102169482865195, 
        beta_wood=0.26331553815787545, 
        T_0=1.9560345980471245, 
        E=6.951145421284498, 
        KM=12.73895376887386, 
        r_C_leaf_litter_rh=0.0012830039575098323, 
        r_C_wood_litter_rh=0.0010536416454437036, 
        r_C_root_litter_rh=0.00022271755326847413, 
        r_C_soil_fast_rh=0.00013707839288781872, 
        r_C_soil_slow_rh=3.228645064482276e-05, 
        r_C_soil_passive_rh=3.8079656062059605e-06, 
        r_C_leaf_2_C_leaf_litter=0.011755309034589333, 
        r_C_wood_2_C_wood_litter=0.00012990716959685548, 
        r_C_root_2_C_root_litter=8.243205281709114e-05, 
        r_C_leaf_litter_2_C_soil_fast=0.0014521759026031634, 
        r_C_leaf_litter_2_C_soil_slow=0.000200225210999593, 
        r_C_leaf_litter_2_C_soil_passive=8.380707345301035e-05, 
        r_C_wood_litter_2_C_soil_fast=3.19128931685429e-05, 
        r_C_wood_litter_2_C_soil_slow=7.278721749448471e-05, 
        r_C_wood_litter_2_C_soil_passive=3.275165103336979e-06, 
        r_C_root_litter_2_C_soil_fast=0.00044055481426693227, 
        r_C_root_litter_2_C_soil_slow=3.1019188662910444e-05, 
        r_C_root_litter_2_C_soil_passive=0.00012243099600679218, 
        C_leaf_0=0.042596582017273114, 
        C_wood_0=2.4493874052342517, 
        C_leaf_litter_0=0.16251047110808622, 
        C_wood_litter_0=0.17601405444541945, 
        C_soil_fast_0=1.8746682268250323, 
        C_soil_slow_0=1.9807766341505468
    )

    cpa = msh.Constants(
        cVeg_0=svs_0.cVeg,
        cLitter_0=svs_0.cLitter,
        cSoil_0=svs_0.cSoil,
        npp_0=dvs.npp[0],   # kg/m2/s kg/m2/day
        rh_0=svs_0.rh,   # kg/m2/s kg/m2/day
        #ra_0=svs_0.ra,   # kg/m2/s kg/m2/day
        #r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
        #r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
        #number_of_months=len(svs.rh)
        number_of_months=24 # for testing and tuning mcmc
    )

    ds=nc.Dataset(Path(conf_dict["dataPath"]).joinpath("VISIT_S2_cLitter.nc"))    
    return TestArgs(
        par_dict=par_dict,
        func_dict=msh.make_func_dict(dvs),
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
        start_date=msh.start_date()
    )


def confDict():
    p = Path(__file__).parent.joinpath("config.json")
    with p.open(
            "r"
        ) as f:
        conf_dict=frozendict(json.load(f))
    
    return conf_dict    


def lats_lons():
    conf_dict = confDict() 
    ds=nc.Dataset(Path(conf_dict["dataPath"]).joinpath("VISIT_S2_cLitter.nc"))    
    lats=ds.variables["lat"][:]
    lons=ds.variables["lon"][:]
    return lats.data, lons.data
