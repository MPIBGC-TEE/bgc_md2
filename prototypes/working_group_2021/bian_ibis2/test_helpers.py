# this file provides some variables for testing.
from collections import namedtuple
from sympy import Symbol, Function
from pathlib import Path
import netCDF4 as nc
import json 
from functools import lru_cache

# +
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
    par_dict={
        Symbol(k):v for k,v in 
        {
             "beta_leaf": 0.6,
             "beta_wood": 0.25,
             "r_C_leaf_rh": 0.000833333333333334,
             "r_C_wood_rh": 9.13242009132421e-6,
             "r_C_root_rh": 2.49066002490660e-5,
             "r_C_mll_rh": 0.0900000000000000,
             "r_C_mwl_rh": 0.000600000000000000,
             "r_C_mrl_rh": 0.0600000000000000,
             "r_C_sll_rh": 0.00700000000000000,
             "r_C_swl_rh": 0.000700000000000000,
             "r_C_srl_rh": 0.00350000000000000,
             "r_C_lll_rh": 0,
             "r_C_lwl_rh": 0,
             "r_C_lrl_rh": 0,
             "r_C_mic_rh": 6.08828006088280e-5,
             "r_C_prot_rh": 1.24533001245330e-5,
             "r_C_nonprot_rh": 2.24159402241594e-5,
             "r_C_pass_rh": 6.18395303326810e-6,
             "r_C_leaf_2_C_mll": 0.00583333333333333,
             "r_C_leaf_2_C_sll": 0.000833333333333333,
             "r_C_leaf_2_C_lll": 0.000833333333333333,
             "r_C_wood_2_C_mwl": 9.13242009132420e-6,
             "r_C_wood_2_C_swl": 5.47945205479452e-5,
             "r_C_wood_2_C_lwl": 1.82648401826484e-5,
             "r_C_root_2_C_mrl": 2.49066002490660e-5,
             "r_C_root_2_C_srl": 4.98132004981320e-5,
             "r_C_root_2_C_lrl": 2.49066002490660e-5,
             "r_C_mll_2_C_mic": 0.0600000000000000,
             "r_C_mwl_2_C_mic": 0.000400000000000000,
             "r_C_mrl_2_C_mic": 0.0400000000000000,
             "r_C_sll_2_C_mic": 0.00300000000000000,
             "r_C_swl_2_C_mic": 0.000300000000000000,
             "r_C_srl_2_C_mic": 0.00150000000000000,
             "r_C_lll_2_C_prot": 0.00500000000000000,
             "r_C_lll_2_C_nonprot": 0.00500000000000000,
             "r_C_lwl_2_C_prot": 0.000500000000000000,
             "r_C_lwl_2_C_nonprot": 0.000500000000000000,
             "r_C_lrl_2_C_prot": 0.00250000000000000,
             "r_C_lrl_2_C_nonprot": 0.00250000000000000,
             "r_C_mic_2_C_prot": 4.56621004566210e-5,
             "r_C_mic_2_C_nonprot": 4.56621004566210e-5,
             "r_C_prot_2_C_mic": 8.71731008717310e-5,
             "r_C_prot_2_C_pass": 0.000149439601494396,
             "r_C_nonprot_2_C_mic": 0.000102117061021171,
             "r_C_nonprot_2_C_pass": 0.000124533001245330,
             "r_C_pass_2_C_mic": 1.64383561643836e-6
        }.items()
    }
    svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
    svs_0=msh.Observables(*map(lambda v: v[0],svs))
    # create a start parameter tuple for the mcmc.
#     epa_0=msh.EstimatedParameters(
#          beta_leaf=0.3,
#          beta_wood=0.3,
         
#          r_C_leaf_2_C_mll=0.4/365,    
#          r_C_wood_2_C_mwl=0.3/365,     
#          r_C_root_2_C_mrl=0.3/365,     
#          r_C_leaf_2_C_sll=0.35/365,    
#          r_C_wood_2_C_swl=0.35/365,    
#          r_C_root_2_C_srl=0.3/365,    
#          r_C_prot_2_C_mic=0.5/365,     
#          r_C_nonprot_2_C_mic=0.5/365,  
#          r_C_mic_2_C_prot=0.5/365,    
#          r_C_mic_2_C_nonprot=0.5/365,  
#          r_C_prot_2_C_pass=0.5/365,  
#          r_C_nonprot_2_C_pass=0.5/365,

#          k_C_mic=0.001,            
#          k_C_protsom=0.0001,       
#          k_C_nonprotsom=0.0001,      
#          k_C_passsom=0.0001,          
#          C_wood_0=svs_0.cVeg/3, # 4.487521733805875,           
#          C_root_0=svs_0.cVeg/3, # 4.936273907192815,           
#          C_mll_0=svs_0.cLitter/6, # 0.0019124902366312892,            
#          C_mwl_0=svs_0.cLitter/6, # 0.04098193364207623,            
#          C_sll_0=svs_0.cLitter/6, # 0.004098193364198974,            
#          C_swl_0=svs_0.cLitter/6, # 0.24589160185151485,            
#          C_lll_0=svs_0.cLitter/6, # 0.00378294772079453,            
#          C_mrl_0=svs_0.cSoil/7, # 0.0012294580092597053,           
#          C_srl_0=svs_0.cSoil/7, # 0.04917832037053079,          
#          C_lrl_0=svs_0.cSoil/7, # 0.024467280828709734,           
#          C_mic_0=svs_0.cSoil/7, # 4.298575136296554,           
#          C_prot_0=svs_0.cSoil/7, # 1.2726585066414637,          
#          C_nonprot_0=svs_0.cSoil/7, # 1.272658506641464,     
#     )

    epa_0=msh.EstimatedParameters(
         beta_leaf=0.6,
         beta_wood=0.25,
            
         r_C_leaf_2_C_mll=0.00583333333333333,
         r_C_leaf_2_C_sll=0.000833333333333333,
         r_C_leaf_2_C_lll=0.000833333333333333,
         r_C_wood_2_C_mwl=5.13242009132420e-6,
         r_C_wood_2_C_swl=2.47945205479452e-5,
         r_C_wood_2_C_lwl=1.82648401826484e-5,
         r_C_root_2_C_mrl=2.49066002490660e-5,
         r_C_root_2_C_srl=4.98132004981320e-5,
         r_C_root_2_C_lrl=2.49066002490660e-5,  
                
              
         r_C_prot_2_C_mic=8.71731008717310e-5,
         r_C_nonprot_2_C_mic=0.000102117061021171,
         r_C_mic_2_C_prot=4.56621004566210e-5,
         r_C_mic_2_C_nonprot=4.56621004566210e-5,
         r_C_prot_2_C_pass=0.000149439601494396,
         r_C_nonprot_2_C_pass=0.000124533001245330, 
          
         r_C_mll_rh=0.0004,
         r_C_mwl_rh=0.0001,
         r_C_mrl_rh=0.0001,
         r_C_sll_rh=0.0004,
         r_C_swl_rh=0.0001,
         r_C_srl_rh=0.0001,         
         r_C_lll_rh=0.0004,
         r_C_lwl_rh=0.0001,
         r_C_lrl_rh=0.0001,            
         r_C_mic_rh=0.0001,           
         r_C_prot_rh=0.00015,
         r_C_nonprot_rh=3.7e-05,        
         r_C_pass_rh=8.8e-06,           
                
         C_wood_0=2.3870336362171884,           
         C_root_0=1.1815816499293954,          
         C_mll_0=0.003999340544806679,            
         C_mwl_0=0.03570839772146303,            
         C_sll_0=0.008570015453134308,           
         C_swl_0=0.21425038632795695,          
         C_lll_0=0.008570015453134308,            
         C_mrl_0=0.00048206336923888824,       
         C_srl_0=0.019282534769611434,          
         C_lrl_0=0.009641267384805717,           
         C_mic_0=0.6843048316422622,           
         C_prot_0=0.18650353669163452,          
         C_nonprot_0=0.18650353669163452,  
    )

    #func_dict={
    #    Function(k):v
    #    for k,v in {
    #        'NPP':msh.make_npp_func(dvs),
    #        'xi':msh.make_xi_func(dvs,epa_0)
    #    }.items()
    #}
    epa_min=msh.EstimatedParameters(
         beta_leaf=0.01,
         beta_wood=0.01,

         r_C_leaf_2_C_mll=epa_0.r_C_leaf_2_C_mll/100,
         r_C_leaf_2_C_sll=epa_0.r_C_leaf_2_C_sll/100,
         r_C_leaf_2_C_lll=epa_0.r_C_leaf_2_C_lll/100,
         r_C_wood_2_C_mwl=epa_0.r_C_wood_2_C_mwl/100,
         r_C_wood_2_C_swl=epa_0.r_C_wood_2_C_swl/100,
         r_C_wood_2_C_lwl=epa_0.r_C_wood_2_C_lwl/100,
         r_C_root_2_C_mrl=epa_0.r_C_root_2_C_mrl/100,
         r_C_root_2_C_srl=epa_0.r_C_root_2_C_srl/100,
         r_C_root_2_C_lrl=epa_0.r_C_root_2_C_lrl/100,  
                
              
         r_C_prot_2_C_mic=epa_0.r_C_prot_2_C_mic/100,
         r_C_nonprot_2_C_mic=epa_0.r_C_nonprot_2_C_mic/100,
         r_C_mic_2_C_prot=epa_0.r_C_mic_2_C_prot/100,
         r_C_mic_2_C_nonprot=epa_0.r_C_mic_2_C_nonprot/100,
         r_C_prot_2_C_pass=epa_0.r_C_prot_2_C_pass/100,
         r_C_nonprot_2_C_pass=epa_0.r_C_nonprot_2_C_pass/100, 
          
         r_C_mll_rh=epa_0.r_C_mll_rh/100,
         r_C_mwl_rh=epa_0.r_C_mwl_rh/100,
         r_C_mrl_rh=epa_0.r_C_mrl_rh/100,
         r_C_sll_rh=epa_0.r_C_sll_rh/100,
         r_C_swl_rh=epa_0.r_C_swl_rh/100,
         r_C_srl_rh=epa_0.r_C_srl_rh/100,        
         r_C_lll_rh=epa_0.r_C_lll_rh/100,
         r_C_lwl_rh=epa_0.r_C_lwl_rh/100,
         r_C_lrl_rh=epa_0.r_C_lrl_rh/100,            
         r_C_mic_rh=epa_0.r_C_mic_rh/100,           
         r_C_prot_rh=epa_0.r_C_prot_rh/100,
         r_C_nonprot_rh=epa_0.r_C_nonprot_rh/100,        
         r_C_pass_rh=epa_0.r_C_pass_rh/100,   
            
         C_wood_0=0,           
         C_root_0=0,           
         C_mll_0=0,            
         C_mwl_0=0,            
         C_sll_0=0,            
         C_swl_0=0,            
         C_lll_0=0,            
         C_mrl_0=0,           
         C_srl_0=0,          
         C_lrl_0=0,           
         C_mic_0=0,           
         C_prot_0=0,          
         C_nonprot_0=0, 
    )


    epa_max=msh.EstimatedParameters(
         beta_leaf=0.99,
         beta_wood=0.99,
         r_C_leaf_2_C_mll=epa_0.r_C_leaf_2_C_mll*100,
         r_C_leaf_2_C_sll=epa_0.r_C_leaf_2_C_sll*100,
         r_C_leaf_2_C_lll=epa_0.r_C_leaf_2_C_lll*100,
         r_C_wood_2_C_mwl=epa_0.r_C_wood_2_C_mwl*100,
         r_C_wood_2_C_swl=epa_0.r_C_wood_2_C_swl*100,
         r_C_wood_2_C_lwl=epa_0.r_C_wood_2_C_lwl*100,
         r_C_root_2_C_mrl=epa_0.r_C_root_2_C_mrl*100,
         r_C_root_2_C_srl=epa_0.r_C_root_2_C_srl*100,
         r_C_root_2_C_lrl=epa_0.r_C_root_2_C_lrl*100,  
                
              
         r_C_prot_2_C_mic=epa_0.r_C_prot_2_C_mic*100,
         r_C_nonprot_2_C_mic=epa_0.r_C_nonprot_2_C_mic*100,
         r_C_mic_2_C_prot=epa_0.r_C_mic_2_C_prot*100,
         r_C_mic_2_C_nonprot=epa_0.r_C_mic_2_C_nonprot*100,
         r_C_prot_2_C_pass=epa_0.r_C_prot_2_C_pass*100,
         r_C_nonprot_2_C_pass=epa_0.r_C_nonprot_2_C_pass*100, 
          
         r_C_mll_rh=epa_0.r_C_mll_rh*100,
         r_C_mwl_rh=epa_0.r_C_mwl_rh*100,
         r_C_mrl_rh=epa_0.r_C_mrl_rh*100,
         r_C_sll_rh=epa_0.r_C_sll_rh*100,
         r_C_swl_rh=epa_0.r_C_swl_rh*100,
         r_C_srl_rh=epa_0.r_C_srl_rh*100,        
         r_C_lll_rh=epa_0.r_C_lll_rh*100,
         r_C_lwl_rh=epa_0.r_C_lwl_rh*100,
         r_C_lrl_rh=epa_0.r_C_lrl_rh*100,            
         r_C_mic_rh=epa_0.r_C_mic_rh*100,           
         r_C_prot_rh=epa_0.r_C_prot_rh*100,
         r_C_nonprot_rh=epa_0.r_C_nonprot_rh*100,        
         r_C_pass_rh=epa_0.r_C_pass_rh*100,   
            
         C_wood_0=svs_0.cVeg,           
         C_root_0=svs_0.cVeg,           
         C_mll_0=svs_0.cLitter,            
         C_mwl_0=svs_0.cLitter,            
         C_sll_0=svs_0.cLitter,            
         C_swl_0=svs_0.cLitter,            
         C_lll_0=svs_0.cLitter,            
         C_mrl_0=svs_0.cSoil,           
         C_srl_0=svs_0.cSoil,          
         C_lrl_0=svs_0.cSoil,           
         C_mic_0=svs_0.cSoil,           
         C_prot_0=svs_0.cSoil,          
         C_nonprot_0=svs_0.cSoil, 
    )

    epa_opt=msh.EstimatedParameters(
        beta_leaf=0.5573756820192807, 
        beta_wood=0.24700794523783445, 
        r_C_leaf_2_C_mll=0.017528600518462494, 
        r_C_wood_2_C_mwl=1.6799097644290409e-06, 
        r_C_root_2_C_mrl=7.844130296528675e-06, 
        r_C_leaf_2_C_sll=0.0038705461605337975, 
        r_C_wood_2_C_swl=9.794679284314383e-05, 
        r_C_root_2_C_srl=7.1768211563668116e-06, 
        r_C_leaf_2_C_lll=0.005764905864736112, 
        r_C_wood_2_C_lwl=5.1804843333305954e-05, 
        r_C_root_2_C_lrl=0.00012753787039784646, 
        r_C_prot_2_C_mic=3.314461112947594e-05, 
        r_C_nonprot_2_C_mic=7.191960979251249e-05, 
        r_C_mic_2_C_prot=5.212347171735692e-05, 
        r_C_mic_2_C_nonprot=8.136318122484778e-05, 
        r_C_prot_2_C_pass=0.00020883969011529447, 
        r_C_nonprot_2_C_pass=0.0002788541490132211, 
        C_wood_0=2.5629353516907303, 
        C_root_0=1.4574759512852167, 
        C_mll_0=0.014534448208301288, 
        C_mwl_0=0.042334122827142925, 
        C_sll_0=0.013945880688603408, 
        C_swl_0=0.2185156564799679, 
        C_lll_0=0.011072695289578287, 
        C_mrl_0=0.32744881372351026, 
        C_srl_0=0.09500687424632417, 
        C_lrl_0=0.3478144969923375, 
        C_mic_0=0.9507532888112129, 
        C_prot_0=0.2592878819753468, 
        C_nonprot_0=0.1895759738112035, 
        r_C_mll_rh=0.00047840517766109655, 
        r_C_mwl_rh=0.00023154691864402867, 
        r_C_mrl_rh=0.00021520341200739518, 
        r_C_sll_rh=0.0017959007188430528, 
        r_C_swl_rh=0.0004692454377863136, 
        r_C_srl_rh=0.00038440762193817134, 
        r_C_lll_rh=0.0006900563753390039, 
        r_C_lwl_rh=0.00022943454556341198, 
        r_C_lrl_rh=0.00044676628181839184, 
        r_C_mic_rh=0.000288705294601301, 
        r_C_prot_rh=0.0008114136119091693, 
        r_C_nonprot_rh=4.761002328529026e-05, 
        r_C_pass_rh=5.492295571043415e-05,
    )

    cpa=msh.Constants(
         cVeg_0=svs_0.cVeg,
         cLitter_0=svs_0.cLitter,
         cSoil_0=svs_0.cSoil,
         npp_0=dvs.npp[0], # kg/m2/day
         rh_0=svs_0.rh,    # kg/m2/day
         ra_0=svs_0.ra,    # kg/m2/day
                       
         r_C_mll_2_C_mic=0.0600000000000000,
         r_C_mwl_2_C_mic=0.000400000000000000,
         r_C_mrl_2_C_mic=0.0400000000000000,
         r_C_sll_2_C_mic=0.00300000000000000,
         r_C_swl_2_C_mic=0.000300000000000000,
         r_C_srl_2_C_mic=0.00150000000000000,    
         r_C_pass_2_C_mic=1.64383561643836e-6,       
            
         r_C_lll_2_C_prot=0.00500000000000000,
         r_C_lll_2_C_nonprot=0.00500000000000000,
         r_C_lwl_2_C_prot=0.000500000000000000,
         r_C_lwl_2_C_nonprot=0.000500000000000000,
         r_C_lrl_2_C_prot=0.00250000000000000,
         r_C_lrl_2_C_nonprot=0.00250000000000000,
       
     number_of_months=len(svs.rh)
    )

    StartVector = msh.make_StartVector(mvs) 
    V_init = StartVector(
        C_leaf=cpa.cVeg_0-(epa_0.C_wood_0 + epa_0.C_root_0),
        C_wood=epa_0.C_wood_0,            
        C_root=epa_0.C_root_0,
            
        C_mll=epa_0.C_mll_0,
        C_mwl=epa_0.C_mwl_0,
        C_mrl=epa_0.C_mrl_0,
        C_sll=epa_0.C_sll_0,
        C_swl=epa_0.C_swl_0,
        C_srl=epa_0.C_srl_0,
        C_lll=epa_0.C_lll_0,
        C_lwl=cpa.cLitter_0-(epa_0.C_mll_0 + epa_0.C_mwl_0 + epa_0.C_sll_0 + epa_0.C_swl_0 + epa_0.C_lll_0),
        C_lrl=epa_0.C_lrl_0,
            
        C_mic=epa_0.C_mic_0,
        C_prot=epa_0.C_prot_0,
        C_nonprot=epa_0.C_nonprot_0,
        C_pass=cpa.cSoil_0-(epa_0.C_mrl_0 + epa_0.C_srl_0 + epa_0.C_lrl_0 + epa_0.C_mic_0 + epa_0.C_prot_0 + epa_0.C_nonprot_0),

        ra=cpa.ra_0,
        rh=cpa.rh_0
    )
    ds=nc.Dataset(Path(conf_dict["dataPath"]).joinpath("IBIS_S2_cLitter.nc"))    
    return TestArgs(
        V_init=V_init,
        par_dict=par_dict,
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
        start_date=msh.start_date()
    )
# -


