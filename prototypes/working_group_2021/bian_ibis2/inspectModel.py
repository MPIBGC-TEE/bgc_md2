# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Use the common infrastructure to inspect the model
#
# This illustrative notebook shows how to use the factored out model description in ```source.py``` the model specific functions in ```model_specific_helpers.py``` and the general infrastructure in ```general_helpers.py``` to inspect the model.
#
# ## Preconditions
# This is the next step after the ```createModel.py``` notebook where everything was defined in one place, the notebook.
# In order to be able to use the symbolic descriptions and functions in other code (e.g tests, scripts or other  notebooks)  
# we have disassembled it, moving the functions into the seperate file ```model_specific_helpers.py``` and the model description to ```source.py```.
#
# ## Applications
# 1. Inspect the structure of the model with symbolic tools
# 1. Run the model forward with a guessed  parameter set
# 1. Optimize the parameter set using data assimilation
# 1. Use the optimized paramerts to run the tracability analysis.
#


# +

# adjust the output to full width
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# make changes to imported files immidiately available 
# avoiding the need to reload (in most cases)
# %load_ext autoreload
# %autoreload 2

from pathlib import Path
import json 
from sympy import  Symbol, Function 
import numpy as np
import matplotlib.pyplot as plt
from ComputabilityGraphs.CMTVS import CMTVS
import CompartmentalSystems.helpers_reservoir as hr
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c
import bgc_md2.helper as h
import bgc_md2.display_helpers as dh

# imports from new files 
import sys
sys.path.insert(0,'..')
from source import mvs 
import model_specific_helpers_2 as msh
from general_helpers import day_2_month_index, make_B_u_funcs_2 
import general_helpers as gh

# -

# we can also print the whole mass balance equation
dh.mass_balance_equation(mvs)

# we can also plot a picture
h.compartmental_graph(mvs)

# +
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

#msh.download_my_TRENDY_output(conf_dict)

# + codehighlighter=[[5, 6], [23, 33], [5, 6], [23, 33]]
#     # Read NetCDF data  ********************************************************
svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
print(dvs.npp.shape)
print(svs.cVeg.shape)
svs_0=msh.Observables(*map(lambda v: v[0],svs))
dvs_0 = msh.Drivers(*map(lambda v: v[0],dvs))
print('svs_0=',svs_0)
print('dvs_0=',dvs_0)
# -
# To be able to run the model forward we not only have to replace parameter symbols by values but symbolic functions by normal python functions.
# In our case the functions for NPP and ξ have to be provided. NPP_fun will interpolate the NPP for the day in question from the data. Which we have to load. 
# We will later store these functions in  `model_specific_helpers.py` which resides in the same folder as this notebook. You will have to adapt them to your data set. 

# +
# ## Data assimilation
# Until now we have used only the initial values of the observations. 
# The next step is to decide which parameters we want to consider fixed and which to be estimated.
# This distinction helps, to keep the to create a function which only takes the estimated parameters and thus can be used by a generalized mcmc as will become clear.
#
# We can change which parameters we fix and which we estimate later or can have several approaches for the same symbolic model.
# The distinction is not model inherent but just a reflection of our choice for data assimilation.
# The more parameter values we can find out from the literature the fewer values we have to estimate. 

# +
cpa=msh.Constants(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 npp_0=dvs.npp[0], # kg/m2/day
 rh_0=svs_0.rh,    # kg/m2/day
 ra_0=svs_0.ra,    # kg/m2/day
        
#  k_C_mll=0.15/365, 
#  k_C_mwl=0.001/365,
#  k_C_mrl=0.1/365,
#  k_C_sll=0.01/365,
#  k_C_swl=0.001/365,
#  k_C_srl=0.005/365,
#  k_C_lll=0.01/365,
#  k_C_lwl=0.001/365,
#  k_C_lrl=0.005/365,
        
#  r_C_mll_2_C_mic=0.4/365, # f13_4
#  r_C_mwl_2_C_mic=0.4/365, # f13_5
#  r_C_mrl_2_C_mic=0.4/365, # f13_6
#  r_C_sll_2_C_mic=0.3/365, # f13_7
#  r_C_swl_2_C_mic=0.3/365, # f13_8
#  r_C_srl_2_C_mic=0.3/365, # f13_9
#  r_C_pass_2_C_mic=0.2/365, # f13_16
    
 r_C_mll_2_C_mic=0.0600000000000000,
 r_C_mwl_2_C_mic=0.000400000000000000,
 r_C_mrl_2_C_mic=0.0400000000000000,
 r_C_sll_2_C_mic=0.00300000000000000,
 r_C_swl_2_C_mic=0.000300000000000000,
 r_C_srl_2_C_mic=0.00150000000000000,    
 r_C_pass_2_C_mic=1.64383561643836e-6,     

#  r_C_lll_2_C_prot=0.5/365, # f14_10
#  r_C_lwl_2_C_prot=0.5/365, # f14_11
#  r_C_lrl_2_C_prot=0.5/365, # f14_12
#  r_C_lll_2_C_nonprot=0.5/365, # f15_10
#  r_C_lwl_2_C_nonprot=0.5/365, # f15_11
#  r_C_lrl_2_C_nonprot=0.5/365, # f15_12     
    
 r_C_lll_2_C_prot=0.00500000000000000,
 r_C_lll_2_C_nonprot=0.00500000000000000,
 r_C_lwl_2_C_prot=0.000500000000000000,
 r_C_lwl_2_C_nonprot=0.000500000000000000,
 r_C_lrl_2_C_prot=0.00250000000000000,
 r_C_lrl_2_C_nonprot=0.00250000000000000,
   

 number_of_months=len(svs.rh)
 #number_of_months=3840 # for testing and tuning mcmc
)


# ### Finding better start values for the data assimilation
# You don't have to do this. It's a heuristic approach to find a better starting position.
# -

from sympy import Symbol
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

# +
import numpy as np 
#svs_0=msh.Observables(*map(lambda v: v[0],svs))

X_0= np.array((
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cLitter/6,
    svs_0.cLitter/6,
    svs_0.cSoil/7,
    svs_0.cLitter/6,
    svs_0.cLitter/6,
    svs_0.cSoil/7,
    svs_0.cLitter/6,
    svs_0.cLitter/6,
    svs_0.cSoil/7,
    svs_0.cSoil/7,
    svs_0.cSoil/7,
    svs_0.cSoil/7,
    svs_0.cSoil/7,
))#.reshape(9,)

# +
# create a start parameter tuple for the mcmc. The order has to be the same as when you created the namedtupl3 
# If you don't you get a "TypeError". 
epa_0=msh.EstimatedParameters(
 beta_leaf=0.3,
 beta_wood=0.3,
 
#  r_C_leaf_2_C_mll=0.4/365,    
#  r_C_wood_2_C_mwl=0.03/365,     
#  r_C_root_2_C_mrl=0.3/365,     
#  r_C_leaf_2_C_sll=0.35/365,    
#  r_C_wood_2_C_swl=0.035/365,    
#  r_C_root_2_C_srl=0.3/365,
    
 r_C_leaf_2_C_mll=0.00583333333333333,
 r_C_leaf_2_C_sll=0.000833333333333333,
 r_C_leaf_2_C_lll=0.000833333333333333,
 r_C_wood_2_C_mwl=9.13242009132420e-6,
 r_C_wood_2_C_swl=5.47945205479452e-5,
 r_C_wood_2_C_lwl=1.82648401826484e-5,
 r_C_root_2_C_mrl=2.49066002490660e-4,
 r_C_root_2_C_srl=4.98132004981320e-4,
 r_C_root_2_C_lrl=2.49066002490660e-4,  
        
#  r_C_prot_2_C_mic=0.5/365,     
#  r_C_nonprot_2_C_mic=0.5/365,  
#  r_C_mic_2_C_prot=0.5/365,    
#  r_C_mic_2_C_nonprot=0.5/365,  
#  r_C_prot_2_C_pass=0.5/365,  
#  r_C_nonprot_2_C_pass=0.5/365,
       
 r_C_prot_2_C_mic=8.71731008717310e-5,
 r_C_nonprot_2_C_mic=0.000102117061021171,
 r_C_mic_2_C_prot=4.56621004566210e-5,
 r_C_mic_2_C_nonprot=4.56621004566210e-5,
 r_C_prot_2_C_pass=0.000149439601494396,
 r_C_nonprot_2_C_pass=0.000124533001245330, 

#  k_C_mic=0.001,            
#  k_C_protsom=0.0001,       
#  k_C_nonprotsom=0.0001,      
#  k_C_passsom=0.0001, 
    
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
        
 C_wood_0=3.739601444837676,           
 C_root_0=1.851102715197606,          
 C_mll_0=0.0038249804732625784,            
 C_mwl_0=0.0341516113683918,            
 C_sll_0=0.008196386728397948,           
 C_swl_0=0.2049096682095654,          
 C_lll_0=0.008196386728397948,            
 C_mrl_0=0.0004610467534724643,       
 C_srl_0=0.01844187013895204,          
 C_lrl_0=0.00922093506947602,           
 C_mic_0=4.418910930305989,           
 C_prot_0=1.2043499895354082,          
 C_nonprot_0=1.2043499895354082, 
)
#-

# +
import sys
sys.path.insert(0,'..') # necessary to import general_helpers
import general_helpers as gh
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2

def make_steady_state_iterator_sym(
        mvs,
        V_init,
        par_dict,
        func_dict
    ):
    B_func, u_func = gh.make_B_u_funcs_2(mvs,par_dict,func_dict)  
    def f(it,tup):
        X,_,_=tup
        b = u_func(it,X)
        B = B_func(it,X)
        return (X,b,B)
  
    return TimeStepIterator2(
        initial_values=V_init,
        f=f)
# calculate steady state
func_dict=msh.make_func_dict(svs,dvs,cpa,epa_0)
B_func, u_func = gh.make_B_u_funcs_2(mvs,par_dict,func_dict)  

# +
it_sym = make_steady_state_iterator_sym(
    mvs,
    V_init=(X_0,u_func(0,X_0),B_func(0,X_0)),
    par_dict=par_dict,
    func_dict=func_dict
)

Bs=[]
bs=[]
max_day=(cpa.number_of_months-1)*30
print(max_day,gh.day_2_month_index(max_day))

for i in range(max_day):
    res=it_sym.__next__()
    bs.append(res[1])
    Bs.append(res[2])
B_mean=np.stack(Bs).mean(axis=0)
b_mean=np.stack(bs).mean(axis=0)
B_mean,b_mean
np.linalg.inv(Bs[0])

# calculate pseudo steady state
X_ss = np.linalg.solve(B_mean, (-b_mean))

steady_state_dict={str(name): X_ss[i,0] for i,name in enumerate(mvs.get_StateVariableTuple())}
# -

#svs
#cpa
print('Steadt_state_dict:',steady_state_dict)

dvs.npp

svs_0.cVeg

# +
# deriving initail pool sizes proportional to steady state 
C_Veg_steady=steady_state_dict.get("C_leaf")+steady_state_dict.get("C_wood")+steady_state_dict.get("C_root")
C_wood_0_steady=svs_0.cVeg*steady_state_dict.get("C_wood")/C_Veg_steady
C_root_0_steady=svs_0.cVeg*steady_state_dict.get("C_root")/C_Veg_steady

C_Litter_steady=steady_state_dict.get("C_mll")+steady_state_dict.get("C_mwl")+steady_state_dict.get("C_mrl")+steady_state_dict.get("C_sll")+steady_state_dict.get("C_swl")+steady_state_dict.get("C_srl")+steady_state_dict.get("C_lll")+steady_state_dict.get("C_lwl")+steady_state_dict.get("C_lrl")
C_mll_0_steady=svs_0.cLitter*steady_state_dict.get("C_mll")/C_Litter_steady          
C_mwl_0_steady=svs_0.cLitter*steady_state_dict.get("C_mwl")/C_Litter_steady             
C_sll_0_steady=svs_0.cLitter*steady_state_dict.get("C_sll")/C_Litter_steady           
C_swl_0_steady=svs_0.cLitter*steady_state_dict.get("C_swl")/C_Litter_steady           
C_lll_0_steady=svs_0.cLitter*steady_state_dict.get("C_lll")/C_Litter_steady              
C_mrl_0_steady=svs_0.cLitter*steady_state_dict.get("C_mrl")/C_Litter_steady         
C_srl_0_steady=svs_0.cLitter*steady_state_dict.get("C_srl")/C_Litter_steady           
C_lrl_0_steady=svs_0.cLitter*steady_state_dict.get("C_lrl")/C_Litter_steady  

C_Soil_steady=steady_state_dict.get("C_mic")+steady_state_dict.get("C_prot")+steady_state_dict.get("C_nonprot")+steady_state_dict.get("C_pass")
         
C_mic_0=svs_0.cSoil*steady_state_dict.get("C_mic")/C_Soil_steady             
C_prot_0=svs_0.cSoil*steady_state_dict.get("C_prot")/C_Soil_steady           
C_nonprot_0=svs_0.cSoil*steady_state_dict.get("C_nonprot")/C_Soil_steady  
print (C_wood_0_steady)
print (C_root_0_steady)
print (C_mll_0_steady)
print (C_mwl_0_steady)
print (C_sll_0_steady)
print (C_swl_0_steady)
print (C_lll_0_steady)
print (C_mrl_0_steady)
print (C_srl_0_steady)
print (C_lrl_0_steady)
print (C_mic_0)
print (C_prot_0)
print (C_nonprot_0)

# +
# now we change epa_0 to incorporate steady state initial pool sizes
epa_1=msh.EstimatedParameters(
 beta_leaf=0.6,
 beta_wood=0.25,
 
#  r_C_leaf_2_C_mll=0.4/365,    
#  r_C_wood_2_C_mwl=0.03/365,     
#  r_C_root_2_C_mrl=0.3/365,     
#  r_C_leaf_2_C_sll=0.35/365,    
#  r_C_wood_2_C_swl=0.035/365,    
#  r_C_root_2_C_srl=0.3/365,
    
 r_C_leaf_2_C_mll=0.00583333333333333,
 r_C_leaf_2_C_sll=0.000833333333333333,
 r_C_leaf_2_C_lll=0.000833333333333333,
 r_C_wood_2_C_mwl=5.13242009132420e-6,
 r_C_wood_2_C_swl=2.47945205479452e-5,
 r_C_wood_2_C_lwl=1.82648401826484e-5,
 r_C_root_2_C_mrl=2.49066002490660e-5,
 r_C_root_2_C_srl=4.98132004981320e-5,
 r_C_root_2_C_lrl=2.49066002490660e-5,  
        
#  r_C_prot_2_C_mic=0.5/365,     
#  r_C_nonprot_2_C_mic=0.5/365,  
#  r_C_mic_2_C_prot=0.5/365,    
#  r_C_mic_2_C_nonprot=0.5/365,  
#  r_C_prot_2_C_pass=0.5/365,  
#  r_C_nonprot_2_C_pass=0.5/365,
       
 r_C_prot_2_C_mic=8.71731008717310e-5,
 r_C_nonprot_2_C_mic=0.000102117061021171,
 r_C_mic_2_C_prot=4.56621004566210e-5,
 r_C_mic_2_C_nonprot=4.56621004566210e-5,
 r_C_prot_2_C_pass=0.000149439601494396,
 r_C_nonprot_2_C_pass=0.000124533001245330, 

#  k_C_mic=0.001,            
#  k_C_protsom=0.0001,       
#  k_C_nonprotsom=0.0001,      
#  k_C_passsom=0.0001, 
    
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
#-

# +
# now test it

import matplotlib.pyplot as plt
from general_helpers import plot_solutions

param2res_sym = msh.make_param2res_sym(mvs,cpa,dvs)

#print(type(param2res_sym))

obs_0 = param2res_sym(epa_1)
#obs=np.column_stack([ np.array(v) for v in svs])
#obs=np.column_stack((np.repeat(svs.cVeg, 12),np.repeat(svs.cLitter, 12),np.repeat(svs.cSoil, 12),svs.rh,svs.ra))
#xs.shape
#obs_0

# +
# obs=obs[0:cpa.number_of_months,:] #cut 
# obs[:,3:4]=obs[:,3:4]
# n=cpa.number_of_months

# convert to yearly output if necessary
# obs_yr=np.zeros(int(cpa.number_of_months/12)*obs.shape[1]).reshape([int(cpa.number_of_months/12),obs.shape[1]])  
# for i in range(obs.shape[1]):
#    obs_yr[:,i]=gh.monthly_to_yearly(obs[:,i])
# obs=obs_yr
# n=int(cpa.number_of_months/12)

out_simu_d=obs_0._asdict()
obs_d=svs._asdict()

print("Forward run with initial parameters (blue) vs TRENDY output (orange)")
obs_d['cVeg'],svs.cVeg

fig = plt.figure(figsize=(10,50))
axs=fig.subplots(5,1)

for ind,f in enumerate(msh.Observables._fields):
    val_sim=out_simu_d[f]
    val_obs=obs_d[f]
    axs[ind].plot(range(len(val_obs)),val_obs,label=f+"_obs")
    axs[ind].plot(range(len(val_sim)),val_sim,label=f+"_sim")    
    axs[ind].legend()

# fig = plt.figure(figsize=(12, 4), dpi=80)
# plot_solutions(
#         fig,
#         times=range(n),
#         var_names=msh.Observables._fields,
#         tup=(xs,obs)
#         #tup=(obs,)
# )
fig.savefig('solutions.pdf')
# -

    epa_min=msh.EstimatedParameters(
         beta_leaf=0.01,
         beta_wood=0.01,

         r_C_leaf_2_C_mll=epa_1.r_C_leaf_2_C_mll/100,
         r_C_leaf_2_C_sll=epa_1.r_C_leaf_2_C_sll/100,
         r_C_leaf_2_C_lll=epa_1.r_C_leaf_2_C_lll/100,
         r_C_wood_2_C_mwl=epa_1.r_C_wood_2_C_mwl/100,
         r_C_wood_2_C_swl=epa_1.r_C_wood_2_C_swl/100,
         r_C_wood_2_C_lwl=epa_1.r_C_wood_2_C_lwl/100,
         r_C_root_2_C_mrl=epa_1.r_C_root_2_C_mrl/100,
         r_C_root_2_C_srl=epa_1.r_C_root_2_C_srl/100,
         r_C_root_2_C_lrl=epa_1.r_C_root_2_C_lrl/100,  
                            
         r_C_prot_2_C_mic=epa_1.r_C_prot_2_C_mic/100,
         r_C_nonprot_2_C_mic=epa_1.r_C_nonprot_2_C_mic/100,
         r_C_mic_2_C_prot=epa_1.r_C_mic_2_C_prot/100,
         r_C_mic_2_C_nonprot=epa_1.r_C_mic_2_C_nonprot/100,
         r_C_prot_2_C_pass=epa_1.r_C_prot_2_C_pass/100,
         r_C_nonprot_2_C_pass=epa_1.r_C_nonprot_2_C_pass/100, 
          
         r_C_mll_rh=epa_1.r_C_mll_rh/100,
         r_C_mwl_rh=epa_1.r_C_mwl_rh/100,
         r_C_mrl_rh=epa_1.r_C_mrl_rh/100,
         r_C_sll_rh=epa_1.r_C_sll_rh/100,
         r_C_swl_rh=epa_1.r_C_swl_rh/100,
         r_C_srl_rh=epa_1.r_C_srl_rh/100,        
         r_C_lll_rh=epa_1.r_C_lll_rh/100,
         r_C_lwl_rh=epa_1.r_C_lwl_rh/100,
         r_C_lrl_rh=epa_1.r_C_lrl_rh/100,            
         r_C_mic_rh=epa_1.r_C_mic_rh/100,           
         r_C_prot_rh=epa_1.r_C_prot_rh/100,
         r_C_nonprot_rh=epa_1.r_C_nonprot_rh/100,        
         r_C_pass_rh=epa_1.r_C_pass_rh/100,   
            
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
         r_C_leaf_2_C_mll=epa_1.r_C_leaf_2_C_mll*100,
         r_C_leaf_2_C_sll=epa_1.r_C_leaf_2_C_sll*100,
         r_C_leaf_2_C_lll=epa_1.r_C_leaf_2_C_lll*100,
         r_C_wood_2_C_mwl=epa_1.r_C_wood_2_C_mwl*100,
         r_C_wood_2_C_swl=epa_1.r_C_wood_2_C_swl*100,
         r_C_wood_2_C_lwl=epa_1.r_C_wood_2_C_lwl*100,
         r_C_root_2_C_mrl=epa_1.r_C_root_2_C_mrl*100,
         r_C_root_2_C_srl=epa_1.r_C_root_2_C_srl*100,
         r_C_root_2_C_lrl=epa_1.r_C_root_2_C_lrl*100,                
              
         r_C_prot_2_C_mic=epa_1.r_C_prot_2_C_mic*100,
         r_C_nonprot_2_C_mic=epa_1.r_C_nonprot_2_C_mic*100,
         r_C_mic_2_C_prot=epa_1.r_C_mic_2_C_prot*100,
         r_C_mic_2_C_nonprot=epa_1.r_C_mic_2_C_nonprot*100,
         r_C_prot_2_C_pass=epa_1.r_C_prot_2_C_pass*100,
         r_C_nonprot_2_C_pass=epa_1.r_C_nonprot_2_C_pass*100, 
          
         r_C_mll_rh=epa_1.r_C_mll_rh*100,
         r_C_mwl_rh=epa_1.r_C_mwl_rh*100,
         r_C_mrl_rh=epa_1.r_C_mrl_rh*100,
         r_C_sll_rh=epa_1.r_C_sll_rh*100,
         r_C_swl_rh=epa_1.r_C_swl_rh*100,
         r_C_srl_rh=epa_1.r_C_srl_rh*100,        
         r_C_lll_rh=epa_1.r_C_lll_rh*100,
         r_C_lwl_rh=epa_1.r_C_lwl_rh*100,
         r_C_lrl_rh=epa_1.r_C_lrl_rh*100,            
         r_C_mic_rh=epa_1.r_C_mic_rh*100,           
         r_C_prot_rh=epa_1.r_C_prot_rh*100,
         r_C_nonprot_rh=epa_1.r_C_nonprot_rh*100,        
         r_C_pass_rh=epa_1.r_C_pass_rh*100,   
            
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

#import test_helpers as th
#ta=th.make_test_args(conf_dict,msh,mvs)
costfunction=msh.make_weighted_cost_func(svs)
print(costfunction(obs_0))
print(np.array(epa_1)-np.array(epa_min))
print(np.array(epa_max)-np.array(epa_1))

# ### mcmc to optimize parameters 
#

# +

from general_helpers import autostep_mcmc, make_param_filter_func, make_feng_cost_func

isQualified = make_param_filter_func(epa_max, epa_min)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa_1,
    filter_func=isQualified,
    param2res=param2res,
    #costfunction=make_feng_cost_func(obs),
    costfunction=msh.make_weighted_cost_func(svs),
    nsimu=200, # for testing and tuning mcmc
    #nsimu=20000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=15,   # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=1,   # default value | increase value to reduce initial step size
    K=2 # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

# +
# optimized parameter set (lowest cost function)
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields),1),axis=1)
epa_opt=msh.EstimatedParameters(*par_opt)
mod_opt = param2res(epa_opt)  
DA = mod_opt._asdict()
obs_d=svs._asdict()
print("Forward run with optimized parameters (orange) vs TRENDY output (blue)")
fig = plt.figure(figsize=(4, 10), dpi=80)
axs_DA = fig.subplots(5,1)

for a, b in enumerate(msh.Observables._fields):
    val_DA=DA[b]
    val_obs=obs_d[b]
    axs_DA[a].plot(range(len(val_obs)),val_obs, label=b+"_obs")
    axs_DA[a].plot(range(len(val_DA)), val_DA, label=b+"_DA")    
    axs_DA[a].legend()

fig.savefig('solution_DA.pdf')

fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

import pandas as pd
pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('IBIS_da_aa.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('IBIS_da_j_aa.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('IBIS_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('IBIS_optimized_solutions.csv'), sep=',')
# -

print("Optimized parameters: ", epa_opt)


# +
# import model_specific_helpers_2 as msh
# # Traceability Analysis
# it_sym_trace = msh.make_traceability_iterator(mvs,dvs,cpa,epa_opt)
# ns= 1500 #10*360
# StartVectorTrace=gh.make_StartVectorTrace(mvs)
# nv=len(StartVectorTrace._fields)
# res_trace= np.zeros((ns,nv))
# for i in range(ns):
#     res_trace[i,:]=it_sym_trace.__next__().reshape(nv)
# #res_trace

# import matplotlib.pyplot as plt
# n=len(mvs.get_StateVariableTuple())
# fig=plt.figure(figsize=(20,(n+1)*10), dpi=80)
# axs=fig.subplots(n+1,2)
# days=list(range(ns))


# for i in range(n):
    
#     ax = axs[i,0]
#     #  the solution
#     pos=i
#     ax.plot(
#         days,
#         res_trace[:,i],
#         label=StartVectorTrace._fields[pos],
#         color='blue'
#     )
#     # X_p
#     pos=i+n
#     ax.plot(
#         days,
#         res_trace[:,pos],
#         label=StartVectorTrace._fields[pos],
#         color='red'
#     )
#     # X_c
#     pos=i+2*n
#     ax.plot(
#         days,
#         res_trace[:,pos],
#         label=StartVectorTrace._fields[pos],
#         color='yellow'
#     )
#     ax.legend()
    
#     ax = axs[i,1]
#     # RT
#     pos=i+3*n
#     ax.plot(
#         days,
#         res_trace[:,pos],
#         label=StartVectorTrace._fields[pos],
#         color='black'
#     )
#     ax.legend()
    
# axs[n,0].plot(
#     days,
#     [msh.make_npp_func(dvs)(d) for d in days],
#     label='NPP',
#     color='green'
# )
# axs[n,0].legend()

# fig.savefig('solution_Traceability.pdf')
