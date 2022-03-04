# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# This illustrative notebook shows how to use the factored out model description in ```source.py``` the model specific functions in ```model_specific_helpers.py``` and the general infrastructure in ```general_helpers.py``` to inspect the model.
#
# This is the next step after the ```createModel.py``` notebook where everything was defined in one place, the notebook.
# In order to be able to use the symbolic descriptions and functions from more than one notebook  
# we now disassemble it, moving the functions into the seperate file ```model_specific_helpers.py``` and the model description to ```source.py```.
#
# This also enables us to test model descriptions and model specific helper functions  in a centralized way.

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
import model_specific_helpers as msh
from general_helpers import day_2_month_index, make_B_u_funcs_2, day_2_month_index
# -

# The last statement in the code imports the variable `mvs` from the module source (defined in ```source.py``` which is 
# an instance of CMTVS which stands for `C`onnected`M`ulti`T`ype`V`ariable`S`et".
# It contains information in two forms. 

mvs.get_StateVariableTuple()

mvs.get_CompartmentalMatrix()

mvs.get_InputTuple()

# we can also print the whole mass balance equation
dh.mass_balance_equation(mvs)

# we can also plot a picture
h.compartmental_graph(mvs)

# #### Intermediate summary:
# We have achieved the symbolic formulation of the model. We can use it to check the structure and compute derived diagnostic variables including the matrices used for the traceability analysis, but up to now only in symbolic form.
#
# ## Connecting symbolic description and data
# The next goal is to connect the symbolic formulation to the data.
# Since the data comes in several shapes this involves several steps. 
# We also want to be able to make this step portable across different models and computers.
# The endproduct is a collection of models that everybody can run who installes the package and executes the code we provide
#
# We will have to: 
# 1. Find as many model parameters as possible in the model description (in the literature or in communication with the modeling group) so that we do not have to estimate them from the model output. 
# 1. provide code to download the output for your model.
# 1. implement functions for the drivers (using the data)
# 1. run the model forward with a possible set of parameters.
# 1. infer unknown parameters by data assimilation.
#
# ### downloading the data
# This assumes that you have **created a small site specific config file**  
# This file specifies:
# - a username and password to download the data 
# - the location where you want to download the data to 
#   which will differ depending on the machine you are using (your laptop or a supercomputer) and     also accross users. You will have to have one everywhere you want to work with the model.
# Here comes a template from my laptop (content of file `../config.json`):
# `{"username": "trendy-v9", "password": "gcb-2020", "dataPath": "/home/data/VISIT_data_CMIP6"}`
#
#
# #### use a small model specific function to download the data
# This function is stored in the file `model_specific_helpers.py` 
# This function **saves us a lot of time** when we want to reproduce your results and run your model since finding the correct dataset can be a time consuming task.
#

# +
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

#msh.download_my_TRENDY_output(conf_dict)
# -

#     # Read NetCDF data  ******************************************************************************************************************************
svs,dvs=msh.get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))


svs,dvs

# ### Forward run
# The next goal is to run the model forward with a given set of parameters.
# So we need:
# 1. values for all the parameters 
# 1. implementations for the symbolic functions 
# 1. start values for all the pool contents

# In this example we have the initial values for the elements of the K and A matrices ($k_i$ and $f_{i,j}$ ) but we want the values for the individual flux rates. 
# So we first create the symbolic $k_{i}$ and $f_{i,j}$ which gives us an alternative description 
# of the product $M = K A=K_{sym} A_{sym}$ from this matrix we can compute the outfluxes and internal fluxes. If we assume that the same fluxes could be expressed as $flux=fluxrate * donorpoolcontent$ we can compute the fluxrates 
#
# **The next few steps until the next heading are hopefully not necessarry for your model since you should find parameters directly**
# Since we only had parameters for a different parameterisation, we have to convert  
#

# ### Providing dictionaries  for the parameters and functions.

# be able to refer to the symbols and functions in the symbolic description
# we recreate them here.
BI=mvs.get_BibInfo()
for k in BI.sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)
for k in BI.func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

# To be able to run the model forward we not only have to replace parameter symbols by values but symbolic functions by normal python functions.
# In our case the functions for $NPP$ and $\xi$ have to be provided. NPP_fun will interpolate the NPP for the day in question from the data. Which we have to load. 
# We will later store these functions in  `model_specific_helpers.py` which resides in the same folder as this notebook. You will have to adapt them to your data set. 

# +
par_dict={
 beta_leaf: 0.6,
 beta_wood: 0.25,
 T_0: 2,
 E: 4,
 KM: 10,
 r_C_leaf_litter_rh: 0.000415110004151100,
 r_C_wood_litter_rh: 0.000124533001245330,
 r_C_root_litter_rh: 0.000122042341220423,
 r_C_soil_fast_rh: 0.000152207001522070,
 r_C_soil_slow_rh: 2.73972602739726e-5,
 r_C_soil_passive_rh: 7.82778864970646e-6,
 r_C_leaf_2_C_leaf_litter: 0.00833333333333333,
 r_C_wood_2_C_wood_litter: 9.13242009132420e-5,
 r_C_root_2_C_root_litter: 0.000124533001245330,
 r_C_leaf_litter_2_C_soil_fast: 0.000340390203403902,
 r_C_leaf_litter_2_C_soil_slow: 5.81154005811540e-5,
 r_C_leaf_litter_2_C_soil_passive: 1.66044001660440e-5,
 r_C_wood_litter_2_C_soil_fast: 7.47198007471980e-5,
 r_C_wood_litter_2_C_soil_slow: 2.98879202988792e-5,
 r_C_wood_litter_2_C_soil_passive: 1.99252801992528e-5,
 r_C_root_litter_2_C_soil_fast: 7.47198007471980e-5,
 r_C_root_litter_2_C_soil_slow: 3.48692403486924e-5,
 r_C_root_litter_2_C_soil_passive: 1.74346201743462e-5
}


# import from modelspecific helpers
func_dict={
    'NPP':msh.make_npp_func(dvs,svs),
     'xi':msh.make_xi_func(dvs,svs)
}


# check the numeric functions for B and u

# for the next line to work the 
# two dictionaries par_dict and func_dict have to be complete.
# In sympy terms this means that the expressions for 
B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
# we create a numeric startvector for the compartmental system
# 
svs_0=msh.Observables(*map(lambda v: v[0],svs))

X_0= np.array((
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
))#.reshape(9,)
# in general B and u are nonlinear and can depend on X, thats why we have to test them with index and X arguments
u_func(0,X_0),B_func(0,X_0)
# -

V_init= msh.StartVector(
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
V_init.__getattribute__("C_leaf")

V_init

np.array(V_init).shape

# We now create another iterator that does run with a different timestep and check that this does not change the results.
# If we can increase the timestep the data assimilation will be much faster.
#

# +

# test the different iterators
    
delta_t_1=1 #daily iterator
it_sym = msh.make_iterator_sym(
    mvs,
    V_init=V_init,
    par_dict=par_dict,
    func_dict=func_dict,
    delta_t_val=delta_t_1
)
delta_t_2=30 #30_day iterator
it_sym_2 = msh.make_iterator_sym(
    mvs,
    V_init=V_init,
    par_dict=par_dict,
    func_dict=func_dict,
    delta_t_val=delta_t_2
)
# we will run the dayly for 15 steps
ns=10000
res= np.zeros((ns,len(V_init)))
times = np.arange(0,ns)
res[0,:]=V_init 
for i in range(1,ns-1):
    res[i,:]=it_sym.__next__().reshape(len(V_init),)

times_2= np.arange(0,ns,delta_t_2)
res_2= np.zeros((len(times_2),len(V_init)))
res_2[0,:]=V_init 
for i in range(1,len(times_2)-1):
    res_2[i,:]=it_sym_2.__next__().reshape(len(V_init),)
# -

times


times_2

res[:10,10]

from matplotlib import pyplot as plt
fig = plt.figure(figsize=(5,50))
n_coord=len(V_init)
axs=fig.subplots(n_coord,1)
for coord in range(n_coord):
    axs[coord].plot(times,res[:,coord])
    axs[coord].plot(times_2,res_2[:,coord])

# ## Data assimilation
# Until now we have used only the initial values of the observations. 
# The next step is to decide which parameters we want to consider fixed and which to be estimated.
# This distinction helps, to keep the to create a function which only takes the estimated parameters and thus can be used by a generalized mcmc as will become clear.
#
# We can change which parameters we fix and which we estimate later or can have several approaches for the same symbolic model.
# The distinction is not model inherent but just a reflection of our choice for data assimilation.
# The more parameter values we can find out from the literature the fewer values we have to estimate.  

msh.EstimatedParameters._fields



cpa = msh.UnEstimatedParameters(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 gpp_0=dvs.gpp[0] * 86400,   # kg/m2/s kg/m2/day
 rh_0=svs_0.rh * 86400,   # kg/m2/s kg/m2/day
 ra_0=svs_0.ra * 86400,   # kg/m2/s kg/m2/day
 r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
 r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
 #number_of_months=len(svs.rh)
 number_of_months=120 # for testing and tuning mcmc
)

# +
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
def make_steady_state_iterator_sym(
        mvs,
        V_init,
        par_dict,
        func_dict
    ):
    B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
    def f(it,tup):
        X,_,_=tup
        b = u_func(it,X)
        B = B_func(it,X)
        return (X,b,B)
  
    return TimeStepIterator2(
        initial_values=V_init,
        f=f)
# calculate steady state
  
it_sym = make_steady_state_iterator_sym(
    mvs,
    V_init=(X_0,u_func(0,X_0),B_func(0,X_0)),
    par_dict=par_dict,
    func_dict=func_dict
)
Bs=[]
bs=[]
for i in range(cpa.number_of_months*30):
    bs.append(it_sym.__next__()[1])
    Bs.append(it_sym.__next__()[2])
B_mean=np.stack(Bs).mean(axis=0)
b_mean=np.stack(bs).mean(axis=0)
B_mean,b_mean
np.linalg.inv(Bs[0])


# +
# calculate pseudo steady state
X_ss = np.linalg.solve(B_mean, (-b_mean))

steady_state_dict={str(name): X_ss[i,0] for i,name in enumerate(mvs.get_StateVariableTuple())}
# -

# create a start parameter tuple for the mcmc. 
epa_0=msh.EstimatedParameters(
    beta_leaf=0.6,
    beta_wood=0.25,
    T_0=2,
    E=4,
    KM=10,
    r_C_leaf_litter_rh=0.000415110004151100,
    r_C_wood_litter_rh=0.000124533001245330,
    r_C_root_litter_rh=0.000122042341220423,
    r_C_soil_fast_rh=0.000152207001522070,
    r_C_soil_slow_rh=2.73972602739726e-5,
    r_C_soil_passive_rh=7.82778864970646e-6,
    r_C_leaf_2_C_leaf_litter=0.00833333333333333,
    r_C_wood_2_C_wood_litter=9.13242009132420e-5,
    r_C_root_2_C_root_litter=0.000124533001245330,
    r_C_leaf_litter_2_C_soil_fast=0.000340390203403902,
    r_C_leaf_litter_2_C_soil_slow=5.81154005811540e-5,
    r_C_leaf_litter_2_C_soil_passive=1.66044001660440e-5,
    r_C_wood_litter_2_C_soil_fast=7.47198007471980e-5,
    r_C_wood_litter_2_C_soil_slow=2.98879202988792e-5,
    r_C_wood_litter_2_C_soil_passive=1.99252801992528e-5,
    r_C_root_litter_2_C_soil_fast=7.47198007471980e-5,
    r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
    r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
    C_leaf_0=steady_state_dict["C_leaf"],
    C_wood_0=steady_state_dict["C_wood"],
    C_leaf_litter_0=steady_state_dict["C_leaf_litter"],
    C_wood_litter_0=steady_state_dict["C_wood_litter"],
    C_soil_fast_0=steady_state_dict["C_soil_fast"],
    C_soil_slow_0=steady_state_dict["C_soil_slow"]
)

# The function `param2res` (which will be used by a general model independent mcmc) only takes the estimated parameters as arguments and produce data in the same shape as the observations.
# We will taylor make it by another function `make_param2res` which depends on the parameters that we decide to fix.
# This function is the main effort to make the data assimilation possible. **Although we give the full definition here it's suggested that you incrementally build and check the code inside it before you make it a function that returns a function...** 
# - You can start with sets of  `UnEstimatedParameters` and `EstimatedParameter` (e.g.  `epa0` from the test above) and try to produce a set of observations by running the model. 
# - then produce `param2res` automatically by code
# - them wrap this code into `make_param2res` so that you can change the unestimated parameters easily.
#
# This function will be later moved to a model specific helper file, and loaded from there by the data assimilation functions.

# +
# now test it 
import matplotlib.pyplot as plt
from general_helpers import plot_solutions

param2res_sym = msh.make_param2res_sym(cpa,dvs,svs)
xs= param2res_sym(epa_0)
obs=np.column_stack([ np.array(v) for v in svs])

obs=obs[0:cpa.number_of_months,:] #cut 
#obs[:,3:4]=obs[:,3:4]
n=cpa.number_of_months

# convert to yearly output if necessary
#obs_yr=np.zeros(int(cpa.number_of_months/12)*obs.shape[1]).reshape([int(cpa.number_of_months/12),obs.shape[1]])  
#for i in range(obs.shape[1]):
#    obs_yr[:,i]=monthly_to_yearly(obs[:,i])
#obs=obs_yr
#n=int(cpa.number_of_months/12)

print("Forward run with initial parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(12, 4), dpi=80)
plot_solutions(
        fig,
        times=range(n),
        var_names=msh.Observables._fields,
        tup=(xs,obs)
        #tup=(obs,)
)
fig.savefig('solutions.pdf')

# +
epa_min=msh.EstimatedParameters(
    beta_leaf=0,
    beta_wood=0,
    T_0=-20,
    E=.1,
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
    T_0=10,
    E=100,
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
# -

# ### mcmc to optimize parameters 
#

np.array(epa_max)

# +
from general_helpers import autostep_mcmc, make_param_filter_func, make_feng_cost_func

isQualified = make_param_filter_func(epa_max, epa_min)
param2res = msh.make_param2res_sym(cpa,dvs,svs)

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa_0,
    filter_func=isQualified,
    param2res=param2res,
    costfunction=make_feng_cost_func(obs),
    nsimu=2000, # for testing and tuning mcmc
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

print("Forward run with optimized parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(12, 4), dpi=80)
plot_solutions(
        fig,
        #times=range(cpa.number_of_months),
        times=range(int(cpa.number_of_months)), # for yearly output
        var_names=msh.Observables._fields,
        tup=(mod_opt,obs)
)

fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

import pandas as pd
pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('visit_da_aa.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('visit_da_j_aa.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('visit_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('visit_optimized_solutions.csv'), sep=',')
# -

print("Optimized parameters: ", epa_opt)
par_dict_opt={
    beta_leaf: epa_opt.beta_leaf,
    beta_wood: epa_opt.beta_wood,
    T_0: epa_opt.T_0,
    E: epa_opt.E,
    KM: epa_opt.KM,
    #r_C_leaf_rh: 0,
    #r_C_wood_rh: 0,
    #r_C_root_rh: 0,
    r_C_leaf_litter_rh: epa_opt.r_C_leaf_litter_rh,
    r_C_wood_litter_rh: epa_opt.r_C_wood_litter_rh,
    r_C_root_litter_rh: epa_opt.r_C_root_litter_rh,
    r_C_soil_fast_rh: epa_opt.r_C_soil_fast_rh,
    r_C_soil_slow_rh: epa_opt.r_C_soil_slow_rh,
    r_C_soil_passive_rh: epa_opt.r_C_soil_passive_rh,
    r_C_leaf_2_C_leaf_litter: epa_opt.r_C_leaf_2_C_leaf_litter,
    r_C_wood_2_C_wood_litter: epa_opt.r_C_wood_2_C_wood_litter,
    r_C_root_2_C_root_litter: epa_opt.r_C_root_2_C_root_litter,
    r_C_leaf_litter_2_C_soil_fast: epa_opt.r_C_leaf_litter_2_C_soil_fast,
    r_C_leaf_litter_2_C_soil_slow: epa_opt.r_C_leaf_litter_2_C_soil_slow,
    r_C_leaf_litter_2_C_soil_passive: epa_opt.r_C_leaf_litter_2_C_soil_passive,
    r_C_wood_litter_2_C_soil_fast: epa_opt.r_C_wood_litter_2_C_soil_fast,
    r_C_wood_litter_2_C_soil_slow: epa_opt.r_C_wood_litter_2_C_soil_slow,
    r_C_wood_litter_2_C_soil_passive: epa_opt.r_C_wood_litter_2_C_soil_passive,
    r_C_root_litter_2_C_soil_fast: epa_opt.r_C_root_litter_2_C_soil_fast,
    r_C_root_litter_2_C_soil_slow: epa_opt.r_C_root_litter_2_C_soil_slow,
    r_C_root_litter_2_C_soil_passive: epa_opt.r_C_root_litter_2_C_soil_passive 
}
print("Optimized parameters dictionary: ", par_dict_opt)

# ### Traceability analysis  
#
# #### Outline
# The traceability analysis defines several diagnostic variables using as much algebraic structure of the mass balance equation as is available.
# Not all diagnostic variables are possible for all compartmental models. 
#
# We chose here to introduce the diagnostic variables not all at once but rather in the order of decreasing generality.
#
# The first diagnostic variables are available for all compartmental models and need no additional assumptions. 
# In the later parts of this section we then assume to be able to identify more and more specific terms in the mass balance equation and use those to derive and trace ever more specific diagnostics.
# Thus the very first part is valid for all models but how many of the later parts are applicable to a specific model  depends on how much we know about it.  
#
#
# #### Derivation of the matrix decomposition 
# Compartmental models (well mixed mass balanced) can be written in as an ordinary differential equation in matrix form that relates the momentary value of the (time) derivative $\frac{d X}{d t}$ of an yet unknown function $X$ to the momentary value of $X$ itself.   
# $$
# \frac{d X}{d t}= M(X,t) X + I(X,t) \quad (1)   
# $$ 
# where $X$ is the statevector representing the pool contents, $M$ the "Compartmental matrix" and $I$ the input vector.
# Together with a startvalue $X_0$ it constitutes an "initial value problem" (ivp) which can be solved numerically by moving step by step forward in time.
#
# Note: 
#
# It is mathematical standard notation to use $X$ in the *formulation* of the ivp (representing the momentary value) althoug *after we have solved it* the solution is expressed as function of time $X(t)$. This avoids confusion since everything appering with arguments is recognizable as explicitly calculable *before* we have solved the ivp.
#
# The system "nonautonomous" (if they depend on time $t$) and "nonlinear" if the dependent on $X$.
# It is always possible to factorize $M(X,t)$ into a product $M=A(X,t) K(X,t)$ where $K$ is a  diagonal matrix.
# and $I=B(t)*u(t)$ where $u$ is a scalar.
# Using these we arrive at 
# $$
# \frac{d X}{d t}= A(X,t) K(X,t) X + B(X,t) u(X,t)  
# $$
# ##### Linearity assumption
# If we assume the model to be linear and nonautonomous the dependency on $X$ vanishes and we have
#
# $$
# \frac{d X}{d t}= A(t) K(t) X + B(t) u(t) . 
# $$
#
# ##### Factorizability  assumption
# Although this is not possible in general in many published models the nonautonous part  can be further localized into a diagonal matrix $\xi(t)$ so that we can achieve constant $A$ and $K$ which allows more specific interpretation.
#
# $$
# \frac{d X}{d t}= A \xi(t) K X + B(t)u(t)
# $$
#
# ##### Factorizability of $\xi$ assumption 
# In some cases we can resolve $\xi$ further.
# $$
# \frac{d X}{d t}= A \xi_temp(t) \xi_mois(t) K X + B(t)u(t)
# $$
#
# #### Definition of diagnostic variables
#
# ##### Storage capacity $X_c$ and storage potential $X_p$
# These variables can be defined for any compartmental system and do not require either linearity nor factorizability. 
# We can rearrange eq. $(1)$ and give names to the two summands. 
# $$
# X = M^{-1}(X,t) \left( \frac{d X}{d t}-I(X,t) \right) \\ 
#   = \underbrace{M^{-1}(X,t) \frac{d X}{d t}}_{X_c} - \underbrace{M^{-1}(X,t)I(X,t)}_{X_p} \\
#   = X_c - X_p
# $$
# Note:
# This is not to be read as a recipe to compute $X$.
# The equation becomes a bit clearer if we adapt the nomenclature to express that we *have solved the ivp* and know its solution $X(t)$  
# <!-- and therefore also  the derivative $\frac{d X}{d t}=M(X(t),t) X(t) + I(X(t),t)=\prime{X}(t)$ -->
# By substituting the solution $X(t)$ we get the recipes to compute:
# $$
# X_p(t) = M^{-1}(X(t),t)I(X(t),t)  \\ 
# X_c(t) = X(t)-X_p(t) \\ 
# $$
# we see that all the ingredients become explicit functions of time.   
# Since all values are taken at the same time $t$ we can drop the time dependence
# in the notation and write an equation we can use in the iterator.
# $$
# X_p = M^{-1}I(X,t)  \\ 
# X_c = X + X_p \\ 
# $$
#
# ##### Residence time
# The influx $I$ can always be written as $I=b u$ where the scalar $u=\sum_{k=1\dots n} I_k$  and the dimensionless vector $b=I/u$ where $\sum_{k=1\dots n} b_k =1$.
# Assumimg that the pool contents (the components of $X$)  have dimension $mass$ we can infer from eq. (1) that $M$ has dimension $\frac{1}{time}$.
# The components of the (inverse) matrix $M^{-1}$ have therefore dimension $time$. Accordingly the product $RT= M^{-1} b$ is a vector of the same shape as $X$  whose components have dimesion $time$.
# In the context of the Traceability Framework $RT$ is therefore called *residence time*.
#
# Notes on nomenclature: 
# 1. The term *residence time* is not universally used with the same connotation outside the context of the *Traceability Analysis*.
#
# 1. It is not *the time of residence* of the particles in the system for the following reasons:
#     1. In well mixed systems particles can reside in a pool for different times from zero to infinity.
#     1. You could compute the mean of these times over all particles exiting a pool, but even then the result is in general not equal to the above mentioned $rt$.
#     1. The mean residence time would only coincide with the definition above if the system was in equilibrium (which it clearly is not as e.g $NPP(t)$ shows.)
#     1. The origin of the term is probably most easily understood as the generalization of a one dimensional rate equation $\frac{d}{dt} x = m x + u$ 
#        If $r$ and $u$ are constant then the mean residence time is $rt= m^{-1}$. If we start with the rate as property of the model the *residence time* 
#        can be defined as the inverse of this rate. The above definition is the generalization of this simple relationship to matrices and vectors.
#        The matrix $M^{-1}$ takes the role of the number $\frac{1}{m}$ . In the context of the *Traceability Analysis* $M^{-1}$ is called *Chasing Time*. 
#

# +
from collections import namedtuple
# lets build an iterator to trace  X, X_c and X_p
# we will extend it further later (for more )
# but we can also make it faster because we are not interested in
# the respiration this time

# build a new template for the StartVector 
# at the moment 3 times the length of the vector of pool contents.
# we will add more components later
svt=mvs.get_StateVariableTuple()

StartVectorTrace=namedtuple(
    "StartVectorTrace",
    [str(v) for v in svt]+
    [str(v)+"_p" for v in svt]+
    [str(v)+"_c" for v in svt]+
    [str(v)+"_RT" for v in svt]
)


# -

# now build the iterator to deal with such vectors
def make_daily_iterator_sym_trace(
        mvs,
        V_init: StartVectorTrace,
        par_dict,
        func_dict
    ):
    B_func, I_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
    V_arr=np.array(V_init).reshape(-1,1) #reshaping for matmul which expects a one column vector (nr,1) 
    
    n=len(mvs.get_StateVariableTuple())
    def f(it,V):
        #the pools are the first n values
        X = V[0:n] 
        I = I_func(it,X) 
        # we decompose I
        u=I.sum()
        b=I/u
        B = B_func(it,X)
        B_inf = np.linalg.inv(B)
        X_new = X + I + B @ X
        X_p = B_inf @ I
        X_c = X_new+X_p
        RT = B_inf @ b
        V_new = np.concatenate(
            (
                X_new.reshape(n,1),
                X_p.reshape(n,1),
                X_c.reshape(n,1),
                RT.reshape(n,1),
            ),
            axis=0
        )
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


# +
# test the new iterator

# first build a new s
# actually we realy can choose only the startvalues for the pools
# but the iterator now produces a 3 times longer vector 
X_0= np.array((
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
)).reshape(9,1)
# we make the X_p and X_c  parts compatible with the ones computed by the iterator 
# for following timesteps (to keep it) 
# As you can see in the definition of the iterator these values have no impact on further results  
I = u_func(0,X_0)
u=I.sum()
b=I/u
B = B_func(0,X_0)
B_inf = np.linalg.inv(B)
X_p_0 = B_inf@I
X_c_0 = X_0+X_p_0
RT_0 = B_inf@b
# combine the three 
#here we rely on order to be consistent 
#(although we could use also the names of the namedtuple)
V_arr=np.concatenate((X_0,X_p_0,X_c_0,RT_0),axis=0 )
V_init=StartVectorTrace(*V_arr)

# -

it_sym_trace = make_daily_iterator_sym_trace(
    mvs,
    V_init=V_init,
    par_dict=par_dict,
    func_dict=func_dict
)
ns=1500
nv=len(V_init)
res_trace= np.zeros((ns,nv))
for i in range(ns):
    res_trace[i,:]=it_sym_trace.__next__().reshape(len(V_init),)
res_trace

# +
import matplotlib.pyplot as plt
n=len(mvs.get_StateVariableTuple())
fig=plt.figure(figsize=(20,(n+1)*10), dpi=80)
axs=fig.subplots(n+1,2)
days=list(range(ns))


for i in range(n):
    
    ax = axs[i,0]
    #  the solution
    pos=i
    ax.plot(
        days,
        res_trace[:,i],
        label=StartVectorTrace._fields[pos],
        color='blue'
    )
    # X_p
    pos=i+n
    ax.plot(
        days,
        res_trace[:,pos],
        label=StartVectorTrace._fields[pos],
        color='red'
    )
    # X_c
    pos=i+2*n
    ax.plot(
        days,
        res_trace[:,pos],
        label=StartVectorTrace._fields[pos],
        color='yellow'
    )
    ax.legend()
    
    ax = axs[i,1]
    # RT
    pos=i+3*n
    ax.plot(
        days,
        res_trace[:,pos],
        label=StartVectorTrace._fields[pos],
        color='black'
    )
    ax.legend()
    
axs[n,0].plot(
    days,
    [msh.make_npp_func(dvs,svs)(d) for d in days],
    label='NPP',
    color='green'
)
axs[n,0].legend()

# -

# ###### Remark:
# # For simple matrices it is possible to compute the inverse M^{-1} 
# # symbolically ONCE
# mvs.get_CompartmentalMatrix().inv()
# # and then just evaluate it for the X and t along the solution.
# # This could potentially be MUCH faster that inverting 
# # the numerical version of the matrix in every timestep.
# # However it could be very difficult to compute the symbolic inverse in 
# # some cases at all (why we will demonstrate the numeric approach) first. 
