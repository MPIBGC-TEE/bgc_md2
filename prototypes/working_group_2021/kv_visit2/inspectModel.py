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
msh.download_my_TRENDY_output(conf_dict)

#     # Read NetCDF data  ******************************************************************************************************************************
svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
svs_0=msh.Observables(*map(lambda v: v[0],svs))
# -

dvs

# ## Data assimilation
#

# +
#msh.EstimatedParameters._fields
# -

cpa = msh.Constants(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 npp_0=dvs.npp[0],# * 86400,   # kg/m2/s kg/m2/day
 #xi_0=dvs.xi[0],
 rh_0=svs_0.rh,# * 86400,   # kg/m2/s kg/m2/day
 #ra_0=svs_0.ra,# * 86400,   # kg/m2/s kg/m2/day
 #r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
 #r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
 number_of_months=len(svs.rh)
 #number_of_months=120 # for testing and tuning mcmc
)

svs

# provide an inital guess for the paramters to be estimated by the data assimilation
# this is the result of the elaborate procedures 
epa_0=msh.EstimatedParameters(
    beta_leaf=0.6,
    beta_wood=0.25,
    T_0=2,
    E=6.5,
    KM=10,
    #env_modifier=1,
    r_C_leaf_litter_rh=0.0004151100041511,
    r_C_wood_litter_rh=0.00012453300124533,
    r_C_root_litter_rh=0.000122042341220423,
    r_C_soil_fast_rh=0.00015220700152207,
    #r_C_soil_slow_rh=2.73972602739726e-05,
    r_C_soil_slow_rh=3.73972602739726e-05,
    #r_C_soil_passive_rh=7.82778864970646e-06,
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


svs_0.cLitter

# +
## now test it 
import matplotlib.pyplot as plt
from general_helpers import plot_solutions

param2res_sym = msh.make_param2res_sym(mvs,cpa,dvs)
obs_simu= param2res_sym(epa_0)
#obs=np.column_stack([ np.array(v) for v in svs])
#obs=obs[0:cpa.number_of_months,:] #cut 
#obs[:,3:4]=obs[:,3:4]
#### clipping ###
n=len(svs.rh)
obs_arr=np.stack([ arr for arr in svs],axis=1); obs=obs_arr[0:n,:5]
obs_T=msh.Observables(
    cVeg=obs[:,0],
    cLitter=obs[:,1],
    cSoil=obs[:,2],
    rh=obs[:,3],
    #ra=obs[:,4]

)
simu_arr=np.stack([ arr for arr in obs_simu],axis=1); simu=simu_arr[0:n,:5]
simu_T=msh.Observables(
    cVeg=simu[:,0],
    cLitter=simu[:,1],
    cSoil=simu[:,2],
    rh=simu[:,3],
    #ra=simu[:,4]
)

####

fig = plt.figure(figsize=(12, 4), dpi=80)

#fig = plt.figure()
from general_helpers import plot_observations_vs_simulations
plot_observations_vs_simulations(fig,obs_T,simu_T)
# gh.plot_solutions(
#         fig,
#         times=range(n),
#         var_names=msh.Observables._fields,
#         tup=(obs_simu,obs)
#         #tup=(obs,)
# )
fig.savefig('solutions_initial.pdf')
# -
# test cost function
feng_cost_function_2=msh.make_feng_cost_func_2(svs)
feng_cost_function_2(obs_simu)

# +
epa_min=msh.EstimatedParameters(
    beta_leaf=0,
    beta_wood=0,
    T_0=-10,
    E=1,
    KM=1,
    #env_modifier=0,
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
    #env_modifier=10,
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

# ### Initial MCMC run to roughly optimize parameters

# +
from general_helpers import autostep_mcmc, make_feng_cost_func 

#obs=np.column_stack([ np.array(v) for v in svs])
#
#obs=obs[0:cpa.number_of_months,:] #cut 
isQualified = msh.make_param_filter_func(epa_max, epa_min)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = gh.autostep_mcmc(
    initial_parameters=epa_0,
    filter_func=isQualified,
    param2res=param2res,
    costfunction=msh.make_feng_cost_func_2(svs),
    nsimu=2000,# for testing and tuning mcmc
    #nsimu=20000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=10,   # target acceptance rate in %
    chunk_size=100, # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=1,   # default value | increase value to reduce initial step size
    K=1.5 # increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

# optimized parameter set (lowest cost function)
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields),1),axis=1)
epa_opt_1=msh.EstimatedParameters(*par_opt)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)
mod_opt = param2res(epa_opt_1)

# +
# optimized parameter set (lowest cost function)
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields),1),axis=1)
epa_opt_1=msh.EstimatedParameters(*par_opt)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)
mod_opt = param2res(epa_opt_1)

# print optimized parameters
print (epa_opt_1)

# +
# full duration plot

fig = plt.figure(figsize=(12, 4), dpi=80)
plot_observations_vs_simulations(
        fig,
        svs,
        mod_opt
    )
fig.savefig('solutions_mcmc1.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

#import pandas as pd
#pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('visit_da_aa.csv'), sep=',')
#pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('visit_da_j_aa.csv'), sep=',')
#pd.DataFrame(epa_opt_1).to_csv(outputPath.joinpath('visit_optimized_pars.csv'), sep=',')
#pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('visit_optimized_solutions.csv'), sep=',')
# -

# ### Deriving initial pool sizes from the results of 1st mcmc run 

# +
# Get initial pool sizes from the optimized model
from model_specific_helpers_2 import make_param2res_full_output
param2res_full_output = make_param2res_full_output(mvs,cpa,dvs)
mod_opt_full = param2res_full_output(epa_opt_1)
# print optimized model output 
leaf_size=np.mean(mod_opt_full.C_leaf[1200:1750])
wood_size=np.mean(mod_opt_full.C_wood[1200:1750])
root_size=np.mean(mod_opt_full.C_root[1200:1750])
leaf_frac=leaf_size / (leaf_size+wood_size+root_size)
wood_frac=wood_size / (leaf_size+wood_size+root_size)

leaf_litter_size=np.mean(mod_opt_full.C_leaf_litter[1200:1750])
wood_litter_size=np.mean(mod_opt_full.C_wood_litter[1200:1750])
root_litter_size=np.mean(mod_opt_full.C_root_litter[1200:1750])
leaf_litter_frac=leaf_litter_size / (leaf_litter_size+wood_litter_size+root_litter_size)
wood_litter_frac=wood_litter_size / (leaf_litter_size+wood_litter_size+root_litter_size)

soil_fast_size=np.mean(mod_opt_full.C_soil_fast[1200:1750])
soil_slow_size=np.mean(mod_opt_full.C_soil_slow[1200:1750])
soil_passive_size=np.mean(mod_opt_full.C_soil_passive[1200:1750])
soil_fast_frac=soil_fast_size / (soil_fast_size+soil_slow_size+soil_passive_size)
soil_slow_frac=soil_slow_size / (soil_fast_size+soil_slow_size+soil_passive_size)

epa_1=msh.EstimatedParameters(
    beta_leaf=epa_opt_1.beta_leaf,
    beta_wood=epa_opt_1.beta_wood,
    T_0=epa_opt_1.T_0,
    E=epa_opt_1.E,
    KM=epa_opt_1.KM,
    r_C_leaf_litter_rh=epa_opt_1.r_C_leaf_litter_rh,
    r_C_wood_litter_rh=epa_opt_1.r_C_wood_litter_rh,
    r_C_root_litter_rh=epa_opt_1.r_C_root_litter_rh,
    r_C_soil_fast_rh=epa_opt_1.r_C_soil_fast_rh,
    r_C_soil_slow_rh=epa_opt_1.r_C_soil_slow_rh,
    r_C_soil_passive_rh=epa_opt_1.r_C_soil_passive_rh,
    r_C_leaf_2_C_leaf_litter=epa_opt_1.r_C_leaf_2_C_leaf_litter,
    r_C_wood_2_C_wood_litter=epa_opt_1.r_C_wood_2_C_wood_litter,
    r_C_root_2_C_root_litter=epa_opt_1.r_C_root_2_C_root_litter,
    r_C_leaf_litter_2_C_soil_fast=epa_opt_1.r_C_leaf_litter_2_C_soil_fast,
    r_C_leaf_litter_2_C_soil_slow=epa_opt_1.r_C_leaf_litter_2_C_soil_slow,
    r_C_leaf_litter_2_C_soil_passive=epa_opt_1.r_C_leaf_litter_2_C_soil_passive,
    r_C_wood_litter_2_C_soil_fast=epa_opt_1.r_C_wood_litter_2_C_soil_fast,
    r_C_wood_litter_2_C_soil_slow=epa_opt_1.r_C_wood_litter_2_C_soil_slow,
    r_C_wood_litter_2_C_soil_passive=epa_opt_1.r_C_wood_litter_2_C_soil_passive,
    r_C_root_litter_2_C_soil_fast=epa_opt_1.r_C_root_litter_2_C_soil_fast,
    r_C_root_litter_2_C_soil_slow=epa_opt_1.r_C_root_litter_2_C_soil_slow,
    r_C_root_litter_2_C_soil_passive=epa_opt_1.r_C_root_litter_2_C_soil_passive,
    C_leaf_0=leaf_frac * svs_0.cVeg,
    C_wood_0=wood_frac * svs_0.cVeg,
    C_leaf_litter_0=leaf_litter_frac * svs_0.cLitter,
    C_wood_litter_0=wood_litter_frac * svs_0.cLitter,
    C_soil_fast_0=soil_fast_frac * svs_0.cSoil,
    C_soil_slow_0=soil_slow_frac * svs_0.cSoil
)

print(epa_1)
# -

# ### Final MCMC run to optimize parameters 

# +
isQualified = msh.make_param_filter_func(epa_max, epa_min)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = gh.autostep_mcmc(
    initial_parameters=epa_1,
    filter_func=isQualified,
    param2res=param2res,
    costfunction=msh.make_feng_cost_func_2(svs),
    nsimu=2000,# for testing and tuning mcmc
    #nsimu=20000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=10,   # target acceptance rate in %
    chunk_size=100, # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=2,   # increase value to reduce initial step size
    K=1 # increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

# +
# optimized parameter set (lowest cost function)
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields),1),axis=1)
epa_opt=msh.EstimatedParameters(*par_opt)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)
mod_opt = param2res(epa_opt)

# print optimized parameters
print (epa_opt)

# +
# full duration plot

fig = plt.figure(figsize=(12, 4), dpi=80)
plot_observations_vs_simulations(
        fig,
        svs,
        mod_opt
    )
fig.savefig('solutions_full.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

import pandas as pd
pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('visit_da_aa.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('visit_da_j_aa.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('visit_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('visit_optimized_solutions.csv'), sep=',')

# +
# close-up plots (1st 20 years)
n=240
obs_arr=np.stack([ arr for arr in svs],axis=1); obs=obs_arr[0:n,:5]
obs_T=msh.Observables(
    cVeg=obs[:,0],
    cLitter=obs[:,1],
    cSoil=obs[:,2],
    rh=obs[:,3],
    #ra=obs[:,4]

)
simu_arr=np.stack([ arr for arr in mod_opt],axis=1); simu=simu_arr[0:n,:5]
simu_T=msh.Observables(
    cVeg=simu[:,0],
    cLitter=simu[:,1],
    cSoil=simu[:,2],
    rh=simu[:,3],
    #ra=simu[:,4]
)

####

fig = plt.figure(figsize=(12, 4), dpi=80)

#fig = plt.figure()
#from general_helpers import plot_observations_vs_simulations
plot_observations_vs_simulations(fig,obs_T,simu_T)
# gh.plot_solutions(
#         fig,
#         times=range(n),
#         var_names=msh.Observables._fields,
#         tup=(obs_simu,obs)
#         #tup=(obs,)
# )
fig.savefig('solutions_closeup.pdf')


# fig = plt.figure(figsize=(12, 4), dpi=80)
# plot_observations_vs_simulations(
#         fig,
#         svs,
#         mod_opt
#     )

# gh.plot_solutions(
#         fig,
#         #times=range(cpa.number_of_months),
#         times=range(int(cpa.number_of_months)), # for yearly output
#         var_names=msh.Observables._fields,
#         tup=(mod_opt,obs)
# )

# fig.savefig('solutions_opt.pdf')

# -

# ### Traceability analysis  
#
#

# +
it_sym_trace = msh.make_traceability_iterator(mvs,dvs,cpa,epa_opt)
ns=10*360 #1500
StartVectorTrace=gh.make_StartVectorTrace(mvs)
nv=len(StartVectorTrace._fields)
res_trace= np.zeros((ns,nv))
for i in range(ns):
    res_trace[i,:]=it_sym_trace.__next__().reshape(nv)
#res_trace

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
    [msh.make_npp_func(dvs)(d) for d in days],
    label='NPP',
    color='green'
)
axs[n,0].legend()
# -


