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

#msh.download_my_TRENDY_output(conf_dict)

#     # Read NetCDF data  ******************************************************************************************************************************
svs,dvs=msh.get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))
svs_0=msh.Observables(*map(lambda v: v[0],svs))
# -

# ## Data assimilation
#

# +
#msh.EstimatedParameters._fields
# -

cpa = msh.Constants(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 npp_0=dvs.npp[0] * 86400,   # kg/m2/s kg/m2/day
 rh_0=svs_0.rh * 86400,   # kg/m2/s kg/m2/day
 ra_0=svs_0.ra * 86400,   # kg/m2/s kg/m2/day
 r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
 r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
 #number_of_months=len(svs.rh)
 number_of_months=120 # for testing and tuning mcmc
)

# provide an inital guess for the paramters to be estimated by the data assimilation
# this is the result of the elaborate procedures 
epa_0=msh.EstimatedParameters(
    beta_leaf=0.6,
    beta_wood=0.25,
    T_0=2,
    E=4,
    KM=10,
    r_C_leaf_litter_rh=0.0004151100041511,
    r_C_wood_litter_rh=0.00012453300124533,
    r_C_root_litter_rh=0.000122042341220423,
    r_C_soil_fast_rh=0.00015220700152207,
    r_C_soil_slow_rh=2.73972602739726e-05,
    r_C_soil_passive_rh=7.82778864970646e-06,
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
    C_leaf_litter_0=0.5202311902470766,
    C_wood_litter_0=0.7225433197876749,
    C_soil_fast_0=1.7309510511856925,
    C_soil_slow_0=2.4435101360092473
)


# +
## now test it 
#import matplotlib.pyplot as plt
#from general_helpers import plot_solutions
#
#param2res_sym = msh.make_param2res_sym(mvs,cpa,dvs)
#xs= param2res_sym(epa_0)
#obs=np.column_stack([ np.array(v) for v in svs])
#obs=obs[0:cpa.number_of_months,:] #cut 
##obs[:,3:4]=obs[:,3:4]
#n=cpa.number_of_months
#
#
#fig = plt.figure(figsize=(12, 4), dpi=80)
#gh.plot_solutions(
#        fig,
#        times=range(n),
#        var_names=msh.Observables._fields,
#        tup=(xs,obs)
#        #tup=(obs,)
#)
#fig.savefig('solutions.pdf')

# +
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
# -

# ### mcmc to optimize parameters 
#

# +
from general_helpers import autostep_mcmc, make_feng_cost_func 

obs=np.column_stack([ np.array(v) for v in svs])
#
obs=obs[0:cpa.number_of_months,:] #cut 
isQualified = msh.make_param_filter_func(epa_max, epa_min)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)

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
gh.plot_solutions(
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

# ### Traceability analysis  
#
#

it_sym_trace = msh.make_traceability_iterator(mvs,dvs,cpa,epa_opt)
ns=1500
StartVectorTrace=gh.make_StartVectorTrace(mvs)
nv=len(StartVectorTrace._fields)
res_trace= np.zeros((ns,nv))
for i in range(ns):
    res_trace[i,:]=it_sym_trace.__next__().reshape(nv)
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
    [msh.make_npp_func(dvs)(d) for d in days],
    label='NPP',
    color='green'
)
axs[n,0].legend()

# -




