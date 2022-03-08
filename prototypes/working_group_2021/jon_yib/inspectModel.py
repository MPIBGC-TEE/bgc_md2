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

# # Minimal Code Example: Symbolic YIBs Model
# ## Python/Jupyter Setup (No edits)
# Jupyter Settings:

# +
#load HTML to adjust jupyter settings
from IPython.display import HTML

#adjust jupyter display to full screen width
display(HTML("<style>.container { width:100% !important; }</style>"))

#set auto reload for notebook
# %load_ext autoreload
# %autoreload 2
# -

# Python Packages:

# +
# Packages for symbolic code: 
from sympy import Symbol, Function, diag, ImmutableMatrix 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple
)
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
import CompartmentalSystems.helpers_reservoir as hr
import bgc_md2.resolve.computers as bgc_c
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h
from collections import namedtuple

# Other packages
import sys
sys.path.insert(0,'..') # necessary to import general_helpers
from general_helpers import (
    download_TRENDY_output,
    day_2_month_index,
    month_2_day_index,
    make_B_u_funcs_2,
    monthly_to_yearly,
    plot_solutions,
    autostep_mcmc,  
    make_jon_cost_func,
    make_param_filter_func
)
from pathlib import Path
from copy import copy, deepcopy
from functools import reduce
from typing import Callable
from pprint import pprint
import json 
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from source import mvs
import model_specific_helpers_2 as msh
# -

# ## Model Figure and Matrix Equations
# #### Model Figure:

# +

h.compartmental_graph(mvs)
# -

# #### Matrix equations:

dh.mass_balance_equation(mvs)

# + [markdown] codehighlighter=[[0, 0]]
# ## Download Data (Must Edit)
# #### TRENDY Data
# Make sure you have a config.json file in your model folder: <br>
# Config.jon file contents: `{"username": "trendy-v9", "password": "gcb-2020", "dataPath": "/path/to/data/folder"}`

# + codehighlighter=[[11, 12], [16, 17], [8, 28], [41, 43], [8, 24], [42, 44]]
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

#msh.download_my_TRENDY_output(conf_dict)
# -

# ## Connect Data and Symbols (Must Edit)
# Define function to subset netCDF files and link to data symbols:

# + codehighlighter=[[5, 6], [23, 33], [5, 6], [23, 33]]
svs,dvs=msh.get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))

# + codehighlighter=[[5, 6], [23, 33], [5, 6], [23, 33]]
#look at data
dvs.npp
# -

msh.EstimatedParameters._fields

svs_0=msh.Observables(*map(lambda v: v[0],svs))
dvs.npp

# ## Define Forward Model
# #### Create constants for forward sim:

# + codehighlighter=[[1, 9], [1, 8]]
cpa = msh.Constants(             #use Constants namedtuple to define constant values
    npp_0 = dvs.npp[0],
    rh_0 = svs.rh[0],   
    c_veg_0 = svs.cVeg[0],
    c_soil_0 = svs.cSoil[0],
    clay = 0.2028,
    silt = 0.2808,
    nyears = 320
)
cpa._asdict()    #print - everything should have a numeric value
# -

# #### Create start values for parameters to be optimized during data assimilation:

# how we transform given startvalues for the f and k to these is shown in createModel
# but once we have them, we can print them out and use them from now on directly
epa0 =msh.EstimatedParameters(
    c_leaf_0=svs_0.cVeg/3,          #set inital pool values to svs values 
    c_root_0=svs_0.cVeg/3,          #you can set numerical values here directly as well
    c_lit_cwd_0=svs_0.cSoil/35,
    c_lit_met_0=svs_0.cSoil/35,
    c_lit_str_0=svs_0.cSoil/35,
    c_lit_mic_0=svs_0.cSoil/35,
    c_soil_met_0=svs_0.cSoil/20,
    c_soil_str_0=svs_0.cSoil/15,
    c_soil_mic_0=svs_0.cSoil/10,
    c_soil_slow_0=svs_0.cSoil/3,
    beta_leaf=0.3,
    beta_root=0.3,
    r_c_leaf_rh=0,
    r_c_root_rh=0,
    r_c_wood_rh=0,
    r_c_lit_cwd_rh=0.00150000000000000,
    r_c_lit_met_rh=0.0400000000000000,
    r_c_lit_str_rh=0.0300000000000000,
    r_c_lit_mic_rh=0.0400000000000000,
    r_c_soil_met_rh=0.0320000000000000,
    r_c_soil_str_rh=0.0234000000000000,
    r_c_soil_mic_rh=0.0300000000000000,
    r_c_soil_slow_rh=7.99034960000000e-5,
    r_c_soil_passive_rh=0.000261600000000000,
    r_c_leaf_2_c_lit_met=0.0100000000000000,
    r_c_leaf_2_c_lit_str=0.0100000000000000,
    r_c_root_2_c_soil_met=0.00500000000000000,
    r_c_root_2_c_soil_str=0.00500000000000000,
    r_c_wood_2_c_lit_cwd=0.00700000000000000,
    r_c_lit_cwd_2_c_lit_mic=0.00350000000000000,
    r_c_lit_cwd_2_c_soil_slow=0.00500000000000000,
    r_c_lit_met_2_c_lit_mic=0.0100000000000000,
    r_c_lit_str_2_c_lit_mic=0.0100000000000000,
    r_c_lit_str_2_c_soil_slow=0.0100000000000000,
    r_c_lit_mic_2_c_soil_slow=0.0100000000000000,
    r_c_soil_met_2_c_soil_mic=0.00800000000000000,
    r_c_soil_str_2_c_soil_mic=0.00780000000000000,
    r_c_soil_str_2_c_soil_slow=0.00780000000000000,
    r_c_soil_mic_2_c_soil_slow=0.0100000000000000,
    r_c_soil_mic_2_c_soil_passive=0.0100000000000000,
    r_c_soil_slow_2_c_soil_mic=2.00000000000000e-5,
    r_c_soil_slow_2_c_soil_passive=9.65040000000000e-8,
    r_c_soil_passive_2_c_soil_mic=6.54000000000000e-5,
)    


# +
def npp_func(day):
        month=day_2_month_index(day)
        return dvs.npp[month]

n=cpa.nyears*12
npp_obs = np.array([npp_func(d) for d in range(n)])

# Plot simulation output for observables
fig = plt.figure()
plot_solutions(
        fig,
        times=range(n),
        var_names=msh.Observables._fields,
        tup=(npp_obs,)
)
fig.savefig('solutions.pdf')
# -

# #### Create forward model function:

# ## Forward Model Run
# #### Run model forward:

param2res_sym = msh.make_param2res_sym(mvs,cpa,dvs) # Define forward model
obs_simu = param2res_sym(epa0)                # Run forward model from initial conditions

# #### Plot data-model fit:

fig = plt.figure()
from general_helpers import plot_observations_vs_simulations
plot_observations_vs_simulations(fig,svs,obs_simu)


# ## Data Assimilation
# #### Define parameter min/max values:

# +
# set min/max parameters to +- 100 times initial values
epa_min=msh.EstimatedParameters._make(tuple(np.array(epa0)*0.01))
epa_max=msh.EstimatedParameters._make(tuple(np.array(epa0)*100))

# fix values that are problematic from calculation
epa_max = epa_max._replace(beta_leaf = 0.99)
epa_max = epa_max._replace(beta_root = 0.99)
epa_max = epa_max._replace(c_leaf_0 = svs_0.cVeg)
epa_max = epa_max._replace(c_root_0 = svs_0.cVeg)
epa_max = epa_max._replace(c_lit_cwd_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_lit_met_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_lit_str_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_lit_mic_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_soil_met_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_soil_str_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_soil_mic_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_soil_slow_0 = svs_0.cSoil)

#print - all names should have numerical values
#epa_max._asdict()
# -

# #### Conduct data assimilation:

param2res=msh.make_param2res_sym(mvs,cpa,dvs)
print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa0,
    filter_func=msh.make_param_filter_func(epa_max, epa_min),
    param2res=param2res,
    costfunction=msh.make_weighted_cost_func(svs),
    #nsimu=200, # for testing and tuning mcmc
    nsimu=200,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=15,   # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=1,   # default value | increase value to reduce initial step size
    K=2 # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

# #### Graph data assimilation results:

# +
# optimized parameter set (lowest cost function)
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields),1),axis=1)
epa_opt=msh.EstimatedParameters(*par_opt)
mod_opt = param2res(epa_opt)  

print("Forward run with optimized parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(12, 4), dpi=80)
plot_observations_vs_simulations(
        fig,
        svs,
        mod_opt
    )
fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

import pandas as pd
pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('YIBs_da_pars.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('YIBS_da_cost.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('YIBs_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('YIBs_optimized_solutions.csv'), sep=',')
# +
import model_specific_helpers_2 as msh
import general_helpers as gh
it_sym_trace = msh.make_traceability_iterator(mvs,dvs,cpa,epa_opt)
ns=1500
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

