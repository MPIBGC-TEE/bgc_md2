# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Minimal Code Example: Symbolic JULES Model
# ## Python/Jupyter Setup (No edits)
# Jupyter Settings:
#
# Notes:
# 2020-02-18: copied from jon_yib/createModel3.py
#

# +
# load HTML to adjust jupyter settings
from IPython.display import HTML

# adjust jupyter display to full screen width
display(HTML("<style>.container { width:100% !important; }</style>"))

# set auto reload for notebook
# %load_ext autoreload
# %autoreload 2

# +
# Packages for symbolic code:
from sympy import Symbol, Function
from ComputabilityGraphs.CMTVS import CMTVS
from sympy import Symbol, Function, diag, ImmutableMatrix
from pathlib import Path
from copy import copy, deepcopy
from functools import reduce
from typing import Callable
from pprint import pprint
import json
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
#from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
#import CompartmentalSystems.helpers_reservoir as hr

from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

import bgc_md2.resolve.computers as bgc_c
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h
from collections import namedtuple

# Other packages
import sys

sys.path.insert(0, '..')  # necessary to import general_helpers
#from general_helpers import (
#    download_TRENDY_output,
#    day_2_month_index,
#    month_2_day_index,
#    make_B_u_funcs_2,
#    monthly_to_yearly,
#    plot_solutions
#)
import general_helpers as gh
import model_specific_helpers_2 as msh
import test_helpers as th
from source import mvs
# -

# ## Model Figure and Matrix Equations
# #### Model Figure:

h.compartmental_graph(mvs)

mvs.get_CompartmentalMatrix()

mvs.get_InternalFluxesBySymbol()

# + [markdown] codehighlighter=[[0, 0]]
# ## Download Data
# #### TRENDY Data
# Make sure you have a config.json file in your model folder: <br>
# Config.jon file contents: `{"username": "trendy-v9", "password": "gcb-2020", "dataPath": "/path/to/data/folder"}`

# +
import json

# Read username, password, dataPath from config.json file
with Path('config.json').open(mode='r') as f:
    conf_dict = json.load(f)
    
# -

#msh.download_my_TRENDY_output(conf_dict)
ta=th.make_test_args(conf_dict,msh,mvs)

# ## Connect Data and Symbols (Must Edit)
# Define function to subset netCDF files and link to data symbols:

# + codehighlighter=[[5, 6], [23, 33], [5, 6], [23, 33]]
## call function to link symbols and data
svs, dvs = msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
# -

Path(conf_dict["dataPath"])

# +
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms


fig, ax = plt.subplots() #[(12*10-1):12*30]
plt.plot(list(range(0, 12*320)), np.array(svs.cVeg), label = "observation")
# -

# ## Create Symbols for $\xi$, $K$, and $A$ (No Edits)
# Setup Yiqi matrix format:

# + codehighlighter=[[3, 22], [26, 45], [48, 79], [83, 85], [3, 22], [26, 45], [48, 79], [83, 85]]
par_dict = {
    Symbol(k): v
    for k, v in {
        'beta_leaf': 0.35,
        'beta_wood': 0.30,
        # 'r_c_leaf_rh': 0,
        # 'r_c_wood_rh': 0,
        # 'r_c_root_rh': 0,
        'r_c_DPM_rh':0.0218855218855219,
        'r_c_RPM_rh':0.000866666666666667,
        'r_c_BIO_rh':0.00174841269841270,
        'r_c_HUM_rh':5.87450980392157e-5,
        'r_c_leaf_2_c_DPM':0.000152777777777778,
        'r_c_leaf_2_c_RPM':0.000541666666666667,
        'r_c_wood_2_c_DPM':2.00364298724954e-5,
        'r_c_wood_2_c_RPM':7.10382513661202e-5,
        'r_c_root_2_c_DPM':0.000152777777777778,
        'r_c_root_2_c_RPM':0.000541666666666667,
        'r_c_DPM_2_c_BIO':0.00283950617283951,
        'r_c_DPM_2_c_HUM':0.00333333333333333,
        'r_c_RPM_2_c_BIO':0.000112444444444444,
        'r_c_RPM_2_c_HUM':0.000132000000000000,
        'r_c_BIO_2_c_HUM':0.000235714285714286,
        'r_c_HUM_2_c_BIO':6.61437908496732e-6
    }.items()
}
#par_dict

# -


svs_0 = msh.Observables(*map(lambda v: v[0], svs))
svs_0.cSoil

# ## Assign Initial Values for the iterator

# + codehighlighter=[[5, 17], [5, 17]]
# Create vector of initial pool values
svs_0 = msh.Observables(*map(lambda v: v[0], svs))

StartVector = msh.make_StartVector(mvs)
StartVector._fields

# + codehighlighter=[[5, 17], [5, 17]]
# Assign values to initial pools using InitialPools named tuple
V_init = StartVector(
    c_leaf=svs_0.cVeg * 0.12,  # set inital pool values to svs values
    c_root=svs_0.cVeg * 0.12,  # you can set numerical values here directly as well
    c_wood=svs_0.cVeg * 0.76,
    c_DPM=svs_0.cSoil * 0.0025,
    c_RPM=svs_0.cSoil * 0.248,
    c_BIO=svs_0.cSoil * 0.022,
    c_HUM=svs_0.cSoil * 0.7275,
    rh=svs_0.rh,
    fVegSoil=svs_0.fVegSoil
    # f_veg2soil=svs_0.f_veg2soil# add the fraction
)
V_init._asdict()  # print - everything should have a numeric value
# -

# ## Define Forward Model
# #### Create constants for forward sim:

# + codehighlighter=[[1, 9], [1, 8]]
cpa = msh.Constants(
    # use Constants namedtuple to define constant values #Define the constant values of parameters NOT affected by data assimilation
    npp_0=dvs.npp[0],
    rh_0=svs_0.rh,
    c_veg_0=svs_0.cVeg,
    c_soil_0=svs_0.cSoil,
    fVegSoil_0=svs_0.fVegSoil,  # add the fraction
    nyears=320
)
cpa._asdict()  # print - everything should have a numeric value
# -

# #### Create list of parameters to be optimized during data assimilation:

epa_0 = msh.EstimatedParameters(
    **{
        'c_leaf_0': svs_0.cVeg * 0.12,  # set inital pool values to svs values
        'c_wood_0': svs_0.cVeg * 0.76,  # you can set numerical values here directly as well
        'c_DPM_0': svs_0.cSoil * 0.0025,  # set according to QY's single site results: 0.0025 DPM, 0.22 RPM, 0.02 BIO, 0.7575 HUM
        'c_RPM_0': svs_0.cSoil * 0.248,
        'c_BIO_0': svs_0.cSoil * 0.022
    },
    **{
        'beta_leaf': 0.35,
        'beta_wood': 0.3,
        'Mw': 0.1,
        'Ms': np.max(dvs.mrsos) + 500, #, may need add a condition here ## ASK MARKUS
        'Topt': 18.32,
        'Tcons': 47.91,
        # 'r_c_leaf_rh': 0,
        # 'r_c_wood_rh': 0,
        # 'r_c_root_rh': 0,
        'r_c_DPM_rh':0.0218855218855219,
        'r_c_RPM_rh':0.000866666666666667,
        'r_c_BIO_rh':0.00174841269841270,
        'r_c_HUM_rh':5.87450980392157e-5,
        'r_c_leaf_2_c_DPM':0.000152777777777778,
        'r_c_leaf_2_c_RPM':0.000541666666666667,
        'r_c_wood_2_c_DPM':2.00364298724954e-5,
        'r_c_wood_2_c_RPM':7.10382513661202e-5,
        'r_c_root_2_c_DPM':0.000152777777777778,
        'r_c_root_2_c_RPM':0.000541666666666667,
        'r_c_DPM_2_c_BIO':0.00283950617283951,
        'r_c_DPM_2_c_HUM':0.00333333333333333,
        'r_c_RPM_2_c_BIO':0.000112444444444444,
        'r_c_RPM_2_c_HUM':0.000132000000000000,
        'r_c_BIO_2_c_HUM':0.000235714285714286,
        'r_c_HUM_2_c_BIO':6.61437908496732e-6
    }
    # **{str(key): value for key,value in  par_dict.items() if}
)

# #### Create forward model function:

# + codehighlighter=[[37, 51], [67, 69], [64, 65], [137, 139], [133, 135], [32, 45], [112, 113], [117, 118], [120, 123]]
func_dict = msh.make_func_dict(mvs, dvs, cpa, epa_0)

# +
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms


fig, ax = plt.subplots() #[(12*10-1):12*30]
plt.plot(list(range(0, 12*320)), np.array(svs.cVeg), label = "observation")
# -

# ## Forward Model Run
# #### Run model forward:

# +
epa_0 = ta.epa_0
cpa = ta.cpa
dvs = ta.dvs
mvs = ta.mvs
svs = ta.svs
param2res = msh.make_param2res_sym(mvs, cpa, dvs)  # Define forward model
obs_simu = param2res(epa_0)  # Run forward model from initial conditions

obs_simu

svs_cut = msh.Observables(
    cVeg = svs.cVeg[:cpa.nyears * 12],
    cSoil = svs.cSoil[:cpa.nyears * 12],
    rh = svs.rh[:cpa.nyears * 12],
    fVegSoil = svs.fVegSoil[:cpa.nyears*12]
)
svs_cut.rh
cost_func = msh.make_weighted_cost_func_2(svs_cut)
cost_func(obs_simu)
# -

# #### Create array of yearly observation data:

# +
#fig = plt.figure()
#gh.plot_observations_vs_simulations(fig, svs_cut, obs_simu)
#fig.savefig("test.pdf")
# -

# #### Plot data-model fit:

# +
# Plot simulation output for observables
n_plots = len(svs_cut)
fig = plt.figure(figsize=(10, n_plots * 5))
axs = fig.subplots(n_plots)
for i, name in enumerate(msh.Observables._fields):
    var = svs_cut[i]
    var_simu = obs_simu[i]
    axs[i].plot(range(len(var_simu)), var_simu, label="simulation")
    axs[i].plot(range(len(var)), var, label='observation')
    axs[i].legend()
    axs[i].set_title(name)

sum((np.array(svs.rh) - np.array(obs_simu.rh))**2)
# -

# ## Data Assimilation
# #### Define parameter min/max values:

# #### Conduct data assimilation:

# +

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = gh.autostep_mcmc(
    initial_parameters=ta.epa_0,
    filter_func=gh.make_param_filter_func_2(ta.epa_max, ta.epa_min,["beta_leaf","beta_wood"]),
    param2res=msh.make_param2res_sym(mvs, cpa, dvs),
    costfunction=msh.make_weighted_cost_func_2(svs),
    nsimu=5000,  # for testing and tuning mcmc
    c_max=np.array(ta.epa_max),
    c_min=np.array(ta.epa_min),
    acceptance_rate=15,  # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=1,  # default value | increase value to reduce initial step size
    K=2  # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")
# -

# #### Graph data assimilation results:

# +
# optimized parameter set (lowest cost function)
par_opt = np.min(
    C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields), 1),
    axis=1)
epa_opt = msh.EstimatedParameters(*par_opt)
mod_opt = param2res(epa_opt)

n_plots = len(svs)
print("Forward run with optimized parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(10, n_plots * 5))
# Plot simulation output for observables
# n_plots=len(svs_cut)
fig = plt.figure(figsize=(10, 10 * n_plots))
axs = fig.subplots(n_plots)
plt.rcParams['font.size'] = 18
for i, name in enumerate(msh.Observables._fields):
    var_simu = mod_opt[i]
    var = svs[i]
    axs[i].plot(range(len(var_simu)), var_simu, label='simulation(opt)')
    axs[i].plot(range(len(var)), var, label='observation')
    axs[i].set_title(name)
    axs[i].legend()

fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath = Path(conf_dict["dataPath"])  # save output to data directory (or change it)

import pandas as pd

pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('JULES_da_pars.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('JULES_da_cost.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('JULES_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('JULES_optimized_solutions.csv'), sep=',')
# -
ta.epa_min

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
fig=plt.figure(figsize=(20,(n+1)*10), dpi = 400)
axs=fig.subplots(n+1,2)
plt.rcParams['font.size'] = 20

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

