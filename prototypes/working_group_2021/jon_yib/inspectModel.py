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

# ## Create Symbols for $\xi$, $K$, and $A$ (No Edits)
# Setup Yiqi matrix format:

# +
sv=mvs.get_StateVariableTuple()                            # Get state variables from symbolic matrix code
n=len(sv)                                                  # Count number of pools
srm = mvs.get_SmoothReservoirModel()                       # Create smooth resevoir model
_,A,N,_,_=srm.xi_T_N_u_representation(factor_out_xi=False) # Create A and N symbols
BI=mvs.get_BibInfo()
for k in BI.sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")
beta_wood = 1.0-(beta_leaf+beta_root)
beta_wood = 1.0-(beta_leaf+beta_root)

#create symbols for scaler and input functions
func_dict={
    'xi': 'Environmental scaler as a function of time',
    'NPP': 'Inputs as a function of time',
}
for k in BI.func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

# -

# $\xi$ Matrix:

# Create environmental scaler matrix
xi_d=diag([1,1,1]+[xi(t) for i in range(n-3)],unpack=True)
xi_d

# $K$ Matrix:

# +
# Create empty K matrix
K=xi_d.inv()*N
# Create new symbols for the k_{i}
for i in range(n):
    if K[i,i]!=0:
        name="k_{0}".format(sv[i])
        code="{0}=Symbol('{0}')".format(name)
        #print(code)
        exec(code)

# Create $K$ matrix      
K_sym=ImmutableMatrix(
    n,n,
    lambda i,j: Symbol("k_" + str(sv[i])) if i==j else 0
)
K_sym
# -

# $f$ symbols in $A$ Matrix:

# +
# Create new symbols for the f_{i,j}
for i in range(n):
    for j in range(n):
        if A[i,j]!=0 and i!=j:
            name="f_" + str(sv[j]) + "_2_" + str(sv[i])
            code="{0}=Symbol('{0}')".format(name)
            #print(code)
            exec(code)

# Place $f$ values in $A$ matrix            
A_sym=ImmutableMatrix(
    n,n,
    lambda i,j:  -1 if i==j else (
        0 if A[i,j]==0 else Symbol("f_" + str(sv[j]) + "_2_" + str(sv[i]))
    )
)
A_sym
# -

# $A$ matrix:

# Create A matrix
M_sym=A_sym*K_sym
M_sym

# ## Create Dictionary of All Fluxes (No Edits)

# +
# Create a dictionary for the external and internal fluxes (flux divided by dono pool content)
outflux_rates = {"r_"+str(key)+"_rh":value/key for key,value in hr.out_fluxes_by_symbol(sv,M_sym).items() }
internal_flux_rates = {"r_"+str(key[0])+"_2_"+str(key[1]):value/key[0] for key,value in hr.internal_fluxes_by_symbol(sv,M_sym).items()}

# Create dictionary of all flux rates
all_rates=deepcopy(outflux_rates)
all_rates.update(internal_flux_rates)
all_rates
# -

# ## Calculate Rates from $f$ and $k$ values (Must Edit)
# I have $k$ and $f$ values describing my model. we can define them here and use them to assign values to $r$s

# + codehighlighter=[[3, 22], [26, 45], [48, 79], [83, 85], [3, 22], [26, 45], [48, 79], [83, 85]]
# fixme mm
# The followiwng namedtuple is only used once.
# It would be much simpler to just add the values directly to old_par_dict
# Is this a reference to the old implementation?
ParameterValues = namedtuple(
    "ParameterValues",
    [
        "beta_leaf",
        "beta_root",
        "clay",
        "silt",
        "k_leaf",
        "k_root",
        "k_wood",
        "k_cwd",
        "k_samet",
        "k_sastr",
        "k_samic",
        "k_slmet",
        "k_slstr",
        "k_slmic",
        "k_slow",
        "k_arm",
        "f_samet_leaf",
        "f_slmet_root",
        "f_samic_cwd",
    ]
)

epa0 = ParameterValues(
    beta_leaf=0.3,
    beta_root=0.3,
    clay=0.2028,
    silt=0.2808,
    k_leaf=0.020,
    k_root=0.010,
    k_wood=0.007,
    k_cwd=0.01,
    k_samet=0.05,
    k_sastr=0.05,
    k_samic=0.05,
    k_slmet=0.040,
    k_slstr=0.039,
    k_slmic=0.05,
    k_slow=0.0001,
    k_arm=3.27E-04,
    f_samet_leaf=0.50,
    f_slmet_root=0.50,
    f_samic_cwd=0.50,
)

old_par_dict = {
    'k_c_leaf': epa0.k_leaf, # define all k values
    'k_c_root': epa0.k_root,
    'k_c_wood': epa0.k_wood,
    'k_c_lit_cwd': epa0.k_cwd,
    'k_c_lit_met': epa0.k_samet,
    'k_c_lit_str': epa0.k_sastr,
    'k_c_lit_mic': epa0.k_samic,
    'k_c_soil_met': epa0.k_slmet,
    'k_c_soil_str': epa0.k_slstr,
    'k_c_soil_mic': epa0.k_slmic,
    'k_c_soil_slow': epa0.k_slow,
    'k_c_soil_passive': epa0.k_arm,
    'f_c_leaf_2_c_lit_met': epa0.f_samet_leaf,    #define all f values
    'f_c_root_2_c_soil_met': epa0.f_slmet_root,
    'f_c_lit_cwd_2_c_lit_mic': epa0.f_samic_cwd*0.7,
    'f_c_leaf_2_c_lit_str': (1-epa0.f_samet_leaf),
    'f_c_root_2_c_soil_str': (1-epa0.f_slmet_root),
    'f_c_wood_2_c_lit_cwd': 1,
    'f_c_lit_cwd_2_c_soil_slow': (1-epa0.f_samic_cwd),
    'f_c_lit_met_2_c_lit_mic': 0.2,
    'f_c_lit_str_2_c_lit_mic': 0.2,
    'f_c_lit_str_2_c_soil_slow': 0.2,
    'f_c_lit_mic_2_c_soil_slow': 0.2,
    'f_c_soil_met_2_c_soil_mic': 0.2,
    'f_c_soil_str_2_c_soil_mic': 0.2,
    'f_c_soil_str_2_c_soil_slow': 0.2,
    'f_c_soil_mic_2_c_soil_slow': 0.2,
    'f_c_soil_mic_2_c_soil_passive': 0.2,
    'f_c_soil_slow_2_c_soil_mic': 0.2,
    'f_c_soil_slow_2_c_soil_passive': 0.2*(0.003+0.009*epa0.clay),
    'f_c_soil_passive_2_c_soil_mic': 0.2, 
}

# Define allocation parameters to be optimized
par_dict = {
    'beta_leaf': epa0.beta_leaf,
    'beta_root': epa0.beta_root,
}

# translate rates from previous parameters to create dictionary of rates to optimize
par_dict.update(
    {str(k):v.subs(old_par_dict) for k,v in all_rates.items()}
)

# Create namedtuple of parameters to optimize and their translated values
#makeTuple = namedtuple('makeTuple', par_dict)
#parameters = makeTuple(**par_dict)

#If symbols remain in output below then set them to numerical values in old_par_dict.
#parameters._asdict() # print - everything below should have a numeric value
# -

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

epa0 = msh.EstimatedParameters(
    **{
        "c_leaf_0": svs_0.cVeg/3,          #set inital pool values to svs values 
        "c_root_0": svs_0.cVeg/3,          #you can set numerical values here directly as well
        "c_lit_cwd_0": svs_0.cSoil/35,
        "c_lit_met_0": svs_0.cSoil/35,
        "c_lit_str_0": svs_0.cSoil/35,
        "c_lit_mic_0": svs_0.cSoil/35,
        "c_soil_met_0": svs_0.cSoil/20,
        "c_soil_str_0": svs_0.cSoil/15,
        "c_soil_mic_0": svs_0.cSoil/10,
        "c_soil_slow_0": svs_0.cSoil/3
    },
    **par_dict
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
epa_max._asdict()
# -

# #### Conduct data assimilation:

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa0,
    filter_func=make_param_filter_func(epa_max, epa_min),
    param2res=msh.make_param2res_sym(mvs,cpa,dvs),
    costfunction=msh.make_weighted_cost_func(svs),
    #nsimu=200, # for testing and tuning mcmc
    nsimu=2000,
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
plot_solutions(
        fig,
        #times=range(cpa.number_of_months),
        times=range(n), # for yearly output
        var_names=msh.observables._fields,
        tup=(mod_opt,obs)
        #tup=(obs,)
)

fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

import pandas as pd
pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('YIBs_da_pars.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('YIBS_da_cost.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('YIBs_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('YIBs_optimized_solutions.csv'), sep=',')
# -



