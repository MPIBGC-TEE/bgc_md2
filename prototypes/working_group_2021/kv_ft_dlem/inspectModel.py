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
    autostep_mcmc_2,  
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

h.compartmental_graph(mvs)

# #### Matrix equations:

import bgc_md2.display_helpers as dh
dh.mass_balance_equation(mvs)

# ## Download Data
# #### TRENDY Data
# Make sure you have a config.json file in your model folder: <br>
# Config.jon file contents: `{"username": "trendy-v9", "password": "gcb-2020", "dataPath": "/path/to/data/folder"}`

with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 
#msh.download_my_TRENDY_output(conf_dict)

# Subset data to site for simulation:

svs,dvs=msh.get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))
# Look at data
svs

# +
# Subset inital conditions:
# -
svs_0=msh.Observables(*map(lambda v: v[0],svs))

# ## Define Forward Model
# #### Create constants for forward sim:

cpa=msh.Constants(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 npp_0=dvs.npp[0],   # kg/m2/s kg/m2/day
 rh_0=svs_0.rh,   # kg/m2/s kg/m2/day
 #mrso_0=dvs.mrso[0],
 #tsl_0=dvs.tsl[0],
 number_of_years=int(320) # for testing and tuning mcmc
 #number_of_months=len(svs.rh)
)
cpa._asdict()

# #### Create start values for parameters to be optimized during data assimilation:

epa0 = msh.EstimatedParameters(
    beta_leaf=0.3028851224272647,
    beta_wood=0.4564030561461778,
    Theta_sat=0.09469315267122164,
    Theta_fc=0.0795237425542056,
    r_C_leaf_rh=0,
    r_C_wood_rh=0,
    r_C_root_rh=0,
    r_C_aom1_rh=0.0003037126928450523,
    r_C_aom2_rh=0.00029780017612134036,
    r_C_smb1_rh=1.0689106224513612e-05,
    r_C_smb2_rh=0.00017382741384128404,
    r_C_smr_rh=7.29448370877164e-06,
    r_C_nom_rh=1.0164098602070474e-05,
    r_C_dom_rh=1.1500691016736846e-05,
    r_C_psom_rh=1.7644513592821266e-05,
    r_C_leaf_2_C_aom1=0.01970013741245839,
    r_C_leaf_2_C_aom2=0.015317429471496852,
    r_C_wood_2_C_aom1=0.0004097261277630772,
    r_C_wood_2_C_aom2=0.0003689152669226263,
    r_C_root_2_C_aom1=0.0006603079638088481,
    r_C_root_2_C_aom2=0.00043562236684542056,
    r_C_aom1_2_C_smb1=0.00023145662814734428,
    r_C_aom1_2_C_smb2=0.0002977928058255672,
    r_C_aom1_2_C_nom=0.00015287835816890072,
    r_C_aom1_2_C_dom=0.0001781924058628594,
    r_C_aom2_2_C_smb1=2.0810976139531716e-05,
    r_C_aom2_2_C_smb2=2.8349756784116157e-05,
    r_C_aom2_2_C_dom=2.14750752962314e-05,
    r_C_smb1_2_C_nom=2.3404172674296687e-05,
    r_C_smb1_2_C_psom=3.3408243666981594e-05,
    r_C_smb2_2_C_smr=1.2338527499527314e-05,
    r_C_smr_2_C_smb1=4.801573811093888e-05,
    r_C_nom_2_C_smb1=1.0587791011076786e-05,
    r_C_nom_2_C_dom=2.551756974855787e-05,
    r_C_nom_2_C_psom=1.307830769309431e-05,
    r_C_dom_2_C_smb1=1.0950958379980324e-05,
    r_C_dom_2_C_nom=1.4109489107131771e-05,
    r_C_psom_2_C_smb1=2.222172994302753e-06,
    C_leaf_0=0.0034735199394882624,
    C_wood_0=0.23539047059980747,
    C_aom1_0=0.2312532697753018,
    C_smb1_0=2.267901320473235,
    C_smb2_0=0.5198209480949317,
    C_smr_0=0.11596093546336099,
    C_nom_0=2.3560315471577784,
    C_dom_0=3.309888037557665
)

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
epa_max = epa_max._replace(beta_leaf = 0.9)
epa_max = epa_max._replace(beta_wood = 0.9)
epa_max = epa_max._replace(C_leaf_0 = svs_0.cVeg)
epa_max = epa_max._replace(C_wood_0 = svs_0.cVeg)
epa_max = epa_max._replace(C_aom1_0 = svs_0.cLitter)
epa_max = epa_max._replace(C_smb1_0 = svs_0.cSoil)
epa_max = epa_max._replace(C_smb2_0 = svs_0.cSoil)
epa_max = epa_max._replace(C_smr_0 = svs_0.cSoil)
epa_max = epa_max._replace(C_nom_0 = svs_0.cSoil)
epa_max = epa_max._replace(C_dom_0 = svs_0.cSoil)

#print - all names should have numerical values
#epa_max._asdict()
# -

# #### Conduct data assimilation:

param2res=msh.make_param2res_sym(mvs,cpa,dvs)
print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc_2(
    initial_parameters=epa0,
    filter_func=msh.make_param_filter_func(epa_max, epa_min),
    param2res=param2res,
    costfunction=msh.make_weighted_cost_func(svs),
    #nsimu=200, # for testing and tuning mcmc
    nsimu=4000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=0.23,   # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=0.25,   # default value | increase value to reduce initial step size
    K=2 # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

# #### Graph data assimilation results:

# optimized parameter set (lowest cost function)
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields),1),axis=1)
epa_opt=msh.EstimatedParameters(*par_opt)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)
mod_opt = param2res(epa_opt)  
#obs = msh.Observables(cVeg=svs.cVeg,cSoil=svs.cSoil,rh=svs.rh)
print("Forward run with optimized parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(12, 4), dpi=80)
plot_observations_vs_simulations(
        fig,
        svs,
        mod_opt
    )
fig.savefig('solutions_opt_DLEM.pdf')
# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)
