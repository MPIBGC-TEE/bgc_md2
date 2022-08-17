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

# -

# ## Symbolic Setup (Must Edit)
# Define Symbols using named tuples for allocation, pools, fluxes, constants:

# + codehighlighter=[[4, 7], [14, 26], [33, 45], [52, 80], [93, 101], [117, 118], [4, 7], [14, 26], [33, 45], [52, 80], [93, 100], [116, 117]]
#Create namedtuple of allocation coefficients
Allocation = namedtuple(
    "Allocation",
    [
        'beta_leaf',              #Names: beta_poolname
        'beta_root',
        'beta_wood'
    ]
)

#Create namedtuple of pools
Pools = namedtuple(               
    "Pools",
    [
        'c_leaf',                 #Names: c_poolname
        'c_root',
        'c_wood',
        'c_lit_cwd',
        'c_lit_met',
        'c_lit_str',
        'c_lit_mic',
        'c_soil_met',
        'c_soil_str',
        'c_soil_mic',
        'c_soil_slow',
        'c_soil_passive'
    ]
)

#Create namedtuple of initial pools
InitialPools = namedtuple(               
    "InitialPools",
    [
        'c_leaf_0',               #Names: c_poolname_0
        'c_root_0',               #Only initial pools that are estimated
        'c_lit_cwd_0',
        'c_lit_met_0',
        'c_lit_str_0',
        'c_lit_mic_0',
        'c_soil_met_0',
        'c_soil_str_0',
        'c_soil_mic_0',
        'c_soil_slow_0',
    ]
)

#Create namedtuple of flux rates leaving the system
FluxRates = namedtuple(
    "FluxRates",
    [
        'r_c_lit_cwd_rh',         #Pools with loss from system will be listed here
        'r_c_lit_met_rh',         #Names: r_c_poolname_rh
        'r_c_lit_str_rh',
        'r_c_lit_mic_rh',
        'r_c_soil_met_rh',
        'r_c_soil_str_rh',
        'r_c_soil_mic_rh',
        'r_c_soil_slow_rh',
        'r_c_soil_passive_rh', 
        'r_c_leaf_2_c_lit_met',   #Pool transfer paths
        'r_c_leaf_2_c_lit_str',   #Names: r_c_donorPool_2_recievingPool
        'r_c_root_2_c_soil_met',    
        'r_c_root_2_c_soil_str',
        'r_c_wood_2_c_lit_cwd',
        'r_c_lit_cwd_2_c_lit_mic',
        'r_c_lit_cwd_2_c_soil_slow',
        'r_c_lit_met_2_c_lit_mic',
        'r_c_lit_str_2_c_lit_mic',
        'r_c_lit_str_2_c_soil_slow',
        'r_c_lit_mic_2_c_soil_slow',
        'r_c_soil_met_2_c_soil_mic',
        'r_c_soil_str_2_c_soil_mic',
        'r_c_soil_str_2_c_soil_slow',
        'r_c_soil_mic_2_c_soil_slow',
        'r_c_soil_mic_2_c_soil_passive',
        'r_c_soil_slow_2_c_soil_mic',
        'r_c_soil_slow_2_c_soil_passive',
        'r_c_soil_passive_2_c_soil_mic'
    ]
)

#define time
Time = namedtuple(
    "Time",
    ['t']
)

#Create namedtuple of constants used in model
Constants = namedtuple(
    "Constants",
    [
        'npp_0',       #Initial input/pools
        'rh_0',
        'c_veg_0',
        'c_soil_0',
        'clay',        #Constants like clay
        'silt',
        'nyears'       #Run time (years for my model)
    ]
)

#Combine all named tuples to create symbols
Symbols = namedtuple(
    "Symbols",
    Allocation._fields + Pools._fields + InitialPools._fields + \
    FluxRates._fields + Time._fields + Constants._fields
)

#Create symbols
for k in Symbols._fields:
    code=k+" = Symbol('{0}')".format(k)
    exec(code)
    
#define beta wood from other allocation values
beta_wood = 1.0-(beta_leaf+beta_root)
# -

# Create functions for environmental scaler and input:

# + codehighlighter=[[2, 4]]
#create symbols for scaler and input functions
func_dict={
    'xi': 'Environmental scaler as a function of time',
    'NPP': 'Inputs as a function of time',
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)
    
#define t as a symbol for time
t=TimeSymbol("t")
# -

# ## Symbolic Model Description (Must Edit)
# Define your model using sympy:

# + codehighlighter=[[11, 14], [19, 28], [33, 52], [11, 14], [19, 28], [33, 52]]
#define model in sympy format
mvs = CMTVS(
    {
        t,
        StateVariableTuple(
            #Define State variables
            Pools._fields
        ),
        InFluxesBySymbol( 
            {   #define input/allocation
                #RecievingPool: Input * Allocation
                c_leaf: NPP(t) * beta_leaf, 
                c_root: NPP(t) * beta_root, 
                c_wood: NPP(t) * beta_wood
            }
        ),
        OutFluxesBySymbol( 
            {   #define fluxes leaving the system
                #Fluxes leaving the system: FluRate * DonorPool * EnvironmentalScaler
                c_lit_cwd: r_c_lit_cwd_rh * c_lit_cwd * xi(t),
                c_lit_met: r_c_lit_met_rh * c_lit_met * xi(t),
                c_lit_str: r_c_lit_str_rh * c_lit_str * xi(t),
                c_lit_mic: r_c_lit_mic_rh * c_lit_mic * xi(t),
                c_soil_met: r_c_soil_met_rh * c_soil_met * xi(t),
                c_soil_str: r_c_soil_str_rh * c_soil_str * xi(t),
                c_soil_mic: r_c_soil_mic_rh * c_soil_mic * xi(t),
                c_soil_slow: r_c_soil_slow_rh * c_soil_slow * xi(t),
                c_soil_passive: r_c_soil_passive_rh * c_soil_passive * xi(t),
            }
        ),
        InternalFluxesBySymbol( 
            {   #define fluxes between pools
                #(Donor pool, recieving pool): FluxRate * DonorPool
                (c_leaf, c_lit_met): r_c_leaf_2_c_lit_met * c_leaf,
                (c_leaf, c_lit_str): r_c_leaf_2_c_lit_str * c_leaf,
                (c_root, c_soil_met): r_c_root_2_c_soil_met * c_root,
                (c_root, c_soil_str): r_c_root_2_c_soil_str * c_root,
                (c_wood, c_lit_cwd): r_c_wood_2_c_lit_cwd * c_wood,
                (c_lit_cwd, c_lit_mic): r_c_lit_cwd_2_c_lit_mic * c_lit_cwd * xi(t),
                (c_lit_cwd, c_soil_slow): r_c_lit_cwd_2_c_soil_slow * c_lit_cwd * xi(t),
                (c_lit_met, c_lit_mic): r_c_lit_met_2_c_lit_mic * c_lit_met * xi(t),
                (c_lit_str, c_lit_mic): r_c_lit_str_2_c_lit_mic * c_lit_str * xi(t),
                (c_lit_str, c_soil_slow): r_c_lit_str_2_c_soil_slow * c_lit_str * xi(t),
                (c_lit_mic, c_soil_slow): r_c_lit_mic_2_c_soil_slow * c_lit_mic * xi(t),
                (c_soil_met, c_soil_mic): r_c_soil_met_2_c_soil_mic * c_soil_met * xi(t),
                (c_soil_str, c_soil_mic): r_c_soil_str_2_c_soil_mic * c_soil_str * xi(t),
                (c_soil_str, c_soil_slow): r_c_soil_str_2_c_soil_slow * c_soil_str* xi(t),
                (c_soil_mic, c_soil_slow): r_c_soil_mic_2_c_soil_slow * c_soil_mic* xi(t),
                (c_soil_mic, c_soil_passive): r_c_soil_mic_2_c_soil_passive * c_soil_mic * xi(t),
                (c_soil_slow, c_soil_mic): r_c_soil_slow_2_c_soil_mic * c_soil_slow * xi(t),
                (c_soil_slow, c_soil_passive): r_c_soil_slow_2_c_soil_passive * c_soil_slow * xi(t),
                (c_soil_passive, c_soil_mic): r_c_soil_passive_2_c_soil_mic * c_soil_passive * xi(t),               
            }
        ),
    },
    computers=module_computers(bgc_c)
)
# -

# ## Model Figure and Matrix Equations
# #### Model Figure:

h.compartmental_graph(mvs)

# #### Matrix equations:

dh.mass_balance_equation(mvs)

# + [markdown] codehighlighter=[[0, 0]]
# ## Download Data (Must Edit)
# #### TRENDY Data
# Make sure you have a config.json file in your model folder: <br>
# Config.jon file contents: `{"username": "trendy-v9", "password": "gcb-2020", "dataPath": "/path/to/data/folder"}`

# + codehighlighter=[[11, 12], [16, 17], [8, 28], [41, 43], [8, 24], [42, 44]]
# Read username, password, dataPath from config.json file
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

# Create tuples of data streams
# For YIBs I have annual and monthly file names so I defined them separately
# If all your data is similarly named this would be a single "observables" tuple

# Annual data streams on TRENDY server
observables_annual = namedtuple(
    'observables_annual',
    ["cVeg", "cSoil"]
)
# Monthly data streams on TRENDY server
observables_monthly = namedtuple(
    'observables_monthly', 
    ["rh"]

)
# Combine tuples to request data from TRENDY server
Observables = namedtuple(
    "Observables",
    observables_annual._fields+observables_monthly._fields
)
# Driver data streams on TRENDY server
drivers=namedtuple(
    "Drivers",
    ["npp"]
) 

#when downloading data make sure model names match TRENDY server names:
#"CABLE-POP","CLASSIC","CLM5","DLEM","IBIS","ISAM","ISBA_CTRIP",
#"JSBACH","JULES-ES","LPJ-GUESS","LPJwsl","LPX-Bern","OCN",
#"ORCHIDEE","ORCHIDEE-CNP","ORCHIDEEv3","ORCHIDEEv3_0.5deg"
#"SDGVM","VISIT","YIBs"

# Define function to download data from TRENDY server
def download_my_TRENDY_output():
    download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['YIBs'],
        variables =observables._fields + drivers._fields
    )

# call above function to download TRENDY data
# This can takes a minute to connect to the server
# The download itself can take hours depending on internet speed.
# If "can't connect to TRENDY server" error occurs, try again later.
#download_my_TRENDY_output()
# -

# ## Connect Data and Symbols (Must Edit)
# Define function to subset netCDF files and link to data symbols:

# + codehighlighter=[[5, 6], [23, 33], [5, 6], [23, 33]]
# Define function to download, scale, map TRENDY data
def get_example_site_vars(dataPath):

    # Define single geospatial cell 
    s = slice(None, None, None)  # this is the same as :
    t = s, 74, 118  # [t] = [:,49,325]

    # Define function to select geospatial cell and scale data
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them 
        ds = nc.Dataset(str(path))
        #check for npp/gpp/rh/ra to convert from kg/m2/s to kg/m2/day
        if vn in ["npp","gpp","rh","ra"]:
            for name, variable in ds.variables.items():            
                for attrname in variable.ncattrs():
                    print("{} -- {}".format(attrname, getattr(variable, attrname)))
            return ds.variables[vn][t]
        else:
            for name, variable in ds.variables.items():            
                for attrname in variable.ncattrs():
                    print("{} -- {}".format(attrname, getattr(variable, attrname)))
            return ds.variables[vn][t]

    # Link symbols and data:
    # YIBS has annual vs monthly file names so they are linked separately
    # If all your data is similarly named you can do this in one step

    # Create annual file names (single step if files similarly named)
    o_names=[(f,"YIBs_S2_Annual_{}.nc".format(f)) for f in observables_annual._fields]

    # Create monthly file names (can remove if done in one step above)
    monthly_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in observables_monthly._fields]
    # Extend name list with monthly names
    o_names.extend(monthly_names)

    # create file names for drivers
    d_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in drivers._fields]

    # Link symbols and data for observables/drivers
    return (Observables(*map(f, o_names)),drivers(*map(f,d_names)))

#call function to link symbols and data
svs,dvs=get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))

# + codehighlighter=[[5, 6], [23, 33], [5, 6], [23, 33]]
#look at data
dvs.npp.shape,svs.cSoil.shape


# -

# ## Create Symbols for $\xi$, $K$, and $A$ (No Edits)
# Setup Yiqi matrix format:

sv=mvs.get_StateVariableTuple()                            # Get state variables from symbolic matrix code
n=len(sv)                                                  # Count number of pools
srm = mvs.get_SmoothReservoirModel()                       # Create smooth resevoir model
_,A,N,_,_=srm.xi_T_N_u_representation(factor_out_xi=False) # Create A and N symbols

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
par_dict
# -

import test_helpers as th
import model_specific_helpers_2 as msh
par_dict=th.make_test_args(conf_dict,msh,mvs)

svs_0=Observables(*map(lambda v: v[0],svs))
svs_0

# ## Define Forward Model
# #### Create constants for forward sim:

# + codehighlighter=[[1, 9], [1, 8]]
cpa = Constants(             #Define the constant values of parameters NOT affected by data assimilation
    npp_0 = dvs.npp[0],
    rh_0 = svs.rh[0],   
    c_veg_0 = svs_0.cVeg,
    c_soil_0 = svs_0.cSoil,
    clay = 0.2028,
    silt = 0.2808,
    nyears = 320
    #nyears = 10
)
cpa._asdict()    #print - everything should have a numeric value
# -

# #### Create list of parameters to be optimized during data assimilation:

EstimatedParameters = namedtuple(
    'EstimatedParameters', 
    [
        'c_leaf_0',               #Names: c_poolname_0
        'c_root_0',               #Only initial pools that are estimated
        'c_lit_cwd_0',
        'c_lit_met_0',
        'c_lit_str_0',
        'c_lit_mic_0',
        'c_soil_met_0',
        'c_soil_str_0',
        'c_soil_mic_0',
        'c_soil_slow_0',
        'beta_leaf',
        'beta_root',
        'r_c_leaf_rh',
        'r_c_root_rh',
        'r_c_wood_rh',
        'r_c_lit_cwd_rh',
        'r_c_lit_met_rh',
        'r_c_lit_str_rh',
        'r_c_lit_mic_rh',
        'r_c_soil_met_rh',
        'r_c_soil_str_rh',
        'r_c_soil_mic_rh',
        'r_c_soil_slow_rh',
        'r_c_soil_passive_rh',
        'r_c_leaf_2_c_lit_met',
        'r_c_leaf_2_c_lit_str',
        'r_c_root_2_c_soil_met',
        'r_c_root_2_c_soil_str',
        'r_c_wood_2_c_lit_cwd',
        'r_c_lit_cwd_2_c_lit_mic',
        'r_c_lit_cwd_2_c_soil_slow',
        'r_c_lit_met_2_c_lit_mic',
        'r_c_lit_str_2_c_lit_mic',
        'r_c_lit_str_2_c_soil_slow',
        'r_c_lit_mic_2_c_soil_slow',
        'r_c_soil_met_2_c_soil_mic',
        'r_c_soil_str_2_c_soil_mic',
        'r_c_soil_str_2_c_soil_slow',
        'r_c_soil_mic_2_c_soil_slow',
        'r_c_soil_mic_2_c_soil_passive',
        'r_c_soil_slow_2_c_soil_mic',
        'r_c_soil_slow_2_c_soil_passive',
        'r_c_soil_passive_2_c_soil_mic'
    ]
)


epa_0 = EstimatedParameters(
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

n_days=cpa.nyears*36#0
days=range(n_days)
npp_obs = np.array([npp_func(d) for d in days])

# Plot simulation output for observables
fig = plt.figure(figsize=(12,4))
plot_solutions(
        fig,
        times=days,
        var_names=["npp"],
        tup=(npp_obs,)
)
fig.savefig('solutions.pdf')
# -

# #### Create forward model function:

# + codehighlighter=[[37, 51], [67, 69], [64, 65], [137, 139], [133, 135], [32, 45], [112, 113], [117, 118], [120, 123]]
from general_helpers import numfunc
# Create namedtuple function for initial values
StartVector=namedtuple(
    "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()]+["rh"]
)
def make_iterator_sym(
        mvs,
        V_init: StartVector,
        par_dict,
        func_dict,
        delta_t_val=1 # defaults to 1day timestep
    ):
    B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict,delta_t_val)  
    sv=mvs.get_StateVariableTuple()
    n=len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct 
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr=np.array(
        [V_init.__getattribute__(str(v)) for v in sv]+
        [V_init.rh]
    ).reshape(n+1,1) #reshaping is neccessary for matmul (the @ in B @ X)

    #numOutFluxesBySymbol={
    #    k:numfunc(expr_cont,delta_t_val=delta_t_val) 
    #    for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
    #} 
    numOutFluxesBySymbol={
        k: numfunc(expr_cont, mvs, delta_t_val=delta_t_val, par_dict=par_dict, func_dict=func_dict) 
    
        for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
    } 

    def f(it,V):
        X = V[0:n]
        b = u_func(it,X)
        B = B_func(it,X)
        X_new = X + b + B @ X
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
        
        l=[
                numOutFluxesBySymbol[Symbol(k)](it,*X.reshape(n,))
                for k in ["c_lit_cwd","c_lit_met","c_lit_str","c_lit_mic","c_soil_met","c_soil_str","c_soil_mic","c_soil_slow","c_soil_passive"] 
                if Symbol(k) in numOutFluxesBySymbol.keys()
        ]
        rh = np.array(l).sum()
        V_new = np.concatenate(
            (
                X_new.reshape(n,1),
                np.array([rh]).reshape(1,1)
            )
            , axis=0
        )
        return V_new
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


# + codehighlighter=[[37, 51], [67, 69], [64, 65], [137, 139], [133, 135], [32, 45], [112, 113], [117, 118], [120, 123]]
def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: drivers,
    ) -> Callable[[np.ndarray], np.ndarray]:
    
    # Build iterator 
    # Need dictionary of numeric values for all parameters that are not state variables/time
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    
    # Time dependent driver function does not change with the estimated parameters
    # Defined once outside param2res function
    #seconds_per_day = 86400
    
    def npp_func(day):
        month=day_2_month_index(day)
        return dvs.npp[month]*86400
    
    # Build environmental scaler function
    def xi_func(day):
        return 1.0     # Set to 1 if no scaler implemented 

    # Define function dictionary
    func_dict={
        'NPP':npp_func,
        'xi':xi_func
    }
    
    
    # Define actual forward simulation function
    def param2res(pa):
        
        # Parameter vector
        epa=EstimatedParameters(*pa)
        
        # Create a startvector for the iterator 
        V_init = StartVector(
            c_leaf=epa.c_leaf_0,
            c_root=epa.c_root_0,
            c_wood=cpa.c_veg_0-(epa.c_leaf_0 + epa.c_root_0),
            c_lit_cwd=epa.c_lit_cwd_0,
            c_lit_met=epa.c_lit_met_0,
            c_lit_str=epa.c_lit_str_0,
            c_lit_mic=epa.c_lit_mic_0,
            c_soil_met=epa.c_soil_met_0,
            c_soil_str=epa.c_soil_str_0,
            c_soil_mic=epa.c_soil_mic_0,
            c_soil_slow=epa.c_soil_slow_0,
            c_soil_passive=cpa.c_soil_0-(epa.c_lit_cwd_0+epa.c_lit_met_0+epa.c_lit_str_0+epa.c_lit_mic_0+epa.c_soil_met_0+epa.c_soil_str_0+epa.c_soil_mic_0+epa.c_soil_slow_0),
            rh=cpa.rh_0
        )
        
        # Parameter dictionary for the iterator
        apa = {**cpa._asdict(),**epa._asdict()}
        model_par_dict = {
            Symbol(k):v for k,v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }
        
        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2) 
        delta_t_val=10 #time step length in days
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
        )
        ######################################################################
        # calling the iterator and projecting the results
        def cVegF(V):
            return float(V.c_leaf+V.c_wood+V.c_root)
        
        def cSoilF(V): 
            return float(
                V.c_lit_cwd+
                V.c_lit_met+
                V.c_lit_str+
                V.c_lit_mic+ 
                V.c_soil_met+
                V.c_soil_str+
                V.c_soil_mic+
                V.c_soil_slow+
                V.c_soil_passive
            )
        
        #empty arrays for saving data
        cVeg_arr=np.zeros(cpa.nyears)
        cSoil_arr=np.zeros(cpa.nyears)
        rh_arr=np.zeros(cpa.nyears*12)
        
        dpm = 30                             # Set days per month
        im = 0
        steps_per_month=int(dpm/delta_t_val)
        steps_per_year=int(dpm/delta_t_val)*12
        for y in range(0,cpa.nyears):
            # reset the yearly average values to 0
            cVeg_avg= 0  
            cSoil_avg = 0
            for m in range(12):
                rh_avg=0
                for d in range(steps_per_month):    
                    V = StartVector(*it_sym.__next__())                  
                    rh_avg=V.rh/steps_per_month
                    cVeg_avg += cVegF(V)/steps_per_year
                    cSoil_avg += cSoilF(V)/steps_per_year
                rh_arr[im] = rh_avg/(24*60*60) #convert to kg/s from kg/day
                im+=1
                
            cVeg_arr[y] = cVeg_avg
            cSoil_arr[y] = cSoil_avg
        return Observables(cVeg=cVeg_arr,cSoil=cSoil_arr,rh=rh_arr)
    return param2res
# -

# ## Forward Model Run
# #### Run model forward:

# +
param2res= make_param2res_sym(mvs,cpa,dvs) # Define forward model
obs_simu = param2res(epa_0)                # Run forward model from initial conditions

obs_simu.cVeg.shape,obs_simu.rh.shape
# -

svs_cut=Observables(
    cVeg=svs.cVeg[:cpa.nyears],
    cSoil=svs.cSoil[:cpa.nyears],
    rh=svs.rh[:cpa.nyears*12]
)
svs_cut.rh


# +
def make_weighted_cost_func(
        obs: Observables
    ) -> Callable[[Observables],np.float64]:
    # first unpack the observation array into its parts
    #cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
    def costfunction(out_simu: Observables) -> np.float64:
        # fixme 
        #   as indicated by the fact that the function lives in this  
        #   model-specific module it is NOT apropriate for (all) other models.
        #   There are model specific properties:
        #   1.) The weight for the different observation streams
        #   
        number_of_ys=out_simu.cVeg.shape[0]
        number_of_ms=out_simu.rh.shape[0]

        J_obj1 = np.mean (( out_simu.cVeg - obs.cVeg )**2)/(2*np.var(obs.cVeg))
        J_obj2 = np.mean (( out_simu.cSoil -  obs.cSoil )**2)/(2*np.var(obs.cSoil))
        
        J_obj3 = np.mean (( out_simu.rh - obs.rh )**2)/(2*np.var(obs.rh))
        
        J_new = (J_obj1 + J_obj2)+ J_obj3/12 #the 12 is purely conjectural
        return J_new
    return costfunction     

#test it
cost_func = make_weighted_cost_func(svs_cut)
cost_func(obs_simu)

#for f in Observables._fields:
#    print(obs_0._asdict()[f])
# -

# #### Create array of yearly observation data:

fig = plt.figure()
from general_helpers import plot_observations_vs_simulations
plot_observations_vs_simulations(fig,svs_cut,obs_simu)

len(svs_cut)

# #### Plot data-model fit:

# Plot simulation output for observables
n_plots=len(svs_cut)
fig = plt.figure(figsize=(10, n_plots*5))
axs=fig.subplots(n_plots)
for i,name in enumerate(Observables._fields):
    var=svs_cut[i]
    var_simu=obs_simu[i]
    axs[i].plot(range(len(var_simu)),var_simu,label="simulation")
    axs[i].plot(range(len(var)),var,label='observation')
    axs[i].legend()
    axs[i].set_title(name)


# ## Data Assimilation
# #### Define parameter min/max values:

# +
# set min/max parameters to +- 100 times initial values
epa_min=EstimatedParameters._make(tuple(np.array(epa_0)*0.01))
epa_max=EstimatedParameters._make(tuple(np.array(epa_0)*100))

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
    initial_parameters=epa_0,
    filter_func=make_param_filter_func(epa_max, epa_min),
    param2res=make_param2res_sym(mvs,cpa,dvs),
    costfunction=make_weighted_cost_func(svs),
    nsimu=200, # for testing and tuning mcmc
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
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(EstimatedParameters._fields),1),axis=1)
epa_opt=EstimatedParameters(*par_opt)
mod_opt = param2res(epa_opt)  

print("Forward run with optimized parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(10, n_plots*5))
# Plot simulation output for observables
#n_plots=len(svs_cut)
n_plots=len(svs)
fig = plt.figure(figsize=(10,10*n_plots))
axs=fig.subplots(n_plots)
for i,name in enumerate(Observables._fields):
    var_simu=mod_opt[i]
    var=svs[i]
    axs[i].plot(range(len(var_simu)),var_simu,label='simulation(opt)')
    axs[i].plot(range(len(var)),var,label='observation')
    #axs[i].set_title(name)
    #axs[i].legend()

fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

import pandas as pd
pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('YIBs_da_pars.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('YIBS_da_cost.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('YIBs_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('YIBs_optimized_solutions.csv'), sep=',')
# -
mod_opt.rh
