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
    make_jon_cost_func
)
from model_specific_helpers import (
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
        'c_root_0',
        'c_wood_0',
        'c_lit_cwd_0',
        'c_lit_met_0',
        'c_lit_str_0',
        'c_lit_mic_0',
        'c_soil_met_0',
        'c_soil_str_0',
        'c_soil_mic_0',
        'c_soil_slow_0',
        'c_soil_passive_0'
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
observables = namedtuple(
    "observables",
    observables_annual._fields+observables_monthly._fields
)
# Driver data streams on TRENDY server
drivers=namedtuple(
    "drivers",
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
download_my_TRENDY_output()
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
        if vn in ["rh","ra"]:
            return ds.variables[vn][t]*86400
        else:
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
    return (observables(*map(f, o_names)),drivers(*map(f,d_names)))

#call function to link symbols and data
svs,dvs=get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))

#look at data
svs,dvs
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
outflux_rates = {"r_"+str(key)+"_rh":value/key for key,value in hr.out_fluxes_by_symbol(sv,M_sym).items()}
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
    beta_leaf=0.21,
    beta_root=0.27,
    clay=0.2028,
    silt=0.2808,
    k_leaf=0.014,
    k_root=0.022,
    k_wood=0.003,
    k_cwd=0.005,
    k_samet=0.01,
    k_sastr=0.001,
    k_samic=0.05,
    k_slmet=0.040,
    k_slstr=0.0039,
    k_slmic=0.005,
    k_slow=0.00001,
    k_arm=3.27E-06,
    f_samet_leaf=0.28,
    f_slmet_root=0.34,
    f_samic_cwd=0.29,
)

old_par_dict = {
    'k_c_leaf': 0.014, # define all k values
    'k_c_wood': 0.003,
    'k_c_root': 0.022,
    'k_c_lit_cwd': 0.005,
    'k_c_lit_met': 0.01,
    'k_c_lit_str': 0.001,
    'k_c_lit_mic': 0.05,
    'k_c_soil_met': 0.040,
    'k_c_soil_str': 0.0039,
    'k_c_soil_mic': 0.005,
    'k_c_soil_slow': 0.00001,
    'k_c_soil_passive': 3.27E-06,
    'f_c_leaf_2_c_lit_met': epa0.f_samet_leaf,    #define all f values
    'f_c_root_2_c_soil_met': epa0.f_slmet_root, 
    'f_c_lit_cwd_2_c_lit_mic': epa0.f_samic_cwd,
    'f_c_leaf_2_c_lit_str': (1-epa0.f_samet_leaf) * epa0.k_leaf,
    'f_c_root_2_c_soil_str': (1-epa0.f_slmet_root) * epa0.k_root,
    'f_c_wood_2_c_lit_cwd': 0.4 * epa0.k_wood,
    'f_c_lit_cwd_2_c_soil_slow': (1-epa0.f_samic_cwd) * epa0.k_cwd,
    'f_c_lit_met_2_c_lit_mic': 0.3 * epa0.k_samet,
    'f_c_lit_str_2_c_lit_mic': 0.1 * epa0.k_sastr,
    'f_c_lit_str_2_c_soil_slow': 0.1 * epa0.k_sastr,
    'f_c_lit_mic_2_c_soil_slow': 0.1 * epa0.k_samic,
    'f_c_soil_met_2_c_soil_mic': 0.4 * epa0.k_slmet,
    'f_c_soil_str_2_c_soil_mic': 0.3 * epa0.k_slstr,
    'f_c_soil_str_2_c_soil_slow': 0.2 * epa0.k_slstr,
    'f_c_soil_mic_2_c_soil_slow': 0.4 * epa0.k_slmic,
    'f_c_soil_mic_2_c_soil_passive': 0.4 * epa0.k_slmic,
    'f_c_soil_slow_2_c_soil_mic': 0.10 * epa0.k_slow,
    'f_c_soil_slow_2_c_soil_passive': 0.45*(0.003+0.009*epa0.clay) * epa0.k_slow,
    'f_c_soil_passive_2_c_soil_mic': 0.10 * epa0.k_arm, 
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
makeTuple = namedtuple('makeTuple', par_dict)
parameters = makeTuple(**par_dict)

#If symbols remain in output below then set them to numerical values in old_par_dict.
parameters._asdict() # print - everything below should have a numeric value
# -

# ## Assign Initial Pool Values

# + codehighlighter=[[5, 17], [5, 17]]
# Create vector of initial pool values
svs_0=observables(*map(lambda v: v[0],svs))

# Assign values to initial pools using InitialPools named tuple
V_init = InitialPools(
    c_leaf_0 = svs_0.cVeg/3,          #set inital pool values to svs values 
    c_root_0 = svs_0.cVeg/3,          #you can set numerical values here directly as well
    c_wood_0 = svs_0.cVeg/3,
    c_lit_cwd_0 = svs_0.cSoil/9,
    c_lit_met_0 = svs_0.cSoil/9,
    c_lit_str_0 = svs_0.cSoil/9,
    c_lit_mic_0 = svs_0.cSoil/9,
    c_soil_met_0 = svs_0.cSoil/9,
    c_soil_str_0 = svs_0.cSoil/9,
    c_soil_mic_0 = svs_0.cSoil/9,
    c_soil_slow_0 = svs_0.cSoil/9,
    c_soil_passive_0 = svs_0.cSoil/9
)
V_init._asdict()   #print - everything should have a numeric value
# -

# ## Define Forward Model
# #### Create constants for forward sim:

# + codehighlighter=[[1, 9], [1, 8]]
cpa = Constants(             #use Constants namedtuple to define constant values
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

# #### Create list of parameters to be optimized during data assimilation:

# +
EstimatedInitPools = namedtuple(
    "EstimatedInitPools",
        [
            'c_leaf_0',               
            'c_root_0',
            'c_lit_cwd_0',
            'c_lit_met_0',
            'c_lit_str_0',
            'c_lit_mic_0',
            'c_soil_met_0',
            'c_soil_str_0',
            'c_soil_mic_0',
            'c_soil_slow_0'
        ]
)

X_init = EstimatedInitPools(
    c_leaf_0 = svs_0.cVeg/3,        
    c_root_0 = svs_0.cVeg/3,          
    c_lit_cwd_0 = svs_0.cSoil/9,
    c_lit_met_0 = svs_0.cSoil/9,
    c_lit_str_0 = svs_0.cSoil/9,
    c_lit_mic_0 = svs_0.cSoil/9,
    c_soil_met_0 = svs_0.cSoil/9,
    c_soil_str_0 = svs_0.cSoil/9,
    c_soil_mic_0 = svs_0.cSoil/9,
    c_soil_slow_0 = svs_0.cSoil/9
)

estimated = {**parameters._asdict(),**V_init._asdict()}            # Create dictionary of parameters and initial pools
EstimatedParameters = namedtuple('EstimatedParameters', estimated) # Create function to convert dictionary to namedtuple
epa0 = EstimatedParameters(**estimated)                            # Create namedtuple of all parameters optimized an initial values
epa0._asdict()   #print


# -

# #### Create forward model function:

# + codehighlighter=[[37, 51], [67, 69], [64, 65], [137, 139], [133, 135], [32, 45], [112, 113], [117, 118], [120, 123]]
def make_param2res_sym(
        cpa: Constants
    ) -> Callable[[np.ndarray], np.ndarray]:
    
    # Build iterator 
    # Need dictionary of numeric values for all parameters that are not state variables/time
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    # Create namedtuple function for initial values
    StartVector=namedtuple(
        "StartVector",
            [str(v) for v in mvs.get_StateVariableTuple()]+["rh"]
    )
    
    # Time dependent driver function does not change with the estimated parameters
    # Defined once outside param2res function
    seconds_per_day = 86400
    
    def npp_func(day):
        month=day_2_month_index(day)
        return dvs.npp[month]
    
    # Build environmental scaler function
    def xi_func(day):
        return 1.0     # Set to 1 if no scaler implemented 

    # Define function dictionary
    func_dict={
        'NPP':npp_func,
        'xi':xi_func
    }
    
    def numfunc(expr_cont,delta_t_val):
        # build the discrete expression (which depends on it,delta_t instead of
        # the continius one that depends on t (TimeSymbol))
        it=Symbol("it")           #arbitrary symbol for the step index )
        t=mvs.get_TimeSymbol()
        delta_t=Symbol('delta_t')
        expr_disc = expr_cont.subs({t:delta_t*it})
        return hr.numerical_function_from_expression(
            expr=expr_disc.subs({delta_t:delta_t_val}),
            tup=(it, *mvs.get_StateVariableTuple()),
            parameter_dict=par_dict,
            func_set=func_dict
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

        numOutFluxesBySymbol={
            k:numfunc(expr_cont,delta_t_val=delta_t_val) 
            for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
        } 
        def f(it,V):
            X = V[0:n]
            b = u_func(it,X)
            B = B_func(it,X)
            X_new = X + b + B @ X
            # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
            rh = np.sum(
                [
                    numOutFluxesBySymbol[k](it,*X)
                    for k in [c_lit_cwd,c_lit_met,c_lit_str,c_lit_mic,c_soil_met,c_soil_str,c_soil_mic,c_soil_slow,c_soil_passive] 
                    if k in numOutFluxesBySymbol.keys()
                ]
            )
            V_new = np.concatenate((X_new.reshape(n,1),np.array([rh]).reshape(1,1)), axis=0)
            return V_new
        return TimeStepIterator2(
            initial_values=V_arr,
            f=f,
        )
    
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
            k:v for k,v in apa.items()
            if k in model_par_dict_keys
        }
    
        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2) 
        delta_t_val=30 
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
        )
        
        # Now that we have the iterator we can start to compute.
        # Note: check if TRENDY months are like this...
        #days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        
        #empty array for saving data
        sols=[]
        for m in range(cpa.nyears*12):           # Loop through months
            dpm = 30                             # Set days for each month
            mrh = 0                              # Start respiration sum at zero each month
            for d in range(int(dpm/delta_t_val)):    # Loop through days in month
                v = it_sym.__next__()                # Update initial vector each day
            V = StartVector(*v)                    
            o = observables(                   
                cVeg = float(V.c_leaf+V.c_wood+V.c_root),
                cSoil = float(V.c_lit_cwd+V.c_lit_met+V.c_lit_str+V.c_lit_mic+ \
                            V.c_soil_met+V.c_soil_str+V.c_soil_mic+V.c_soil_slow+V.c_soil_passive),
                rh = v[12,0]
            )
            sols.append(o) # Append monthly value to results
        sol=np.stack(sols) # Stack arrays 
        # convert to yearly output
        sol_yr=np.zeros(cpa.nyears*sol.shape[1]).reshape([cpa.nyears,sol.shape[1]])  
        for i in range(sol.shape[1]):
           sol_yr[:,i]=monthly_to_yearly(sol[:,i])
        sol=sol_yr
        return sol
    return param2res


# -

# ## Forward Model Run
# #### Run model forward:

param2res_sym = make_param2res_sym(cpa) # Define forward model
xs = param2res_sym(epa0)                # Run forward model from initial conditions
xs

# #### Create array of yearly observation data:

n = cpa.nyears                                   # define number of years
obs = np.zeros(n*3).reshape([n,3])               # create empty yearly dataframe 
obs[:,0] = svs.cVeg                              # add yearly cVeg data to obs data
obs[:,1] = svs.cSoil                             # add yearly cSoil data to obs data
obs[:,2] = monthly_to_yearly(svs.rh)             # convert rh to yearly data 
obs

# #### Plot data-model fit:

# Plot simulation output for observables
fig = plt.figure()
plot_solutions(
        fig,
        times=range(n),
        var_names=observables._fields,
        tup=(xs,obs)
)
fig.savefig('solutions.pdf')

# ## Data Assimilation
# #### Define parameter min/max values:

# +
# set min/max parameters to +- 100 times initial values
epa_min=EstimatedParameters._make(tuple(np.array(epa0)*0.01))
epa_max=EstimatedParameters._make(tuple(np.array(epa0)*100))

# fix values that are problematic from calculation
epa_max = epa_max._replace(beta_leaf = 0.99)
epa_max = epa_max._replace(beta_root = 0.99)
epa_max = epa_max._replace(c_leaf_0 = svs_0.cVeg)
epa_max = epa_max._replace(c_root_0 = svs_0.cVeg)
epa_max = epa_max._replace(c_wood_0 = svs_0.cVeg)
epa_max = epa_max._replace(c_lit_cwd_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_lit_met_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_lit_str_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_lit_mic_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_soil_met_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_soil_str_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_soil_mic_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_soil_slow_0 = svs_0.cSoil)
epa_max = epa_max._replace(c_soil_passive_0 = svs_0.cSoil)

#print - all names should have numerical values
epa_max._asdict()
# -

# #### Conduct data assimilation:

# +
isQualified = make_param_filter_func(epa_max, epa_min, obs[0,0], obs[0,1])
param2res = make_param2res_sym(cpa)
costfunction = make_jon_cost_func(obs)

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa0,
    filter_func=isQualified,
    param2res=param2res,
    costfunction=costfunction,
    nsimu=20, # for testing and tuning mcmc
    #nsimu=20000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=15,   # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=1,   # default value | increase value to reduce initial step size
    K=2 # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")
# -

# #### Graph data assimilation results:

# +
# optimized parameter set (lowest cost function)
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(EstimatedParameters._fields),1),axis=1)
epa_opt=EstimatedParameters(*par_opt)
mod_opt = param2res(epa_opt)  

print("Forward run with optimized parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(12, 4), dpi=80)
plot_solutions(
        fig,
        #times=range(cpa.number_of_months),
        times=range(int(cpa.nyears)), # for yearly output
        var_names=observables._fields,
        tup=(mod_opt,obs)
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


