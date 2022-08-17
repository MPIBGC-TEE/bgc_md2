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
#display(HTML("<style>.container { width:100% !important; }</style>"))

# set auto reload for notebook
# %load_ext autoreload
# %autoreload 2

# +
# Packages for symbolic code:
from sympy import Symbol, Function
from ComputabilityGraphs.CMTVS import CMTVS
from sympy import Symbol, Function, diag, ImmutableMatrix

from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
import CompartmentalSystems.helpers_reservoir as hr
import bgc_md2.resolve.computers as bgc_c
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h
from collections import namedtuple

# Other packages
import sys

sys.path.insert(0, '..')  # necessary to import general_helpers
from general_helpers import (
    download_TRENDY_output,
    day_2_month_index,
    # month_2_day_index,
    make_B_u_funcs_2,
    monthly_to_yearly,
    plot_solutions
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

# Make a small dictionary for the variables we will use
# sym_dict={
#    # "vl": "vegegation leaf pool",
#    # "vl": "vegegation leaf pool",
#    # "vw": "vegetation wood pool content",
#    # "I_vw": "Influx into vegetation wood pool",
#    # "k_vw_o": "out flux rate of wood pool",
#    # "k_vl_2_vw": "internal flux rate from leaf to wood",
#    # "k_vw_2_vl": "internal flux rate from wood to leaf",
#    "vl": "vegegation leaf pool",
#    "vw": "vegetation wood pool content",
#    "r_vl_2_vw": "internal flux rate from leaf to wood",
#    "r_vw_2_vl": "internal flux rate from wood to leaf",
#    'C_soil_fast': '',
#    'C_soil_slow': '',
#    'C_soil_passive': '',
#    'C_leaf': '',
#    'C_root': '',
#    'C_wood': '',
#    'C_leaf_litter': '',
#    'C_root_litter': '',
#    'C_wood_litter': '',
#    'r_C_leaf_2_C_leaf_litter': '',
#    'r_C_root_2_C_root_litter': '',
#    'r_C_wood_2_C_wood_litter': '',
#    'r_C_leaf_litter_rh': '',
#    'r_C_root_litter_rh': '',
#    'r_C_wood_litter_rh': '',
#    'r_C_soil_fast_rh': '',
#    'r_C_soil_slow_rh': '',
#    'r_C_soil_passive_rh': '',
#    'r_C_leaf_litter_2_C_soil_fast': '',
#    'r_C_leaf_litter_2_C_soil_slow': '',
#    'r_C_leaf_litter_2_C_soil_passive': '',
#    'r_C_wood_litter_2_C_soil_fast': '',
#    'r_C_wood_litter_2_C_soil_slow': '',
#    'r_C_wood_litter_2_C_soil_passive': '',
#    'r_C_root_litter_2_C_soil_fast': '',
#    'r_C_root_litter_2_C_soil_slow': '',
#    'r_C_root_litter_2_C_soil_passive': '',
#
#    'mrso': 'Total Soil Moisture Content, in kg m-2',
#    'tsl': 'Temperature of Soil - layer, four layers, in K',
#    'Mw': 'soil moisture contetnt at the wilting point',
#    'beta_leaf': 'NPP allocation fraction to leaf',
#    'beta_wood': 'NPP allocation fraction to wood',
#    'beta_root': 'NPP allocation fraction to root',
# }
#
## Create namedtuple of allocation coefficients
# Allocation = namedtuple(
#    "Allocation",
#    [
#        'beta_leaf',  # Names: beta_poolname
#        'beta_wood',
#        'beta_root'
#    ]
# )
#
## Create namedtuple of pools
# Pools = namedtuple(
#    "Pools",
#    [
#        'c_leaf',  # Names: c_poolname
#        'c_wood',
#        'c_root',
#        'c_DPM',  # full name: decomposable plant material
#        'c_RPM',  # full name: resistant plant material
#        'c_BIO',  # full name: microbial biomass
#        'c_HUM'   # full name: long-lived humified
#    ]
# )
#
## Create namedtuple of initial pools
# InitialPools = namedtuple(
#    "InitialPools",
#    [
#        'c_leaf_0',  # Names: c_poolname_0
#        'c_root_0',
#        'c_wood_0',
#        'c_DPM_0',
#        'c_RPM_0',
#        'c_BIO_0',
#        'c_HUM_0'
#    ]
# )
#
## Create namedtuple of flux rates leaving the system
# FluxRates = namedtuple(
#    "FluxRates",
#    [
#        'r_c_DPM_rh',  # Pools with loss from system will be listed here
#        'r_c_RPM_rh',  # Names: r_c_poolname_rh
#        'r_c_BIO_rh',
#        'r_c_HUM_rh',
#        'r_c_leaf_2_c_DPM',  # Pool transfer paths
#        'r_c_leaf_2_c_RPM',  # Names: r_c_donorPool_2_recievingPool
#        'r_c_wood_2_c_DPM',
#        'r_c_wood_2_c_RPM',
#        'r_c_root_2_c_DPM',
#        'r_c_root_2_c_RPM',
#        'r_c_DPM_2_c_BIO',
#        'r_c_DPM_2_c_HUM',
#        'r_c_RPM_2_c_BIO',
#        'r_c_RPM_2_c_HUM',
#        'r_c_BIO_2_c_HUM',
#        'r_c_HUM_2_c_BIO'
#    ]
# )
#
## define time
# Time = namedtuple(
#    "Time",
#    ['t']
# )
#
## Create namedtuple of constants used in model
# Constants = namedtuple(
#    "Constants",
#    [
#        'npp_0',  # Initial input/pools
#        'rh_0',
#        'c_veg_0',
#        'c_soil_0',
#        'f_veg2soil_0',  # Total carbon mass from vegetation directly into the soil
#        'nyears'  # Run time (years for my model)
#    ]
# )
#
## Combine all named tuples to create symbols
# Symbols = namedtuple(
#    "Symbols",
#    Allocation._fields + Pools._fields + InitialPools._fields + \
#    FluxRates._fields + Time._fields + Constants._fields
# )
#
## Create symbols
# for k in Symbols._fields:
#    code = k + " = Symbol('{0}')".format(k)
#    exec(code)
#
## define beta wood from other allocation values
# beta_wood = 1.0 - (beta_leaf + beta_root)
#
## create symbols for scaler and input functions
# func_dict = {
#    'xi': 'Environmental scaler as a function of time',
#    'NPP': 'Inputs as a function of time',
# }
# for k in func_dict.keys():
#    code = k + " = Function('{0}')".format(k)
#    exec(code)
#
## define t as a symbol for time
# t = TimeSymbol("t")
#
### Symbolic Model Description (Must Edit)
# Define your model using sympy:
#
## define model in sympy format
# mvs = CMTVS(
#    {
#        t,
#        StateVariableTuple(
#            # Define State variables
#            Pools._fields
#        ),
#        InFluxesBySymbol(
#            {  # define input/allocation
#                # RecievingPool: Input * Allocation
#                c_leaf: NPP(t) * beta_leaf,
#                c_root: NPP(t) * beta_root,
#                c_wood: NPP(t) * beta_wood
#            }
#        ),
#        OutFluxesBySymbol(
#            {  # define fluxes leaving the system
#                # Fluxes leaving the system: FluRate * DonorPool * EnvironmentalScaler
#                c_DPM: r_c_DPM_rh * c_DPM * xi(t),
#                c_RPM: r_c_RPM_rh * c_RPM * xi(t),
#                c_BIO: r_c_BIO_rh * c_BIO * xi(t),
#                c_HUM: r_c_HUM_rh * c_HUM * xi(t)
#            }
#        ),
#        InternalFluxesBySymbol(
#            {  # define fluxes between pools
#                # (Donor pool, recieving pool): FluxRate * DonorPool
#                (c_leaf, c_DPM): r_c_leaf_2_c_DPM * c_leaf,
#                (c_leaf, c_RPM): r_c_leaf_2_c_RPM * c_leaf,
#                (c_wood, c_DPM): r_c_wood_2_c_DPM * c_wood,
#                (c_wood, c_RPM): r_c_wood_2_c_RPM * c_wood,
#                (c_root, c_DPM): r_c_root_2_c_DPM * c_root,
#                (c_root, c_RPM): r_c_root_2_c_RPM * c_root,
#                (c_DPM,  c_BIO): r_c_DPM_2_c_BIO  * c_DPM  * xi(t),
#                (c_DPM,  c_HUM): r_c_DPM_2_c_HUM  * c_DPM  * xi(t),
#                (c_RPM,  c_BIO): r_c_RPM_2_c_BIO  * c_RPM  * xi(t),
#                (c_RPM,  c_HUM): r_c_RPM_2_c_HUM  * c_RPM  * xi(t),
#                (c_BIO,  c_HUM): r_c_BIO_2_c_HUM  * c_BIO  * xi(t),
#                (c_HUM,  c_BIO): r_c_HUM_2_c_BIO  * c_HUM  * xi(t)
#            }
#        ),
#    },
#    computers=module_computers(bgc_c)
# )
#
# -

from source import mvs

# ## Model Figure and Matrix Equations
# #### Model Figure:

h.compartmental_graph(mvs)

# #### Matrix equations:

# + active=""
# # #### Intermediate summary:
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
# 2. provide code to download the output for your model.
# 3. implement functions for the drivers (using the data)
# 4. run the model forward with a possible set of parameters.
# 5. infer unknown parameters by data assimilation.

# + [markdown] codehighlighter=[[0, 0]]
# ## Download Data (Must Edit)
# #### TRENDY Data
# Make sure you have a config.json file in your model folder: <br>
# Config.jon file contents: `{"username": "trendy-v9", "password": "gcb-2020", "dataPath": "/path/to/data/folder"}`

# + codehighlighter=[[11, 12], [16, 17], [8, 28], [41, 43], [8, 24], [42, 44]]
# import sys
# sys.path.insert(0,'..') # necessary to import general_helpers
# from general_helpers import download_TRENDY_output
# import json
# from pathlib import Path
# from collections import namedtuple
#
## Read username, password, dataPath from config.json file
# with Path('config.json').open(mode='r') as f:
#    conf_dict = json.load(f)
#
## Create tuples of data streams
## For YIBs I have annual and monthly file names so I defined them separately
## If all your data is similarly named this would be a single "observables" tuple
#
#
## Monthly data streams on TRENDY server
# Observables = namedtuple(
#    'observables_monthly',
#    ["cVeg", "cSoil", "rh", "fVegSoil"]
# )
#
## Driver data streams on TRENDY server
# Drivers = namedtuple(
#    "drivers",
#    ["npp","mrso", "tsl","landCoverFrac"]
# )
#
## when downloading data make sure model names match TRENDY server names:
## "CABLE-POP","CLASSIC","CLM5","DLEM","IBIS","ISAM","ISBA_CTRIP",
## "JSBACH","JULES-ES","LPJ-GUESS","LPJwsl","LPX-Bern","OCN",
## "ORCHIDEE","ORCHIDEE-CNP","ORCHIDEEv3","ORCHIDEEv3_0.5deg"
## "SDGVM","VISIT","YIBs"
#
## Define function to download data from TRENDY server
# def download_my_TRENDY_output():
#    download_TRENDY_output(
#        username=conf_dict["username"],
#        password=conf_dict["password"],
#        dataPath=Path(conf_dict["dataPath"]),  # platform independent path desc. (Windows vs. linux)
#        models=["JULES-ES"],
#        variables=Observables._fields + Drivers._fields
#    )
#
## call above function to download TRENDY data
## This can takes a minute to connect to the server
## The download itself can take hours depending on internet speed.
## If "can't connect to TRENDY server" error occurs, try again later.
# download_my_TRENDY_output()

# +
import json

# Read username, password, dataPath from config.json file
with Path('config.json').open(mode='r') as f:
    conf_dict = json.load(f)
import model_specific_helpers_2 as msh
# -

msh.download_my_TRENDY_output(conf_dict)
# #%pwd

# ## Connect Data and Symbols (Must Edit)
# Define function to subset netCDF files and link to data symbols:

# + codehighlighter=[[5, 6], [23, 33], [5, 6], [23, 33]]
# import netCDF4 as nc
# import numpy as np
# from pathlib import Path
# import json
## Define function to download, scale, map TRENDY data
# def get_example_site_vars(dataPath):
#    # Define single geospatial cell from (3840, 144, 192)
#    s = slice(None, None, None)  # this is the same as :
#    t = s, 70, 160  # a site in South America
#
#    # Define function to select geospatial cell and scale data
#    def f(tup):
#        vn, fn = tup
#        path = dataPath.joinpath(fn)
#        # Read NetCDF data but only at the point where we want them
#        ds = nc.Dataset(str(path))
#
#        if vn in ["gpp", "rh", "ra","f_veg2soil"]: # (3840, 144, 192), kg m-2 s-1
#                                                   # f_veg2soil: Total carbon mass from vegetation directly into the soil
#            print("reading ", vn, ", size is ", ds.variables[vn].shape)
#            return ds.variables[vn][t] * 86400 # convert from kg/m2/s to kg/m2/day
#        elif vn in ["npp"]:
#            print("reading ", vn, ", size is ", ds.variables["npp_nlim"].shape)
#            return ds.variables["npp_nlim"][t] * 86400
#        elif vn in ["tsl"]: # Temperature of Soil - layer, 192x144x4 (depth) x3840, 'K'
#            print("reading ", vn, ", size is ", ds.variables[vn].shape)  ## (3840, 4, 144, 192)
#            tmp = np.mean(ds.variables[vn], axis = 1)
#            return  tmp[t] # average soil temperature at different depth
#            print("converted size is ", tmp.shape)
#        elif vn in ["landCoverFrac"]: # Plant Functional Type Grid Fraction, 192x144x17 (vegtype) x3840
#            print("reading ", vn, ", size is ", ds.variables[vn].shape)  ## (3840, 17, 144, 192)
#            tmp = np.sum(ds.variables[vn][:, 0:13, :, :], axis = 1)
#            print("converted size is ", tmp.shape)
#            return tmp[t]  # sum the vegetation coverages
#        else:
#            print("reading ", vn, ", size is ", ds.variables[vn].shape)
#            return ds.variables[vn][t]
#
#    # Link symbols and data:
#
#    # Create file names (single step if files similarly named)
#    o_names = [(f, "JULES-ES-1p0_S2_{}.nc".format(f)) for f in Observables._fields]
#    d_names = [(f, "JULES-ES-1p0_S2_{}.nc".format(f)) for f in Drivers._fields]
#
#    # Link symbols and data for observables/drivers
#    return (Observables(*map(f, o_names)), Drivers(*map(f, d_names)))
#
#
## call function to link symbols and data
#svs, dvs = msh.get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))

svs, dvs = msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
# look at data
# svs, dvs
# -

Path(conf_dict["dataPath"])

# +
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms


fig, ax = plt.subplots() #[(12*10-1):12*30]
plt.plot(list(range(0, 12*320)), np.array(svs.cVeg), label = "observation")

# +
# Plot simulation output for observables
c_obs = np.array(svs.rh) # dvs.npp
c_obs.shape

n_days = 320 * 12  # 0
days = range(n_days)
fig = plt.figure(figsize=(12, 4))
plot_solutions(
    fig,
    times=days,
    var_names=["rh"],
    tup=(c_obs,)
)
# -

# ## Create Symbols for $\xi$, $K$, and $A$ (No Edits)
# Setup Yiqi matrix format:

# +
# sv = mvs.get_StateVariableTuple()  # Get state variables from symbolic matrix code
# n = len(sv)  # Count number of pools
# srm = mvs.get_SmoothReservoirModel()  # Create smooth resevoir model
# _, A, N, _, _ = srm.xi_T_N_u_representation(factor_out_xi=False)  # Create A and N symbols
# t=Symbol("t")

# BI=mvs.get_BibInfo()
# for k in BI.func_dict.keys():
#    code = k + " = Function('{0}')".format(k)
#    exec(code)
#

## Create environmental scaler matrix
## Set as for YIBs model; JULES will add a xi(t) function (Rh_scaler) in the forward modeling
# xi_d=diag([1,1,1]+[xi(t) for i in range(n-3)],unpack=True)
# xi_d

## Create empty K matrix
# K = xi_d.inv() * N
## Create new symbols for the k_{i}
# for i in range(n):
#    if K[i, i] != 0:
#        name = "k_{0}".format(sv[i])
#        code = "{0}=Symbol('{0}')".format(name)
#        # print(code)
#        exec(code)
#
## Create $K$ matrix
# K_sym = ImmutableMatrix(
#    n, n,
#    lambda i, j: Symbol("k_" + str(sv[i])) if i == j else 0
# )
# K_sym

## Create new symbols for the f_{i,j}
# for i in range(n):
#    for j in range(n):
#        if A[i, j] != 0 and i != j:
#            name = "f_" + str(sv[j]) + "_2_" + str(sv[i])
#            code = "{0}=Symbol('{0}')".format(name)
#            # print(code)
#            exec(code)
#
## Place $f$ values in $A$ matrix
# A_sym = ImmutableMatrix(
#    n, n,
#    lambda i, j: -1 if i == j else (
#        0 if A[i, j] == 0 else Symbol("f_" + str(sv[j]) + "_2_" + str(sv[i]))
#    )
# )
# A_sym

## Create A matrix
# M_sym = A_sym * K_sym
# M_sym

## Create a dictionary for the external and internal fluxes (flux divided by dono pool content)
# outflux_rates = {"r_" + str(key) + "_rh": value / key for key, value in hr.out_fluxes_by_symbol(sv, M_sym).items()}
# internal_flux_rates = {"r_" + str(key[0]) + "_2_" + str(key[1]): value / key[0] for key, value in
#                       hr.internal_fluxes_by_symbol(sv, M_sym).items()}
#
# from copy import  deepcopy
## Create dictionary of all flux rates
# all_rates = deepcopy(outflux_rates)
# all_rates.update(internal_flux_rates)
# all_rates

# ParameterValues = namedtuple(
#    "ParameterValues",
#    [
#        "beta_leaf",    # 0 (indices uses in original code)
#        "beta_wood",    # 1
#        "f_leaf2DPM",   # 2, f41
#        # "f_leaf2RPM",   #  f51
#        "f_wood2DPM",   # 3, f42
#        # "f_wood2RPM",   #  f52
#        "f_root2DPM",   # 4, f43
#        # "f_root2RPM",   #  f53
#        "f_DPM2BIO",    # 5, f64
#        "f_DPM2HUM",    # 6, f74
#        "f_RPM2BIO",    # 7, f65
#        "f_RPM2HUM",    # 8, f75
#        "f_BIO2HUM",    # 9, f76
#        "f_HUM2BIO",    # 10, f67
#        "k_leaf",       # 11
#        "k_wood",       # 12
#        "k_root",       # 13
#        "k_DPM",        # 14
#        "k_RPM",        # 15
#        "k_BIO",        # 16
#        "k_HUM",        # 17
#        "c_leaf0",      # 18
#        "c_wood0",      # 19
#        "c_RPM0",       # 20
#        "c_DPM0",       # 21
#        "c_BIO0",       # 22
#        'Mw'            # 23, soil moisture content at the wilting point
#    ]
# )
#
## save a good hand-tuned parameter in case it's needed in the future
## only for save, if run here, got error
epa0 = ParameterValues( ## need to modify!!
    beta_leaf = 0.35,
    beta_wood = 0.3,
    f_leaf2DPM = 0.22,
    f_wood2DPM = 0.22,
    f_root2DPM = 0.22,
    f_DPM2BIO = 0.46*betaR,
    f_DPM2HUM = 0.54*betaR,
    f_RPM2BIO = 0.46*betaR,
    f_RPM2HUM = 0.54*betaR,
    f_BIO2HUM = 0.54*betaR,
    f_HUM2BIO = 0.46*betaR,
    k_leaf = (1 / 4) / (360),
    k_wood = (1 / 30.5)/ (360),
    k_root = (1 / 4) /360,
    k_DPM  = (1 / 0.065) / 360,
    k_RPM = (1 / 2.5) / 360,
    k_BIO = (1 / 0.85) / 360,
    k_HUM = (1 / 30) / 360,


    c_leaf0 = svs.cVeg[0] *0.12, #/ 3,  # set inital pool values to svs values
    c_wood0 = svs.cVeg[0] *0.76, # / 3,
    c_DPM0 = svs.cSoil[0] *0.005,
    c_RPM0 = svs.cSoil[0] *0.22,
    c_BIO0 = svs.cSoil[0] *0.02,

    Mw=0.01
)

betaR = 0.22   # was derived from clay content(RothC model): 1/(x + 1) , where x = 1.67*(1.85 + 1.60exp(-0.0786%clay))
epa0 = ParameterValues( ## need to modify!!
    beta_leaf = 0.35,
    beta_wood = 0.3,
    
    f_leaf2DPM = 0.22,
    f_wood2DPM = 0.22,
    f_root2DPM = 0.22,
    f_DPM2BIO = 0.46*betaR,
    f_DPM2HUM = 0.54*betaR,
    f_RPM2BIO = 0.46*betaR,
    f_RPM2HUM = 0.54*betaR,
    f_BIO2HUM = 0.54*betaR,
    f_HUM2BIO = 0.46*betaR,
    k_leaf = (1 / 4) / (360),
    k_wood = (1 / 30.5)/ (360),
    k_root = (1 / 4) /360,
    k_DPM  = (1 / 0.099) / 360,
    k_RPM = (1 / 2.5) / 360,
    k_BIO = (1 / 1.4) / 360,
    k_HUM = (1 / 42) / 360,


    c_leaf0 = svs.cVeg[0] *0.12, #/ 3,  # set inital pool values to svs values
    c_wood0 = svs.cVeg[0] *0.76, # / 3,
    c_DPM0 = svs.cSoil[0] *0.002,
    c_RPM0 = svs.cSoil[0] *0.22,
    c_BIO0 = svs.cSoil[0] *0.02,

    Mw=0.1,
    Ms=np.max(dvs.mrsos) + 500, #, may need add a condition here ## ASK MARKUS
    Topt=18.32,
    Tcons=47.91
)
#
## need to check here!!
# old_par_dict = {
#     # define all k values
#    'k_c_leaf': epa0.k_leaf,
#    'k_c_wood': epa0.k_wood,
#    'k_c_root': epa0.k_root,
#    'k_c_DPM':  epa0.k_DPM,
#    'k_c_RPM':  epa0.k_RPM,
#    'k_c_BIO':  epa0.k_BIO,
#    'k_c_HUM':  epa0.k_HUM,
#     # define all f values
#    'f_c_leaf_2_c_DPM': epa0.f_leaf2DPM,
#    'f_c_leaf_2_c_RPM': 1 - epa0.f_leaf2DPM,
#    'f_c_wood_2_c_DPM': epa0.f_wood2DPM,
#    'f_c_wood_2_c_RPM': 1 - epa0.f_wood2DPM,
#    'f_c_root_2_c_DPM': epa0.f_root2DPM,
#    'f_c_root_2_c_RPM': 1 - epa0.f_root2DPM,
#    'f_c_DPM_2_c_BIO': epa0.f_DPM2BIO,
#    'f_c_DPM_2_c_HUM': (1 - epa0.f_DPM2BIO),
#    'f_c_RPM_2_c_BIO': epa0.f_RPM2BIO ,
#    'f_c_RPM_2_c_HUM': (1 - epa0.f_RPM2BIO),
#    'f_c_BIO_2_c_HUM': epa0.f_BIO2HUM ,
#    'f_c_HUM_2_c_BIO': epa0.f_HUM2BIO
# }
#
## Define allocation parameters to be optimized
# par_dict = {
#    'beta_leaf': epa0.beta_leaf,
#    'beta_root': epa0.beta_wood,
# }
#
## translate rates from previous parameters to create dictionary of rates to optimize
# par_dict.update(
#    {str(k):v.subs(old_par_dict) for k,v in all_rates.items()}
# )

# + codehighlighter=[[3, 22], [26, 45], [48, 79], [83, 85], [3, 22], [26, 45], [48, 79], [83, 85]]
par_dict = {
    Symbol(k): v
    for k, v in {
        'beta_leaf': 0.33,
        'beta_wood': 0.34,
        'r_c_DPM_rh': 0.00036551,
        'r_c_RPM_rh': 0.000518173,
        'r_c_BIO_rh': 0.000246165,
        'r_c_HUM_rh': 7.76e-05,
        'r_c_leaf_2_c_DPM': 0.000259686,
        'r_c_leaf_2_c_RPM': 0.006112809,
        'r_c_wood_2_c_DPM': 3.95e-07,
        'r_c_wood_2_c_RPM': 0.000114536,
        'r_c_root_2_c_DPM': 1.88e-07,
        'r_c_root_2_c_RPM': 0.000147164,
        'r_c_DPM_2_c_BIO': 0.00778322,
        'r_c_DPM_2_c_HUM': 0.040097969,
        'r_c_RPM_2_c_BIO': 2.78e-06,
        'r_c_RPM_2_c_HUM': 0.000116877,
        'r_c_BIO_2_c_HUM': 7.68e-05,
        'r_c_HUM_2_c_BIO': 1.98e-07
    }.items()
}
par_dict

# -


for k, v in mvs.get_OutFluxesBySymbol().items():
    display(v)

par_dict

# +
# import model_specific_helpers_2 as msh
# import test_helpers as th
# par_dict=th.make_test_args(conf_dict,msh,mvs)
# -

svs_0 = msh.Observables(*map(lambda v: v[0], svs))
svs_0.fVegSoil

# ## Assign Initial Values for the iterator

# + codehighlighter=[[5, 17], [5, 17]]
# Create vector of initial pool values
svs_0 = msh.Observables(*map(lambda v: v[0], svs))

StartVector = msh.make_StartVector(mvs)
StartVector._fields

# + codehighlighter=[[5, 17], [5, 17]]
# Assign values to initial pools using InitialPools named tuple
V_init = StartVector(
    c_leaf=svs_0.cVeg *0.2, #/ 3,  # set inital pool values to svs values
    c_root=svs_0.cVeg *0.2, #/ 3,  # you can set numerical values here directly as well
    c_wood=svs_0.cVeg *0.6, # / 3,
    c_DPM=svs_0.cSoil *0.0025,
    c_RPM=svs_0.cSoil *0.22,
    c_BIO=svs_0.cSoil *0.02,
    c_HUM=svs_0.cSoil *0.7575,
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

# +
# estimated = {**parameters._asdict(), **V_init._asdict()}  # Create dictionary of parameters and initial pools
# OptimizedParameters = namedtuple('OptimizedParameters',
#                                 estimated)  # Create function to convert dictionary to namedtuple
# epa0 = OptimizedParameters(**estimated)  # Create namedtuple of all parameters optimized an initial values
# epa0._asdict()  # print
# -

svs_0.cSoil

# +
# 'r_c_DPM_rh': k_c_DPM*(-f_c_DPM_2_c_BIO - f_c_DPM_2_c_HUM + 1),
# 'r_c_RPM_rh': k_c_RPM*(-f_c_RPM_2_c_BIO - f_c_RPM_2_c_HUM + 1),
# 'r_c_BIO_rh': k_c_BIO*(1 - f_c_BIO_2_c_HUM),
# 'r_c_HUM_rh': k_c_HUM*(1 - f_c_HUM_2_c_BIO),
# 'r_c_leaf_2_c_DPM': f_c_leaf_2_c_DPM*k_c_leaf,
# 'r_c_leaf_2_c_RPM': f_c_leaf_2_c_RPM*k_c_leaf,
# 'r_c_wood_2_c_DPM': f_c_wood_2_c_DPM*k_c_wood,
# 'r_c_wood_2_c_RPM': f_c_wood_2_c_RPM*k_c_wood,
# 'r_c_root_2_c_DPM': f_c_root_2_c_DPM*k_c_root,
# 'r_c_root_2_c_RPM': f_c_root_2_c_RPM*k_c_root,
# 'r_c_DPM_2_c_BIO': f_c_DPM_2_c_BIO*k_c_DPM,
# 'r_c_DPM_2_c_HUM': f_c_DPM_2_c_HUM*k_c_DPM,
# 'r_c_RPM_2_c_BIO': f_c_RPM_2_c_BIO*k_c_RPM,
# 'r_c_RPM_2_c_HUM': f_c_RPM_2_c_HUM*k_c_RPM,
# 'r_c_BIO_2_c_HUM': f_c_BIO_2_c_HUM*k_c_BIO,
# 'r_c_HUM_2_c_BIO': f_c_HUM_2_c_BIO*k_c_HUM


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
func_dict = msh.make_func_dict(mvs, dvs, cpa, epa_0)
import general_helpers as gh
mvs.get_InternalFluxesBySymbol()

import model_specific_helpers_2 as msh
param2res = msh.make_param2res_sym(mvs, cpa, dvs)  # Define forward model
obs_simu = param2res(epa_0)  # Run forward model from initial conditions

obs_simu

# +
fig = plt.figure()
from general_helpers import plot_observations_vs_simulations
plot_observations_vs_simulations(fig,svs,obs_simu)

svs_0.cSoil * 0.0011

# +
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms


fig, ax = plt.subplots() #[(12*10-1):12*30]
ax.scatter(np.array(svs.rh)*1000, np.array(obs_simu.rh)*1000, c='black') #[(12*10-1):12*30]
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
plt.show()

# +
plt.plot(list(range(0, 12*20)), np.array(svs.rh[(12*300):12*320])*1000, label = "observation")
plt.plot(list(range(0, 12*20)), np.array(obs_simu.rh[(12*300):12*320])*1000, label = "simulation")
#plt.plot(list(range(0, 12*20)), np.array(dvs.npp[(12*300):12*320])*1000, label = "NPP")

#plt.plot(list(range(0, 12*320)), np.array(dvs.npp) - np.array(svs.fVegSoil), label = "net change in veg C")
#plt.plot(list(range(0, 12*320)), np.array(svs.cVeg), label = "net veg C")
plt.legend()
plt.plot()
plt.show()
sum((np.array(svs.rh) - np.array(obs_simu.rh))**2)
# -

plt.plot(list(range(0, 12*320)), np.array(svs.rh)*1000, label = "observation")
#plt.plot(list(range(0, 12*320)), np.array(obs_simu.rh)*1000, label = "simulation")
#plt.plot(list(range(0, 12*20)), np.array(dvs.npp[(12*300):12*320])*1000, label = "NPP")


FT = 47.9 / (1 + np.exp(106/(np.array(dvs.tsl) - 254.85)))
plt.plot(np.array(dvs.tsl), 47.9 / (1 + np.exp(106/(np.array(dvs.tsl) - 254.85))), label = "254.85")
plt.plot(np.array(dvs.tsl), 47.9 / (1 + np.exp(106/(np.array(dvs.tsl) - 264.85))), label = "264.85")


# +
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(list(range(0, 12*20)), np.array(obs_simu.rh)[(12*300):12*320], 'r-')
#FT = 2.0 ** ((np.array(dvs.tsl) - 298.15) / 10)
FT = 45.9 / (1 + np.exp(106/(np.array(dvs.tsl) - 254.85)))
FV = 0.6 + 0.4 * (1 - np.array(dvs.landCoverFrac) / 100)

Mw = 0.01
Ms = np.max(dvs.mrsos) + 50
S0 = 0.5 * (1 + Mw)  # optimum soil moisture
Smin = 1.7 * Mw  # lower threshold soil moisture for soil respiration
mois = np.array(dvs.mrsos)
FS = np.zeros(shape=mois.shape)

for mi in range(len(mois)):
    if S0 < mois[mi]/Ms:
        FS[mi] = 1 - 0.8 * (mois[mi]/Ms - S0)  # effect of soil moisture
    if (Smin < mois[mi]/Ms) and (mois[mi]/Ms <= S0):
        FS[mi] = 0.2 + 0.8 * (mois[mi]/Ms - Smin) / (S0 - Smin)
    if mois[mi]/Ms <= Smin:
        FS[mi] = 0.2


ax2.plot(list(range(0, 12*20)), (FT * FV * FS )[(12*300):12*320], 'b-',)
ax1.set_ylabel('Rh simulation', color='r')
ax2.set_ylabel('Environmental scaler i', color='b')
plt.show()
Ms

# +
Mw = 0.01
Ms = np.max(dvs.mrsos) + 50
S0 = 0.5 * (1 + Mw)  # optimum soil moisture
Smin = 1.7 * Mw  # lower threshold soil moisture for soil respiration
mois =  np.arange(1, 75)
FS = np.zeros(shape= np.array(mois.shape))
for mi in range(len(mois)):
    if S0 < mois[mi]/Ms:
        FS[mi] = 1 - 0.8 * (mois[mi]/Ms - S0)  # effect of soil moisture
    if (Smin < mois[mi]/Ms) and (mois[mi]/Ms <= S0):
        FS[mi] = 0.2 + 0.8 * (mois[mi]/Ms - Smin) / (S0 - Smin)
    if mois[mi]/Ms <= Smin:
        FS[mi] = 0.2
        
        
plt.plot(list(range(1, 75)), F, label = "observation")

# -

lst = np.arange(1,4)  
print(lst)

# +
## fixme: 

#### semi-analytical solution of steady state, haven't implemented this part to proposer yet
#### may need work with Marks and/or Feng
beta_leaf = 0.35
beta_wood = 0.3
f_leaf2DPM = 0.22
f_wood2DPM = 0.22
f_root2DPM = 0.22
betaR = 0.295
f_DPM2BIO = 0.46*betaR
f_DPM2HUM = 0.54*betaR
f_RPM2BIO = 0.46*betaR
f_RPM2HUM = 0.54*betaR
f_BIO2HUM = 0.54*betaR
f_HUM2BIO = 0.46*betaR
k_leaf = (1 / 4) / (360)
k_wood = (1 / 30.5)/ (360)
k_root = (1 / 4) /360
k_DPM  = (1 / 0.065) / 360
k_RPM = (1 / 2.5) / 360
k_BIO = (1 / 0.85) / 360
k_HUM = (1 / 30) / 360

Bmat = np.array([[beta_leaf,         0,                      0,  0,  0,  0,   0],
                  [0,        beta_wood,                      0,  0,  0,  0,   0],
                  [0,                0,  1-beta_leaf-beta_wood,  0,  0,  0,   0],
                  [0,  0,  0, 0,  0,  0,   0],
                  [0,  0,  0,  0, 0,  0,   0],
                  [0,  0,  0,  0,  0, 0,   0],
                  [0,  0,  0,  0,  0,  0,  0]])
k = [k_leaf, k_wood, k_root, k_DPM, k_RPM, k_BIO, k_HUM]
Kmat = np.diag(k)
xi = np.mean(FT * FV * FS)
ximat = np.array([[1,  0,  0,  0,  0,  0,   0],
                  [0,  1,  0,  0,  0,  0,   0],
                  [0,  0,  1,  0,  0,  0,   0],
                  [0,  0,  0, xi,  0,  0,   0],
                  [0,  0,  0,  0, xi,  0,   0],
                  [0,  0,  0,  0,  0, xi,   0],
                  [0,  0,  0,  0,  0,  0,  xi]])
Amat = np.array([[-1,  0,  0,  0,  0,  0,   0], 
              [0,  -1,  0,  0,  0,  0,   0],
              [0,  0,  -1,  0,  0,  0,   0],
              [f_leaf2DPM,      f_wood2DPM,    f_root2DPM,         -1,          0,          0,           0],
              [1-f_leaf2DPM,  1-f_wood2DPM,  1-f_root2DPM,          0,         -1,          0,           0],
              [0,                        0,             0,  f_DPM2BIO,  f_RPM2BIO,         -1,   f_HUM2BIO],
              [0,                        0,             0,  f_DPM2HUM,  f_RPM2HUM,  f_BIO2HUM,         -1]])
U = np.mean(dvs.npp)
Umat = np.diag([U, U, U, U, U, U, U])
# -

ximat * Kmat * Amat * Xcarbonpools

# +
## fixme: 
########### SASU-inferred steady state 

for i in range(10):
    X = - np.linalg.inv(ximat * A * Kmat) * Bmat * Umat
    U = X
   
Xcarbonpools = np.zeros(7)
carbonpools = np.zeros(shape = [14, 400])


for i_spinup in range(0, 1):
    if i_spinup == 1:
        for i_year in range(0, 80):
            Xcarbonpools = Xcarbonpools +  Bmat * U + ximat * Kmat * Amat * Xcarbonpools;
            # carbonpools(:,i_spinup) = Xcarbonpools;
    else:
        for i_year in range(0, 8):
            Xcarbonpools = - np.linalg.inv(ximat * A * Kmat) * Bmat * Xcarbonpools
    
Xcarbonpools

# +
fig, ax = plt.subplots()
ax.scatter(np.array(svs.fVegSoil), np.array(obs_simu.fVegSoil), c='black')
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
plt.show()


# +
def npp_func(day):
    month = day_2_month_index(day)
    return dvs.npp[month]


n_days = cpa.nyears * 12 * 30  # 0
days = range(n_days)
npp_obs = np.array([npp_func(d) for d in days])
npp_obs.shape
# -

# #### Create forward model function:

# + codehighlighter=[[37, 51], [67, 69], [64, 65], [137, 139], [133, 135], [32, 45], [112, 113], [117, 118], [120, 123]]
func_dict = msh.make_func_dict(mvs, dvs, cpa, epa_0)

# import general_helpers as gh
## Create namedtuple function for initial values
# StartVector=namedtuple(
#    "StartVector",
#        [str(v) for v in mvs.get_StateVariableTuple()]+["rh"]
# )
#
# def make_iterator_sym(
#        mvs,
#        V_init: StartVector,
#        par_dict,
#        func_dict,
#        delta_t_val=1 # defaults to 1day timestep
#    ):
#    B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict,delta_t_val)
#    sv=mvs.get_StateVariableTuple()
#    n=len(sv)
#    # build an array in the correct order of the StateVariables which in our case is already correct
#    # since V_init was built automatically but in case somebody handcrafted it and changed
#    # the order later in the symbolic formulation....
#    V_arr=np.array(
#        [V_init.__getattribute__(str(v)) for v in sv]+
#        [V_init.rh]
#    ).reshape(n+1,1) #reshaping is neccessary for matmul (the @ in B @ X)
#
#    #numOutFluxesBySymbol={
#    #    k:numfunc(expr_cont,delta_t_val=delta_t_val)
#    #    for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
#    #}
#    numOutFluxesBySymbol={
#        k: gh.numfunc(expr_cont, mvs, delta_t_val=delta_t_val, par_dict=par_dict, func_dict=func_dict)
#
#        for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
#    }
#
#    def f(it,V):
#        X = V[0:n]
#        b = u_func(it,X)
#        B = B_func(it,X)
#        X_new = X + b + B @ X
#        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
#
#        l=[
#                numOutFluxesBySymbol[Symbol(k)](it,*X.reshape(n,))
#                for k in ["c_DPM","c_RPM","c_BIO","c_HUM"]
#                if Symbol(k) in numOutFluxesBySymbol.keys()
#        ]
#        rh = np.array(l).sum()
#        V_new = np.concatenate(
#            (
#                X_new.reshape(n,1),
#                np.array([rh]).reshape(1,1)
#            )
#            , axis=0
#        )
#        return V_new
#    return TimeStepIterator2(
#        initial_values=V_arr,
#        f=f,
#    )

# + codehighlighter=[[37, 51], [67, 69], [64, 65], [137, 139], [133, 135], [32, 45], [112, 113], [117, 118], [120, 123]]
import general_helpers as gh
mvs.get_InternalFluxesBySymbol()
#v2sfl={k:v for k,v in ifl.items() if k[0] in map(Symbol,["c_leaf","c_wood","c_root"])}
#delta_t_val=1
#numv2sfl = {
#    k: gh.numfunc(
#        expr_cont,
#        mvs,delta_t_val=delta_t_val,
#        par_dict=par_dict,
#        func_dict=func_dict
#    )
#    for k, expr_cont in v2sfl.items()
#}
#it=1
#X=np.ones(7)
#f=numv2sfl[(Symbol("c_leaf"),Symbol("c_DPM"))]
#f(it, *X.reshape(7, )) 
#np.array([
#    f(it, *X.reshape(7, )) 
#    for f in numv2sfl.values()
#]).sum()

# + codehighlighter=[[37, 51], [67, 69], [64, 65], [137, 139], [133, 135], [32, 45], [112, 113], [117, 118], [120, 123]]
import model_specific_helpers_2 as msh
it_sym = msh.make_iterator_sym(mvs, V_init, par_dict, func_dict, delta_t_val=2)
V_init
#
#for i in range(10):
#    print(it_sym.__next__())

# +
# from sympy import lambdify,symbols, Function, sin
# import numpy as np
# x,y,s,t=symbols("x y,s,t")
# z=Function("z")
# f_ex=x**2+y**2+z(t,s)
# f_ex
# def z_num(a1,a2):
#    return a1**2+a2
#
# func_dict={"z":z_num}
##par_dict={x:2,y:3}
##f_ex.subs(par_dict)
# f=lambdify((x,y,s,t),f_ex,[func_dict,'numpy'])
# f(2,3,3,4)

# +
#from general_helpers import day_2_month_index
#
#
#def xi_maker(tsl, Mw, Ms, mrso, landCoverFrac):  # check the mrso unit here!!, in kg m-2
#    def xi_func(day):
#        mi = day_2_month_index(day)
#        FT = 2.0 ** ((tsl[mi] - 298.15) / 10)  # temperature rate modifier
#        FV = 0.6 + 0.4 * (1 - landCoverFrac[mi] / 100)  # effect of vegetation cover
#        # Mw is soil moisture at wilting point as a fraction of saturation
#        # Ms is soil moisture content at satration 
#        S0 = 0.5 * (1 + Mw)  # optimum soil moisture
#        Smin = 1.7 * Mw  # lower threshold soil moisture for soil respiration
#        if S0 > Smin:
#            FS = 1 - 0.8 * (mrso[mi]/Ms - S0)  # effect of soil moisture
#        if (Smin < mrso[mi]/Ms) and (mrso[mi]/Ms <= S0):
#            FS = 0.2 + 0.8 * (mrso[mi]/Ms - Smin) / (S0 - Smin)
#        if mrso[mi]/Ms <= Smin:
#            FS = 0.2
#        rh_factor = FT * FV * FS
#        return rh_factor  # 1.0     # Set to 1 if no scaler implemented
#
#    return xi_func
#
#
#xi = xi_maker(
#    tsl=np.array([1, 2, 3, 4]),
#    Mw=0.5,
#    Ms= 2000,
#    mrso=np.array([0.1, 0.2, 0.3, 0.5]),
#    landCoverFrac=np.array([10, 20, 40, 80])
#)
#
#cpa._fields

# +
# def make_param2res_sym(
#        mvs,
#        cpa: msh.Constants,
#        dvs,# : drivers,
#    ) -> Callable[[np.ndarray], np.ndarray]:
#
#    # Build iterator
#    # Need dictionary of numeric values for all parameters that are not state variables/time
#    srm=mvs.get_SmoothReservoirModel()
#    model_par_dict_keys=srm.free_symbols.difference(
#        [Symbol(str(mvs.get_TimeSymbol()))]+
#        list(mvs.get_StateVariableTuple())
#    )
## Time dependent driver function does not change with the estimated parameters
## Defined once outside param2res function
##seconds_per_day = 86400
#
#    def npp_func(day):
#        month=day_2_month_index(day)
#        return dvs.npp[month]   # kg/m2/s kg/m2/day;
#
############ Ask Markus here
## Build environmental scaler function  ############### day or monthly, monthly inputs here, Mw is the to-be-estimated parameter
#
#
#    # Define actual forward simulation function
#    def param2res(pa):
#        # Define function dictionary
#        func_dict={
#            'NPP':npp_func,
#            'xi': xi_maker(dvs.tsl, pa.Mw, dvs.mrso, dvs.landCoverFrac)
#        }
#
#        # Parameter vector
#        epa=EstimatedParameters(*pa)
#
#        # Create a startvector for the iterator
#        V_init = StartVector(
#            c_leaf = epa.c_leaf_0,
#            c_wood = epa.c_wood_0,
#            c_root = cpa.c_veg_0 - (epa.c_leaf_0 + epa.c_wood_0),
#
#            c_DPM = epa.c_DPM_0,
#            c_RPM = epa.c_RPM_0,
#            c_BIO = epa.c_BIO_0,
#            c_HUM = cpa.c_soil_0-(epa.c_DPM_0 + epa.c_RPM_0 + epa.c_BIO_0) ,
#
#            rh=cpa.rh_0
#        )
#
#        # Parameter dictionary for the iterator
#        apa = {**cpa._asdict(),**epa._asdict()}
#        model_par_dict = {
#            Symbol(k):v for k,v in apa.items()
#            if Symbol(k) in model_par_dict_keys
#        }
#
#        # size of the timestep in days
#        # We could set it to 30 o
#        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2)
#        delta_t_val=15
#        it_sym = make_iterator_sym(
#            mvs,
#            V_init=V_init,
#            par_dict=model_par_dict,
#            func_dict=func_dict,
#            delta_t_val=delta_t_val
#        )
#
#        # Now that we have the iterator we can start to compute.
#        # the first thing to notice is that we don't need to store all values (daili yi delta_t_val=1)
#        # since the observations are recorded monthly while our iterator possibly has a smaller timestep.
#        # - for the pool contents we only need to record the last day of the month
#        # - for the respiration(s) ra and rh we want an accumulated value (unit mass)
#        #   have to sum up the products of the values*delta_t over a month
#        #
#        # Note: check if TRENDY months are like this...
#        # days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
#        sols=[]
#        dpm=30 #
#        n=len(V_init)
#        for m in range(cpa.number_of_months):
#            #dpm = days_per_month[ m % 12]
#            mrh=0
#            for d in range(int(dpm/delta_t_val)):
#                v = it_sym.__next__().reshape(n,)
#                # actually the original datastream seems to be a flux per area (kg m^-2 ^-s )
#                # at the moment the iterator also computes a flux but in kg^-2 ^day
#            V=StartVector(*v)
#            #from IPython import embed;embed()
#            o=Observables(
#                cVeg=float(V.c_leaf+V.c_wood+V.c_root),
#                cSoil = float(V.c_DPM + V.c_RPM + V.c_BIO + V.c_HUM),
#                rh = V.rh/seconds_per_day # the data is per second while the time units for the iterator refer to days
#                ########### Ask Markus here
#                #fVegSoil =  #Total carbon mass from vegetation directly into the soil
#            )
#            sols.append(o)
#
#        sol=np.stack(sols)
#        #convert to yearly output if necessary (the monthly pool data looks very "yearly")
#        #sol_yr=np.zeros(int(cpa.number_of_months/12)*sol.shape[1]).reshape([int(cpa.number_of_months/12),sol.shape[1]])
#        #for i in range(sol.shape[1]):
#        #   sol_yr[:,i]=monthly_to_yearly(sol[:,i])
#        #sol=sol_yr
#        return sol
#
#    return param2res
#
## Observables = namedtuple('observables_monthly', ["cVeg", "cSoil", "rh", "fVegSoil"]) # all monthly
# -

# ## Forward Model Run
# #### Run model forward:

svs_cut = msh.Observables(
    cVeg=svs.cVeg[:cpa.nyears * 12],
    cSoil=svs.cSoil[:cpa.nyears * 12],
    rh=svs.rh[:cpa.nyears * 12],
    fVegSoil =svs.fVegSoil[:cpa.nyears*12]
)
svs_cut.rh

# +
import model_specific_helpers_2 as msh
param2res = msh.make_param2res_sym(mvs, cpa, dvs)  # Define forward model
obs_simu = param2res(epa_0)  # Run forward model from initial conditions

obs_simu

# +
#def make_weighted_cost_func(
#        obs: Observables
#) -> Callable[[Observables], np.float64]:
#    # first unpack the observation array into its parts
#    # cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
#    def costfunction(out_simu: Observables) -> np.float64:
#        # fixme
#        #   as indicated by the fact that the function lives in this
#        #   model-specific module it is NOT apropriate for (all) other models.
#        #   There are model specific properties:
#        #   1.) The weight for the different observation streams
#        #
#        number_of_ys = out_simu.cVeg.shape[0]
#        number_of_ms = out_simu.rh.shape[0]
#
#        J_obj1 = np.mean((out_simu.cVeg - obs.cVeg) ** 2) / (2 * np.var(obs.cVeg))
#        J_obj2 = np.mean((out_simu.cSoil - obs.cSoil) ** 2) / (2 * np.var(obs.cSoil))
#
#        J_obj3 = np.mean((out_simu.rh - obs.rh) ** 2) / (2 * np.var(obs.rh))
#        J_obj3 = np.mean((out_simu.rh - obs.rh) ** 2) / (2 * np.var(obs.rh))
#
#        J_new = (J_obj1 + J_obj2) + J_obj3 #/ 12  # the 12 is purely conjectural
#        return J_new
#
#    return costfunction
#
#
## test it
cost_func = msh.make_weighted_cost_func(svs_cut)
cost_func(obs_simu)

# for f in Observables._fields:
#    print(obs_0._asdict()[f])
# -

sols=[msh.Observables(cVeg=1,cSoil=2,rh=2,fVegSoil=3),msh.Observables(cVeg=1,cSoil=2,rh=2,fVegSoil=3)]
sols
ref=msh.Observables(cVeg=np.array([1,1]),cSoil=np.array([2,2]),rh=np.array([2,2]),fVegSoil=np.array([3,3]))
def arr(field):
    n=len(sols)
    a=np.zeros(n)
    for i in range(n):
        a[i]=sols[i].__getattribute__(field)
    return a
msh.Observables(*( arr(field) for field in msh.Observables._fields))

# #### Create array of yearly observation data:

# +
fig = plt.figure()
from general_helpers import plot_observations_vs_simulations
plot_observations_vs_simulations(fig,svs,obs_simu)

#fig.savefig("test.pdf")
# -

gh.plot_observations_vs_simulations(fig, svs_cut, obs_simu)
fig.savefig("test.pdf")

# #### Plot data-model fit:

# +
# Plot simulation output for observables
# n_plots = len(svs_cut)
# fig = plt.figure(figsize=(10, n_plots * 5))
# axs = fig.subplots(n_plots)
# for i, name in enumerate(msh.Observables._fields):
#     var = svs_cut[i]
#     var_simu = obs_simu[i]
#     axs[i].plot(range(len(var_simu)), var_simu, label="simulation")
#     axs[i].plot(range(len(var)), var, label='observation')
#     axs[i].legend()
#     axs[i].set_title(name)
# -

# ## Data Assimilation
# #### Define parameter min/max values:

Test=namedtuple("Test",["a","b"])
t1=Test(a=1,b=3)
t1*2
t1.__class__._fields

# +
# set min/max parameters to +- 100 times initial values
epa_min = msh.EstimatedParameters(* tuple(np.array(epa_0) * 0.001)) # for fluxes 
epa_max = msh.EstimatedParameters(* tuple(np.array(epa_0) * 200))

# fix values that are problematic from calculation
epa_max = epa_max._replace(beta_leaf = epa_0.beta_leaf * 1.5)
epa_max = epa_max._replace(beta_wood = epa_0.beta_wood * 1.5)
epa_max = epa_max._replace(c_leaf_0 = svs_0.cVeg * epa_0.c_leaf_0 * 1.5)
epa_max = epa_max._replace(c_wood_0 = svs_0.cVeg * epa_0.c_wood_0 * 1.5)
epa_max = epa_max._replace(c_DPM_0 = svs_0.cSoil * epa_0.c_DPM_0 * 10)
epa_max = epa_max._replace(c_RPM_0 = svs_0.cSoil * epa_0.c_RPM_0 * 10)
epa_max = epa_max._replace(c_BIO_0 = svs_0.cSoil * epa_0.c_BIO_0 * 10)
epa_max = epa_max._replace(Mw = 0.5)
epa_max = epa_max._replace(Ms = (max(dvs.mrsos) + 500) * 2) # max(dvs.mrso) * 2
# -

epa_0

# print - all names should have numerical values
epa_min._asdict()

epa_max._asdict()

# #### Conduct data assimilation:

import general_helpers as gh
print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = gh.autostep_mcmc(
    initial_parameters = epa_0,
    filter_func = gh.make_param_filter_func_2(epa_max, epa_min,["beta_leaf","beta_wood"]),
    param2res = msh.make_param2res_sym(mvs, cpa, dvs),
    costfunction = msh.make_weighted_cost_func(svs),
    nsimu = 20000,  # for testing and tuning mcmc
    c_max = np.array(epa_max),
    c_min = np.array(epa_min),
    acceptance_rate = 15,  # default value | target acceptance rate in %
    chunk_size = 100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init = 1,  # default value | increase value to reduce initial step size
    K = 2  # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

# #### Graph data assimilation results:

# +
# optimized parameter set (lowest cost function)
par_opt = np.min(
    C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields), 1),
    axis=1)
epa_opt = msh.EstimatedParameters(*par_opt)
mod_opt = param2res(epa_opt)

n_plots = len(svs)
fig = plt.figure(figsize=(12, 4), dpi=80)
plt.rcParams['font.size'] = 18
plot_observations_vs_simulations(
        fig,
        svs,
        mod_opt
    )
#plt.rcParams.update({'font.size': 22})
fig.savefig('solutions_opt.pdf')

# # n_plots=len(svs_cut)
#fig = plt.figure(figsize=(10, 10 * n_plots))
#axs = fig.subplots(n_plots)
#for i, name in enumerate(msh.Observables._fields):
#    var_simu = mod_opt[i]
#    var = svs[i]
#    axs[i].plot(range(len(var_simu)), var_simu, label='simulation(opt)')
#    axs[i].plot(range(len(var)), var, label='observation')
#    axs[i].set_title(name)
#    axs[i].legend()

#fig.savefig('solutions_opt.pdf')


# +
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms

fig = plt.figure(figsize=(12, 4), dpi=80)
plt.rcParams['font.size'] = 15
fig, ax = plt.subplots() #[(12*10-1):12*30]
ax.scatter(np.array(svs.rh)*1000, np.array(mod_opt.rh)*1000, c='black')
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
plt.show()

# +
plt.plot(list(range(0, 12*20)), np.array(svs.rh[(12*300):12*320])*1000, label = "observation")
plt.plot(list(range(0, 12*20)), np.array(mod_opt.rh[(12*300):12*320])*1000, label = "simulation")
#plt.plot(list(range(0, 12*20)), np.array(dvs.npp[(12*300):12*320])*1000, label = "NPP")

#plt.plot(list(range(0, 12*320)), np.array(dvs.npp) - np.array(svs.fVegSoil), label = "net change in veg C")
#plt.plot(list(range(0, 12*320)), np.array(svs.cVeg), label = "net veg C")
plt.legend()
plt.plot()
plt.show()

# +
# save the parameters and cost function values for postprocessing
outputPath = Path(conf_dict["dataPath"])  # save output to data directory (or change it)

import pandas as pd

pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('JULES_da_pars.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('JULES_da_cost.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('JULES_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('JULES_optimized_solutions.csv'), sep=',')
# -

epa_opt
