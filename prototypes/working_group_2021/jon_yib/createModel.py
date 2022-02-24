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

# ## Jupyter Settings
# Settings for display: 

# +
# adjust jupyter display to full width
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# make changes to imported files immidiately available 
# avoiding the need to reload (in most cases)
# %load_ext autoreload
# %autoreload 2
# -

# ## Symbolic Model Description
# Define your model using symbolic python code:

# +
from sympy import Symbol, Function 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

# Make a dictionary for the variables we will use
sym_dict={
    'c_leaf': '',
    'c_root': '',
    'c_wood': '',
    'c_lit_cwd': '',
    'c_lit_met': '',
    'c_lit_str': '',
    'c_lit_mic': '',
    'c_soil_met': '',
    'c_soil_str': '',
    'c_soil_mic': '',
    'c_soil_slow': '',
    'c_soil_passive': '',
    'r_c_lit_cwd_rh': '',
    'r_c_lit_met_rh': '',
    'r_c_lit_str_rh': '',
    'r_c_lit_mic_rh': '',
    'r_c_soil_met_rh': '',
    'r_c_soil_str_rh': '',
    'r_c_soil_mic_rh': '',
    'r_c_soil_slow_rh': '',
    'r_c_soil_passive_rh': '',
    'r_c_leaf_2_c_lit_met': '',
    'r_c_leaf_2_c_lit_str': '',
    'r_c_root_2_c_soil_met': '',    
    'r_c_root_2_c_soil_str': '',
    'r_c_wood_2_c_lit_cwd': '',
    'r_c_lit_cwd_2_c_lit_mic': '',
    'r_c_lit_cwd_2_c_soil_slow': '',
    'r_c_lit_met_2_c_lit_mic': '',
    'r_c_lit_str_2_c_lit_mic': '',
    'r_c_lit_str_2_c_soil_slow': '',
    'r_c_lit_mic_2_c_soil_slow': '',
    'r_c_soil_met_2_c_soil_mic': '',
    'r_c_soil_str_2_c_soil_mic': '',
    'r_c_soil_str_2_c_soil_slow': '',
    'r_c_soil_mic_2_c_soil_slow': '',
    'r_c_soil_mic_2_c_soil_passive': '',
    'r_c_soil_slow_2_c_soil_mic': '',
    'r_c_soil_slow_2_c_soil_passive': '',
    'r_c_soil_passive_2_c_soil_mic': '',
    't': '',
    'silt': '',
    'clay': '',
    'beta_leaf': '',
    'beta_root': '',
}
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    'xi': 'a scalar function of temperature and moisture and thereby ultimately of time',
    'NPP': '',
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")
beta_wood = 1.0-(beta_leaf+beta_root)
mvs = CMTVS(
    {
        t,
        StateVariableTuple(( 
            c_leaf,
            c_root,
            c_wood,
            c_lit_cwd,
            c_lit_met,
            c_lit_str,
            c_lit_mic,
            c_soil_met,
            c_soil_str,
            c_soil_mic,
            c_soil_slow,
            c_soil_passive,
        )),
        InFluxesBySymbol(
            {
                c_leaf: NPP(t) * beta_leaf, 
                c_root: NPP(t) * beta_root, 
                c_wood: NPP(t) * beta_wood
            }
        ),
        OutFluxesBySymbol(
            {
                c_lit_cwd: r_c_lit_cwd_rh*c_lit_cwd*xi(t),
                c_lit_met: r_c_lit_met_rh*c_lit_met*xi(t),
                c_lit_str: r_c_lit_str_rh*c_lit_str*xi(t),
                c_lit_mic: r_c_lit_mic_rh*c_lit_mic*xi(t),
                c_soil_met: r_c_soil_met_rh*c_soil_met*xi(t),
                c_soil_str: r_c_soil_str_rh*c_soil_str*xi(t),
                c_soil_mic: r_c_soil_mic_rh*c_soil_mic*xi(t),
                c_soil_slow: r_c_soil_slow_rh*c_soil_slow*xi(t),
                c_soil_passive: r_c_soil_passive_rh*c_soil_passive*xi(t),
            }
        ),
        InternalFluxesBySymbol(
            {
                (c_leaf, c_lit_met): r_c_leaf_2_c_lit_met*c_leaf,
                (c_leaf, c_lit_str): r_c_leaf_2_c_lit_str*c_leaf,
                (c_root, c_soil_met): r_c_root_2_c_soil_met*c_root,
                (c_root, c_soil_str): r_c_root_2_c_soil_str*c_root,
                (c_wood, c_lit_cwd): r_c_wood_2_c_lit_cwd*c_wood,
                (c_lit_cwd, c_lit_mic): r_c_lit_cwd_2_c_lit_mic*c_lit_cwd*xi(t),
                (c_lit_cwd, c_soil_slow): r_c_lit_cwd_2_c_soil_slow*c_lit_cwd*xi(t),
                (c_lit_met, c_lit_mic): r_c_lit_met_2_c_lit_mic*c_lit_met*xi(t),
                (c_lit_str, c_lit_mic): r_c_lit_str_2_c_lit_mic*c_lit_str*xi(t),
                (c_lit_str, c_soil_slow): r_c_lit_str_2_c_soil_slow*c_lit_str*xi(t),
                (c_lit_mic, c_soil_slow): r_c_lit_mic_2_c_soil_slow*c_lit_mic*xi(t),
                (c_soil_met, c_soil_mic): r_c_soil_met_2_c_soil_mic*c_soil_met*xi(t),
                (c_soil_str, c_soil_mic): r_c_soil_str_2_c_soil_mic*c_soil_str*xi(t),
                (c_soil_str, c_soil_slow): r_c_soil_str_2_c_soil_slow*c_soil_str*xi(t),
                (c_soil_mic, c_soil_slow): r_c_soil_mic_2_c_soil_slow*c_soil_mic*xi(t),
                (c_soil_mic, c_soil_passive): r_c_soil_mic_2_c_soil_passive*c_soil_mic*xi(t),
                (c_soil_slow, c_soil_mic): r_c_soil_slow_2_c_soil_mic*c_soil_slow*xi(t),
                (c_soil_slow, c_soil_passive): r_c_soil_slow_2_c_soil_passive*c_soil_slow*xi(t),
                (c_soil_passive, c_soil_mic): r_c_soil_passive_2_c_soil_mic*c_soil_passive*xi(t),               
            }
        ),
    },
    computers=module_computers(bgc_c)
)
# -
import bgc_md2.helper as h
h.compartmental_graph(mvs)

import bgc_md2.display_helpers as dh
dh.mass_balance_equation(mvs)

mvs.get_CompartmentalMatrix().free_symbols

# ### Add information specific to the tracability analysis.
#
# We are aiming for a decomposition of the compartmental matrix $B$ into three factors 
# $B(t)= \xi(t)  A K $ 
# with $ \xi$ and $K$ diagonal. 
# `bgc_md2` can automatically compute a decomposition into $ B=A N(t)$ where $N(t)=\xi(t)K$ but
# which part of $N(t)$ should go into $\xi(t)$ is easier to specify manually. 
#
# We will first compute the $B=A N$ decomposition and then specify $\xi$.
#
#

# +
srm = mvs.get_SmoothReservoirModel()

# The matrices T and N in the function refer to what in Yiqi Luo's group is usually called A and K
# and xi is a single scalar (not the diagonal matrix we are looking for here)
# The function has documentation which you can see by typing the following
# # ?srm.xi_T_N_u_representation
_,A,N,_,_=srm.xi_T_N_u_representation(factor_out_xi=False) 
# -

mvs.computable_mvar_types()

# We can see that we can easily decompose N manually into $N= \xi K $.
# We also see that we have a symbol $\xi$ alreasince the model has actually been specified with $\xi in mind

from sympy import diag
xi_d=diag([1,1,1]+[xi(t) for i in range(9)],unpack=True)
xi_d

# We can go on and decompose N =\xi K -> K=\xi^{-1}N
K=xi_d.inv()*N
K
# we now have the (symbolic) ingredients for the tracebility analysis.
#xi_d,K,A

# We can store the symbolic representation in a file that we can use in for several purposes later.
# Copy your version of the following cell into a text editor and save the file in the current directory as `source.py`)



mvs.get_StateVariableTuple()

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
# #### create a small site specific config file 
# This file specifies:
# - a username and password to download the data 
# - the location where you want to download the data to 
#   which will differ depending on the machine you are using (your laptop or a supercomputer) and     also accross users. You will have to have one everywhere you want to work with the model.
# Here comes a template from my laptop (content of file `../config.json`):
# `{"username": "trendy-v9", "password": "gcb-2020", "dataPath": "/home/data/VISIT_data_CMIP6"}`
#
# Note that 
# - the file resides one level above your current folder since it is not modelspecific
#   (This is a change from the first round of model gathering)
#
# #### create a small model specific function to download the data
# This function  will later be stored in the file `model_specific_helpers.py` 
# This function **saves us a lot of time** when we want to reproduce your results and run your model since finding the correct dataset can be a time consuming task.
# There is a helper function `download_TRENDY_output` in `../general_helpers.py` that you can probably use within your `download_my_TRENDY_output`.
# If it's not general enough to be called it should at least give you an idea. 
#

# +
import sys
sys.path.insert(0,'..') # necessary to import general_helpers
from general_helpers import download_TRENDY_output
import json 
from pathlib import Path
from collections import namedtuple 

with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

# we will use the trendy output names directly in other parts of the output
observables_annual = namedtuple(
    'observables_annual',
    ["cSoil", "cVeg"]
)
observables_monthly = namedtuple(
    'observables_monthly',
    ["rh", "ra"]
)
observables = namedtuple(
    "observables",
    observables_annual._fields+observables_monthly._fields
)

# +
Drivers=namedtuple(
    "Drivers",
    ["npp"]
) 

#when downloading data make sure model names match TRENDY server names:
#"CABLE-POP","CLASSIC","CLM5","DLEM","IBIS","ISAM","ISBA_CTRIP",
#"JSBACH","JULES-ES-1.0","LPJ-GUESS","LPJwsl","LPX-Bern","OCN",
#"ORCHIDEE","ORCHIDEE-CNP","ORCHIDEEv3","ORCHIDEEv3_0.5deg"
#"SDGVM","VISIT","YIBs"

#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output():
    download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['YIBs'],
        variables =observables._fields + Drivers._fields
    )
#call it to test that the download works the data
download_my_TRENDY_output()
# -

# #### provide functions to read the data
#



# +
import netCDF4 as nc
import numpy as np
from pathlib import Path
import json 
def get_example_site_vars(dataPath):
    # pick up 1 site   
    s = slice(None, None, None)  # this is the same as :
    t = s, 74, 118  # [t] = [:,49,325]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them 
        ds = nc.Dataset(str(path))
        #check for npp/gpp/rh/ra to convert from kg/m2/s to kg/m2/day
        if vn in ["npp","gpp","rh","ra"]:
            return ds.variables[vn][t]*86400
        else:
            return ds.variables[vn][t]
    o_names=[(f,"YIBs_S2_Annual_{}.nc".format(f)) for f in observables_annual._fields]
    #YIBs has unique naming convention (Annual vs Monthly), create monthly names and extend 0_names vector
    monthly_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in observables_monthly._fields]
    o_names.extend(monthly_names)
    print(o_names)
    d_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Drivers._fields]
    return (observables(*map(f, o_names)),Drivers(*map(f,d_names)))


#     # Read NetCDF data  ******************************************************************************************************************************
# -

svs,dvs=get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))


svs,dvs

# +
sys.path.insert(0,'..')
from general_helpers import day_2_month_index
def NPP_fun(day ):
    return npp[day_2_month_index(day)] 

func_dict={NPP: NPP_fun}
# -



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

# +
from sympy import ImmutableMatrix
sv=mvs.get_StateVariableTuple()
n=len(sv)
# create new symbols for the f_{i,j}
for i in range(n):
    for j in range(n):
        if A[i,j]!=0 and i!=j:
            name="f_" + str(sv[j]) + "_2_" + str(sv[i])
            code="{0}=Symbol('{0}')".format(name)
            print(code)
            exec(code)
            
A_sym=ImmutableMatrix(
    n,n,
    lambda i,j:  -1 if i==j else (
        0 if A[i,j]==0 else Symbol("f_" + str(sv[j]) + "_2_" + str(sv[i]))
    )
)
A_sym

# +
# create new symbols for the k_{i}
for i in range(n):
    if K[i,i]!=0:
        name="k_{0}".format(sv[i])
        code="{0}=Symbol('{0}')".format(name)
        print(code)
        exec(code)
        
K_sym=ImmutableMatrix(
    n,n,
    lambda i,j: Symbol("k_" + str(sv[i])) if i==j else 0
)
K_sym
# -

M_sym=A_sym*K_sym
M_sym

import  CompartmentalSystems.helpers_reservoir as hr
hr.out_fluxes_by_symbol(sv,M_sym)

# we create a dictionary for the outfluxrates (flux divided by dono pool content)
outflux_rates = {"r_"+str(key)+"_rh":value/key for key,value in hr.out_fluxes_by_symbol(sv,M_sym).items()}
internal_flux_rates = {"r_"+str(key[0])+"_2_"+str(key[1]):value/key[0] for key,value in hr.internal_fluxes_by_symbol(sv,M_sym).items()}
from copy import  deepcopy
all_rates=deepcopy(outflux_rates)
all_rates.update(internal_flux_rates)
all_rates

# +
#I have estimates from previous model main.py file
#we can define them here and use them to assign values to fs
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
        "C_leaf_0",
        "C_root_0",
        "C_cwd_0",
        "C_samet_0",
        "C_sastr_0",
        "C_samic_0",
        "C_slmet_0",
        "C_slstr_0",
        "C_slmic_0",
        "C_slow_0"
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
    C_leaf_0=0.06,
    C_root_0=0.08,
    C_cwd_0=0.5,
    C_samet_0=0.4,
    C_sastr_0=1.07,
    C_samic_0=0.49,
    C_slmet_0=0.83,
    C_slstr_0=2.07,
    C_slmic_0=1.04,
    C_slow_0=5.66,
)
# -

old_par_dict = {
    beta_leaf: epa0.beta_leaf,
    beta_root: epa0.beta_root,
    clay: epa0.clay,
    silt: epa0.silt,
    k_c_leaf: 0.014,
    k_c_wood: 0.003,
    k_c_root: 0.022,
    k_c_lit_cwd: 0.005,
    k_c_lit_met: 0.01,
    k_c_lit_str: 0.001,
    k_c_lit_mic: 0.05,
    k_c_soil_met: 0.040,
    k_c_soil_str: 0.0039,
    k_c_soil_mic: 0.005,
    k_c_soil_slow: 0.00001,
    k_c_soil_passive: 3.27E-06,
    #these first three fs were set for DA
    f_c_leaf_2_c_lit_met: epa0.f_samet_leaf,  
    f_c_root_2_c_soil_met: epa0.f_slmet_root, 
    f_c_lit_cwd_2_c_lit_mic: epa0.f_samic_cwd,
    #the rest of the values were defined in the forward sim so are calculated here
    f_c_leaf_2_c_lit_str: (1-epa0.f_samet_leaf) * epa0.k_leaf,
    f_c_root_2_c_soil_str: (1-epa0.f_slmet_root) * epa0.k_root,
    f_c_wood_2_c_lit_cwd: 0.4 * epa0.k_wood,
    f_c_lit_cwd_2_c_soil_slow: (1-epa0.f_samic_cwd) * epa0.k_cwd,
    f_c_lit_met_2_c_lit_mic: 0.3 * epa0.k_samet,
    f_c_lit_str_2_c_lit_mic: 0.1 * epa0.k_sastr,
    f_c_lit_str_2_c_soil_slow: 0.1 * epa0.k_sastr,
    f_c_lit_mic_2_c_soil_slow: 0.1 * epa0.k_samic,
    f_c_soil_met_2_c_soil_mic: 0.4 * epa0.k_slmet,
    f_c_soil_str_2_c_soil_mic: 0.3 * epa0.k_slstr,
    f_c_soil_str_2_c_soil_slow: 0.2 * epa0.k_slstr,
    f_c_soil_mic_2_c_soil_slow: 0.4 * epa0.k_slmic,
    f_c_soil_mic_2_c_soil_passive: 0.4 * epa0.k_slmic,
    f_c_soil_slow_2_c_soil_mic: 0.10 * epa0.k_slow,
    f_c_soil_slow_2_c_soil_passive: 0.45*(0.003+0.009*epa0.clay) * epa0.k_slow,
    f_c_soil_passive_2_c_soil_mic: 0.10 * epa0.k_arm, 
}


# Now we can translate the old paramterisation to the new one.
#
# ### Providing dictionaries  for the parameters and functions.
#

par_dict={
    beta_leaf: 0.21,
    beta_root: 0.27,
    clay: 0.2028,
    silt: 0.2808,
}
par_dict.update(
    {Symbol(k):v.subs(old_par_dict) for k,v in all_rates.items()}
)
par_dict

all_rates

for ex in par_dict.values():
    try:
        print(ex.free_symbols)
    except AttributeError:
        pass
        #print(ex)

# To be able to run the model forward we not only have to replace parameter symbols by values but symbolic functions by normal python functions.
# In our case the functions for $NPP$ and $\xi$ have to be provided. NPP_fun will interpolate the NPP for the day in question from the data. Which we have to load. 
# We will later store these functions in  `model_specific_helpers.py` which resides in the same folder as this notebook. You will have to adapt them to your data set. 

# +
from general_helpers import make_B_u_funcs_2, day_2_month_index
# check the numeric functions for B and u

def npp_func(day):
    month=day_2_month_index(day)
    return dvs.npp[month]#-(svs.rh[month]+svs.ra[month])

def xi_func(day):
    return 1.0 # preliminary fake for lack of better data... 

func_dict={
    'NPP':npp_func,
     'xi':xi_func
}
# for the next line to work the 
# two dictionaries par_dict and func_dict have to be complete.
# In sympy terms this means that the expressions for 
B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
# we create a numeric startvector for the compartmental system
# 
svs_0=observables(*map(lambda v: v[0],svs))

X_0= np.array((
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
))#.reshape(9,)
# in general B and u are nonlinear and can depend on X, thats why we have to test them with index and X arguments
u_func(0,X_0),B_func(0,X_0)
# -

# We now build the essential object to run the model forward. Technically it supports the `iterator` interface which means that we can later call its `__next__()` method to move our system one step forward in time. If iterators had not been invented yet we would invent them now, because they capture exactly the mathematical concept of an initial value system, where we have a startvector $V_0$ and a function $f$ to compute the next value: $V_{i+1} =f(V_{i})$ without all the nonessential technical details of e.g. where to store the results and so on.
# If we were only interested in the timeseries of the pool contents `bgc_md2` could compute the solution automatically wihtout the need to build an iterator.
#
# In our case we are also interested in tracking the autotrophic and heterotrophic respiration.
# So we will let `bgc_md2` derive numeric functions for the Compartmental matrix $B$ and the input $u$ from our symbolic description but use them to build our own iterator. (which we do by a function)    
#
# We will start by creating $V_0$ and then build the function $f$
#
#

# to guard agains accidentally changed order we use a namedtuple again. Since B_func and u_func rely 
# on the correct ordering of the statevariables we build V dependent on this order 
StartVector=namedtuple(
    "StartVector",
    [str(v) for v in mvs.get_StateVariableTuple()]+
    ["ra","rh"]
)
StartVector._fields

V_init= StartVector(
    c_leaf=svs_0.cVeg/3,
    c_root=svs_0.cVeg/3,
    c_wood=svs_0.cVeg/3,
    c_lit_cwd=svs_0.cSoil/9,
    c_lit_met=svs_0.cSoil/9,
    c_lit_str=svs_0.cSoil/9,
    c_lit_mic=svs_0.cSoil/9,
    c_soil_met=svs_0.cSoil/9,
    c_soil_str=svs_0.cSoil/9,
    c_soil_mic=svs_0.cSoil/9,
    c_soil_slow=svs_0.cSoil/9,
    c_soil_passive=svs_0.cSoil/9,
    ra=svs_0.ra,
    rh=svs_0.rh  
)
V_init.__getattribute__("c_leaf")

# +
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
        TimeStepIterator2,
)
from copy import copy

def make_daily_iterator_sym(
        mvs,
        V_init: StartVector,
        par_dict,
        func_dict
    ):
    B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
    sv=mvs.get_StateVariableTuple()
    n=len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct 
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr=np.array(
        [V_init.__getattribute__(str(v)) for v in sv]+
        [V_init.ra,V_init.rh]
    ).reshape(n+2,1) #reshaping is neccessary for matmux
    
    def f(it,V):
        X = V[0:n]
        b = u_func(it,X)
        B = B_func(it,X)
        outfluxes = B @ X
        X_new = X + b + outfluxes
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
        #ra= -np.sum(outfluxes[0:3]) # check with  mvs.get_StateVariableTuple()[0:3]
        #rh= -np.sum(outfluxes[3:n]) # check with  mvs.get_StateVariableTuple()[3:9]
        ra=0
        rh=0
        V_new = np.concatenate((X_new.reshape(n,1),np.array([ra,rh]).reshape(2,1)), axis=0)
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


# -

V_init

np.array(V_init).shape

# +
# test the daily iterator
    
it_sym = make_daily_iterator_sym(
    mvs,
    V_init=V_init,
    par_dict=par_dict,
    func_dict=func_dict
)
# we will run the model for 15 steps
ns=15
res= np.zeros((ns,len(V_init)))
res_sym = copy(res)
for i in range(ns):
    res_sym[i,:]=it_sym.__next__().reshape(len(V_init),)
res_sym
# -

# ## Data assimilation
# Until now we have used only the initial values of the observations. 
# The next step is to decide which parameters we want to consider fixed and which to be estimated.
# This distinction helps, to keep the to create a function which only takes the estimated parameters and thus can be used by a generalized mcmc as will become clear.
#
# We can change which parameters we fix and which we estimate later or can have several approaches for the same symbolic model.
# The distinction is not model inherent but just a reflection of our choice for data assimilation.
# The more parameter values we can find out from the literature the fewer values we have to estimate.  

# +
# As a safety measure we specify those parameters again as 'namedtuples', which are like a mixture of dictionaries and tuples
# They preserve order as numpy arrays which is great (and fast) for the numeric computations
# and they allow to access values by keys (like dictionaries) which makes it difficult to accidentally mix up values.

UnEstimatedParameters = namedtuple(
    "UnEstimatedParameters",
    [
        'npp_0',
        'rh_0',
        'ra_0',
        'c_veg_0',
        'c_soil_0',
        'clay',
        'silt',
        'nyears'
    ]
)

# Note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues and their length (number of months)
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 
# It is better to start with only a few

RateParameters = namedtuple(
    "RateParameters",
    par_dict.keys(),
)
InitialPools = namedtuple(
    "InitialPools",
    [
        'c_leaf_0',
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

#pool_dict = {
#    c_leaf_0: X_0[0],
#    c_root_0: X_0[1],
#    c_wood_0: c_veg_0 - np.sum(X_0[0:2]),
#    c_lit_cwd_0: X_0[3],
#    c_lit_met_0: X_0[4],
#    c_lit_str_0: X_0[5],
#    c_lit_mic_0: X_0[6],
#    c_soil_met_0: X_0[7],
#    c_soil_str_0: X_0[8],
#    c_soil_mic_0: X_0[9],
#    c_soil_slow_0: X_0[10],
#    c_soil_passive_0: c_soil_0 - np.sum(X_0[3:10])
#}

EstimatedParameters = namedtuple(
    "EstimatedParameters",
    RateParameters._fields + InitialPools._fields
)
# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues. 
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 


# -

EstimatedParameters._fields

from general_helpers import monthly_to_yearly
cpa=UnEstimatedParameters(
    npp_0 = monthly_to_yearly(dvs.npp)[0],
    rh_0 = monthly_to_yearly(svs.rh)[0],
    ra_0 = monthly_to_yearly(svs.ra)[0],    
    c_veg_0 = svs.cVeg[0],
    c_soil_0 = svs.cSoil[0],
    clay = 0.2028,
    silt = 0.2808,
    nyears = 320
)

# create a start parameter tuple for the mcmc. The order has to be the same as when you created the namedtupl3 
# If you don't you get a "TypeError". 
epa_0=EstimatedParameters(
    beta_leaf=0.21,
    beta_root=0.27,
    clay=0.2028,
    silt=0.2808,
    r_c_leaf_rh=0.00993888000000000,
    r_c_root_rh=0.0142005600000000,
    r_c_wood_rh=0.00299640000000000,
    r_c_lit_cwd_rh=0.00353225000000000,
    r_c_lit_met_rh=0.00997000000000000,
    r_c_lit_str_rh=0.000999800000000000,
    r_c_lit_mic_rh=0.0497500000000000,
    r_c_soil_met_rh=0.0393600000000000,
    r_c_soil_str_rh=0.00389239500000000,
    r_c_soil_mic_rh=0.00498000000000000,
    r_c_soil_slow_rh=9.99998978286600e-6,
    r_c_soil_passive_rh=3.26999893071000e-6,
    r_c_leaf_2_c_lit_met=0.00392000000000000,
    r_c_leaf_2_c_lit_str=0.000141120000000000,
    r_c_root_2_c_soil_met=0.00748000000000000,
    r_c_root_2_c_soil_str=0.000319440000000000,
    r_c_wood_2_c_lit_cwd=3.60000000000000e-6,
    r_c_lit_cwd_2_c_lit_mic=0.00145000000000000,
    r_c_lit_cwd_2_c_soil_slow=1.77500000000000e-5,
    r_c_lit_met_2_c_lit_mic=3.00000000000000e-5,
    r_c_lit_str_2_c_lit_mic=1.00000000000000e-7,
    r_c_lit_str_2_c_soil_slow=1.00000000000000e-7,
    r_c_lit_mic_2_c_soil_slow=0.000250000000000000,
    r_c_soil_met_2_c_soil_mic=0.000640000000000000,
    r_c_soil_str_2_c_soil_mic=4.56300000000000e-6,
    r_c_soil_str_2_c_soil_slow=3.04200000000000e-6,
    r_c_soil_mic_2_c_soil_slow=1.00000000000000e-5,
    r_c_soil_mic_2_c_soil_passive=1.00000000000000e-5,
    r_c_soil_slow_2_c_soil_mic=1.00000000000000e-11,
    r_c_soil_slow_2_c_soil_passive=2.17134000000000e-13,
    r_c_soil_passive_2_c_soil_mic=1.06929000000000e-12,
    c_leaf_0=V_init[0],
    c_root_0=V_init[1],
    c_wood_0=V_init[2],
    c_lit_cwd_0=V_init[3],
    c_lit_met_0=V_init[4],
    c_lit_str_0=V_init[5],
    c_lit_mic_0=V_init[6],
    c_soil_met_0=V_init[7],
    c_soil_str_0=V_init[8],
    c_soil_mic_0=V_init[9],
    c_soil_slow_0=V_init[10],
    c_soil_passive_0=V_init[11]
)

epa_0

# The function `param2res` (which will be used by a general model independent mcmc) only takes the estimated parameters as arguments and produce data in the same shape as the observations.
# We will taylor make it by another function `make_param2res` which depends on the parameters that we decide to fix.
# This function is the main effort to make the data assimilation possible. **Although we give the full definition here it's suggested that you incrementally build and check the code inside it before you make it a function that returns a function...** 
# - You can start with sets of  `UnEstimatedParameters` and `EstimatedParameter` (e.g.  `epa0` from the test above) and try to produce a set of observations by running the model. 
# - then produce `param2res` automatically by code
# - them wrap this code into `make_param2res` so that you can change the unestimated parameters easily.
#
# This function will be later moved to a model specific helper file, and loaded from there by the data assimilation functions.

# +
from typing import Callable
from general_helpers import month_2_day_index
from functools import reduce

def make_param2res_sym(
        cpa: UnEstimatedParameters
    ) -> Callable[[np.ndarray], np.ndarray]: 
    # To compute numeric solutions we will need to build and iterator 
    # as we did before. As it will need numeric values for all the parameters 
    # we will have to create a complete dictionary for all the symbols
    # exept those for the statevariables and time.
    # This set of symbols does not change during the mcmc procedure, since it only
    # depends on the symbolic model.
    # Therefore we create it outside the mcmc loop and bake the result into the 
    # param2res function.
    # The iterator does not care if they are estimated or not it only 
    # wants a dictionary containing key: value pairs for all
    # parameters that are not state variables or the symbol for time
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test
    def npp_func(day):
        month=day_2_month_index(day)
        return dvs.npp[month]#-(svs.rh[month]+svs.ra[month])
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
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
            ra=cpa.ra_0,
            rh=cpa.rh_0
        )
        # next we create the parameter dict for the iterator
        # The iterator does not care if they are estimated or not so we look for them
        # in the combination
        apa = {**cpa._asdict(),**epa._asdict()}
        model_par_dict = {
            k:v for k,v in apa.items()
            if k in model_par_dict_keys
        }
        # Beside the par_dict the iterator also needs the python functions to replace the symbolic ones with
        # our fake xi_func could of course be defined outside of param2res but in general it
        # could be build from estimated parameters and would have to live here...
        def xi_func(day):
            return 1.0 # preliminary fake for lack of better data... 
    
        func_dict={
            'NPP':npp_func,
             'xi':xi_func
        }
    
        it_sym = make_daily_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=par_dict,
            func_dict=func_dict
        )
        
        # Now that we have the iterator we can start to compute.
        # the first thing to notice is that we don't need to store all daily values,
        # since the observations are recorded monthly while our iterator has a
        # daily timestep.
        # - for the pool contents we only need to record the last day of the month
        # - for the respiration(s) ra and rh we have to sum up the daily values 
        #   over a month
        # 
        # Note: check if TRENDY months are like this...
        days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        sols=[]
        for m in range(cpa.nyears*12):
            dpm = days_per_month[ m % 12]  
            mra=0
            mrh=0
            for d in range(dpm):
                v = it_sym.__next__()
                mra +=v[9,0]
                mrh +=v[10,0]
            V=StartVector(*v)
            o=observables(
                cVeg=float(V.c_leaf+V.c_wood+V.c_root),
                cSoil=float(V.c_lit_cwd+V.c_lit_met+V.c_lit_str+V.c_lit_mic+V.c_soil_met+V.c_soil_str+V.c_soil_mic+V.c_soil_slow+V.c_soil_passive),
                ra=mra,
                rh=mrh,
            )
            # equivalent
            #o=np.array([
            #    np.sum(v[0:3]),
            #    np.sum(v[3:6]),
            #    np.sum(v[6:9]),
            #    mra,
            #    mrh,
            #])
            sols.append(o)
            
        sol=np.stack(sols)       
        return sol
        
    return param2res



# +
# now test it 
import matplotlib.pyplot as plt
from general_helpers import plot_solutions
const_params = cpa

param2res_sym = make_param2res_sym(const_params)
xs= param2res_sym(epa_0)

day_indices=month_2_day_index(range(cpa.nyears*12)),

fig = plt.figure()
plot_solutions(
        fig,
        #times=day_indices,
        times=range(cpa.nyears*12),
        var_names=observables._fields,
        tup=(xs,)
)
fig.savefig('solutions.pdf')

# -

# ### mcmc to optimize parameters 
# coming soon


