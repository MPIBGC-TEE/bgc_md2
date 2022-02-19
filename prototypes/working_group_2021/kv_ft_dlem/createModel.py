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

# ## DLEM model for TRENDY
#
# This is a simplified version of DLEM model with 3 vegetation pools: leaf, root and wood; 2 litter pools: aom1 (added organic matter) and aom2; 6 soil pools: smb1 (soil microbial pool), smb2, smr (soil microbial residues), dom (dissolved organcic matter), nom (native organcic matter), psom (passive organcic matter)

# +
from sympy import var, Symbol, Function  
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

# Make a small dictionary for the variables we will use
sym_dict={
    'C_leaf': '',
    'C_wood': '',
    'C_root': '',
    'C_aom1': '',
    'C_aom2': '',
    'C_smb1': '',
    'C_smb2': '',
    'C_smr': '',
    'C_nom': '',
    'C_dom': '',
    'C_psom': '',
    'r_C_leaf_2_C_aom1': '',
    'r_C_leaf_2_C_aom2': '',
    'r_C_wood_2_C_aom1': '',
    'r_C_wood_2_C_aom2': '',
    'r_C_root_2_C_aom1': '',
    'r_C_root_2_C_aom2': '',
    'r_C_aom1_2_C_smb1': '',
    'r_C_aom1_2_C_smb2': '',
    'r_C_aom1_2_C_dom': '',
    'r_C_aom1_2_C_nom': '',
    'r_C_aom2_2_C_smb1': '',
    'r_C_aom2_2_C_smb2': '',
    'r_C_aom2_2_C_dom': '',
    'r_C_smb2_2_C_smr': '',
    'r_C_smr_2_C_smb1': '',
    'r_C_dom_2_C_smb1': '',
    'r_C_dom_2_C_nom': '',
    'r_C_nom_2_C_dom': '',
    'r_C_nom_2_C_psom': '',
    'r_C_smb1_2_C_psom': '',
    'r_C_smb1_2_C_nom': '',
    'r_C_nom_2_C_smb1': '',
    'r_C_psom_2_C_smb1': '',
    'r_C_aom1_rh': '',
    'r_C_aom2_rh': '',
    'r_C_smb1_rh': '',
    'r_C_smb2_rh': '',
    'r_C_smr_rh': '',
    'r_C_dom_rh': '',
    'r_C_nom_rh': '',
    'r_C_psom_rh': '',
    'tsl': '',
    'mrso': '',
    't': '',
    'Theta_sat': '',
    'Theta_fc': '',
    'beta_leaf': '',
    'beta_wood': '',
}
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    #"I_vl": "Influx into vegetation leaf pool",
    #"k_vl_o": "out flux rate of leaf pool",
    'xi': 'a scalar function of temperature and moisture and thereby ultimately of time',
    'NPP': '',
   # 'tsl': '',
   # 'mrso': '',
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")
beta_root = 1.0- (beta_leaf+beta_wood)
mvs = CMTVS(
    {
        t,
        StateVariableTuple((
            #vl, 
            #vw, 
            C_leaf,
            C_wood,
            C_root,
            C_aom1,
            C_aom2,
            C_smb1,
            C_smb2,
            C_smr,
            C_nom,
            C_dom,
            C_psom,
        )),
        InFluxesBySymbol(
            {
                #vl: I_vl, vw: I_vw
                C_leaf: NPP(t) * beta_leaf, 
                C_root: NPP(t) * beta_root, 
                C_wood: NPP(t) * beta_wood
            }
        ),
        OutFluxesBySymbol(
            {
                #vl: k_vl_o * vl, vw: k_vw_o * vw
                C_aom1: r_C_aom1_rh*C_aom1*xi(t),
                C_aom2: r_C_aom2_rh*C_aom2*xi(t),
                C_smb1: r_C_smb1_rh*C_smb1*xi(t),
                C_smb2: r_C_smb2_rh*C_smb2*xi(t),
                C_smr: r_C_smr_rh*C_smr*xi(t),
                C_nom: r_C_nom_rh*C_nom*xi(t),
                C_dom: r_C_dom_rh*C_dom*xi(t),
                C_psom: r_C_psom_rh*C_psom*xi(t),
            }
        ),
        InternalFluxesBySymbol(
            {
                #(vl, vw): k_vl_2_vw * vl, (vw, vl): k_vw_2_vl * vw
                (C_leaf, C_aom1): r_C_leaf_2_C_aom1*C_leaf,
                (C_leaf, C_aom2): r_C_leaf_2_C_aom2*C_leaf,
                (C_wood, C_aom1): r_C_wood_2_C_aom1*C_wood,
                (C_wood, C_aom2): r_C_wood_2_C_aom2*C_wood,
                (C_root, C_aom1): r_C_root_2_C_aom1*C_root,
                (C_root, C_aom2): r_C_root_2_C_aom2*C_root,
                (C_aom1, C_smb1): r_C_aom1_2_C_smb1 * C_aom1*xi(t),
                (C_aom1, C_smb2): r_C_aom1_2_C_smb2 * C_aom1*xi(t),
                (C_aom1, C_dom): r_C_aom1_2_C_dom * C_aom1*xi(t),
                (C_aom1, C_nom): r_C_aom1_2_C_nom * C_aom1*xi(t),
                (C_aom2, C_smb1): r_C_aom2_2_C_smb1 * C_aom2*xi(t),
                (C_aom2, C_smb2): r_C_aom2_2_C_smb2 * C_aom2*xi(t),
                (C_aom2, C_dom): r_C_aom2_2_C_dom * C_aom2*xi(t),
                (C_smb2, C_smr): r_C_smb2_2_C_smr * C_smb2*xi(t),                
                (C_smr, C_smb1): r_C_smr_2_C_smb1 * C_smr*xi(t), 
                (C_dom, C_smb1): r_C_dom_2_C_smb1 * C_dom*xi(t), 
                (C_dom, C_nom): r_C_dom_2_C_nom * C_dom*xi(t),                
                (C_nom, C_dom): r_C_nom_2_C_dom * C_nom*xi(t),  
                (C_nom, C_psom): r_C_nom_2_C_psom * C_nom*xi(t),   
                (C_smb1, C_psom): r_C_smb1_2_C_psom * C_smb1*xi(t),
                (C_smb1, C_nom): r_C_smb1_2_C_nom * C_smb1*xi(t),
                (C_nom, C_smb1): r_C_nom_2_C_smb1 * C_nom*xi(t), 
                (C_psom, C_smb1): r_C_psom_2_C_smb1 * C_psom*xi(t),
            }
        ),
    },


    computers=module_computers(bgc_c)
)
# -

mvs.get_StateVariableTuple()

# we can also plot a picture
import bgc_md2.helper as h
h.compartmental_graph(mvs)

# we can also print the whole mass balance equation
import bgc_md2.display_helpers as dh
dh.mass_balance_equation(mvs)

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
xi_d=diag([1,1,1]+[xi(t) for i in range(8)],unpack=True)
xi_d

# We can go on and decompose N =\xi K -> K=\xi^{-1}N
K=xi_d.inv()*N
K
# we now have the (symbolic) ingredients for the tracebility analysis.
#xi_d,K,A

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
Observables = namedtuple(
    'Observables',
    ["cVeg", "cLitter", "cSoil", "rh" ]
)
Drivers=namedtuple(
    "Drivers",
    ["npp","mrso"]
)    
    
#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output():
    download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['DLEM'],
        variables = Observables._fields + Drivers._fields
    )
#call it to test that the download works the data
#download_my_TRENDY_output()


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
    t = s, 50, 33  # [t] = [:,49,325]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them 
        ds = nc.Dataset(str(path))
        return ds.variables[vn][t]

    o_names=[(f,"DLEM_S2_{}.nc".format(f)) for f in Observables._fields]
    d_names=[(f,"DLEM_S2_{}.nc".format(f)) for f in Drivers._fields]
    return (Observables(*map(f, o_names)),Drivers(*map(f,d_names)))


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
# and one for the internal fluxrates
# -

old_par_dict = {
    beta_leaf: 0.3,
    beta_wood: 0.4,
    Theta_sat: 0.1,
    Theta_fc:0.2,
    f_C_leaf_2_C_aom1: 0.5,
    f_C_leaf_2_C_aom2: 0.5,
    f_C_wood_2_C_aom1: 0.5,
    f_C_wood_2_C_aom2: 0.5,
    f_C_root_2_C_aom1: 0.5,
    f_C_root_2_C_aom2: 0.5,
           
    f_C_aom1_2_C_smb1: 0.2, 
    f_C_aom1_2_C_smb2: 0.2,
    f_C_aom1_2_C_dom: 0.2,
    f_C_aom1_2_C_nom: 0.2,
    f_C_aom2_2_C_smb1: 0.2,
    f_C_aom2_2_C_smb2: 0.2,
    f_C_aom2_2_C_dom: 0.2,
    f_C_smb2_2_C_smr: 0.7,
    f_C_smr_2_C_smb1: 0.7,
    f_C_dom_2_C_smb1: 0.4,
    f_C_dom_2_C_nom: 0.4,
    f_C_nom_2_C_dom: 0.2,
    f_C_nom_2_C_psom: 0.2,
    f_C_smb1_2_C_psom: 0.4,
    f_C_smb1_2_C_nom: 0.4,
    f_C_nom_2_C_smb1: 0.2,
    f_C_psom_2_C_smb1: 0.2,
    
    k_C_leaf: 1 / (60 * 2),
    k_C_wood: 1 / (365 * 30),
    k_C_root: 1 / (365 * 22),
    k_C_aom1: 1 / (365 * 4.5),
    k_C_aom2: 1 / (365 * 18),
    k_C_smb1: 1 / (365 * 1.5),
    k_C_smb2: 1 / (365 * 6),
    k_C_smr: 1 / (365 * 7.3),
    k_C_dom: 1 / (365 * 0.003),
    k_C_nom: 1 / (365 * 0.25),
    k_C_psom: 1 / (365 * 1.25),
}


# Now we can translate the old paramterisation to the new one.
#
# ### Providing dictionaries  for the parameters and functions.
#

par_dict={
    beta_leaf: 0.3,
    beta_wood: 0.4,
    Theta_sat: 0.1,
    Theta_fc: 0.2, 
}
par_dict.update(
    {Symbol(k):v.subs(old_par_dict) for k,v in all_rates.items()}
)
list(par_dict.keys())[4]
par_dict

# To be able to run the model forward we not only have to replace parameter symbols by values but symbolic functions by normal python functions.
# In our case the functions for $NPP$ and $\xi$ have to be provided. NPP_fun will interpolate the NPP for the day in question from the data. Which we have to load. 
# We will later store these functions in  `model_specific_helpers.py` which resides in the same folder as this notebook. You will have to adapt them to your data set. 

# +
from general_helpers import make_B_u_funcs_2, day_2_month_index
# check the numeric functions for B and u

def npp_func(day):
    month=day_2_month_index(day)
    return dvs.npp[month] * 86400   # kg/m2/s kg/m2/day

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
svs_0=Observables(*map(lambda v: v[0],svs))

X_0= np.array((
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cLitter/2,
    svs_0.cLitter/2,
    svs_0.cSoil/6,
    svs_0.cSoil/6,
    svs_0.cSoil/6,
    svs_0.cSoil/6,
    svs_0.cSoil/6,
    svs_0.cSoil/6,
))#.reshape(11,)
# in general B and u are nonlinear and can depend on X, thats why we have to test them with index and X arguments
u_func(0,X_0),B_func(0,X_0)

# +
import CompartmentalSystems.helpers_reservoir as hr

# To compute the ra and rh we have to some up the values for autotrophic and heterotrophic respiration we have to sum up outfluxes.
# We first create numerical functions for all the outfluxes from our symbolic description.
# lets do it for the outflux of one pool first (and then for the whole lot of them)
expr_cont=mvs.get_OutFluxesBySymbol()[C_aom1]

# this is the momentary outflux as a function of t,C_leaf_litter as it would occure in the differential equation
expr_cont
# -

# what we want is acutally the accumulated flux over one timestep in the simplest approximation (euler forward)
# the framework has a helper function to create an euler forward discretisation of a flux
delta_t=Symbol("delta_t")
it=Symbol("it") #arbitrary symbol (in our case it=days_since_start )
expr_disc=hr.euler_forward_net_flux_sym(
    flux_sym_cont=expr_cont,
    cont_time=t,
    delta_t=delta_t,
    iteration=it #arbitrary
)
# If you look at the result you see that the euler forwar approximation assumes the flux at time t to be constant for the timestep  
expr_disc

# if we assume tat delta_t is 1 day and it counts days 
# it becomes even simpler
expr_disc.subs({delta_t:1})
#Which is the same as if we had 't' replaced in the above formula wiht it

# +
# this expression we turn now into a numeric funciton of it
# although in our example it depends only on t and C_leaf_litter we make it a function of ALL state variables to be able to call it in the same way as u_func and B_func 
argtup=(mvs.get_TimeSymbol(), *mvs.get_StateVariableTuple())

C_leaf_litter_func = hr.numerical_function_from_expression(
    expr=expr_cont,
    tup=argtup, 
    parameter_dict=par_dict,
    func_set=func_dict
)
# call it for testing
C_leaf_litter_func(0,*X_0)


# +
#mass production of output functions
def numfunc(expr):
    return hr.numerical_function_from_expression(
    expr=expr,
    tup=(mvs.get_TimeSymbol(), *mvs.get_StateVariableTuple()),
    parameter_dict=par_dict,
    func_set=func_dict
)
    
numOutFluxesBySymbol={k:numfunc(expr) for k,expr in mvs.get_OutFluxesBySymbol().items()} 

# now we can compute ra 
# apply all the outfluxfunctions of the veg pools and add up the result
ra_0=np.sum(
    [
        numOutFluxesBySymbol[k](0,*X_0) 
        for k in [C_leaf,C_wood,C_root] 
        if k in numOutFluxesBySymbol.keys()
    ]
)
rh_0=np.sum(
    [
        numOutFluxesBySymbol[k](0,*X_0) 
        for k in [C_aom1,C_aom2,C_smb1,C_smb2,C_smr,C_nom,C_dom,C_psom] 
        if k in numOutFluxesBySymbol.keys()
    ]
)
ra_0,rh_0
# -

OutFluxesBySymbol

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
    ["rh"]
)
StartVector._fields

V_init= StartVector(
    C_leaf=svs_0.cVeg/3,
    C_wood=svs_0.cVeg/3,
    C_root=svs_0.cVeg/3,
    C_aom1=svs_0.cLitter/2,
    C_aom2=svs_0.cLitter/2,
    C_smb1=svs_0.cSoil/6,
    C_smb2=svs_0.cSoil/6,
    C_smr=svs_0.cSoil/6,
    C_dom=svs_0.cSoil/6,
    C_nom=svs_0.cSoil/6,
    C_psom=svs_0.cSoil/6,
    rh=svs_0.rh * 86400   # kg/m2/s kg/m2/day        
)
#V_init.__getattribute__("rh")
V_init

mvs.get_StateVariableTuple()

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
        [V_init.rh]
    ).reshape(n+1,1) #reshaping is neccessary for matmux

    def f(it,V):
        X = V[0:n]
        b = u_func(it,X)
        B = B_func(it,X)
        #outfluxes = B @ X
        X_new = X + b + B @ X
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
        rh = np.sum(
            [
                numOutFluxesBySymbol[k](0,*X_0) 
                for k in [C_aom1,C_aom2,C_smb1,C_smb2,C_smr,C_nom,C_dom,C_psom] 
                if k in numOutFluxesBySymbol.keys()
            ]
        )  
        #print(numOutFluxesBySymbol.keys())
#         for k in [C_aom1,C_aom2,C_smb1,C_smb2,C_smr,C_nom,C_dom,C_psom]:
#             print(numOutFluxesBySymbol[k](0,*X_0))
        V_new = np.concatenate((X_new.reshape(n,1),rh.reshape(1,1)), axis=0)
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )

# def make_steady_state_iterator_sym(
#         mvs,
#         V_init: StartVector,
#         par_dict,
#         func_dict
#     ):
#     B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
#     sv=mvs.get_StateVariableTuple()
#     n=len(sv)
#     print(n)
#     # build an array in the correct order of the StateVariables which in our case is already correct 
#     # since V_init was built automatically but in case somebody handcrafted it and changed
#     # the order later in the symbolic formulation....
#     V_arr=np.array(
#         [V_init.__getattribute__(str(v)) for v in sv]+
#         [V_init.rh]
#     ).reshape(n+1,1) #reshaping is neccessary for matmul
    
#     def f(it,V):
#         #V=V_arr    
    
#         X = V[0:n]
#         #print(X)
#         b = u_func(it,X)
#         #print(b)
#         B = B_func(it,X)
        
#         V_new = np.concatenate((b.reshape(n,1),B.reshape(n,n)), axis=1)
#         #print("V_new shape: ", V_new.shape)
#         return V_new
    
#     return TimeStepIterator2(
#         initial_values=V_arr,
#         f=f)


# -

V_init

np.array(V_init).shape

# +
# # calculate steady state
    
# it_sym = make_steady_state_iterator_sym(
#     mvs,
#     V_init=V_init,
#     par_dict=par_dict,
#     func_dict=func_dict
# )
# # we will run the model for 15 steps
# ns=10
# #res= np.zeros((ns,len(V_init)))
# res=np.zeros((ns,11,12))
# #print(res)
# res_sym = copy(res)
# it_sym.__next__()
# for i in range(ns):
#     #print("it_sym... ")
#     #print(it_sym.__next__().shape)
#     res_sym[i,:,:]=it_sym.__next__()#.reshape(len(V_init),)
# #print("res_sym: ",res_sym[:,:,0])
# #print(res_sym.shape)

# # calculate steady state

# b=res_sym[:,:,0]
# B=res_sym[:,:,1:12]

# B_mean = np.mean(B, axis = 0)
# b_mean = np.mean(b, axis = 0).reshape(11,1)

# #print("B_mean: ", B_mean.shape)
# #print("b_mean: ", b_mean.shape)

# X_ss = np.linalg.solve(B_mean, (b_mean))

# print("steady: ", X_ss)
print("NPP: ", dvs.npp[0:12])
b=u_func(180,X_0)
print("b: ",b)
B=B_func(180,X_0)
print ("B: ",B)
X_ss = np.linalg.solve(B, (-b))
print("steady: ", X_ss)

# +
# deriving initial pool proportions
print("Total pools: ", svs_0)
fraction_leaf=X_ss[0]/(X_ss[0]+X_ss[1]+X_ss[2])
fraction_wood=X_ss[1]/(X_ss[0]+X_ss[1]+X_ss[2])
fraction_root=X_ss[2]/(X_ss[0]+X_ss[1]+X_ss[2])
fraction_aom1=X_ss[3]/(X_ss[3]+X_ss[4])
fraction_aom2=X_ss[4]/(X_ss[3]+X_ss[4])
fraction_smb1=X_ss[5]/(X_ss[5]+X_ss[6]+X_ss[7]+X_ss[8]+X_ss[9]+X_ss[10])
fraction_smb2=X_ss[6]/(X_ss[5]+X_ss[6]+X_ss[7]+X_ss[8]+X_ss[9]+X_ss[10])
fraction_smr=X_ss[7]/(X_ss[5]+X_ss[6]+X_ss[7]+X_ss[8]+X_ss[9]+X_ss[10])
fraction_nom=X_ss[8]/(X_ss[5]+X_ss[6]+X_ss[7]+X_ss[8]+X_ss[9]+X_ss[10])
fraction_dom=X_ss[9]/(X_ss[5]+X_ss[6]+X_ss[7]+X_ss[8]+X_ss[9]+X_ss[10])
fraction_psom=X_ss[10]/(X_ss[5]+X_ss[6]+X_ss[7]+X_ss[8]+X_ss[9]+X_ss[10])

# initial pools
C_leaf_0=svs_0.cVeg*fraction_leaf
C_wood_0=svs_0.cVeg*fraction_wood
C_root_0=svs_0.cVeg*fraction_root
C_aom1_0=svs_0.cLitter*fraction_aom1
C_aom2_0=svs_0.cLitter*fraction_aom2
C_smb1_0=svs_0.cSoil*fraction_smb1
C_smb2_0=svs_0.cSoil*fraction_smb2
C_smr_0=svs_0.cSoil*fraction_smr
C_nom_0=svs_0.cSoil*fraction_nom
C_dom_0=svs_0.cSoil*fraction_dom
C_psom_0=svs_0.cSoil*fraction_psom
C_leaf_0
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
        "cVeg_0",
        "cLitter_0",
        "cSoil_0",
        "npp_0",
        "rh_0",
       # "mrso",
       # "tsl",
        "number_of_months" # necessary to prepare the output in the correct lenght 
    ]
)
# Note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues and their length (number of months)
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 
# It is better to start with only a few

# deriving Estimated Parameters from par_dict
temp_list=list()
for name,dict_ in par_dict.items():
    temp_list.append(name)
# adding initial pool values that need to be estimated for TRENDY
temp_list.append("C_leaf_0")
temp_list.append("C_wood_0")
temp_list.append("C_aom1_0")
temp_list.append("C_smb1_0")
temp_list.append("C_smb2_0")
temp_list.append("C_smr_0")
temp_list.append("C_nom_0")
temp_list.append("C_dom_0")

# fix non-string elements
temp_list = [str(element) for element in temp_list]
EstimatedParameters=namedtuple("EstimatedParameters",field_names=temp_list)

# EstimatedParameters = namedtuple(
#     "EstimatedParameters",[ 
#         "beta_leaf",
#         "beta_wood",
#         "Theta_sat",
#         "Theta_fc",
#         "r_C_leaf_rh",
#         "r_C_wood_rh",
#         "r_C_root_rh",
#         "r_C_leaf_litter_rh",
#         "r_C_wood_litter_rh",
#         "r_C_root_litter_rh",
#         "r_C_soil_fast_rh",
#         "r_C_soil_slow_rh",
#         "r_C_soil_passive_rh",
#         "r_C_leaf_2_C_leaf_litter",
#         "r_C_wood_2_C_wood_litter",
#         "r_C_root_2_C_root_litter",
#         "r_C_leaf_litter_2_C_soil_fast",
#         "r_C_leaf_litter_2_C_soil_slow",
#         "r_C_leaf_litter_2_C_soil_passive",
#         "r_C_wood_litter_2_C_soil_fast",
#         "r_C_wood_litter_2_C_soil_slow",
#         "r_C_wood_litter_2_C_soil_passive",
#         "r_C_root_litter_2_C_soil_fast",
#         'C_leaf_0',#for the trendy data also the startvalues have to be estimated but 
#         'C_wood_0',
#         #C_root_0 can be inferred as cVeg_0-(C_leaf_0+C_wood_0)
#         'C_leaf_litter_0',
#         'C_wood_litter_0',
#         #C_root_litter_0 can be inferred
#         'C_soil_fast_0',
#         'C_soil_slow_0',
#         #C_soil_passive_0 can be inferred 
#     ]
# )

# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues. 
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 


# -

[ par_dict[key] for key in [beta_leaf,beta_wood]]
#par_dict.keys()

UnEstimatedParameters._fields

cpa=UnEstimatedParameters(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 npp_0=dvs.npp[0] * 86400,   # kg/m2/s kg/m2/day
 rh_0=svs_0.rh * 86400,   # kg/m2/s kg/m2/day
 #mrso_0=dvs.mrso[0],
 #tsl_0=dvs.tsl[0],
 number_of_months=120#len(svs.rh)
)
print(cpa)

# check initial values for the pool sizes
V_init

# +
# create a start parameter tuple for the mcmc. The order has to be the same as when you created the namedtupl3 
# If you don't you get a "TypeError". 

# deriving epa_0 (initial Estimated Parameter values) from par_dict
temp_list=list()
for name,dict_ in par_dict.items():
    temp_list.append(dict_)

# adding initial pool values that need to be estimated for TRENDY
temp_list.append(C_leaf_0) # C_leaf_0
temp_list.append(C_wood_0) # C_wood_0
temp_list.append(C_aom1_0) # C_aom1_0
temp_list.append(C_smb1_0) # C_smb1_0
temp_list.append(C_smb2_0) # C_smb2_0
temp_list.append(C_smr_0) # C_smr_0
temp_list.append(C_nom_0) # C_nom_0
temp_list.append(C_dom_0) # C_dom_0

epa_0=EstimatedParameters(*temp_list)

# epa_0=EstimatedParameters(
#  beta_leaf=0.6,
#  beta_wood=0.25,
#  T_0=2,
#  E=4,
#  KM=10,
#  r_C_leaf_rh=0,
#  r_C_wood_rh=0,
#  r_C_root_rh=0,
#  r_C_leaf_litter_rh=0.000415110004151100,
#  r_C_wood_litter_rh=0.000124533001245330,
#  r_C_root_litter_rh=0.000122042341220423,
#  r_C_soil_fast_rh=0.000152207001522070,
#  r_C_soil_slow_rh=2.73972602739726e-5,
#  r_C_soil_passive_rh=7.82778864970646e-6,
#  r_C_leaf_2_C_leaf_litter=0.00833333333333333,
#  r_C_wood_2_C_wood_litter=9.13242009132420e-5,
#  r_C_root_2_C_root_litter=0.000124533001245330,
#  r_C_leaf_litter_2_C_soil_fast=0.000340390203403902,
#  r_C_leaf_litter_2_C_soil_slow=5.81154005811540e-5,
#  r_C_leaf_litter_2_C_soil_passive=1.66044001660440e-5,
#  r_C_wood_litter_2_C_soil_fast=7.47198007471980e-5,
#  r_C_wood_litter_2_C_soil_slow=2.98879202988792e-5,
#  r_C_wood_litter_2_C_soil_passive=1.99252801992528e-5,
#  r_C_root_litter_2_C_soil_fast=7.47198007471980e-5,
#  C_leaf_0=svs_0.cVeg/3,
#  C_wood_0=svs_0.cVeg/3,
#  C_leaf_litter_0=svs_0.cLitter/3,
#  C_wood_litter_0=svs_0.cLitter/3,
#  C_soil_fast_0=svs_0.cSoil/3,
#  C_soil_slow_0=svs_0.cSoil/3,
# )
# -

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
from general_helpers import month_2_day_index, monthly_to_yearly # convert all observables to yearly timescale
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
        return dvs.npp[month] * 86400   # kg/m2/s kg/m2/day
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
        V_init = StartVector(
            C_leaf=epa.C_leaf_0,
            C_wood=epa.C_wood_0,
            C_root=cpa.cVeg_0-(epa.C_leaf_0 + epa.C_wood_0), #X_ss[2],
            C_aom1=epa.C_aom1_0,
            C_aom2=cpa.cLitter_0-epa.C_aom1_0, #X_ss[4],
            C_smb1=epa.C_smb1_0,
            C_smb2=epa.C_smb2_0,
            C_smr=epa.C_smr_0,
            C_nom=epa.C_nom_0,
            C_dom=epa.C_dom_0,
            C_psom=cpa.cSoil_0-(epa.C_smb1_0+epa.C_smb2_0+epa.C_smr_0+epa.C_nom_0+epa.C_dom_0), #X_ss[10],
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
        for m in range(cpa.number_of_months):
            dpm = days_per_month[ m % 12]  
            #mra=0
            mrh=0
            for d in range(dpm):
                v = it_sym.__next__()
                #mra +=v[9,0]
                mrh +=v[11,0]
            mrh=mrh/dpm
            V=StartVector(*v)
            o=Observables(
                cVeg=float(V.C_leaf+V.C_wood+V.C_root),
                cLitter=float(V.C_aom1+V.C_aom2),
                cSoil=float(V.C_smb1+V.C_smb2+V.C_smr+V.C_nom+V.C_dom+V.C_psom),
                #ra=mra,
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
        mod=np.zeros(int(cpa.number_of_months/12)*sol.shape[1]).reshape([int(cpa.number_of_months/12),sol.shape[1]])  
        for i in range(sol.shape[1]):
            mod[:,i]=monthly_to_yearly(sol[:,i])
        return mod

    return param2res


# -

param2res_sym = make_param2res_sym(cpa)
xs= param2res_sym(epa_0)
print(xs.shape)
print(xs)
#print(cpa)

obs=np.column_stack((np.array(svs.cVeg),np.array(svs.cLitter),np.array(svs.cSoil),monthly_to_yearly(np.array(svs.rh))))
obs=obs[0:int(cpa.number_of_months/12),:]
#print(obs.shape)
#print(obs[:,3])
obs[:,3]=obs[:,3]*86400
#print(obs[:,3])

# +
# now test it 
import matplotlib.pyplot as plt
from general_helpers import plot_solutions

param2res_sym = make_param2res_sym(cpa)
mod= param2res_sym(epa_0)

# convert model output to yearly
# mod=np.zeros(obs.shape[0]*obs.shape[1]).reshape([obs.shape[0],obs.shape[1]])  
# for i in range(obs.shape[1]):
#     mod[:,i]=monthly_to_yearly(xs[:,i])

day_indices=month_2_day_index(range(cpa.number_of_months)),

fig = plt.figure()
plot_solutions(
        fig,
        #times=day_indices,
        times=range(int(cpa.number_of_months/12)),
        var_names=Observables._fields,
        tup=(mod,obs)
)
fig.savefig('solutions.pdf')

# -

# ### mcmc to optimize parameters 
# coming soon

# +
c_max = np.array(
    EstimatedParameters(
        beta_leaf=0.3, 
        beta_wood=0.4, 
        Theta_sat=0.1, 
        Theta_fc=0.2, 
        r_C_leaf_rh=0, 
        r_C_wood_rh=0, 
        r_C_root_rh=0, 
        r_C_aom1_rh=0.000121765601217656, 
        r_C_aom2_rh=6.08828006088280e-5, 
        r_C_smb1_rh=0.000365296803652968, 
        r_C_smb2_rh=0.000136986301369863, 
        r_C_smr_rh=0.000112591480577970, 
        r_C_nom_rh=0.00438356164383562, 
        r_C_dom_rh=0.182648401826484, 
        r_C_psom_rh=0.00175342465753425, 
        r_C_leaf_2_C_aom1=0.00416666666666667, 
        r_C_leaf_2_C_aom2=0.00416666666666667, 
        r_C_wood_2_C_aom1=4.56621004566210e-5, 
        r_C_wood_2_C_aom2=4.56621004566210e-5, 
        r_C_root_2_C_aom1=6.22665006226650e-5, 
        r_C_root_2_C_aom2=6.22665006226650e-5, 
        r_C_aom1_2_C_smb1=0.000121765601217656, 
        r_C_aom1_2_C_smb2=0.000121765601217656, 
        r_C_aom1_2_C_nom=0.000121765601217656, 
        r_C_aom1_2_C_dom=0.000121765601217656, 
        r_C_aom2_2_C_smb1=3.04414003044140e-5, 
        r_C_aom2_2_C_smb2=3.04414003044140e-5, 
        r_C_aom2_2_C_dom=3.04414003044140e-5, 
        r_C_smb1_2_C_nom=0.000730593607305936, 
        r_C_smb1_2_C_psom=0.000730593607305936, 
        r_C_smb2_2_C_smr=0.000319634703196347, 
        r_C_smr_2_C_smb1=0.000262713454681929, 
        r_C_nom_2_C_smb1=0.00219178082191781, 
        r_C_nom_2_C_dom=0.00219178082191781, 
        r_C_nom_2_C_psom=0.00219178082191781, 
        r_C_dom_2_C_smb1=0.365296803652968, 
        r_C_dom_2_C_nom=0.365296803652968, 
        r_C_psom_2_C_smb1=0.000438356164383562, 
        C_leaf_0= 0.0017252, 
        C_wood_0=0.2098991, 
        C_aom1_0= 0.16313342, 
        C_smb1_0= 3.30185076, 
        C_smb2_0= 4.68136884, 
        C_smr_0= 3.98696579, 
        C_nom_0= 0.43008121, 
        C_dom_0= 0.00337288
    )
)

c_min = np.array(
    EstimatedParameters(
        beta_leaf=0.3, 
        beta_wood=0.4, 
        Theta_sat=0.1, 
        Theta_fc=0.2, 
        r_C_leaf_rh=0, 
        r_C_wood_rh=0, 
        r_C_root_rh=0, 
        r_C_aom1_rh=0.000121765601217656, 
        r_C_aom2_rh=6.08828006088280e-5, 
        r_C_smb1_rh=0.000365296803652968, 
        r_C_smb2_rh=0.000136986301369863, 
        r_C_smr_rh=0.000112591480577970, 
        r_C_nom_rh=0.00438356164383562, 
        r_C_dom_rh=0.182648401826484, 
        r_C_psom_rh=0.00175342465753425, 
        r_C_leaf_2_C_aom1=0.00416666666666667, 
        r_C_leaf_2_C_aom2=0.00416666666666667, 
        r_C_wood_2_C_aom1=4.56621004566210e-5, 
        r_C_wood_2_C_aom2=4.56621004566210e-5, 
        r_C_root_2_C_aom1=6.22665006226650e-5, 
        r_C_root_2_C_aom2=6.22665006226650e-5, 
        r_C_aom1_2_C_smb1=0.000121765601217656, 
        r_C_aom1_2_C_smb2=0.000121765601217656, 
        r_C_aom1_2_C_nom=0.000121765601217656, 
        r_C_aom1_2_C_dom=0.000121765601217656, 
        r_C_aom2_2_C_smb1=3.04414003044140e-5, 
        r_C_aom2_2_C_smb2=3.04414003044140e-5, 
        r_C_aom2_2_C_dom=3.04414003044140e-5, 
        r_C_smb1_2_C_nom=0.000730593607305936, 
        r_C_smb1_2_C_psom=0.000730593607305936, 
        r_C_smb2_2_C_smr=0.000319634703196347, 
        r_C_smr_2_C_smb1=0.000262713454681929, 
        r_C_nom_2_C_smb1=0.00219178082191781, 
        r_C_nom_2_C_dom=0.00219178082191781, 
        r_C_nom_2_C_psom=0.00219178082191781, 
        r_C_dom_2_C_smb1=0.365296803652968, 
        r_C_dom_2_C_nom=0.365296803652968, 
        r_C_psom_2_C_smb1=0.000438356164383562, 
        C_leaf_0= 0.0017252, 
        C_wood_0=0.2098991, 
        C_aom1_0= 0.16313342, 
        C_smb1_0= 3.30185076, 
        C_smb2_0= 4.68136884, 
        C_smr_0= 3.98696579, 
        C_nom_0= 0.43008121, 
        C_dom_0= 0.00337288
    )
)


# +
from general_helpers import autostep_mcmc
from model_specific_helpers import make_param_filter_func

isQualified = make_param_filter_func(c_max, c_min)
param2res = make_param2res_sym(cpa)


# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa_0,
    filter_func=isQualified,
    param2res=param2res,
    costfunction=make_feng_cost_func(obs),
    nsimu=100,
    c_max=c_max,
    c_min=c_min,
    acceptance_rate=10,   # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=1   # default value | increase value to reduce initial step size
)
# save the parameters and cost function values for postprocessing
#pd.DataFrame(C_autostep).to_csv(dataPath.joinpath('visit_autostep_da_aa.csv'), sep=',')
#pd.DataFrame(J_autostep).to_csv(dataPath.joinpath('visit_autostep_da_j_aa.csv'), sep=',')
