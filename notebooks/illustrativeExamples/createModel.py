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

# This illustrative notebook shows how to create a representation for a new model.
#
# For the sake of simplicity we assume here that we start with a description of pools and fluxes.
# This is the form the models are usually described in the literature.
# We will see that the framework can derive a matrix representation automatically. 
#
# It would also be possible to start from the matrices or even
# mix the two approaches. 
# We will point to some more advanced examples where you can see this in action. If you get familiar with the framework you will find many more combinations than we can cover. 
#
# ## Inspect a minimal model
#
# We will start with an extremely simple model.
# whose original you can find in 
# ```bash
# bgc_md2/src/bgc_md2/models/testVectorFree/source.py
# ```
# or via the following code:

import inspect
import bgc_md2.models.testVectorFree.source 
print(inspect.getsource(bgc_md2.models.testVectorFree.source))

# # copy its contents into a new cell. (we can later save it under a new name) 

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
from importlib import import_module

# Make a small dictionary for the variables we will use
sym_dict={
    "vl": "vegegation leaf pool",
    "vw": "vegetation wood pool content",
    "I_vw": "Influx into vegetation wood pool",
    "k_vw_o": "out flux rate of wood pool",
    "k_vl_2_vw": "internal flux rate from leaf to wood", 
    "k_vw_2_vl": "internal flux rate from wood to leaf", 
}
# Make symbols from  the strings that we can later use in expressions  
# vl, vw,...
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# We will also use some symbolic functions ("symbols" with an argument) 
func_dict={
    "I_vl": "Influx into vegetation leaf pool",
    "k_vl_o": "out flux rate of leaf pool",
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")

mvs = CMTVS(
    {
        StateVariableTuple((vl, vw)),
        t,
        InFluxesBySymbol({vl: I_vl(t), vw: I_vw}),
        OutFluxesBySymbol({vl: k_vl_o(t) * vl, vw: k_vw_o * vw}),
        InternalFluxesBySymbol({(vl, vw): k_vl_2_vw * vl, (vw, vl): k_vw_2_vl * vw}),
    },

    computers=module_computers(bgc_md2.resolve.computers)
)
# -

# The last statement in the code defines a variable `mvs` which is 
# an instance of CMTVS which stands for `C`onnected`M`ulti`T`ype`V`ariable`S`et".
# It contains information in two forms. 
# 1. Variables of certain types (like InFluxesBySymbol)
# 2. Computers, Functions that "connect" these Variables and to other results we did not specify but which can be computed.
#    
#
#
# Taks:
# To see what it can do with this information add a new cell an type `mvs.` and press the `tab` key `->`. This will show you the available methods, in other words what can be computed from the provided information.
#
# I wanted to see the compartmental the pools, the matrix and the inputs.

mvs.get_StateVariableTuple()

mvs.get_CompartmentalMatrix()

mvs.get_InputTuple()

# we can also print the whole mass balance equation
import bgc_md2.display_helpers as dh
dh.mass_balance_equation(mvs)

# we can also plot a picture
import bgc_md2.helper as h
h.compartmental_graph(mvs)

# +
# adjust the output to full width
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# make changes to imported files immidiately available 
# avoiding the need to reload (in most cases)
# %load_ext autoreload
# %autoreload 2
# -

# ## Extend the minimal model
#
# ### add more state variables and add a small description
#
# In a first step I will 
# - add new symbols to desricbe pools and parameters for the VISIT  model (from Kostias code) but
#   leave the old ones there since they are presently used in the fluxes.
#   We will replace them one bye one later 
#   but this way we are not forced to replace all the fluxes at
#   once.It's always good to be able to make small steps...
#

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

# Make a small dictionary for the variables we will use
sym_dict={
    "vl": "vegegation leaf pool",
    "vw": "vegetation wood pool content",
    "k_vl_2_vw": "internal flux rate from leaf to wood", 
    "k_vw_2_vl": "internal flux rate from wood to leaf", 
    'C_soil_fast': '',
    'C_soil_slow': '',
    'C_soil_passive': '',
    'C_leaf': '',
    'C_root': '',
    'C_wood': '',
    'C_leaf_litter': '',
    'C_root_litter': '',
    'C_wood_litter': '',
    'k_leaf_2_leaf_litter': '',
    'k_root_2_root_litter': '',
    'k_wood_2_wood_litter': '',
    'k_leaf_litter_rh': '',
    'k_root_litter_rh': '',
    'k_wood_litter_rh': '',
    'k_soil_fast_rh': '',
    'k_soil_slow_rh': '',
    'k_soil_passive_rh': '',
    'k_leaf_litter_2_soil_fast': '',
    'k_leaf_litter_2_soil_slow': '',
    'k_leaf_litter_2_soil_passive': '',
    'k_wood_litter_2_soil_fast': '',
    'k_wood_litter_2_soil_slow': '',
    'k_wood_litter_2_soil_passive': '',
    'k_root_litter_2_soil_fast': '',
    'k_root_litter_2_soil_slow': '',
    'k_root_litter_2_soil_passive': '',
    'tsl': '',
    'mrso': '',
    't': '',
    'T_0': '',
    'E': '',
    'KM': '',
    'beta_leaf': '',
    'beta_wood': '',
}
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    "I_vl": "Influx into vegetation leaf pool",
    "k_vl_o": "out flux rate of leaf pool",
    'xi': '',
    'NPP': '',
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)



t=TimeSymbol("t")

mvs = CMTVS(
    {
        StateVariableTuple((vl, vw)),
        t,
        InFluxesBySymbol({vl: I_vl(t), vw: I_vw}),
        OutFluxesBySymbol({vl: k_vl_o(t) * vl, vw: k_vw_o * vw}),
        InternalFluxesBySymbol({(vl, vw): k_vl_2_vw * vl, (vw, vl): k_vw_2_vl * vw}),
    },

    computers=module_computers(bgc_c)
)

# -

# Nothing has changed in the model description but be have some more symbols to work with.
# We can now use the new symbols in expression without errors since they are now defined.
# In the next step we add the new state variables.
# For the moment we leave the old ones (vl and vw) in to have a template for how to formulate fluxes)

# +
from sympy import var 
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
    "vl": "vegegation leaf pool",
    "vl": "vegegation leaf pool",
    "vw": "vegetation wood pool content",
    "I_vw": "Influx into vegetation wood pool",
    "k_vw_o": "out flux rate of wood pool",
    "k_vl_2_vw": "internal flux rate from leaf to wood", 
    "k_vw_2_vl": "internal flux rate from wood to leaf", 
    'C_soil_fast': '',
    'C_soil_slow': '',
    'C_soil_passive': '',
    'C_leaf': '',
    'C_root': '',
    'C_wood': '',
    'C_leaf_litter': '',
    'C_root_litter': '',
    'C_wood_litter': '',
    'k_leaf_2_leaf_litter': '',
    'k_root_2_root_litter': '',
    'k_wood_2_wood_litter': '',
    'k_leaf_litter_rh': '',
    'k_root_litter_rh': '',
    'k_wood_litter_rh': '',
    'k_soil_fast_rh': '',
    'k_soil_slow_rh': '',
    'k_soil_passive_rh': '',
    'k_leaf_litter_2_soil_fast': '',
    'k_leaf_litter_2_soil_slow': '',
    'k_leaf_litter_2_soil_passive': '',
    'k_wood_litter_2_soil_fast': '',
    'k_wood_litter_2_soil_slow': '',
    'k_wood_litter_2_soil_passive': '',
    'k_root_litter_2_soil_fast': '',
    'k_root_litter_2_soil_slow': '',
    'k_root_litter_2_soil_passive': '',
    'tsl': '',
    'mrso': '',
    't': '',
    'T_0': '',
    'E': '',
    'KM': '',
    'beta_leaf': '',
    'beta_wood': '',
}
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    "I_vl": "Influx into vegetation leaf pool",
    "k_vl_o": "out flux rate of leaf pool",
    'xi': '',
    'NPP': '',
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t =  TimeSymbol("t")
mvs = CMTVS(
    {
        t,
        StateVariableTuple((
            vl, 
            vw, 
            C_leaf,
	        C_wood,
	        C_root,
	        C_leaf_litter,
	        C_wood_litter,
	        C_root_litter,
	        C_soil_fast,
	        C_soil_slow,
	        C_soil_passive,
        )),
        InFluxesBySymbol({vl: I_vl(t), vw: I_vw}),
        OutFluxesBySymbol({vl: k_vl_o(t) * vl, vw: k_vw_o * vw}),
        InternalFluxesBySymbol({(vl, vw): k_vl_2_vw * vl, (vw, vl): k_vw_2_vl * vw}),
    },

    computers=module_computers(bgc_c)
)
# -

h.compartmental_graph(mvs)

dh.mass_balance_equation(mvs)

# We can see only the old fluxes.
# Now we add more fluxes according to the description of the model and remove the old ones, along with the symbols and fucntions that we do no longer need
#
#

# +
from sympy import var 
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
    # "vl": "vegegation leaf pool",
    # "vl": "vegegation leaf pool",
    # "vw": "vegetation wood pool content",
    # "I_vw": "Influx into vegetation wood pool",
    # "k_vw_o": "out flux rate of wood pool",
    # "k_vl_2_vw": "internal flux rate from leaf to wood", 
    # "k_vw_2_vl": "internal flux rate from wood to leaf", 
    'C_soil_fast': '',
    'C_soil_slow': '',
    'C_soil_passive': '',
    'C_leaf': '',
    'C_root': '',
    'C_wood': '',
    'C_leaf_litter': '',
    'C_root_litter': '',
    'C_wood_litter': '',
    'k_leaf_2_leaf_litter': '',
    'k_root_2_root_litter': '',
    'k_wood_2_wood_litter': '',
    'k_leaf_litter_rh': '',
    'k_root_litter_rh': '',
    'k_wood_litter_rh': '',
    'k_soil_fast_rh': '',
    'k_soil_slow_rh': '',
    'k_soil_passive_rh': '',
    'k_leaf_litter_2_soil_fast': '',
    'k_leaf_litter_2_soil_slow': '',
    'k_leaf_litter_2_soil_passive': '',
    'k_wood_litter_2_soil_fast': '',
    'k_wood_litter_2_soil_slow': '',
    'k_wood_litter_2_soil_passive': '',
    'k_root_litter_2_soil_fast': '',
    'k_root_litter_2_soil_slow': '',
    'k_root_litter_2_soil_passive': '',
    'tsl': '',
    'mrso': '',
    't': '',
    'T_0': '',
    'E': '',
    'KM': '',
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
	        C_leaf_litter,
	        C_wood_litter,
	        C_root_litter,
	        C_soil_fast,
	        C_soil_slow,
	        C_soil_passive,
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
                C_leaf_litter: k_leaf_litter_rh*C_leaf_litter*xi(t),
                C_wood_litter: k_wood_litter_rh*C_wood_litter*xi(t),
                C_root_litter: k_root_litter_rh*C_root_litter*xi(t),
                C_soil_fast: k_soil_fast_rh*C_soil_fast*xi(t),
                C_soil_slow: k_soil_slow_rh*C_soil_slow*xi(t),
                C_soil_passive: k_soil_passive_rh*C_soil_passive*xi(t),
            }
        ),
        InternalFluxesBySymbol(
            {
                #(vl, vw): k_vl_2_vw * vl, (vw, vl): k_vw_2_vl * vw
                (C_leaf, C_leaf_litter): k_leaf_2_leaf_litter*C_leaf, 
                (C_wood, C_wood_litter): k_wood_2_wood_litter*C_wood, 
                (C_root, C_root_litter): k_root_2_root_litter*C_root, 
                (C_leaf_litter, C_soil_fast)    : k_leaf_litter_2_soil_fast * C_leaf_litter*xi(t),
                (C_leaf_litter, C_soil_slow)    : k_leaf_litter_2_soil_slow * C_leaf_litter*xi(t),
                (C_leaf_litter, C_soil_passive) : k_leaf_litter_2_soil_passive * C_leaf_litter*xi(t),
                (C_wood_litter, C_soil_fast)    : k_wood_litter_2_soil_fast * C_wood_litter*xi(t),
                (C_wood_litter, C_soil_slow)    : k_wood_litter_2_soil_slow * C_wood_litter*xi(t),
                (C_wood_litter, C_soil_passive) : k_wood_litter_2_soil_passive * C_wood_litter*xi(t),
                (C_root_litter, C_soil_fast)    : k_root_litter_2_soil_fast * C_root_litter*xi(t),
                (C_root_litter, C_soil_slow)    : k_root_litter_2_soil_slow * C_root_litter*xi(t),
                (C_root_litter, C_soil_passive) : k_root_litter_2_soil_passive * C_root_litter*xi(t),
            }
        ),
    },


    computers=module_computers(bgc_c)
)
# -
h.compartmental_graph(mvs)

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
xi_d=diag([1,1,1]+[xi(t) for i in range(6)],unpack=True)
xi_d

# We can go on and decompose N =\xi K -> K=\xi^{-1}N
K=xi_d.inv()*N
K
# we now have the (symbolic) ingredients for the tracebility analysis.
#xi_d,K,A

mvs.get_OutFluxesBySymbol().keys()

mvs.get_StateVariableTuple()

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
# **This step is hopefully not necessarry for your model since you should find parameters directly**
# Since we unfortunately started with a different description, we have to convert 
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
outflux_rates = {"k_"+str(key)+"_rh":value/key for key,value in hr.out_fluxes_by_symbol(sv,M_sym).items()}
internal_flux_rates = {"k_"+str(key[0])+"_2_"+str(key[1]):value/key[0] for key,value in hr.internal_fluxes_by_symbol(sv,M_sym).items()}
from copy import  deepcopy
all_rates=deepcopy(outflux_rates)
all_rates.update(internal_flux_rates)
all_rates

# +
# and one for the internal fluxrates
# -

old_par_dict = {
    beta_leaf: 0.6,  # 0 (parameters used in original code)
    beta_wood: 0.25,  # 1
    T_0: 2,  # 21
    E: 4,  # 22
    KM: 10,  # 23
    f_C_leaf_2_C_leaf_litter: 1,
    f_C_wood_2_C_wood_litter: 1, 
    f_C_root_2_C_root_litter: 1,
    f_C_leaf_litter_2_C_soil_fast: 0.41,  # 2
    f_C_leaf_litter_2_C_soil_slow: 0.07,  # 3
    f_C_leaf_litter_2_C_soil_passive: 0.02,  # 4
    f_C_wood_litter_2_C_soil_fast: 0.30,  # 5
    f_C_wood_litter_2_C_soil_slow: 0.12,  # 6
    f_C_wood_litter_2_C_soil_passive: 0.08,  # 7
    f_C_root_litter_2_C_soil_fast: 0.30,  # 8
    f_C_root_litter_2_C_soil_slow: 0.14,  # 9
    f_C_root_litter_2_C_soil_passive: 0.07,  # 10
    k_C_leaf: 1 / (60 * 2),  # 11
    k_C_wood: 1 / (365 * 30),  # 12
    k_C_root: 1 / (365 * 22),  # 13
    k_C_leaf_litter: 1 / (365 * 3.3),  # 14
    k_C_wood_litter: 1 / (365 * 11),  # 15
    k_C_root_litter: 1 / (365 * 11),  # 16
    k_C_soil_fast: 1 / (365 * 18),  # 17
    k_C_soil_slow: 1 / (365 * 100),  # 18
    k_C_soil_passive: 1 / (365 * 350),  # 19
}


# Now we can translate the old paramterisation to the new one.

par_dict={
    beta_leaf: 0.6,
    beta_wood: 0.25,
    T_0: 2, 
    E: 4, 
    KM: 10 
}
par_dict.update(
    {k:v.subs(old_par_dict) for k,v in all_rates.items()}
)
par_dict


# To be able to run the model forward we not only have to replace parameter symbols by values but symbolic functions by normal python functions.
# In our case the functions for $NPP$ and $\xi$ have to be provided. NPP_fun will interpolate the NPP for the day in question from the data. Which we have to load. 

# +
def NPP_fun(t):
    
func_dict={NPP: NPP_fun, ksi: ksi_fun}
# -

To make a 


