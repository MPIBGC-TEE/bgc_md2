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
# and copy its contents into a new cell. (we will later save it under a new name) 

# +
from sympy import  Symbol, Function 
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
    "I_vw": "Influx into vegetation wood pool",
    "r_vw_o": "out flux rate of wood pool",
    "r_vl_2_vw": "internal flux rate from leaf to wood", 
    "r_vw_2_vl": "internal flux rate from wood to leaf", 
}
# Make symbols from  the strings that we can later use in expressions  
# vl, vw,...
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    "I_vl": "Influx into vegetation leaf pool",
    "r_vl_o": "out flux rate of leaf pool",
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
        OutFluxesBySymbol({vl: r_vl_o(t) * vl, vw: r_vw_o * vw}),
        InternalFluxesBySymbol({(vl, vw): r_vl_2_vw * vl, (vw, vl): r_vw_2_vl * vw}),
    },

    computers=module_computers(bgc_c)
)
# -

# The last statement in the code defines a variable `mvs` which is 
# an instance of CMTVS which stands for `C`onnected`M`ulti`T`ype`V`ariable`S`et".
# It contains information in two forms. 
# 1. Variables of certain type (like InFluxesBySymbol)
# 2. a set of functions that "connect" these Variables and to other results we did not specify but which can be computed.
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
    "r_vl_2_vw": "internal flux rate from leaf to wood", 
    "r_vw_2_vl": "internal flux rate from wood to leaf", 
    'C_soil_fast': '',
    'C_soil_slow': '',
    'C_soil_passive': '',
    'C_leaf': '',
    'C_root': '',
    'C_wood': '',
    'C_leaf_litter': '',
    'C_root_litter': '',
    'C_wood_litter': '',
    'r_C_leaf_2_C_leaf_litter': '',
    'r_C_root_2_C_root_litter': '',
    'r_C_wood_2_C_wood_litter': '',
    'r_C_leaf_litter_rh': '',
    'r_C_root_litter_rh': '',
    'r_C_wood_litter_rh': '',
    'r_C_soil_fast_rh': '',
    'r_C_soil_slow_rh': '',
    'r_C_soil_passive_rh': '',
    'r_C_leaf_litter_2_C_soil_fast': '',
    'r_C_leaf_litter_2_C_soil_slow': '',
    'r_C_leaf_litter_2_C_soil_passive': '',
    'r_C_wood_litter_2_C_soil_fast': '',
    'r_C_wood_litter_2_C_soil_slow': '',
    'r_C_wood_litter_2_C_soil_passive': '',
    'r_C_root_litter_2_C_soil_fast': '',
    'r_C_root_litter_2_C_soil_slow': '',
    'r_C_root_litter_2_C_soil_passive': '',
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
    "r_vl_o": "out flux rate of leaf pool",
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
        OutFluxesBySymbol({vl: r_vl_o(t) * vl, vw: r_vw_o * vw}),
        InternalFluxesBySymbol({(vl, vw): r_vl_2_vw * vl, (vw, vl): r_vw_2_vl * vw}),
    },

    computers=module_computers(bgc_c)
)

# -

# Nothing has changed in the model description but be have some more symbols to work with.
# We can now use the new symbols in expression without errors since they are now defined.
# In the next step we add the new state variables.
# For the moment we leave the old ones (vl and vw) in to have a template for how to formulate fluxes)

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
    "r_vl_2_vw": "internal flux rate from leaf to wood", 
    "r_vw_2_vl": "internal flux rate from wood to leaf", 
    'C_soil_fast': '',
    'C_soil_slow': '',
    'C_soil_passive': '',
    'C_leaf': '',
    'C_root': '',
    'C_wood': '',
    'C_leaf_litter': '',
    'C_root_litter': '',
    'C_wood_litter': '',
    'r_C_leaf_2_C_leaf_litter': '',
    'r_C_root_2_C_root_litter': '',
    'r_C_wood_2_C_wood_litter': '',
    'r_C_leaf_litter_rh': '',
    'r_C_root_litter_rh': '',
    'r_C_wood_litter_rh': '',
    'r_C_soil_fast_rh': '',
    'r_C_soil_slow_rh': '',
    'r_C_soil_passive_rh': '',
    'r_C_leaf_litter_2_C_soil_fast': '',
    'r_C_leaf_litter_2_C_soil_slow': '',
    'r_C_leaf_litter_2_C_soil_passive': '',
    'r_C_wood_litter_2_C_soil_fast': '',
    'r_C_wood_litter_2_C_soil_slow': '',
    'r_C_wood_litter_2_C_soil_passive': '',
    'r_C_root_litter_2_C_soil_fast': '',
    'r_C_root_litter_2_C_soil_slow': '',
    'r_C_root_litter_2_C_soil_passive': '',
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
        OutFluxesBySymbol({vl: r_vl_o(t) * vl, vw: r_vw_o * vw}),
        InternalFluxesBySymbol({(vl, vw): r_vl_2_vw * vl, (vw, vl): r_vw_2_vl * vw}),
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
    # "vl": "vegegation leaf pool",
    # "vl": "vegegation leaf pool",
    # "vw": "vegetation wood pool content",
    # "I_vw": "Influx into vegetation wood pool",
    # "k_vw_o": "out flux rate of wood pool",
    # "k_vl_2_vw": "internal flux rate from leaf to wood", 
    # "k_vw_2_vl": "internal flux rate from wood to leaf", 
    "vl": "vegegation leaf pool",
    "vw": "vegetation wood pool content",
    "r_vl_2_vw": "internal flux rate from leaf to wood", 
    "r_vw_2_vl": "internal flux rate from wood to leaf", 
    'C_soil_fast': '',
    'C_soil_slow': '',
    'C_soil_passive': '',
    'C_leaf': '',
    'C_root': '',
    'C_wood': '',
    'C_leaf_litter': '',
    'C_root_litter': '',
    'C_wood_litter': '',
    'r_C_leaf_2_C_leaf_litter': '',
    'r_C_root_2_C_root_litter': '',
    'r_C_wood_2_C_wood_litter': '',
    'r_C_leaf_litter_rh': '',
    'r_C_root_litter_rh': '',
    'r_C_wood_litter_rh': '',
    'r_C_soil_fast_rh': '',
    'r_C_soil_slow_rh': '',
    'r_C_soil_passive_rh': '',
    'r_C_leaf_litter_2_C_soil_fast': '',
    'r_C_leaf_litter_2_C_soil_slow': '',
    'r_C_leaf_litter_2_C_soil_passive': '',
    'r_C_wood_litter_2_C_soil_fast': '',
    'r_C_wood_litter_2_C_soil_slow': '',
    'r_C_wood_litter_2_C_soil_passive': '',
    'r_C_root_litter_2_C_soil_fast': '',
    'r_C_root_litter_2_C_soil_slow': '',
    'r_C_root_litter_2_C_soil_passive': '',
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
                C_leaf_litter: r_C_leaf_litter_rh*C_leaf_litter*xi(t),
                C_wood_litter: r_C_wood_litter_rh*C_wood_litter*xi(t),
                C_root_litter: r_C_root_litter_rh*C_root_litter*xi(t),
                C_soil_fast: r_C_soil_fast_rh*C_soil_fast*xi(t),
                C_soil_slow: r_C_soil_slow_rh*C_soil_slow*xi(t),
                C_soil_passive: r_C_soil_passive_rh*C_soil_passive*xi(t),
            }
        ),
        InternalFluxesBySymbol(
            {
                #(vl, vw): k_vl_2_vw * vl, (vw, vl): k_vw_2_vl * vw
                (C_leaf, C_leaf_litter): r_C_leaf_2_C_leaf_litter*C_leaf, 
                (C_wood, C_wood_litter): r_C_wood_2_C_wood_litter*C_wood, 
                (C_root, C_root_litter): r_C_root_2_C_root_litter*C_root, 
                (C_leaf_litter, C_soil_fast)    : r_C_leaf_litter_2_C_soil_fast * C_leaf_litter*xi(t),
                (C_leaf_litter, C_soil_slow)    : r_C_leaf_litter_2_C_soil_slow * C_leaf_litter*xi(t),
                (C_leaf_litter, C_soil_passive) : r_C_leaf_litter_2_C_soil_passive * C_leaf_litter*xi(t),
                (C_wood_litter, C_soil_fast)    : r_C_wood_litter_2_C_soil_fast * C_wood_litter*xi(t),
                (C_wood_litter, C_soil_slow)    : r_C_wood_litter_2_C_soil_slow * C_wood_litter*xi(t),
                (C_wood_litter, C_soil_passive) : r_C_wood_litter_2_C_soil_passive * C_wood_litter*xi(t),
                (C_root_litter, C_soil_fast)    : r_C_root_litter_2_C_soil_fast * C_root_litter*xi(t),
                (C_root_litter, C_soil_slow)    : r_C_root_litter_2_C_soil_slow * C_root_litter*xi(t),
                (C_root_litter, C_soil_passive) : r_C_root_litter_2_C_soil_passive * C_root_litter*xi(t),
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

# We can store the symbolic representation in a file that we can use in for several purposes later.
# Copy your version of the following cell into a text editor and save the file in the current directory as `source.py`)

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
    "r_vl_2_vw": "internal flux rate from leaf to wood", 
    "r_vw_2_vl": "internal flux rate from wood to leaf", 
    'C_soil_fast': '',
    'C_soil_slow': '',
    'C_soil_passive': '',
    'C_leaf': '',
    'C_root': '',
    'C_wood': '',
    'C_leaf_litter': '',
    'C_root_litter': '',
    'C_wood_litter': '',
    'r_C_leaf_2_C_leaf_litter': '',
    'r_C_root_2_C_root_litter': '',
    'r_C_wood_2_C_wood_litter': '',
    'r_C_leaf_litter_rh': '',
    'r_C_root_litter_rh': '',
    'r_C_wood_litter_rh': '',
    'r_C_soil_fast_rh': '',
    'r_C_soil_slow_rh': '',
    'r_C_soil_passive_rh': '',
    'r_C_leaf_litter_2_C_soil_fast': '',
    'r_C_leaf_litter_2_C_soil_slow': '',
    'r_C_leaf_litter_2_C_soil_passive': '',
    'r_C_wood_litter_2_C_soil_fast': '',
    'r_C_wood_litter_2_C_soil_slow': '',
    'r_C_wood_litter_2_C_soil_passive': '',
    'r_C_root_litter_2_C_soil_fast': '',
    'r_C_root_litter_2_C_soil_slow': '',
    'r_C_root_litter_2_C_soil_passive': '',
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
                C_leaf: NPP(t) * beta_leaf, 
                C_root: NPP(t) * beta_root, 
                C_wood: NPP(t) * beta_wood
            }
        ),
        OutFluxesBySymbol(
            {
                C_leaf_litter: r_C_leaf_litter_rh*C_leaf_litter*xi(t),
                C_wood_litter: r_C_wood_litter_rh*C_wood_litter*xi(t),
                C_root_litter: r_C_root_litter_rh*C_root_litter*xi(t),
                C_soil_fast: r_C_soil_fast_rh*C_soil_fast*xi(t),
                C_soil_slow: r_C_soil_slow_rh*C_soil_slow*xi(t),
                C_soil_passive: r_C_soil_passive_rh*C_soil_passive*xi(t),
            }
        ),
        InternalFluxesBySymbol(
            {
                (C_leaf, C_leaf_litter): r_C_leaf_2_C_leaf_litter*C_leaf, 
                (C_wood, C_wood_litter): r_C_wood_2_C_wood_litter*C_wood, 
                (C_root, C_root_litter): r_C_root_2_C_root_litter*C_root, 
                (C_leaf_litter, C_soil_fast)    : r_C_leaf_litter_2_C_soil_fast * C_leaf_litter*xi(t),
                (C_leaf_litter, C_soil_slow)    : r_C_leaf_litter_2_C_soil_slow * C_leaf_litter*xi(t),
                (C_leaf_litter, C_soil_passive) : r_C_leaf_litter_2_C_soil_passive * C_leaf_litter*xi(t),
                (C_wood_litter, C_soil_fast)    : r_C_wood_litter_2_C_soil_fast * C_wood_litter*xi(t),
                (C_wood_litter, C_soil_slow)    : r_C_wood_litter_2_C_soil_slow * C_wood_litter*xi(t),
                (C_wood_litter, C_soil_passive) : r_C_wood_litter_2_C_soil_passive * C_wood_litter*xi(t),
                (C_root_litter, C_soil_fast)    : r_C_root_litter_2_C_soil_fast * C_root_litter*xi(t),
                (C_root_litter, C_soil_slow)    : r_C_root_litter_2_C_soil_slow * C_root_litter*xi(t),
                (C_root_litter, C_soil_passive) : r_C_root_litter_2_C_soil_passive * C_root_litter*xi(t),
            }
        ),
    },


    computers=module_computers(bgc_c)
)
# -

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
Observables = namedtuple(
    'Observables',
    ["cLitter", "cSoil", "cVeg", "gpp", "rh" ,"ra" ]
)
    
    
#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output():
    download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['VISIT'],
        variables = Observables._fields
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
    names=[(f,"VISIT_S2_{}.nc".format(f)) for f in Observables._fields]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them 
        ds = nc.Dataset(str(path))
        return ds.variables[vn][t]

    return Observables(*map(f, names))


#     # Read NetCDF data  ******************************************************************************************************************************
# -

svs=get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))


svs

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
    beta_leaf: 0.6,
    beta_wood: 0.25,
    T_0: 2,
    E: 4,
    KM: 10,
    f_C_leaf_2_C_leaf_litter: 1,
    f_C_wood_2_C_wood_litter: 1,
    f_C_root_2_C_root_litter: 1,
    f_C_leaf_litter_2_C_soil_fast: 0.41,
    f_C_leaf_litter_2_C_soil_slow: 0.07,
    f_C_leaf_litter_2_C_soil_passive: 0.02,
    f_C_wood_litter_2_C_soil_fast: 0.30,
    f_C_wood_litter_2_C_soil_slow: 0.12,
    f_C_wood_litter_2_C_soil_passive: 0.08,
    f_C_root_litter_2_C_soil_fast: 0.30,
    f_C_root_litter_2_C_soil_slow: 0.14,
    f_C_root_litter_2_C_soil_passive: 0.07,
    k_C_leaf: 1 / (60 * 2),
    k_C_wood: 1 / (365 * 30),
    k_C_root: 1 / (365 * 22),
    k_C_leaf_litter: 1 / (365 * 3.3),
    k_C_wood_litter: 1 / (365 * 11),
    k_C_root_litter: 1 / (365 * 11),
    k_C_soil_fast: 1 / (365 * 18),
    k_C_soil_slow: 1 / (365 * 100),
    k_C_soil_passive: 1 / (365 * 350),
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
    {Symbol(k):v.subs(old_par_dict) for k,v in all_rates.items()}
)
par_dict

# To be able to run the model forward we not only have to replace parameter symbols by values but symbolic functions by normal python functions.
# In our case the functions for $NPP$ and $\xi$ have to be provided. NPP_fun will interpolate the NPP for the day in question from the data. Which we have to load. 
# We will later store these functions in  `model_specific_helpers.py` which resides in the same folder as this notebook. You will have to adapt them to your data set. 

# +
from general_helpers import make_B_u_funcs_2, day_2_month_index
# check the numeric functions for B and u

# test the dayly iterator
def npp_func(day):
    month=day_2_month_index(day)
    return svs.gpp[month]-(svs.rh[month]+svs.ra[month])

def xi_func(day):
    return 1.0 # preliminary fake for lack of better data... 

func_dict={
    'NPP':npp_func,
     'xi':xi_func
}

B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
X_0= np.array((
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
)).reshape(9,1)
# in general B and u are nonlinear and can depend on X, thats why we have to test them with index and X arguments
u_func(0,X_0),B_func(0,X_0)

# +
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
        TimeStepIterator2,
)
from copy import copy

def make_daily_iterator_sym(
        mvs,
        V_init,
        par_dict,
        func_dict
    ):
    B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
    def f(it,V):
        X = V[0:9]
        co2 = V[9]
        b = u_func(it,X)
        B = B_func(it,X)
        X_new = X + b + B@X

        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
        out = B @ X
        ra= -np.sum(out[0:3]) # check with  mvs.get_StateVariableTuple()[0:3]
        rh= -np.sum(out[3:9]) # check with  mvs.get_StateVariableTuple()[3:9]
        
        V_new = np.concatenate((X_new,np.array([ra,rh]).reshape(2,1)), axis=0)
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_init,
        f=f#,
        #max_it=max(day_indices)+1
    )


# -

mvs.get_StateVariableTuple()[3:9]

Observables._fields

# +
# test the dayly iterator

svs_0=Observables(*map(lambda v: v[0],svs))
V_init=np.array((
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.ra,
    svs_0.rh
)).reshape(11,1)
    
it_sym = make_daily_iterator_sym(
    mvs,
    V_init=V_init,
    par_dict=par_dict,
    func_dict=func_dict
)
n=5
res= np.zeros((n,len(V_init)))
res_sym = copy(res)
for i in range(n):
    res_sym[i,:]=it_sym.__next__().reshape(len(V_init),)
res_sym
# -




