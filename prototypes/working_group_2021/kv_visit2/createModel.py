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
# - Note that for the following steps you **do not** have to copy and paste the previous cell
#   as I did to show the incremental progress.

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
    'tas': '',
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
    'tas': '',
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
    'tas': '',
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
# we now have the (symbolic) ingredients for the traceability analysis.
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
    'tas': '',
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
    ["cVeg","cLitter","cSoil","rh","ra"]
)
Drivers=namedtuple(
    "Drivers",
    ["gpp", "mrso", "tas"]
)    
    
#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output():
    download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['VISIT'],
        variables = Observables._fields + Drivers._fields
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
    t = s, 50, 33  # [t] = [:,49,325]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them 
        ds = nc.Dataset(str(path))
        return ds.variables[vn][t]

    o_names=[(f,"VISIT_S2_{}.nc".format(f)) for f in Observables._fields]
    d_names=[(f,"VISIT_S2_{}.nc".format(f)) for f in Drivers._fields]
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
#
# ### Providing dictionaries  for the parameters and functions.
#

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

def npp_func(day):
    month=day_2_month_index(day)
    return (dvs.gpp[month]-svs.ra[month]) * 86400   # kg/m2/s kg/m2/day;

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
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cLitter/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
    svs_0.cSoil/3,
))#.reshape(9,)
# in general B and u are nonlinear and can depend on X, thats why we have to test them with index and X arguments
u_func(0,X_0),B_func(0,X_0)

# +
import CompartmentalSystems.helpers_reservoir as hr

# To compute the ra and rh we have to some up the values for autotrophic and heterotrophic respiration we have to sum up outfluxes.
# We first create numerical functions for all the outfluxes from our symbolic description.
# lets do it for the outflux of one pool first (and then for the whole lot of them)
expr_cont=mvs.get_OutFluxesBySymbol()[C_leaf_litter]

# this is the momentary outflux as a function of t,C_leaf_litter as it would occure in the differential equation
expr_cont

# -

# what we want is acutally the flux in the simplest approximation 
delta_t=Symbol("delta_t") #
it=Symbol("it")           #arbitrary symbol for the step index )
t=mvs.get_TimeSymbol()
expr_disc=expr_cont.subs({t:delta_t*it})
# If you look at the result you see that the euler forwar approximation assumes the flux at time t to be constant for the timestep  
expr_disc

# Remark
# if we assume tat delta_t is 1 day and it counts days 
# it becomes even simpler
expr_disc.subs({delta_t:1})
# Which is the same as if we had 't' replaced in the above formula wiht 'it'
# but this does not generalise to a bigger stepsize

# +
# this expression we turn now into a numeric funciton of it
# although in our example it depends only on the iteration "it" and C_leaf_litter 
# we make it a function of it and ALL state variables to be able to call it in the same way as u_func and B_func 
argtup=(it, *mvs.get_StateVariableTuple())

C_leaf_litter_func = hr.numerical_function_from_expression(
    expr=expr_disc.subs({delta_t:1}),
    tup=argtup, 
    parameter_dict=par_dict,
    func_set=func_dict
)
# call it for testing
C_leaf_litter_func(0,*X_0)


# +
#mass production of output functions
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
    
numOutFluxesBySymbol={
    k:numfunc(expr_cont,1) 
    for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
} 
numOutFluxesBySymbol2={
    k:numfunc(expr_cont,2) 
    for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
} 

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
        for k in [C_leaf_litter,C_wood_litter,C_root_litter,C_soil_fast,C_soil_slow,C_soil_passive] 
        if k in numOutFluxesBySymbol.keys()
    ]
)
ra_0,rh_0
# -

numOutFluxesBySymbol[C_leaf_litter](0,*X_0)

numOutFluxesBySymbol2[C_leaf_litter](0,*X_0)

# We now build the essential object to run the model forward. We have a 
# - startvector $V_0$ and 
# - a function $f$ to compute the next value: $V_{it+1} =f(it,V_{it})$
#   the dependence on the iteration $it$ allows us to represent drivers that
#   are functions of time 
#
# So we can compute all values:
#
# $V_1=f(0,V_0)$
#
# $V_2=f(1,V_1)$
#
# ...
#
# $V_n+1=f(n,V_n)$
#
# Technically this can be implemented as an `iterator` object with a `__next__()` method to move our system one step forward in time. 
#
# What we want to build is an object `it_sym` that can use as follows.
# ```python
# for i in range(10):
#     print(it_sym.__next__())
# ```
# to print out the first 10 values.
#
# If iterators had not been invented yet we would invent them now, because they capture exactly the mathematical concept of an initial value system, 
# without all the nonessential technical details of e.g. how many steps we want to make or which of the results we want to store.
# This is essential if we want to test and use the iterator later in different scenarios but avoid reimplemtation of the basic timestep. 
#
# Remark:
#
# If we were only interested in the timeseries of the pool contents `bgc_md2` could compute the solution automatically without the need to build an iterator ourselves. 
# In our case we are also interested in tracking the autotrophic and heterotrophic respiration and therefore have to recreate and extend the code `bgc_md2` would use in the background.
# We will let `bgc_md2` do part of the work and derive numeric functions for the Compartmental matrix $B$ and the input $u$ and the Outfluxes - from which to compute $ra$ $rh$ - from our symbolic description but we will build our own iterator by combining these functions.    
# We will proceed as follows:
# - create $V_0$ 
# - build the function $f$
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
    C_leaf=svs_0.cVeg/3,
    C_wood=svs_0.cVeg/3,
    C_root=svs_0.cVeg/3,
    C_leaf_litter=svs_0.cLitter/3,
    C_wood_litter=svs_0.cLitter/3,
    C_root_litter=svs_0.cLitter/3,
    C_soil_fast=svs_0.cSoil/3,
    C_soil_slow=svs_0.cSoil/3,
    C_soil_passive=svs_0.cSoil/3,
    ra=svs_0.ra*86400,   # kg/m2/s kg/m2/day;,
    rh=svs_0.rh*86400   # kg/m2/s kg/m2/day;        
)
V_init.__getattribute__("C_leaf")

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
    B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict,delta_t_val=1)  
    sv=mvs.get_StateVariableTuple()
    n=len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct 
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr=np.array(
        [V_init.__getattribute__(str(v)) for v in sv]+
        [V_init.ra,V_init.rh]
    ).reshape(n+2,1) #reshaping is neccessary for matmul
    
    numOutFluxesBySymbol={
        k:numfunc(expr_cont,delta_t_val=1) 
        for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
    } 
    def f(it,V):
        X = V[0:n]
        b = u_func(it,X)
        B = B_func(it,X)
        X_new = X + b + B @ X
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
        
        ra = np.sum(
            [
              numOutFluxesBySymbol[k](it,*X) 
                for k in [C_leaf,C_wood,C_root] 
                if k in numOutFluxesBySymbol.keys()
            ]
        )
        rh = np.sum(
            [
                numOutFluxesBySymbol[k](it,*X) 
                for k in [C_leaf_litter,C_wood_litter,C_root_litter,C_soil_fast,C_soil_slow,C_soil_passive] 
                if k in numOutFluxesBySymbol.keys()
            ]
        )
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

# We now create another iterator that does run with a different timestep and check that this does not change the results.
# If we can increase the timestep the data assimilation will be much faster.
#

# +

from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
        TimeStepIterator2,
)
from copy import copy

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
        [V_init.ra,V_init.rh]
    ).reshape(n+2,1) #reshaping is neccessary for matmul (the @ in B @ X)
    
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
        
        ra = np.sum(
            [
              numOutFluxesBySymbol[k](it,*X)
                for k in [C_leaf,C_wood,C_root] 
                if k in numOutFluxesBySymbol.keys()
            ]
        )
        rh = np.sum(
            [
                numOutFluxesBySymbol[k](it,*X)
                for k in [C_leaf_litter,C_wood_litter,C_root_litter,C_soil_fast,C_soil_slow,C_soil_passive] 
                if k in numOutFluxesBySymbol.keys()
            ]
        )
        V_new = np.concatenate((X_new.reshape(n,1),np.array([ra,rh]).reshape(2,1)), axis=0)
        
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


# +

# test the different iterators
    
it_sym = make_daily_iterator_sym(
    mvs,
    V_init=V_init,
    par_dict=par_dict,
    func_dict=func_dict
)
delta_t_val=30 #n_day iterator
it_sym_2 = make_iterator_sym(
    mvs,
    V_init=V_init,
    par_dict=par_dict,
    func_dict=func_dict,
    delta_t_val=delta_t_val
)
# we will run the dayly for 15 steps
ns=10000
res= np.zeros((ns,len(V_init)))
times = np.arange(0,ns)
res[0,:]=V_init 
for i in range(1,ns-1):
    res[i,:]=it_sym.__next__().reshape(len(V_init),)

times_2= np.arange(0,ns,delta_t_val)
res_2= np.zeros((len(times_2),len(V_init)))
res_2[0,:]=V_init 
for i in range(1,len(times_2)-1):
    res_2[i,:]=it_sym_2.__next__().reshape(len(V_init),)
# -

times


times_2

from matplotlib import pyplot as plt
fig = plt.figure(figsize=(5,50))
n_coord=len(V_init)
axs=fig.subplots(n_coord,1)
for coord in range(n_coord):
    axs[coord].plot(times,res[:,coord])
    axs[coord].plot(times_2,res_2[:,coord])

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
        "cLitter_0",
        "cSoil_0",
        "cVeg_0",
        "gpp_0",
        "rh_0",
        "ra_0",
        "r_C_root_litter_2_C_soil_passive",# here  we pretend to know these two rates 
        "r_C_root_litter_2_C_soil_slow",# it would be much better to know more  
        "number_of_months" # necessary to prepare the output in the correct lenght 
    ]
)
# Note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues and their length (number of months)
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 
# It is better to start with only a few

EstimatedParameters = namedtuple(
    "EstimatedParameters",[ 
        "beta_leaf",
        "beta_wood",
        "T_0",
        "E",
        "KM",
        #"r_C_leaf_rh",
        #"r_C_wood_rh",
        #"r_C_root_rh",
        "r_C_leaf_litter_rh",
        "r_C_wood_litter_rh",
        "r_C_root_litter_rh",
        "r_C_soil_fast_rh",
        "r_C_soil_slow_rh",
        "r_C_soil_passive_rh",
        "r_C_leaf_2_C_leaf_litter",
        "r_C_wood_2_C_wood_litter",
        "r_C_root_2_C_root_litter",
        "r_C_leaf_litter_2_C_soil_fast",
        "r_C_leaf_litter_2_C_soil_slow",
        "r_C_leaf_litter_2_C_soil_passive",
        "r_C_wood_litter_2_C_soil_fast",
        "r_C_wood_litter_2_C_soil_slow",
        "r_C_wood_litter_2_C_soil_passive",
        "r_C_root_litter_2_C_soil_fast",
        "r_C_root_litter_2_C_soil_slow",
        "r_C_root_litter_2_C_soil_passive",
        'C_leaf_0',#for the trendy data also the startvalues have to be estimated but 
        'C_wood_0',
        #C_root_0 can be inferred as cVeg_0-(C_leaf_0+C_wood_0)
        'C_leaf_litter_0',
        'C_wood_litter_0',
        #C_root_litter_0 can be inferred
        'C_soil_fast_0',
        'C_soil_slow_0',
        #C_soil_passive_0 can be inferred 
    ]
)
# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues. 
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 


# -

EstimatedParameters._fields

dvs.gpp[1:100]

cpa=UnEstimatedParameters(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 gpp_0=dvs.gpp[0] * 86400,   # kg/m2/s kg/m2/day
 rh_0=svs_0.rh * 86400,   # kg/m2/s kg/m2/day
 ra_0=svs_0.ra * 86400,   # kg/m2/s kg/m2/day
 r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
 r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
 #number_of_months=len(svs.rh)
 number_of_months=120 # for testing and tuning mcmc
)

len(svs.rh)


# +
def make_steady_state_iterator_sym(
        mvs,
        V_init,
        par_dict,
        func_dict
    ):
    B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
    def f(it,tup):
        X,_,_=tup
        b = u_func(it,X)
        B = B_func(it,X)
        return (X,b,B)
  
    return TimeStepIterator2(
        initial_values=V_init,
        f=f)
# calculate steady state
  
it_sym = make_steady_state_iterator_sym(
    mvs,
    V_init=(X_0,u_func(0,X_0),B_func(0,X_0)),
    par_dict=par_dict,
    func_dict=func_dict
)
Bs=[]
bs=[]
for i in range(cpa.number_of_months*30):
    bs.append(it_sym.__next__()[1])
    Bs.append(it_sym.__next__()[2])
B_mean=np.stack(Bs).mean(axis=0)
b_mean=np.stack(bs).mean(axis=0)
B_mean,b_mean
np.linalg.inv(Bs[0])


# +
# calculate pseudo steady state
X_ss = np.linalg.solve(B_mean, (-b_mean))

steady_state_dict={str(name): X_ss[i,0] for i,name in enumerate(mvs.get_StateVariableTuple())}
# -

# create a start parameter tuple for the mcmc. The order has to be the same as when you created the namedtupl3 
# If you don't you get a "TypeError". 
epa_0=EstimatedParameters(
 beta_leaf=0.6,
 beta_wood=0.25,
 T_0=2,
 E=4,
 KM=10,
 r_C_leaf_litter_rh=0.000415110004151100,
 r_C_wood_litter_rh=0.000124533001245330,
 r_C_root_litter_rh=0.000122042341220423,
 r_C_soil_fast_rh=0.000152207001522070,
 r_C_soil_slow_rh=2.73972602739726e-5,
 r_C_soil_passive_rh=7.82778864970646e-6,
 r_C_leaf_2_C_leaf_litter=0.00833333333333333,
 r_C_wood_2_C_wood_litter=9.13242009132420e-5,
 r_C_root_2_C_root_litter=0.000124533001245330,
 r_C_leaf_litter_2_C_soil_fast=0.000340390203403902,
 r_C_leaf_litter_2_C_soil_slow=5.81154005811540e-5,
 r_C_leaf_litter_2_C_soil_passive=1.66044001660440e-5,
 r_C_wood_litter_2_C_soil_fast=7.47198007471980e-5,
 r_C_wood_litter_2_C_soil_slow=2.98879202988792e-5,
 r_C_wood_litter_2_C_soil_passive=1.99252801992528e-5,
 r_C_root_litter_2_C_soil_fast=7.47198007471980e-5,
 r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
 r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
 C_leaf_0=steady_state_dict["C_leaf"],
 C_wood_0=steady_state_dict["C_wood"],
 C_leaf_litter_0=steady_state_dict["C_leaf_litter"],
 C_wood_litter_0=steady_state_dict["C_wood_litter"],
 C_soil_fast_0=steady_state_dict["C_soil_fast"],
 C_soil_slow_0=steady_state_dict["C_soil_slow"]
)

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
from general_helpers import month_2_day_index, monthly_to_yearly
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
    seconds_per_day = 86400
    def npp_func(day):
        month=day_2_month_index(day)
        return (dvs.gpp[month]-svs.ra[month]) * seconds_per_day   # kg/m2/s kg/m2/day;
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
        V_init = StartVector(
            C_leaf=epa.C_leaf_0,
            C_wood=epa.C_wood_0,
            C_root=cpa.cVeg_0-(epa.C_leaf_0 + epa.C_wood_0),
            C_leaf_litter=epa.C_leaf_litter_0,
            C_wood_litter=epa.C_wood_litter_0,
            C_root_litter=cpa.cLitter_0-(epa.C_leaf_litter_0 + epa.C_wood_litter_0),
            C_soil_fast=epa.C_soil_fast_0,
            C_soil_slow=epa.C_soil_slow_0,
            C_soil_passive=cpa.cSoil_0-(epa.C_soil_fast_0 + epa.C_soil_slow_0),
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
        
        #mass production of output functions
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
        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2) 
        delta_t_val=15 
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
        )
        
        # Now that we have the iterator we can start to compute.
        # the first thing to notice is that we don't need to store all values (daili yi delta_t_val=1)
        # since the observations are recorded monthly while our iterator possibly has a smaller timestep.
        # - for the pool contents we only need to record the last day of the month
        # - for the respiration(s) ra and rh we want an accumulated value (unit mass) 
        #   have to sum up the products of the values*delta_t over a month
        # 
        # Note: check if TRENDY months are like this...
        # days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        sols=[]
        for m in range(cpa.number_of_months):
            #dpm = days_per_month[ m % 12]  
            dpm=30 # 
            mra=0
            mrh=0
            for d in range(int(dpm/delta_t_val)):
                v = it_sym.__next__()
                #mra +=v[9,0] 
                #mrh +=v[10,0]
                # actually the original datastream seems to be a flux per area (kg m^-2 ^-s )
                # at the moment the iterator also computes a flux but in kg^-2 ^day
            V=StartVector(*v)
            o=Observables(
                cVeg=float(V.C_leaf+V.C_wood+V.C_root),
                cLitter=float(V.C_leaf_litter+V.C_wood_litter+V.C_root_litter),
                cSoil=float(V.C_soil_fast+V.C_soil_slow+V.C_soil_passive),
                ra=v[9,0]/seconds_per_day,
                rh=v[10,0]/seconds_per_day # the data is per second while the time units for the iterator refer to days
                #ra=mra/dpm, # monthly respiration back to kg/m2/day units
                #rh=mrh/dpm, # monthly respiration back to kg/m2/day units
            )
            sols.append(o)
            
        sol=np.stack(sols)       
        #convert to yearly output if necessary (the monthly pool data looks very "yearly")
        #sol_yr=np.zeros(int(cpa.number_of_months/12)*sol.shape[1]).reshape([int(cpa.number_of_months/12),sol.shape[1]])  
        #for i in range(sol.shape[1]):
        #   sol_yr[:,i]=monthly_to_yearly(sol[:,i])
        #sol=sol_yr
        return sol 
        
    return param2res

# +
# now test it 
import matplotlib.pyplot as plt
from general_helpers import plot_solutions

param2res_sym = make_param2res_sym(cpa)
xs= param2res_sym(epa_0)
obs=np.column_stack([ np.array(v) for v in svs])
xs.shape

# +
obs=obs[0:cpa.number_of_months,:] #cut 
obs[:,3:4]=obs[:,3:4]
n=cpa.number_of_months

# convert to yearly output if necessary
#obs_yr=np.zeros(int(cpa.number_of_months/12)*obs.shape[1]).reshape([int(cpa.number_of_months/12),obs.shape[1]])  
#for i in range(obs.shape[1]):
#    obs_yr[:,i]=monthly_to_yearly(obs[:,i])
#obs=obs_yr
#n=int(cpa.number_of_months/12)

print("Forward run with initial parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(12, 4), dpi=80)
plot_solutions(
        fig,
        times=range(n),
        var_names=Observables._fields,
        tup=(xs,obs)
        #tup=(obs,)
)
fig.savefig('solutions.pdf')

# +
epa_min=EstimatedParameters(
 beta_leaf=0,
 beta_wood=0,
 T_0=-20,
 E=.1,
 KM=1,
 r_C_leaf_litter_rh=epa_0.r_C_leaf_litter_rh/100,
 r_C_wood_litter_rh=epa_0.r_C_wood_litter_rh/100,
 r_C_root_litter_rh=epa_0.r_C_root_litter_rh/100,
 r_C_soil_fast_rh=epa_0.r_C_soil_fast_rh/100,
 r_C_soil_slow_rh=epa_0.r_C_soil_slow_rh/100,
 r_C_soil_passive_rh=epa_0.r_C_soil_passive_rh/100,
 r_C_leaf_2_C_leaf_litter=epa_0.r_C_leaf_2_C_leaf_litter/100,       
 r_C_wood_2_C_wood_litter=epa_0.r_C_wood_2_C_wood_litter/100,
 r_C_root_2_C_root_litter=epa_0.r_C_root_2_C_root_litter/100,
 r_C_leaf_litter_2_C_soil_fast=epa_0.r_C_leaf_litter_2_C_soil_fast/100,
 r_C_leaf_litter_2_C_soil_slow=epa_0.r_C_leaf_litter_2_C_soil_slow/100,
 r_C_leaf_litter_2_C_soil_passive=epa_0.r_C_leaf_litter_2_C_soil_passive/100,
 r_C_wood_litter_2_C_soil_fast=epa_0.r_C_wood_litter_2_C_soil_fast/100,
 r_C_wood_litter_2_C_soil_slow=epa_0.r_C_wood_litter_2_C_soil_slow/100,
 r_C_wood_litter_2_C_soil_passive=epa_0.r_C_wood_litter_2_C_soil_passive/100,
 r_C_root_litter_2_C_soil_fast=epa_0.r_C_root_litter_2_C_soil_fast/100,
 r_C_root_litter_2_C_soil_slow=epa_0.r_C_root_litter_2_C_soil_slow/100,
 r_C_root_litter_2_C_soil_passive=epa_0.r_C_root_litter_2_C_soil_passive/100,
 C_leaf_0=0,
 C_wood_0=0,
 C_leaf_litter_0=0,
 C_wood_litter_0=0,
 C_soil_fast_0=0,
 C_soil_slow_0=0,
)


epa_max=EstimatedParameters(
 beta_leaf=0.99,
 beta_wood=0.99,
 T_0=10,
 E=100,
 KM=100,
 r_C_leaf_litter_rh=epa_0.r_C_leaf_litter_rh*100,
 r_C_wood_litter_rh=epa_0.r_C_wood_litter_rh*100,
 r_C_root_litter_rh=epa_0.r_C_root_litter_rh*100,
 r_C_soil_fast_rh=epa_0.r_C_soil_fast_rh*100,
 r_C_soil_slow_rh=epa_0.r_C_soil_slow_rh*100,
 r_C_soil_passive_rh=epa_0.r_C_soil_passive_rh*100,
 r_C_leaf_2_C_leaf_litter=epa_0.r_C_leaf_2_C_leaf_litter*100,       
 r_C_wood_2_C_wood_litter=epa_0.r_C_wood_2_C_wood_litter*100,
 r_C_root_2_C_root_litter=epa_0.r_C_root_2_C_root_litter*100,
 r_C_leaf_litter_2_C_soil_fast=epa_0.r_C_leaf_litter_2_C_soil_fast*100,
 r_C_leaf_litter_2_C_soil_slow=epa_0.r_C_leaf_litter_2_C_soil_slow*100,
 r_C_leaf_litter_2_C_soil_passive=epa_0.r_C_leaf_litter_2_C_soil_passive*100,
 r_C_wood_litter_2_C_soil_fast=epa_0.r_C_wood_litter_2_C_soil_fast*100,
 r_C_wood_litter_2_C_soil_slow=epa_0.r_C_wood_litter_2_C_soil_slow*100,
 r_C_wood_litter_2_C_soil_passive=epa_0.r_C_wood_litter_2_C_soil_passive*100,
 r_C_root_litter_2_C_soil_fast=epa_0.r_C_root_litter_2_C_soil_fast*100,
 r_C_root_litter_2_C_soil_slow=epa_0.r_C_root_litter_2_C_soil_slow*100,
 r_C_root_litter_2_C_soil_passive=epa_0.r_C_root_litter_2_C_soil_passive*100,
 C_leaf_0=svs_0.cVeg,
 C_wood_0=svs_0.cVeg,
 C_leaf_litter_0=svs_0.cLitter,
 C_wood_litter_0=svs_0.cLitter,
 C_soil_fast_0=svs_0.cSoil,
 C_soil_slow_0=svs_0.cSoil,
)
# -

# ### mcmc to optimize parameters 
#

np.array(epa_max)

# +
from general_helpers import autostep_mcmc, make_param_filter_func, make_feng_cost_func

isQualified = make_param_filter_func(epa_max, epa_min)
param2res = make_param2res_sym(cpa)

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa_0,
    filter_func=isQualified,
    param2res=param2res,
    costfunction=make_feng_cost_func(obs),
    nsimu=2000, # for testing and tuning mcmc
    #nsimu=20000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=15,   # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=1,   # default value | increase value to reduce initial step size
    K=2 # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

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
        times=range(int(cpa.number_of_months)), # for yearly output
        var_names=Observables._fields,
        tup=(mod_opt,obs)
)

fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

import pandas as pd
pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('visit_da_aa.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('visit_da_j_aa.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('visit_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('visit_optimized_solutions.csv'), sep=',')
# -

print("Optimized parameters: ", epa_opt)
par_dict_opt={
    beta_leaf: epa_opt.beta_leaf,
    beta_wood: epa_opt.beta_wood,
    T_0: epa_opt.T_0,
    E: epa_opt.E,
    KM: epa_opt.KM,
    #r_C_leaf_rh: 0,
    #r_C_wood_rh: 0,
    #r_C_root_rh: 0,
    r_C_leaf_litter_rh: epa_opt.r_C_leaf_litter_rh,
    r_C_wood_litter_rh: epa_opt.r_C_wood_litter_rh,
    r_C_root_litter_rh: epa_opt.r_C_root_litter_rh,
    r_C_soil_fast_rh: epa_opt.r_C_soil_fast_rh,
    r_C_soil_slow_rh: epa_opt.r_C_soil_slow_rh,
    r_C_soil_passive_rh: epa_opt.r_C_soil_passive_rh,
    r_C_leaf_2_C_leaf_litter: epa_opt.r_C_leaf_2_C_leaf_litter,
    r_C_wood_2_C_wood_litter: epa_opt.r_C_wood_2_C_wood_litter,
    r_C_root_2_C_root_litter: epa_opt.r_C_root_2_C_root_litter,
    r_C_leaf_litter_2_C_soil_fast: epa_opt.r_C_leaf_litter_2_C_soil_fast,
    r_C_leaf_litter_2_C_soil_slow: epa_opt.r_C_leaf_litter_2_C_soil_slow,
    r_C_leaf_litter_2_C_soil_passive: epa_opt.r_C_leaf_litter_2_C_soil_passive,
    r_C_wood_litter_2_C_soil_fast: epa_opt.r_C_wood_litter_2_C_soil_fast,
    r_C_wood_litter_2_C_soil_slow: epa_opt.r_C_wood_litter_2_C_soil_slow,
    r_C_wood_litter_2_C_soil_passive: epa_opt.r_C_wood_litter_2_C_soil_passive,
    r_C_root_litter_2_C_soil_fast: epa_opt.r_C_root_litter_2_C_soil_fast,
    r_C_root_litter_2_C_soil_slow: epa_opt.r_C_root_litter_2_C_soil_slow,
    r_C_root_litter_2_C_soil_passive: epa_opt.r_C_root_litter_2_C_soil_passive 
}
print("Optimized parameters dictionary: ", par_dict_opt)

# ### Traceability analysis  
#
# #### Outline
# The traceability analysis defines several diagnostic variables using as much algebraic structure of the mass balance equation as is available.
# Not all diagnostic variables are possible for all compartmental models. 
#
# We chose here to introduce the diagnostic variables not all at once but rather in the order of decreasing generality.
#
# The first diagnostic variables are available for all compartmental models and need no additional assumptions. 
# In the later parts of this section we then assume to be able to identify more and more specific terms in the mass balance equation and use those to derive and trace ever more specific diagnostics.
# Thus the very first part is valid for all models but how many of the later parts are applicable to a specific model  depends on how much we know about it.  
#
#
# #### Derivation of the matrix decomposition 
# Compartmental models (well mixed mass balanced) can be written in as an ordinary differential equation in matrix form that relates the momentary value of the (time) derivative $\frac{d X}{d t}$ of an yet unknown function $X$ to the momentary value of $X$ itself.   
# $$
# \frac{d X}{d t}= M(X,t) X + I(X,t) \quad (1)   
# $$ 
# where $X$ is the statevector representing the pool contents, $M$ the "Compartmental matrix" and $I$ the input vector.
# Together with a startvalue $X_0$ it constitutes an "initial value problem" (ivp) which can be solved numerically by moving step by step forward in time.
#
# Note: 
#
# It is mathematical standard notation to use $X$ in the *formulation* of the ivp (representing the momentary value) althoug *after we have solved it* the solution is expressed as function of time $X(t)$. This avoids confusion since everything appering with arguments is recognizable as explicitly calculable *before* we have solved the ivp.
#
# The system "nonautonomous" (if they depend on time $t$) and "nonlinear" if the dependent on $X$.
# It is always possible to factorize $M(X,t)$ into a product $M=A(X,t) K(X,t)$ where $K$ is a  diagonal matrix.
# and $I=B(t)*u(t)$ where $u$ is a scalar.
# Using these we arrive at 
# $$
# \frac{d X}{d t}= A(X,t) K(X,t) X + B(X,t) u(X,t)  
# $$
# ##### Linearity assumption
# If we assume the model to be linear and nonautonomous the dependency on $X$ vanishes and we have
#
# $$
# \frac{d X}{d t}= A(t) K(t) X + B(t) u(t) . 
# $$
#
# ##### Factorizability  assumption
# Although this is not possible in general in many published models the nonautonous part  can be further localized into a diagonal matrix $\xi(t)$ so that we can achieve constant $A$ and $K$ which allows more specific interpretation.
#
# $$
# \frac{d X}{d t}= A \xi(t) K X + B(t)u(t)
# $$
#
# ##### Factorizability of $\xi$ assumption 
# In some cases we can resolve $\xi$ further.
# $$
# \frac{d X}{d t}= A \xi_temp(t) \xi_mois(t) K X + B(t)u(t)
# $$
#
# #### Definition of diagnostic variables
#
# ##### Storage capacity $X_c$ and storage potential $X_p$
# These variables can be defined for any compartmental system and do not require either linearity nor factorizability. 
# We can rearrange eq. $(1)$ and give names to the two summands. 
# $$
# X = M^{-1}(X,t) \left( \frac{d X}{d t}-I(X,t) \right) \\ 
#   = \underbrace{M^{-1}(X,t) \frac{d X}{d t}}_{X_c} - \underbrace{M^{-1}(X,t)I(X,t)}_{X_p} \\
#   = X_c - X_p
# $$
# Note:
# This is not to be read as a recipe to compute $X$.
# The equation becomes a bit clearer if we adapt the nomenclature to express that we *have solved the ivp* and know its solution $X(t)$  
# <!-- and therefore also  the derivative $\frac{d X}{d t}=M(X(t),t) X(t) + I(X(t),t)=\prime{X}(t)$ -->
# By substituting the solution $X(t)$ we get the recipes to compute:
# $$
# X_p(t) = M^{-1}(X(t),t)I(X(t),t)  \\ 
# X_c(t) = X(t)-X_p(t) \\ 
# $$
# we see that all the ingredients become explicit functions of time.   
# Since all values are taken at the same time $t$ we can drop the time dependence
# in the notation and write an equation we can use in the iterator.
# $$
# X_p = M^{-1}I(X,t)  \\ 
# X_c = X + X_p \\ 
# $$
#
# ##### Residence time
# The influx $I$ can always be written as $I=b u$ where the scalar $u=\sum_{k=1\dots n} I_k$  and the dimensionless vector $b=I/u$ where $\sum_{k=1\dots n} b_k =1$.
# Assumimg that the pool contents (the components of $X$)  have dimension $mass$ we can infer from eq. (1) that $M$ has dimension $\frac{1}{time}$.
# The components of the (inverse) matrix $M^{-1}$ have therefore dimension $time$. Accordingly the product $RT= M^{-1} b$ is a vector of the same shape as $X$  whose components have dimesion $time$.
# In the context of the Traceability Framework $RT$ is therefore called *residence time*.
#
# Notes on nomenclature: 
# 1. The term *residence time* is not universally used with the same connotation outside the context of the *Traceability Analysis*.
#
# 1. It is not *the time of residence* of the particles in the system for the following reasons:
#     1. In well mixed systems particles can reside in a pool for different times from zero to infinity.
#     1. You could compute the mean of these times over all particles exiting a pool, but even then the result is in general not equal to the above mentioned $rt$.
#     1. The mean residence time would only coincide with the definition above if the system was in equilibrium (which it clearly is not as e.g $NPP(t)$ shows.)
#     1. The origin of the term is probably most easily understood as the generalization of a one dimensional rate equation $\frac{d}{dt} x = m x + u$ 
#        If $r$ and $u$ are constant then the mean residence time is $rt= m^{-1}$. If we start with the rate as property of the model the *residence time* 
#        can be defined as the inverse of this rate. The above definition is the generalization of this simple relationship to matrices and vectors.
#        The matrix $M^{-1}$ takes the role of the number $\frac{1}{m}$ . In the context of the *Traceability Analysis* $M^{-1}$ is called *Chasing Time*. 
#

# +
# lets build an iterator to trace  X, X_c and X_p
# we will extend it further later (for more )
# but we can also make it faster because we are not interested in
# the respiration this time

# build a new template for the StartVector 
# at the moment 3 times the length of the vector of pool contents.
# we will add more components later
svt=mvs.get_StateVariableTuple()

StartVectorTrace=namedtuple(
    "StartVectorTrace",
    [str(v) for v in svt]+
    [str(v)+"_p" for v in svt]+
    [str(v)+"_c" for v in svt]+
    [str(v)+"_RT" for v in svt]
)


# -

# now build the iterator to deal with such vectors
def make_daily_iterator_sym_trace(
        mvs,
        V_init: StartVectorTrace,
        par_dict,
        func_dict
    ):
    B_func, I_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
    V_arr=np.array(V_init).reshape(-1,1) #reshaping for matmul which expects a one column vector (nr,1) 
    
    n=len(mvs.get_StateVariableTuple())
    def f(it,V):
        #the pools are the first n values
        X = V[0:n] 
        I = I_func(it,X) 
        # we decompose I
        u=I.sum()
        b=I/u
        B = B_func(it,X)
        B_inf = np.linalg.inv(B)
        X_new = X + I + B @ X
        X_p = B_inf @ I
        X_c = X_new+X_p
        RT = B_inf @ b
        V_new = np.concatenate(
            (
                X_new.reshape(n,1),
                X_p.reshape(n,1),
                X_c.reshape(n,1),
                RT.reshape(n,1),
            ),
            axis=0
        )
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


# +
# test the new iterator

# first build a new s
# actually we realy can choose only the startvalues for the pools
# but the iterator now produces a 3 times longer vector 
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
# we make the X_p and X_c  parts compatible with the ones computed by the iterator 
# for following timesteps (to keep it) 
# As you can see in the definition of the iterator these values have no impact on further results  
I = u_func(0,X_0)
u=I.sum()
b=I/u
B = B_func(0,X_0)
B_inf = np.linalg.inv(B)
X_p_0 = B_inf@I
X_c_0 = X_0+X_p_0
RT_0 = B_inf@b
# combine the three 
#here we rely on order to be consistent 
#(although we could use also the names of the namedtuple)
V_arr=np.concatenate((X_0,X_p_0,X_c_0,RT_0),axis=0 )
V_init=StartVectorTrace(*V_arr)

# -

it_sym_trace = make_daily_iterator_sym_trace(
    mvs,
    V_init=V_init,
    par_dict=par_dict,
    func_dict=func_dict
)
# we will run the model for 15 steps
ns=1500
nv=len(V_init)
res_trace= np.zeros((ns,nv))
for i in range(ns):
    res_trace[i,:]=it_sym_trace.__next__().reshape(len(V_init),)
    #print(it_sym_trace.__next__().shape)
res_trace

# +
import matplotlib.pyplot as plt
n=len(mvs.get_StateVariableTuple())
fig=plt.figure(figsize=(20,(n+1)*10), dpi=80)
axs=fig.subplots(n+1,2)
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
    [npp_func(d) for d in days],
    label='NPP',
    color='green'
)
axs[n,0].legend()

# -

#Remark:
# For simple matrices it is possible to compute the inverse M^{-1} 
# symbolically ONCE
mvs.get_CompartmentalMatrix().inv()
# and then just evaluate it for the X and t along the solution.
# This could potentially be MUCH faster that inverting 
# the numerical version of the matrix in every timestep.
# However it could be very difficult to compute the symbolic inverse in 
# some cases at all (why we will demonstrate the numeric approach) first. 






