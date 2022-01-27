# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# This illustrative notebook shows how to create a representation for a new model.
#
# Technical Note:
# To keep the paths short we assume that you have 
# - created a folder where your model should live, 
# - have copied this notebook there, 
# - changed directory there, and 
# - run the `jupyter notebook` command in this directory. 
#
# This is just convinient since the python interpreter includes it's working directory automatically in its search path for modules or scripts. 
# If your jupyter server was started somewhere else you can use `pwd` and  `cd` directly in this notebook to find out where you are and go to the target directory or deliberatly work from somewhere else and change the paths.
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
# We will start with an extremely simple model, that you can find in 
# ```bash
# bgc_md2/src/bgc_md2/models/testVectorFree/source.py
# ```
# Copy its contents into a new cell. (we will later save it under a new name) 

# +
from sympy import var 
from ComputabilityGraphs.CMTVS import CMTVS
#from bgc_md2.helper import bgc_md2_computers
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

var("I_vl I_vw vl vw k_vl k_vw")

mvs = CMTVS(
    {
        StateVariableTuple((vl, vw)),
        TimeSymbol("t"),
        InFluxesBySymbol({vl: I_vl, vw: I_vw}),
        OutFluxesBySymbol({vl: k_vl * vl, vw: k_vw * vw}),
        InternalFluxesBySymbol({(vl, vw): k_vl * vl, (vw, vl): k_vw * vw}),
    },

    computers=module_computers(bgc_c)
)
# -

# The last statement in the code defines a variable `mvs` which is 
# an instance of CMTVS which stands for `C`onnected`M`ulti`T`ype`V`ariable`S`et".
# It contains information in two forms. 
# 1. Variables of certain type (like InFluxesBySymbol)
# 2. a set of functions that "connect" these Variables and to other results we did not specify 
#    explicitly.
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
# - add new symbols to desricbe pools for the cable model (from Alisons code) but
#   leave the old ones there since they are presently used in the fluxes.
#   We will replace them one bye one later 
#   but this way we are not forced to replace all the fluxes at
#   once.It's always good to be able to make small steps...
# - use a different way to declare symbols 
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

var("I_vl I_vw vl vw k_vl k_vw")
# we organize the new symbols a bit better by putting them in a dictionary
# together with some description that we can use later to display some metainformation
sym_dict = {
        'C_leaf': 'content of leaf pool',
        'C_root': 'content of root pool',
        'C_wood': 'content of wood pool',
        'C_metlit': 'content of metabolic litter pool',
        'C_strlit': 'content of static litter pool',
        'C_cwd': 'content of corse woody debris pool',
        'C_mic': 'content of microbial pool',
        'C_slowsom': 'content of the slow soil pool',
        'C_passsom': 'content of the passive soil pool'
}
# for the moment only use the 
var(list(sym_dict.keys()))


mvs = CMTVS(
    {
        StateVariableTuple((vl, vw)),
        TimeSymbol("t"),
        InFluxesBySymbol({vl: I_vl, vw: I_vw}),
        OutFluxesBySymbol({vl: k_vl * vl, vw: k_vw * vw}),
        InternalFluxesBySymbol({(vl, vw): k_vl * vl, (vw, vl): k_vw * vw}),
    },

    computers=module_computers(bgc_c)
)
# -

# Nothing has changed in the model description but be have some more symbols to work with.
# We could type `C_leaf` somewhere without an error since it is now known as a variable.
# In the next step we replace all occurences of `vl` by `C_leaf` 

from sympy import var
var("a b")
a**2


