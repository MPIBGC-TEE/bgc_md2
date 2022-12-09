# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
# adjust the output to full width
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# make changes to imported files immidiately available 
# avoiding the need to reload (in most cases)
# %load_ext autoreload
# %autoreload 2
# -

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
from ComputabilityGraphs.helpers import module_computers
#from bgc_md2.helper import module_computers
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
    "leaf": "vegegation leaf pool",
    "wood": "vegetation wood pool content",
    "I_wood": "Influx into vegetation wood pool",
    "k_wood_o": "out flux rate of wood pool",
    "k_leaf_2_wood": "constant internal flux rate from leaf to wood", 
    "k_wood_2_leaf": "constant internal flux rate from wood to leaf", 
}
# Make symbols from  the strings that we can later use in expressions  
# leaf, wood,...
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# We will also use some symbolic functions ("symbols" with an argument) 
func_dict={
    "I_leaf": "Influx into vegetation leaf pool",
    "k_leaf_o": "out flux rate of leaf pool",
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")

mvs = CMTVS(
    {
        StateVariableTuple((leaf, wood)),
        t,
        InFluxesBySymbol({leaf: I_leaf(t), wood: I_wood}),
        OutFluxesBySymbol({leaf: k_leaf_o(t) * leaf, wood: k_wood_o * wood}),
        InternalFluxesBySymbol({(leaf, wood): k_leaf_2_wood * leaf, (wood, leaf): k_wood_2_leaf * wood}),
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
# Tasks:
# To see what it can do with this information add a new cell an type `mvs.` and press the `tab` key `->`. This will show you the available methods, in other words what can be computed from the provided information.
#
# I wanted to see the compartmental the pools, the matrix and the inputs.

mvs.get_CompartmentalMatrix()

mvs.get_CompartmentalMatrix()

mvs.get_InputTuple()

# we can also print the whole mass balance equation
import bgc_md2.display_helpers as dh
dh.mass_balance_equation(mvs)

# we can also plot a picture
import bgc_md2.helper as h
h.compartmental_graph(mvs)

# ## Extend the minimal model
#
# ### Add more state variables and fluxes for soil pools
# Task:
# Add new pools for litter (slit), coarse woody debris (scwd) and soil organic matter (som)
# and fluxes from the vegetation pools to the new pools.
#
# To avoid (at the moment) unintelligibleerror messages from the depth of the package, go in small steps:
# - first add only a new symbol (e.g. for a new pool), 
# - then add the pool to StateVariableTuple
# - then add a new symbol for a fluxrate or a symbolic function  
# - then add the flux using the rates 
# - repeat with the next pool,flux
#
# The result will look some what like this:

# +
# Make a small dictionary for the variables we will use
sym_dict={
    "leaf": "vegegation leaf pool",
    "wood": "vegetation wood pool content",
    "lit":"soil litter pool",
    "som": "soil organic litter",
    "cwd": "soil coarse woody debris",
    "I_wood": "Influx into vegetation wood pool",
    "k_wood_o": "out flux rate of wood pool",
    "k_som_o": "out flux rate of som pool",
    "k_leaf_2_wood": "constant internal flux rate from leaf to wood", 
    "k_wood_2_leaf": "constant internal flux rate from wood to leaf", 
    "k_leaf_2_lit": "flux rate from leaf to litter",
    "k_wood_2_cwd": "flux rate from wood to coarse woody debris",
    "k_lit_2_som": "flux rate from litter to som",
    "k_cwd_2_som": "flux rate from coarse woody debris to som",
}
# Make symbols from  the strings that we can later use in expressions  
# leaf, wood,...
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# We will also use some symbolic functions ("symbols" with an argument) 
func_dict={
    "I_leaf": "Influx into vegetation leaf pool",
    "k_leaf_o": "out flux rate of leaf pool",
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")

mvs = CMTVS(
    {
        StateVariableTuple((leaf, wood, lit, cwd, som)),
        t,
        InFluxesBySymbol({leaf: I_leaf(t), wood: I_wood}),
        OutFluxesBySymbol({
            leaf: k_leaf_o(t) * leaf,
            wood: k_wood_o * wood,
            som: k_som_o * som,
        }),
        InternalFluxesBySymbol({
            (leaf, wood): k_leaf_2_wood * leaf, 
            (wood, leaf): k_wood_2_leaf * wood,
            (leaf, lit): k_leaf_2_lit * leaf,
            (wood, cwd): k_wood_2_cwd * wood,
            (lit, som): k_lit_2_som *lit,
            (cwd, som): k_cwd_2_som *cwd,
        }),
    },

    computers=module_computers(bgc_md2.resolve.computers)
)
# -

h.compartmental_graph(mvs)

# ### Task: Compute the aggregated flux from the vegetation to soil compartments.
# We want to decompose the model into a soil and vegetation part specify the pools belonging to  vegetation  and soil. 
# We show how do it by hand first and then use ComputabilityGraphs to show us the shortcuts.
# #### manual Approach

# +
vcsvt = (wood,leaf)
scsvt = (lit, som, cwd)
# find the internal fluxes that have a vegetation pool as source and a soil pool target
v2sfls= {k:v for k,v in mvs.get_InternalFluxesBySymbol().items() if (k[0] in vcsvt and k[1] in scsvt) }

# find the internal fluxes that have a vegetation pool as source and a soil pool target
s2vfls= {k:v for k,v in mvs.get_InternalFluxesBySymbol().items() if (k[0] in scsvt and k[1] in vcsvt) }
v2sfls,s2vfls

# now we can sum those fluxes up
v2sAgg = sum(
    v2sfls.values()
)
s2vAgg = sum (
    s2vfls.values()
)
v2sAgg,s2vAgg

# -

# #### Use ComputatbilityGraphs
#
# 1. To find out what `mvars` are available to describe what we know about the model we first look at all potentially computable properties.
# 1. We pick one with a promising name and check if it is available. 

import ComputabilityGraphs.helpers as cgh
computers=module_computers(bgc_md2.resolve.computers)
types=cgh.all_mvars(computers)

# Since we are interested in the Fluxes from the Vegetation part to the soil Part, *AggregatedVegetation2SoilCarbonFlux* sounds promising.
# Let's look up the documentation:

# +
# #?bgc_md2.resolve.mvars.AggregatedVegetation2SoilCarbonFlux
# -

# We now use the `ComputabilityGraphs` package to find out what information (mvars) we have to provide to compute the `AggregatedVegetation2SoilCarbonFlux`

# +

ca=h.numbered_aliases("f",computers)
ta=h.numbered_aliases("T",types)
ta

# +
import ComputabilityGraphs.or_graph_helpers as ogh
from bgc_md2.resolve.mvars import AggregatedVegetation2SoilCarbonFlux

og=ogh.t_tree(
    root_type=AggregatedVegetation2SoilCarbonFlux,
    available_computers=computers,
    avoid_types=frozenset({}),
    #given_types=mvs.provided_mvar_types
)
#import matplotlib.pyplot as plt
#ax=plt.subplot()
#og.to_networkx_graph(avoid_types=ogh.TypeSet({})).draw_matplotlib(ax)
og.psts
# -

og.jupyter_widget(
    computer_aliases_tup=ca,
    type_aliases_tup=ta,
    given=mvs.provided_mvar_types
)

import matplotlib.pyplot as plt
from ComputabilityGraphs.dep_graph_helpers import ( 
    #duplicated_computer_dict, 
    all_dep_graphs
)
root=bgc_md2.resolve.mvars.AggregatedVegetation2SoilCarbonFlux
cs=computers
gs= set(
    all_dep_graphs(
        root_type=root,
        cs=cs,
        #given=frozenset()
        given=mvs.provided_mvar_types
    )
)
len(gs)


# +
#n=len(gs)
#fig=plt.figure(figsize=(20,20*n))
#if n>1:
#    axs = fig.subplots(n,1)
#    for i,ax in enumerate(axs):
#        dg = list(gs)[i]
#        dg.draw_matplotlib(ax)
#else:
#    dg=list(gs)[0]
#    ax=fig.subplots(1,1)
#    dg.draw_matplotlib(ax)
#fig.savefig("depgraphs_" + root.__name__ + '.pdf')
#dg.required_mvars(given=mvs.provided_mvar_types)


type(mvs.provided_mvar_types)


# +
fig=plt.figure(figsize=(20,10))
dg=list(gs)[0]

ax = fig.add_subplot(1, 1, 1)
B=dg.to_bipartite()
mvs.provided_mvar_types.difference(B.types())

B.draw_matplotlib(
    ax,
    given=mvs.provided_mvar_types,
    #computer_aliases=ca,
    #type_aliases=ta
)
#B.draw_matplotlib(ax,target_node=root)


# +
#[n for n in B.nodes if B.out_degree(n)==0]
#dg.jupyter_widget(
#    given=mvs.provided_mvar_types,
#    #computer_aliases=ca,
#    type_aliases=ta
#)
# -

#from bokeh.io import output_file, show, output_notebook
#from ComputabilityGraphs.rec_graph_helpers import fast_graph
#from ComputabilityGraphs.fast_graph_helpers import project_to_multiDiGraph
#from ComputabilityGraphs.graph_plotting import (
#    draw_ComputerSetMultiDiGraph_matplotlib,
#    bokeh_plot
#)
##from frozenset import frozenset
#
#fg = fast_graph(
#    cs=cs,
#    root_type=root,
#    given=mvs.provided_mvar_types
#    #given=frozenset()
#)
#G=project_to_multiDiGraph(fg)
#output_notebook()
##plot=bokeh_plot(G,frozenset({bgc_md2.resolve.mvars.VegetationCarbonCompartmentalMatrix}))
#plot=bokeh_plot(
#    G
#    ,frozenset({bgc_md2.resolve.mvars.VegetationCarbonInFluxesBySymbol})
#)
#show(plot)
##fig=plt.figure(figsize=(20,20))
##ax1 = fig.add_subplot(1, 1, 1)
##fg.draw_matplotlib(
##    ax1
##
##)



# +

#fig=plt.figure(figsize=(20,20))
#ax1 = fig.add_subplot(1, 1, 1)
#draw_ComputerSetMultiDiGraph_matplotlib(ax=ax1,spsg=project_to_multiDiGraph(fg))