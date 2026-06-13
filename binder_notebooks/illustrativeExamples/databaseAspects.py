# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# This illustrative notebook shows the interplay of all three packages on the collection of predefined models in bgc_md2 and ways to extend this collection by defining a new model that can be compared with respect to diagnostic variables that are computed using LAPM and CompartmentalSystems.
#

# +
# adjust the output
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# %load_ext autoreload
# %autoreload 2
import bgc_md2.helper as h
# -

# # Database aspects of bgc_md
# ## Collection
# First we demonstrate the collection aspect of bgc_md by listing all the models.<p>
#  **Fixme mm**: 
#   * *activate the thumbnail plots again*
#   
#   
#

h.list_target_models()

# ## Inspection with respect to common diagnostics

# ## Queriyng
# One purpose behind the unified (and therefore restricted)  representation of data of a database is the ability to select subsets. In bgc_md a record typically describes what we know about a model or a particular simulations.
# We demonstrate the querying capability by selectiong models for which the variables
# `VegetationCarbonInputPartitioningTuple` and `NumericSolutionArray` are defined or computable which indicates models that include vegetation pools and for which we also have enough example data to run a simulation.<p>
# **Fixme mm** <p>
#
# * There are no examples for the computation of a transit time or age density computation. Not even a variable with that name exists (in bgc_md2) yet.

from bgc_md2.resolve.mvars import(
    CompartmentalMatrix,
    StateVariableTuple,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonInputTuple,
    NumericSolutionArray
)    
import importlib
from ComputabilityGraphs.CMTVS import CMTVS 

li = h.list_target_models(
    target_classes=frozenset(
        {
            VegetationCarbonInputPartitioningTuple,
            NumericSolutionArray
        }
    ),
    # explicit_exclude_models=frozenset({'CARDAMOM'})
)
li    


# +
# we can create a little dictionary of our selected models
model_dict = {
    mn: importlib.import_module(f"bgc_md2.models.{mn}.source").mvs
    for mn in li
}    

model_dict["TECOmm"].computable_mvar_types()
# -

# we can also try to find out with respect to which computable varaibles we can compare them
for k,v in model_dict.items():
    print(k)
    print(v.computable_mvar_types())

from functools import reduce
reduce(
    lambda acc,s: frozenset.intersection(acc,s),
    map(
        lambda val:val.computable_mvar_types(),
        model_dict.values()
    )
)

for k in ['Hilbert1991AnnBot']:
    display(model_dict[k].get_CompartmentalMatrix())


