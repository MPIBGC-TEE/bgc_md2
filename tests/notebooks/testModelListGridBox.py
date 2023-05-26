# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
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
# -

# # Database aspects of bgc_md
# ## Collection
# First we demonstrate the collection aspect of bgc_md by listing all the models.<p>
#  **Fixme mm**: 
#   * *activate the thumbnail plots again*
#   
#   
#

# +
import bgc_md2.helper as h

model_inspection = h.MvarSetInspectionBox()
model_list = h.ModelListGridBox(
    inspection_box=model_inspection,
)
model_list
# -

# ## Inspection with respect to common diagnostics
# The variable in the next cell is connected to the list above and will show a summary of the model you click. This summary contains not only the data the author of the model put into the database but also all the results that can be computed from it. It also contains a button to create a model specific notebook.

model_inspection

# ## Queriyng
# One purpose behind the unified (and therefore restricted)  representation of data of a database is the ability to select subsets. In bgc_md a record typically describes what we know about a model or a particular simulations.
# We demonstrate the querying capability by selectiong models for which the variables
# `VegetationCarbonInputPartitioningTuple` and `NumericSolutionArray` are defined or computable which indicates models that include vegetation pools and for which we also have enough example data to run a simulation.<p>
# **Fixme mm** <p>
#
# * There are no examples for the computation of a transit time or age density computation. Not even a variable with that name exists (in bgc_md2) yet.

from bgc_md2.resolve.mvars import VegetationCarbonInputPartitioningTuple, NumericSolutionArray
from bgc_md2.resolve.MVarSet import MVarSet

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


# From these we chose two models to investigate more thoroughly.
# We load the first one with a helper function and the second by using standard python tools.
#
#

# # ComputabilityGraphs
# Now that we actually have two records (MVarSets) we can exlore what we can comptute from them.
# Just add a dot "." behind the variable in the next cell and press the tab key!
#
