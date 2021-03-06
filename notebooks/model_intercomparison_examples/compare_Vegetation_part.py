# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from IPython.display import HTML

display(HTML("<style>.container { width:100% !important; }</style>"))

# %load_ext autoreload
# %autoreload 2
import bgc_md2.helper as h

model_inspection = h.MvarSetInspectionBox()

from bgc_md2.resolve.mvars import CompartmentalMatrix, StateVariableTuple, VegetationCarbonInputPartitioningTuple,VegetationCarbonInputTuple

# This time we are only interested in Vegetation models and Carbon input partitioning. Therefore we look for models for which the variable
# `VegetationCarbonInputPartitioningTuple` is defined or computable.

li = h.list_target_models(
    target_classes=frozenset(
        {
            CompartmentalMatrix,
            StateVariableTuple,
            VegetationCarbonInputPartitioningTuple,
            VegetationCarbonInputTuple
            
        }
    ),
    # explicit_exclude_models=frozenset({'CARDAMOM'})
)
li    


from bgc_md2.resolve.MVarSet import MVarSet
for mn in li:
    mvs = MVarSet.from_model_name(mn)
    print('######################')
    print(mn)
    display(mvs._get_single_mvar_value(VegetationCarbonInputTuple))
    display(mvs._get_single_mvar_value(VegetationCarbonInputPartitioningTuple))

model_list = h.GeneralMvarSetListGridBox(
    inspection_box=model_inspection,
    target_classes=(CompartmentalMatrix,StateVariableTuple,VegetationCarbonInputPartitioningTuple),
    #explicit_exclude_models=frozenset({'CARDAMOM'})
)
model_list

model_inspection


