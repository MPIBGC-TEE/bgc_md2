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

from IPython.display import HTML

display(HTML("<style>.container { width:100% !important; }</style>"))

# %load_ext autoreload
# %autoreload 2
import bgc_md2.helper as h
from bgc_md2.resolve.mvars import CompartmentalMatrix, StateVariableTuple, VegetationCarbonInputPartitioningTuple,VegetationCarbonInputTuple
model_inspection = h.MvarSetInspectionBox()

# This time we are only interested in Vegetation models and Carbon input partitioning. Therefore we look for models for which the variable
# `VegetationCarbonInputPartitioningTuple` is defined or computable.

# this may take some time since everything is computed
model_list = h.GeneralMvarSetListGridBox(
    inspection_box=model_inspection,
    target_classes=(CompartmentalMatrix,StateVariableTuple,VegetationCarbonInputPartitioningTuple),
    #explicit_exclude_models=frozenset({'CARDAMOM'})
)
model_list


