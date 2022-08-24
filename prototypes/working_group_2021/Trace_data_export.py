# -*- coding: utf-8 -*-
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

# # Exporting yearly traceable components to a .csv

from IPython.display import Markdown, display
display(Markdown("TracebilityText.md"))

# ### Loading required packages  and functions

# %load_ext autoreload
# %autoreload 2
import matplotlib.pyplot as plt
import numpy as np
from functools import lru_cache
import general_helpers as gh
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    StateVariableTuple
)

# +
delta_t_val=30

model_names={
    "yz_jules": "JULES",
    "kv_visit2": "VISIT",
}
model_folders=[(k) for k in model_names]
test_arg_list=gh.get_test_arg_list(model_folders)

model_cols={
    "yz_jules": "blue",
    "kv_visit2": "orange",
    "jon_yib": "green",
    "kv_ft_dlem": "red",
    "Aneesh_SDGVM":"yellow",
    "cj_isam": "purple",
    "bian_ibis2":"magenta",
    "ORCHIDEE-V2":"teal",
}

var_names={
    "x": "X",
    "x_c":"X_c",
    "x_p": "X_p",
    "u": "C input",
    "rt": "Residence time",
}

# +
#var_names
# -

output=gh.write_yearly_components(model_names=model_names,
                        test_arg_list=test_arg_list,   
                        var_names=var_names,
                        delta_t_val=delta_t_val,
                        model_cols=model_cols,
                        part=1,
                        averaging=12*30//delta_t_val # yearly averaging
                       )


output


