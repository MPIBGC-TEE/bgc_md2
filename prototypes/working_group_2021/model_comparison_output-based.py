# -*- coding: utf-8 -*-
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

# # O'Sullivan / Koven method

# +
# from IPython.display import Markdown, display
# display(Markdown("TracebilityText.md"))
# -

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

# ### Selecting models to compare

# define models to compare as a dictionary (folder name : model name)
model_names={
    "yz_jules": "JULES",
    "kv_visit2": "VISIT",
    "jon_yib": "YIBs",
    "kv_ft_dlem": "DLEM",
#     "Aneesh_SDGVM":"SDGVM",
    #"cj_isam": "ISAM",
    "bian_ibis2":"IBIS",
    #"ORCHIDEE-V2":"OCN",
}

# ### Loading TRENDY data and model parameters

# define same step size for each model (in days)
delta_t_val=30

# +
# load data and parameters
model_folders=[(k) for k in model_names]
test_arg_list=gh.get_test_arg_list(model_folders)

# fixme mm 8-12: 
# it would make sense to create a dictionary indexed by the model name 
# so it can be used for model comparisons by name like this one
test_args_dictionary={mf: gh.test_args(mf) for mf in model_folders}
# -

# ### Plots of traceable components

plt.rcParams.update({'font.size': 14})

model_cols={
    "JULES": "blue",
    "VISIT": "orange",
    "YIBs": "green",
    "DLEM": "red",
    "SDGVM":"cyan",
    "ISAM": "purple",
    "IBIS":"magenta",
    "OCN":"teal",
}

all_comp_dict= gh.get_traceable_components(model_names=model_names,
             test_arg_list=test_arg_list,
             delta_t_val=delta_t_val, 
             model_cols=model_cols,
             part=1,
             averaging=12*30//delta_t_val, # yearly averaging
             #averaging=30//delta_t_val, # monthly averaging
             overlap=True
             )

all_comp_dict_out=gh.get_components_from_output(model_names=model_names,
             test_arg_list=test_arg_list,
             delta_t_val=delta_t_val, 
             model_cols=model_cols,
             part=1,
             averaging=12*30//delta_t_val, # yearly averaging
             #averaging=30//delta_t_val, # monthly averaging
             #overlap=True
             )

# +
# all_comp_dict=all_comp_dict_out
# -

x_x_c,sigma_x_x_c=gh.plot_traceable_component(
    all_comp_dict,
    "x_x_c",
    model_cols,
    #delta=True,
)
x_x_c,sigma_x_x_c=gh.plot_traceable_component(
    all_comp_dict_out,
    "x_x_c",
    model_cols,
    #delta=True,
)

x1,sigma_x1=gh.plot_traceable_component(
    all_comp_dict,
    "x",
    model_cols,
    #delta=True,
)
x1,sigma_x1=gh.plot_traceable_component(
    all_comp_dict_out,
    "x",
    model_cols,
    #delta=True,
)
x2,sigma_x2=gh.plot_traceable_component(
    all_comp_dict,
    "x",
    model_cols,
    delta=True,
)
x2,sigma_x2=gh.plot_traceable_component(
    all_comp_dict_out,
    "x",
    model_cols,
    delta=True,
)

x_c1,sigma_x_c1=gh.plot_traceable_component(
    all_comp_dict,
    "x_c",
    model_cols,
    #delta=True,
)
x_c1,sigma_x_c1=gh.plot_traceable_component(
    all_comp_dict_out,
    "x_c",
    model_cols,
    #delta=True,
)
x_c2,sigma_x_c2=gh.plot_traceable_component(
    all_comp_dict,
    "x_c",
    model_cols,
    delta=True,
)
x_c2,sigma_x_c2=gh.plot_traceable_component(
    all_comp_dict_out,
    "x_c",
    model_cols,
    delta=True,
)

x_p1,sigma_x_p1=gh.plot_traceable_component(
    all_comp_dict,
    "x_p",
    model_cols,
    #delta=True,
)
x_p1,sigma_x_p1=gh.plot_traceable_component(
    all_comp_dict_out,
    "x_p",
    model_cols,
    #delta=True,
)

u1, sigma_u1=gh.plot_traceable_component(
    all_comp_dict,
    "u",
    model_cols,
    #delta=True,
)
u1, sigma_u1=gh.plot_traceable_component(
    all_comp_dict_out,
    "u",
    model_cols,
    #delta=True,
)
u2, sigma_u2=gh.plot_traceable_component(
    all_comp_dict,
    "u",
    model_cols,
    delta=True,
)
u2, sigma_u2=gh.plot_traceable_component(
    all_comp_dict_out,
    "u",
    model_cols,
    delta=True,
)

rt1,sigma_rt1=gh.plot_traceable_component(
    all_comp_dict,
    "rt",
    model_cols,
    #delta=True,
)
rt1,sigma_rt1=gh.plot_traceable_component(
    all_comp_dict_out,
    "rt",
    model_cols,
    #delta=True,
)
rt2,sigma_rt2=gh.plot_traceable_component(
    all_comp_dict,
    "rt",
    model_cols,
    delta=True,
)
rt2,sigma_rt2=gh.plot_traceable_component(
    all_comp_dict_out,
    "rt",
    model_cols,
    delta=True,
)

plt.rcParams.update({'font.size': 12})

gh.plot_attribution_sum (
    all_comp_dict=all_comp_dict,
    percent=True,
    part=1,
)
gh.plot_attribution_sum (
    all_comp_dict=all_comp_dict_out,
    percent=True,
    part=1,
)

gh.plot_attribution_per_model(
    all_comp_dict=all_comp_dict_out,
    #percent=True,    
)


