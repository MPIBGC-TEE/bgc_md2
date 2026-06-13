# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.0
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
from trendy9helpers import general_helpers as gh
import MIP_output_helpers as moh
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    StateVariableTuple
)

# ### Selecting models to compare

model_names={
    #"ab_classic":"CLASSIC",  
    #"clm5":"CLM5.0",
    #"kv_ft_dlem": "DLEM", 
    #"bian_ibis2":"IBIS",    
    #"cj_isam": "ISAM",    
    #"isba-ctrip":"ISBA-CTRIP",    
    #"jsbach":"JSBACH",
    #"yz_jules": "JULES-ES-1p0",    
    #"lpj-guess":"LPJ-GUESS",
    #"lpjwsl":"LPJ",
    #"lpx-bern":"LPX-Bern",
    #"ORCHIDEE-V2":"OCN",    
    #"ORCHIDEE":"ORCHIDEE",
    #"ORCHIDEE-CNP":"ORCHIDEE-CNP",    
    #"ORCHIDEEv3":"ORCHIDEEv3",
    #"Aneesh_SDGVM":"SDGVM",
    #"kv_visit2": "VISIT",
    "jon_yib": "YIBs"    
}
model_folders=[(k) for k in model_names]
m_names=list(model_names.values())
experiment_names_S2=list()
for name in m_names:
    experiment_names_S2.append(name + "_S2_") 
experiment_names_S3=list()
for name in m_names:
    experiment_names_S3.append(name + "_S3_") 
#experiment_names_S2

# ### Loading TRENDY data and model parameters

vars_all_list_S2 = moh.get_vars_all_list(model_folders, experiment_names_S2)
#vars_all_list_S3 = moh.get_vars_all_list(model_folders, experiment_names_S3)
#
##np.mean(vars_all_list_S2[2].npp)
#np.mean(vars_all_list_S2[0].npp)
##vars_all_list_S3[2].ra
#
## define same step size for all models (in days)
#delta_t_val=30
#
#all_comp_dict_S2=moh.get_global_components_from_output(model_names=model_names,
#             vars_all_list=vars_all_list_S2,
#             delta_t_val=delta_t_val, 
#             #model_cols=model_cols,
#             part=1,
#             #averaging=12*30//delta_t_val, # yearly averaging
#             #averaging=30//delta_t_val, # monthly averaging
#             #overlap=True,
#             start_shift=109,
#             #end_shift=4
#             )
#all_comp_dict_S3=moh.get_global_components_from_output(model_names=model_names,
#             vars_all_list=vars_all_list_S3,
#             delta_t_val=delta_t_val, 
#             #model_cols=model_cols,
#             part=1,
#             #averaging=12*30//delta_t_val, # yearly averaging
#             #averaging=30//delta_t_val, # monthly averaging
#             #overlap=True,
#             start_shift=109,
#             #end_shift=4
#             )
#
## +
##all_comp_dict_S2["DLEM"]["nep"]
## -
#
## ### Plots of traceable components and their uncertainty
#
#plt.rcParams.update({'font.size': 15})
#import plotly.express as px
#cols = px.colors.qualitative.Light24
#model_cols = {m_names[i]: cols[i] for i in range(len(m_names))}
## import seaborn as sns # requires installation of the "seaborn" package (not part of bgc_md2)
##cols = sns.color_palette("tab20", len(m_names))
## cols = sns.color_palette("cubehelix", len(m_names))
#
## +
## x_x_c1,sigma_x_x_c1=moh.plot_traceable_component(
##     all_comp_dict_S2,
##     "x_x_c",
##     model_cols,
##     #delta=True,
## )
## x_x_c2,sigma_x_x_c2=moh.plot_traceable_component(
##     all_comp_dict_S3,
##     "x_x_c",
##     model_cols,
##     #delta=True,
## )
#
## +
## x1,sigma_x1=moh.plot_traceable_component(
##     all_comp_dict_S2,
##     "x",
##     model_cols,
##     #delta=True,
## )
## x2,sigma_x2=moh.plot_traceable_component(
##     all_comp_dict_S3,
##     "x",
##     model_cols,
##     #delta=True,
## )
## -
#
#all_comp_dict_S3["VISIT"]["cVeg"][38:43]
#
#x1,sigma_x1=moh.plot_traceable_component(
#    all_comp_dict_S2,
#    "x",
#    model_cols,
#    delta=True,
#)
#x2,sigma_x2=moh.plot_traceable_component(
#    all_comp_dict_S3,
#    "x",
#    model_cols,
#    delta=True,
#)
#
#nep1,sigma_nep1=moh.plot_traceable_component(
#    all_comp_dict_S2,
#    "nep",
#    model_cols,
#    #delta=True,
#)
#nep2,sigma_nep2=moh.plot_traceable_component(
#    all_comp_dict_S3,
#    "nep",
#    model_cols,
#    #delta=True,
#)
#
#all_comp_dict_S2_1st_half=moh.get_components_from_output(model_names=model_names,
#             vars_all_list=vars_all_list_S2,
#             delta_t_val=delta_t_val, 
#             part=1,
#             end_shift=52,
#             start_shift=109                                   
#             )
#all_comp_dict_S3_1st_half=moh.get_components_from_output(model_names=model_names,
#             vars_all_list=vars_all_list_S3,
#             delta_t_val=delta_t_val, 
#             part=1,
#             end_shift=52,
#             start_shift=109                                   
#             )
#
#x1,sigma_x1=moh.plot_traceable_component(
#    all_comp_dict_S2_1st_half,
#    "x",
#    model_cols,
#    delta=True,
#)
#x2,sigma_x2=moh.plot_traceable_component(
#    all_comp_dict_S3_1st_half,
#    "x",
#    model_cols,
#    delta=True,
#)
#
#all_comp_dict_S2_2nd_half=moh.get_components_from_output(model_names=model_names,
#             vars_all_list=vars_all_list_S2,
#             delta_t_val=delta_t_val, 
#             part=1,
#             end_shift=2,
#             start_shift=57                                   
#             )
#all_comp_dict_S3_2nd_half=moh.get_components_from_output(model_names=model_names,
#             vars_all_list=vars_all_list_S3,
#             delta_t_val=delta_t_val, 
#             part=1,
#             end_shift=2,
#             start_shift=57                                   
#             )
#
#x1,sigma_x1=moh.plot_traceable_component(
#    all_comp_dict_S2_2nd_half,
#    "x",
#    model_cols,
#    delta=True,
#)
#x2,sigma_x2=moh.plot_traceable_component(
#    all_comp_dict_S3_2nd_half,
#    "x",
#    model_cols,
#    delta=True,
#)
#
#
#
## +
## times=all_comp_dict_S2["Times"]
## var=sigma_x1#/x1*100
## moh.plot_single_trend(var,times,3,"Standard deviation of X over time - S2")
#
## times=all_comp_dict_S3["Times"]
## var=sigma_x2#/x2*100
## moh.plot_single_trend(var,times,3, "Standard deviation of X over time - S3")
#
## times=all_comp_dict_S3["Times"]
## var=sigma_x2/x2*100-sigma_x1/x1*100
## moh.plot_single_trend(var,times,3, "Standard deviation in % of X over time - S3-S2")
#
## +
## x_c1,sigma_x_c1=moh.plot_traceable_component(
##     all_comp_dict_S2,
##     "x_c",
##     model_cols,
##     delta=True,
## )
## x_c2,sigma_x_c2=moh.plot_traceable_component(
##     all_comp_dict_S3,
##     "x_c",
##     model_cols,
##     delta=True,
## )
#
## +
## times=all_comp_dict_S2["Times"]
## var=sigma_x_c1
## moh.plot_single_trend(var,times,3,"Standard deviation of X_c over time - S2")
#
## times=all_comp_dict_S3["Times"]
## var=sigma_x_c2
## moh.plot_single_trend(var,times,3, "Standard deviation of X_c over time - S3")
#
## times=all_comp_dict_S3["Times"]
## var=sigma_x_c2-sigma_x_c1
## moh.plot_single_trend(var,times,3, "Standard deviation of X_c over time - S3-S2")
#
## +
## x_p1,sigma_x_p1=moh.plot_traceable_component(
##     all_comp_dict_S2,
##     "x_p",
##     model_cols,
##     #delta=True,
## )
## x_p2,sigma_x_p2=moh.plot_traceable_component(
##     all_comp_dict_S3,
##     "x_p",
##     model_cols,
##     #delta=True,
## )
#
## +
## times=all_comp_dict_S2["Times"]
## var=x_p1/x1*100
## moh.plot_single_trend(var,times,3,"Mean X_p in % of X over time")
## times=all_comp_dict_S3["Times"]
## var=x_p2/x2*100
## moh.plot_single_trend(var,times,3,"Mean X_p in % of X over time")
## -
#
#u1, sigma_u1=moh.plot_traceable_component(
#    all_comp_dict_S2,
#    "npp",
#    model_cols,
#    delta=True,
#)
#u2, sigma_u2=moh.plot_traceable_component(
#    all_comp_dict_S3,
#    "npp",
#    model_cols,
#    delta=True,
#)
#
## +
## times=all_comp_dict_S2["Times"]
## var=sigma_u1/u1*100
## moh.plot_single_trend(var,times,3,"Standard deviation in % of u over time - S2")
#
## times=all_comp_dict_S3["Times"]
## var=sigma_u2/u2*100
## moh.plot_single_trend(var,times,3, "Standard deviation in % of u over time - S3")
#
## times=all_comp_dict_S3["Times"]
## var=(sigma_u2/u2*100-sigma_u1/u1*100)
## moh.plot_single_trend(var,times,3, "Uncertainty of u increase S3-S2 in % of S2")
## -
#
#rt1,sigma_rt1=moh.plot_traceable_component(
#    all_comp_dict_S2,
#    "rt",
#    model_cols,
#    delta=True,
#)
#rt2,sigma_rt2=moh.plot_traceable_component(
#    all_comp_dict_S3,
#    "rt",
#    model_cols,
#    delta=True,
#)
#
## +
## times=all_comp_dict_S2["Times"]
## var=sigma_rt1/rt1*100
## moh.plot_single_trend(var,times,3,"Standard deviation in % of rt over time - S2")
#
## times=all_comp_dict_S3["Times"]
## var=sigma_rt2/rt2*100
## moh.plot_single_trend(var,times,3, "Standard deviation in % of rt over time - S3")
#
## times=all_comp_dict_S3["Times"]
## var=(sigma_rt2/rt2*100-sigma_rt1/rt1*100)
## moh.plot_single_trend(var,times,3, "Uncertainty of rt increase S3-S2 in % of S2")
## -
#
#
#
## +
## times=all_comp_dict_S2["Times"]
## var=sigma_nep1
## moh.plot_single_trend(var,times,3,"Standard deviation of nep over time - S2")
#
## times=all_comp_dict_S3["Times"]
## var=sigma_nep2
## moh.plot_single_trend(var,times,3, "Standard deviation of nep over time - S3")
#
## times=all_comp_dict_S3["Times"]
## var=(sigma_nep2-sigma_nep1)/sigma_u1*100
## moh.plot_single_trend(var,times,3, "Uncertainty of nep increase S3-S2 in % of S2")
## -
#
#cSoil1,sigma_cSoil1=moh.plot_traceable_component(
#    all_comp_dict_S2,
#    "cSoil",
#    model_cols,
#    delta=True,
#)
#cSoil2,sigma_cSoil2=moh.plot_traceable_component(
#    all_comp_dict_S3,
#    "cSoil",
#    model_cols,
#    delta=True,
#)
#
#rh1,sigma_rh1=moh.plot_traceable_component(
#    all_comp_dict_S2,
#    "rh",
#    model_cols,
#    #delta=True,
#)
#rh2,sigma_rh2=moh.plot_traceable_component(
#    all_comp_dict_S3,
#    "rh",
#    model_cols,
#    #delta=True,
#)
#
#cVeg1,sigma_cVeg1=moh.plot_traceable_component(
#    all_comp_dict_S2,
#    "cVeg",
#    model_cols,
#    delta=True,
#)
#cVeg2,sigma_cVeg2=moh.plot_traceable_component(
#    all_comp_dict_S3,
#    "cVeg",
#    model_cols,
#    delta=True,
#)
#
#ra1,sigma_ra1=moh.plot_traceable_component(
#    all_comp_dict_S2,
#    "ra",
#    model_cols,
#    #delta=True,
#)
#ra2,sigma_ra2=moh.plot_traceable_component(
#    all_comp_dict_S3,
#    "ra",
#    model_cols,
#    #delta=True,
#)
#
#all_comp_dict_S2
#
## ### Uncertainty attribution
#
#plt.rcParams.update({'font.size': 12})
#all_comp_dict_S3['VISIT']['rt']
#
#moh.plot_attribution_sum (
#    all_comp_dict=all_comp_dict_S2,
#    percent=True,
#    part=1,
#)
#moh.plot_attribution_sum (
#    all_comp_dict=all_comp_dict_S3,
#    percent=True,
#    part=1,
#)
#
#
#
## +
#from collections import namedtuple
#
#data_str = namedtuple(
#    'data_str',
#    ["cVeg", "cSoil_total","cVeg_diff", "cSoil_total_diff","nep", "gpp", "ra", "rh",  "dist", "f_v2s", 
#    "X", "X_c", "X_p","RT", "RT_veg", "RT_soil"]
#    )  
#data_str._fields
## -
#
#avd=moh.get_global_mean_uncertainty(
#    dataPath="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",  
#    experiment_name="S2", # string, e.g. "S2"
#    data_str=data_str, # named tuple
#    avd=True
#    )
#mean=moh.get_global_mean_uncertainty(
#    dataPath="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",  
#    experiment_name="S2", # string, e.g. "S2"
#    data_str=data_str, # named tuple
#    avd=False
#    )
#
#avd
#
#mean
#
#mean.gpp*mean.RT/120
#
#moh.plot_attribution_C_storage_global (
#    all_comp_dict=all_comp_dict_S2,
#    #percent=True,
#    #part=1,
#)
#moh.plot_attribution_C_storage_global (
#    all_comp_dict=all_comp_dict_S3,
#    #percent=True,
#    #part=1,
#)
#
#moh.plot_attribution_C_sink_global (
#    all_comp_dict=all_comp_dict_S2,
#    #percent=True,
#    #part=1,
#)
#moh.plot_attribution_C_sink_global (
#    all_comp_dict=all_comp_dict_S3,
#    #percent=True,
#    #part=1,
#)
