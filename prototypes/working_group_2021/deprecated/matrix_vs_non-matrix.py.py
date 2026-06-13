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

# # Comparison with O’Sullivan et al. (2022)
# <a href=" https://doi.org/10.1038/s41467-022-32416-8">https://doi.org/10.1038/s41467-022-32416-8</a>

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

# ### Selecting models to compare

# ### Loading TRENDY data and model parameters

# ### Plots of traceable components

# +
delta_t_val=1
model_names_test={
    "yz_jules": "JULES",
    "kv_visit2": "VISIT",
}
model_folders_test=[(k) for k in model_names_test]
test_arg_list_test=gh.get_test_arg_list(model_folders_test)

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

add_cols={
    #"yz_jules": "red",
    "kv_visit2": "green",
#     "jon_yib": "green",
#     "kv_ft_dlem": "red",
#     "Aneesh_SDGVM":"yellow",
#     "cj_isam": "purple",
#     "bian_ibis2":"magenta",
#     "ORCHIDEE-V2":"teal",
}


# -

def plot_x_xc_trendy(
    model_names,  # dictionary (folder name : model name)
    test_arg_list,  # a list of test_args from all models involved
    delta_t_val,  # model time step
    model_cols,  # dictionary (folder name :color)
    add_cols,
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging,  # number of iterator steps over which to average results. 1 for no averaging
    overlap=True,  # compute overlapping timeframe or plot whole duration for all models
):
    if (part < 0) | (part > 1):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        )
    model_folders = [(k) for k in model_names]
    fig = plt.figure(figsize=(17, 8))
    ax = fig.subplots(1, 1)
    k = 0
    for mf in model_folders:
        itr = gh.traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        if overlap == True:
            start_min, stop_max = gh.min_max_index(
                test_arg_list[k],
                delta_t_val,
                *gh.t_min_tmax_overlap(test_arg_list, delta_t_val)
            )
        else:
            start_min, stop_max = gh.min_max_index(
                test_arg_list[k],
                delta_t_val,
                *gh.t_min_tmax_full(test_arg_list, delta_t_val)
            )
        # if we do not want the whole interval but look at a smaller part to observe the dynamics
        start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        times = (
            gh.times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            / gh.days_per_year()
        )
        vals = itr[start:stop]
        # print("vals.x")
        # print(vals.x[0:240])
        
        ############## method from O'Sallivan et al. (2022)
        cVeg_trendy=test_arg_list[k].svs.cVeg
        if model_names[mf]=="VISIT":
            cSoil_trendy=test_arg_list[k].svs.cLitter+test_arg_list[k].svs.cSoil
        else:
            cSoil_trendy=test_arg_list[k].svs.cSoil
        rh_trendy=test_arg_list[k].svs.rh
        npp_trendy=test_arg_list[k].dvs.npp
        #npp_trendy_yearly=npp_sum(npp_trendy)
        
        delta_cVeg=np.zeros_like(cVeg_trendy)
        for i in range (len(cVeg_trendy)-1):
            delta_cVeg[i]=cVeg_trendy[i+1]-cVeg_trendy[i]
    
        delta_cSoil=np.zeros_like(cSoil_trendy)
        for i in range (len(cSoil_trendy)-1):
            delta_cSoil[i]=cSoil_trendy[i+1]-cSoil_trendy[i]    
    
        tau_v=cVeg_trendy/(npp_trendy-delta_cVeg)
        X_c_v_trendy=tau_v*npp_trendy
        #X_c_v_trendy_yearly=gh.avg_timeline(X_c_v_trendy,12)

        f_vs=delta_cSoil+rh_trendy

        tau_s=cSoil_trendy/rh_trendy
        X_c_s_trendy=f_vs*tau_s
        #X_c_s_trendy_yearly=gh.avg_timeline(X_c_s_trendy,12)

        X_trendy_total=(cVeg_trendy+cSoil_trendy)[start//averaging//delta_t_val:stop//averaging//delta_t_val]#*148940000*1000000*0.000000000001
        X_c_trendy_total=(X_c_v_trendy+X_c_s_trendy)[start//averaging//delta_t_val:stop//averaging//delta_t_val]#*148940000*1000000*0.000000000001 
        u_trendy=npp_trendy[start//averaging//delta_t_val:stop//averaging//delta_t_val]
        tau_trendy=(tau_v+tau_s)[start//averaging//delta_t_val:stop//averaging//delta_t_val]
#         print(str(mf)+": "+str(X_c_trendy_total.shape)+" ; "+str(X_trendy_total.shape))
#         print ("start: "+str(start)+" stop: "+str(stop))
#         print (tau_trendy)
#         print (vals.rt)
        ###################
        ax.plot(
            gh.avg_timeline(times, averaging),#[:-1],
            gh.avg_timeline(vals.x, averaging)#[:-1]
            * 148940000
            * 1000000
            * 0.000000000001,  # convert to global C in Gt
            label=model_names[mf] + " - X-matrix",
            color=model_cols[mf],            
        )
        ax.plot(            
            gh.avg_timeline(times, averaging),#[:-1],
            X_trendy_total
            #gh.avg_timeline(X_trendy_total, averaging)
            * 148940000
            * 1000000
            * 0.000000000001,  # convert to global C in Gt            
            label=model_names[mf] + " - X-trendy",
            color=add_cols[mf],
            #linestyle="dashed",
        )
        ax.plot(
            gh.avg_timeline(times, averaging),#[:-1],
            gh.avg_timeline(vals.x_c, averaging)#[:-1]
            * 148940000
            * 1000000
            * 0.000000000001,  # convert to global C in Gt
            label=model_names[mf] + " - X_c-matrix",
            color=model_cols[mf],
            linestyle="dashed",
        )
        ax.plot(
            gh.avg_timeline(times, averaging),#[:-1],            
            X_c_trendy_total
            #gh.avg_timeline(X_c_trendy_total, averaging)
            * 148940000
            * 1000000
            * 0.000000000001,  # convert to global C in Gt
            label=model_names[mf] + " - X_c-trendy",
            color=add_cols[mf],
            linestyle="dashed",
        )
        
        k += 1
    ax.legend()
    ax.set_title("Total Carbon (X) and Carbon Storage Capacity (X_c)")
    ax.set_ylabel("Gt C")
    ax.grid()
    plt.ylim([-1000, 4000])


# +
def plot_u_trendy(
    model_names,  # dictionary (folder name : model name)
    test_arg_list,  # a list of test_args from all models involved
    delta_t_val,  # model time step
    model_cols,  # dictionary (folder name :color)
    add_cols,
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging,  # number of iterator steps over which to average results. 1 for no averaging
    overlap=True,  # compute overlapping timeframe or plot whole duration for all models
):
    if (part < 0) | (part > 1):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        )
    model_folders = [(k) for k in model_names]
    fig = plt.figure(figsize=(17, 8))
    ax = fig.subplots(1, 1)
    k = 0
    for mf in model_folders:
        itr = gh.traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        if overlap == True:
            start_min, stop_max = gh.min_max_index(
                test_arg_list[k],
                delta_t_val,
                *gh.t_min_tmax_overlap(test_arg_list, delta_t_val)
            )
        else:
            start_min, stop_max = gh.min_max_index(
                test_arg_list[k],
                delta_t_val,
                *gh.t_min_tmax_full(test_arg_list, delta_t_val)
            )
        # if we do not want the whole interval but look at a smaller part to observe the dynamics
        start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        times = (
            gh.times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            / gh.days_per_year()
        )
        vals = itr[start:stop]
        # print("vals.x")
        # print(vals.x[0:240])
        
        ############## method from O'Sallivan et al. (2022)
        cVeg_trendy=test_arg_list[k].svs.cVeg
        if model_names[mf]=="VISIT":
            cSoil_trendy=test_arg_list[k].svs.cLitter+test_arg_list[k].svs.cSoil
        else:
            cSoil_trendy=test_arg_list[k].svs.cSoil
        rh_trendy=test_arg_list[k].svs.rh
        npp_trendy=test_arg_list[k].dvs.npp
        #npp_trendy_yearly=npp_sum(npp_trendy)
        
        delta_cVeg=np.zeros_like(cVeg_trendy)
        for i in range (len(cVeg_trendy)-1):
            delta_cVeg[i]=cVeg_trendy[i+1]-cVeg_trendy[i]
    
        delta_cSoil=np.zeros_like(cSoil_trendy)
        for i in range (len(cSoil_trendy)-1):
            delta_cSoil[i]=cSoil_trendy[i+1]-cSoil_trendy[i]    
    
        tau_v=cVeg_trendy/(npp_trendy-delta_cVeg)
        X_c_v_trendy=tau_v*npp_trendy
        #X_c_v_trendy_yearly=gh.avg_timeline(X_c_v_trendy,12)

        f_vs=delta_cSoil+rh_trendy

        tau_s=cSoil_trendy/rh_trendy
        X_c_s_trendy=f_vs*tau_s
        #X_c_s_trendy_yearly=gh.avg_timeline(X_c_s_trendy,12)

        X_trendy_total=(cVeg_trendy+cSoil_trendy)[start//averaging//delta_t_val:stop//averaging//delta_t_val]#*148940000*1000000*0.000000000001
        X_c_trendy_total=(X_c_v_trendy+X_c_s_trendy)[start//averaging//delta_t_val:stop//averaging//delta_t_val]#*148940000*1000000*0.000000000001 
        u_trendy=npp_trendy[start//averaging//delta_t_val:stop//averaging//delta_t_val]
#         print(str(mf)+": "+str(X_c_trendy_total.shape)+" ; "+str(X_trendy_total.shape))
#         print ("start: "+str(start)+" stop: "+str(stop))
#         print (npp_trendy)
#         print (vals.u)
        ###################
        ax.plot(
            gh.avg_timeline(times, averaging)[:-1],
            gh.avg_timeline(vals.u, averaging)[:-1]
            * 148940000
            * 1000000
            * 0.000000000001,  # convert to global C in Gt
            label=model_names[mf] + " - npp-matrix",
            color=model_cols[mf],            
        )
        ax.plot(            
            gh.avg_timeline(times, averaging)[:-1],
            u_trendy
            #gh.avg_timeline(u_trendy, averaging)
            * 148940000
            * 1000000
            * 0.000000000001,  # convert to global C in Gt            
            label=model_names[mf] + " - npp-trendy",
            color=add_cols[mf],
        )
    ax.legend()
    ax.set_title("Carbon Input (NPP)")
    ax.set_ylabel("Gt C / day")
    ax.grid()
    
#     print(gh.avg_timeline(vals.u, averaging)
#           * 148940000
#           * 1000000
#           * 0.000000000001,
#          )
#     print(gh.avg_timeline(u_trendy, averaging)
#           * 148940000
#           * 1000000
#           * 0.000000000001,
#          )    


# -

def plot_rt_trendy(
    model_names,  # dictionary (folder name : model name)
    test_arg_list,  # a list of test_args from all models involved
    delta_t_val,  # model time step
    model_cols,  # dictionary (folder name :color)
    add_cols,
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging,  # number of iterator steps over which to average results. 1 for no averaging
    overlap=True,  # compute overlapping timeframe or plot whole duration for all models
):
    if (part < 0) | (part > 1):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        )
    model_folders = [(k) for k in model_names]
    fig = plt.figure(figsize=(17, 8))
    ax = fig.subplots(1, 1)
    k = 0
    for mf in model_folders:
        itr = gh.traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        if overlap == True:
            start_min, stop_max = gh.min_max_index(
                test_arg_list[k],
                delta_t_val,
                *gh.t_min_tmax_overlap(test_arg_list, delta_t_val)
            )
        else:
            start_min, stop_max = gh.min_max_index(
                test_arg_list[k],
                delta_t_val,
                *gh.t_min_tmax_full(test_arg_list, delta_t_val)
            )
        # if we do not want the whole interval but look at a smaller part to observe the dynamics
        start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        times = (
            gh.times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            / gh.days_per_year()
        )
        vals = itr[start:stop]
        # print("vals.x")
        # print(vals.x[0:240])
        
        ############## method from O'Sallivan et al. (2022)
        cVeg_trendy=test_arg_list[k].svs.cVeg
        if model_names[mf]=="VISIT":
            cSoil_trendy=test_arg_list[k].svs.cLitter+test_arg_list[k].svs.cSoil
        else:
            cSoil_trendy=test_arg_list[k].svs.cSoil
        rh_trendy=test_arg_list[k].svs.rh
        npp_trendy=test_arg_list[k].dvs.npp
        #npp_trendy_yearly=npp_sum(npp_trendy)
        
        delta_cVeg=np.zeros_like(cVeg_trendy)
        for i in range (len(cVeg_trendy)-1):
            delta_cVeg[i]=cVeg_trendy[i+1]-cVeg_trendy[i]
    
        delta_cSoil=np.zeros_like(cSoil_trendy)
        for i in range (len(cSoil_trendy)-1):
            delta_cSoil[i]=cSoil_trendy[i+1]-cSoil_trendy[i]    
    
        tau_v=cVeg_trendy/(npp_trendy-delta_cVeg)
        X_c_v_trendy=tau_v*npp_trendy
        #X_c_v_trendy_yearly=gh.avg_timeline(X_c_v_trendy,12)

        f_vs=delta_cSoil+rh_trendy

        tau_s=cSoil_trendy/rh_trendy
        X_c_s_trendy=f_vs*tau_s
        #X_c_s_trendy_yearly=gh.avg_timeline(X_c_s_trendy,12)

        X_trendy_total=(cVeg_trendy+cSoil_trendy)[start//averaging//delta_t_val:stop//averaging//delta_t_val]#*148940000*1000000*0.000000000001
        X_c_trendy_total=(X_c_v_trendy+X_c_s_trendy)[start//averaging//delta_t_val:stop//averaging//delta_t_val]#*148940000*1000000*0.000000000001 
        u_trendy=npp_trendy[start//averaging//delta_t_val:stop//averaging//delta_t_val]
        tau_trendy=(tau_v+tau_s)[start//averaging//delta_t_val:stop//averaging//delta_t_val]
#         print(str(mf)+": "+str(X_c_trendy_total.shape)+" ; "+str(X_trendy_total.shape))
#         print ("start: "+str(start)+" stop: "+str(stop))
#         print (tau_trendy)
#         print (vals.rt)
        ###################
        ax.plot(
            gh.avg_timeline(times, averaging)[:-1],
            gh.avg_timeline(vals.rt, averaging)[:-1],  # convert to global C in Gt
            label=model_names[mf] + " - rt-matrix",
            color=model_cols[mf],            
        )
        ax.plot(            
            gh.avg_timeline(times, averaging)[:-1],
            tau_trendy,
            #gh.avg_timeline(u_trendy, averaging),  # convert to global C in Gt            
            label=model_names[mf] + " - tau-trendy",
            color=add_cols[mf],
        )
    ax.legend()
    ax.set_title("Residense time / turnover time")
    ax.set_ylabel("years")
    ax.grid()



plot_x_xc_trendy(model_names=model_names_test,
                     delta_t_val=delta_t_val,
                     test_arg_list=test_arg_list_test,
                     model_cols=model_cols,
                     add_cols=add_cols,
                     part=1,
                     averaging=12*30//delta_t_val,
                     overlap=True
                     )

plot_u_trendy(model_names=model_names_test,
                     delta_t_val=delta_t_val,
                     test_arg_list=test_arg_list_test,
                     model_cols=model_cols,
                     add_cols=add_cols,
                     part=1,
                     averaging=12*30//delta_t_val,
                     overlap=True
                     )

plot_rt_trendy(model_names=model_names_test,
                     delta_t_val=delta_t_val,
                     test_arg_list=test_arg_list_test,
                     model_cols=model_cols,
                     add_cols=add_cols,
                     part=1,
                     averaging=12*30//delta_t_val,
                     overlap=True
                     )


