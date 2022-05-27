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

# +
# %load_ext autoreload
# %autoreload 2
import json
import matplotlib.pyplot as plt
from typing import Tuple
from importlib import import_module
from pathlib import Path
from frozendict import frozendict
import numpy as np
from functools import lru_cache
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    StateVariableTuple
)
import general_helpers as gh


# define some shortcut functions that yield the result depending on the 
# model folder  mf
def sim_day_2_day_aD_func(mf): #->function
    return gh.msh(mf).make_sim_day_2_day_since_a_D(gh.confDict(mf))

def tracebility_iterator(mf,delta_t_val):
    ta=gh.test_args(mf)
    mvs_t=gh.mvs(mf)
    dvs_t=ta.dvs
    cpa_t=ta.cpa
    epa_t=ta.epa_opt
    X_0=gh.msh(mf).numeric_X_0(mvs_t,dvs_t,cpa_t,epa_t)
    func_dict=gh.msh(mf).make_func_dict(mvs_t,dvs_t,cpa_t,epa_t)
    
    return gh.traceability_iterator(
        X_0,
        func_dict,
        mvs=mvs_t,
        dvs=dvs_t,
        cpa=cpa_t,
        epa=epa_t,
        delta_t_val=delta_t_val
    )
    
model_folders=['yz_jules','kv_visit2']#, 'Aneesh_SDGVM']#, 'kv_ft_dlem', 'jon_yib']#,'Aneesh_SDGVM','cable-pop','yz_jules']#,]
mf=model_folders[0]
cd={
        'X':"red",
        'X_c':"orange",
        'X_p':"blue",
        'X_dot':"red",
        'x':"red",
        'x_c':"orange",
        'x_p':"blue",
        'x_dot':"red",
        'I':"green",
        'u':"green",
}

delta_t_val=30 # assuming the same step size for every model (could be made model specific by an additional testarg)

# +
##mf='kv_visit2'
#mf="yz_jules"
#start,stop=gh.min_max_index(mf,delta_t_val,*gh.t_min_tmax_overlap(model_folders,delta_t_val))
#start,stop
# -

# ### Towards attribution of uncertainty
# One of the attractive promises of traceability analysis is the hierachie of attribution to successively finer degree. 
# $$
# \mathbf{X_c}=\mathbf{RT} u
# $$
# $$
# \Delta X_c \approx 
# \left.
# \frac{\partial X_c}{\partial RT} 
# \right|_{RT_1,u_1}
# \Delta RT 
# +
# \left.
# \frac{\partial X_c}{\partial u} 
# \right|_{RT_1,u_1}
# \Delta u
# $$
#
# $$
# =u_1  (RT_2-RT_1)
# + RT_1 (u_2-u_1) 
# $$
#

# +
model_cols={
    "yz_jules": "blue",
    "kv_visit2": "orange",
}
from scipy.interpolate import interp1d, splprep

def plot_diff(mf_1, mf_2, delta_t_val, model_cols):
    
    part=30
    start_min_1,stop_max_1=gh.min_max_index(mf_1,delta_t_val,*gh.t_min_tmax_overlap([mf_1,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_1,stop_1 = start_min_1, int(start_min_1+(stop_max_1-start_min_1)/part)
    itr_1=tracebility_iterator(mf_1,delta_t_val)
    vals_1=itr_1[start_1:stop_1]
    times_1=gh.times_in_days_aD(mf_1,delta_t_val)[start_1:stop_1]/365
    
    start_min_2,stop_max_2=gh.min_max_index(mf_2,delta_t_val,*gh.t_min_tmax_overlap([mf_2,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_2,stop_2 = start_min_2, int(start_min_2+(stop_max_2-start_min_2)/part)
    itr_2=tracebility_iterator(mf_2,delta_t_val)
    vals_2=itr_2[start_2:stop_2]
    times_2=gh.times_in_days_aD(mf_2,delta_t_val)[start_2:stop_2]/365
    
    # Since the two models do not necessarily share the same point in time and not
    # even the same stepsize or number of steps we compute interpolating functions
    # to make them comparable
    # create interpolation functions for common timeline
    x_c_1=interp1d(times_1,vals_1.x_c)
    x_c_2=interp1d(times_2,vals_2.x_c)
    u_1=interp1d(times_1,vals_1.u)
    u_2=interp1d(times_2,vals_2.u)
    rt_1=interp1d(times_1,vals_1.rt)
    rt_2=interp1d(times_2,vals_2.rt)
    
    # common plot times
    start=max(times_1.min(),times_2.min())
    stop=min(times_1.max(),times_2.max())
    nstep=min(len(times_1),len(times_2))
    times=np.linspace(start,stop,nstep)
    
    # values for plots
    delta_u=u_1(times)-u_2(times)
    delta_rt=rt_1(times)-rt_2(times)
    delta_x_c=x_c_1(times)-x_c_2(times)
    
    
    fig=plt.figure(figsize=(2*10,20))
    axs=fig.subplots(3,3)
    ###################################################
    # plot x_c, u and rt for both models 
    ###################################################
    def subp(ax,name):
        def subsubp(mf,vals,times):
            ax.plot(
                times,
                vals.__getattribute__(name),
                color=model_cols[mf],
                label=mf
            )
            
        subsubp( mf_1, vals_1,  times_1 )
        subsubp( mf_2, vals_2,  times_2 )
        ax.legend()
        ax.set_title(name)

        
    subp(axs[0,0],"x_c")
    subp(axs[0,1],"u")
    subp(axs[0,2],"rt")
    #
    ####################################################
    ## plot delta_x_c delta_u and delta_rt  
    ####################################################
    ax=axs[1,0]
    ax.plot(times,delta_x_c,color="green")
    ax.set_title("$\Delta_{X_c}$")
    
    ax=axs[1,1]
    ax.plot(times,delta_u,color="green")
    ax.set_title("$\Delta_{u}$")
    
    ax=axs[1,2]
    ax.plot(times,delta_rt,color="green")
    ax.set_title("$\Delta_{rt}$")
    
    ####################################################
    ## plot x_c_1, x_c_2,  RT_1*delta_u and delta_rt  
    ####################################################
    ## rt1*delta_u
    
    ax=axs[2,1]
    ax.plot(times,delta_u*rt_1(times),color="black")
    ax.plot(times,delta_rt*u_1(times),color="red")
    ax.plot(times,delta_u*rt_1(times)+delta_rt*u_1(times),color="blue")
    ax.plot(times,delta_x_c,color="green")
    ax.set_title("$Approximation$")
    
    
mf_1="yz_jules"
mf_2="kv_visit2"
plot_diff(mf_1, mf_2,delta_t_val=10,model_cols=model_cols)
# +
model_cols={
    "yz_jules": "blue",
    "kv_visit2": "orange",
}
from scipy.interpolate import interp1d, splprep

def plot_diff_1(mf_1, mf_2, delta_t_val, model_cols):
    
    part=30
    start_min_1,stop_max_1=gh.min_max_index(mf_1,delta_t_val,*gh.t_min_tmax_overlap([mf_1,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_1,stop_1 = start_min_1, int(start_min_1+(stop_max_1-start_min_1)/part)
    itr_1=tracebility_iterator(mf_1,delta_t_val)
    vals_1=itr_1[start_1:stop_1]
    times_1=gh.times_in_days_aD(mf_1,delta_t_val)[start_1:stop_1]/365
    
    start_min_2,stop_max_2=gh.min_max_index(mf_2,delta_t_val,*gh.t_min_tmax_overlap([mf_2,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_2,stop_2 = start_min_2, int(start_min_2+(stop_max_2-start_min_2)/part)
    itr_2=tracebility_iterator(mf_2,delta_t_val)
    vals_2=itr_2[start_2:stop_2]
    times_2=gh.times_in_days_aD(mf_2,delta_t_val)[start_2:stop_2]/365
    
    # Since the two models do not necessarily share the same point in time and not
    # even the same stepsize or number of steps we compute interpolating functions
    # to make them comparable
    # create interpolation functions for common timeline
    x_c_1=interp1d(times_1,vals_1.x_c)
    x_c_2=interp1d(times_2,vals_2.x_c)
    u_1=interp1d(times_1,vals_1.u)
    u_2=interp1d(times_2,vals_2.u)
    rt_1=interp1d(times_1,vals_1.rt)
    rt_2=interp1d(times_2,vals_2.rt)
    
    # common plot times
    start=max(times_1.min(),times_2.min())
    stop=min(times_1.max(),times_2.max())
    nstep=min(len(times_1),len(times_2))
    times=np.linspace(start,stop,nstep)
    
    # values for plots
    delta_u=u_1(times)-u_2(times)
    delta_rt=rt_1(times)-rt_2(times)
    delta_x_c=x_c_1(times)-x_c_2(times)
    
    
    fig=plt.figure(figsize=(10,10))
    axs=fig.subplots(1,1)
    ####################################################
    ## plot x_c_1, x_c_2,  RT_1*delta_u and delta_rt  
    ####################################################
    ## rt1*delta_u
    
    ax=axs#[2,1]
    ax.plot(times,delta_u*rt_1(times),label="$\Delta u * rt_1$",color="black")
    ax.plot(times,delta_rt*u_1(times),label="$\Delta rt * u_1$",color="red")
    ax.plot(
        times,
        delta_u*rt_1(times)+delta_rt*u_1(times),
        label="Approximate $\Delta  X_c= \Delta u * rt_1 + \Delta rt * u_1$",
        color="blue"
    )
    ax.plot(times,delta_x_c,label="Exact $\Delta X_c= X_1 - X_2$",color="green")
    ax.set_title("$Approximation$")
    ax.legend()
    
    
mf_1="yz_jules"
mf_2="kv_visit2"
plot_diff_1(mf_1, mf_2,delta_t_val=10,model_cols=model_cols)
# -

# We can see both contributions to $\Delta x_c$. They are of the same order of magnitude although $\Delta rt$ seems a tiny bit bigger.
# We also see that the contributions are sometimes higher for $\Delta u$.
# The error between the real and the linearly approximated $\Delta xc$ is also visible.
# Not too bad...
#


# +
model_cols={
    "yz_jules": "blue",
    "kv_visit2": "orange",
}
from scipy.interpolate import interp1d, splprep

def plot_diff_2(mf_1, mf_2, delta_t_val, model_cols):
    
    part=30
    start_min_1,stop_max_1=gh.min_max_index(mf_1,delta_t_val,*gh.t_min_tmax_overlap([mf_1,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_1,stop_1 = start_min_1, int(start_min_1+(stop_max_1-start_min_1)/part)
    itr_1=tracebility_iterator(mf_1,delta_t_val)
    vals_1=itr_1[start_1:stop_1]
    times_1=gh.times_in_days_aD(mf_1,delta_t_val)[start_1:stop_1]/365
    
    start_min_2,stop_max_2=gh.min_max_index(mf_2,delta_t_val,*gh.t_min_tmax_overlap([mf_2,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_2,stop_2 = start_min_2, int(start_min_2+(stop_max_2-start_min_2)/part)
    itr_2=tracebility_iterator(mf_2,delta_t_val)
    vals_2=itr_2[start_2:stop_2]
    times_2=gh.times_in_days_aD(mf_2,delta_t_val)[start_2:stop_2]/365
    
    # Since the two models do not necessarily share the same point in time and not
    # even the same stepsize or number of steps we compute interpolating functions
    # to make them comparable
    # create interpolation functions for common timeline
    x_c_1=interp1d(times_1,vals_1.x_c)
    x_c_2=interp1d(times_2,vals_2.x_c)
    u_1=interp1d(times_1,vals_1.u)
    u_2=interp1d(times_2,vals_2.u)
    rt_1=interp1d(times_1,vals_1.rt)
    rt_2=interp1d(times_2,vals_2.rt)
    
    # common plot times
    start=max(times_1.min(),times_2.min())
    stop=min(times_1.max(),times_2.max())
    nstep=min(len(times_1),len(times_2))
    times=np.linspace(start,stop,nstep)
    
    # values for plots
    delta_u=u_1(times)-u_2(times)
    delta_rt=rt_1(times)-rt_2(times)
    delta_x_c=x_c_1(times)-x_c_2(times)
    
    
    fig=plt.figure(figsize=(10,10))
    axs=fig.subplots(1,1)
    ####################################################
    ## plot x_c_1, x_c_2,  RT_1*delta_u and delta_rt  
    ####################################################
    ## rt1*delta_u
    
    ax=axs#[2,1]
    #ax.plot(times,delta_u*rt_1(times),label="$\Delta u * rt_1$",color="black")
    #ax.plot(times,delta_rt*u_1(times),label="$\Delta rt * u_1$",color="red")
    ax.stackplot(
        times,
        [
            delta_rt*u_1(times),
            delta_u*rt_1(times)
        ],
        labels=[
            "$\Delta rt * u_1$",
            "$\Delta u * rt_1$"
        ],
        colors=[
            "red","black"
        ]
    )
    ax.plot(
        times,
        delta_u*rt_1(times)+delta_rt*u_1(times),
        label="Approximate $\Delta  X_c= \Delta u * rt_1 + \Delta rt * u_1$",
        color="blue"
    )
    ax.plot(times,delta_x_c,label="Exact $\Delta X_c= X_1 - X_2$",color="green")
    ax.set_title("$Approximation$")
    ax.legend()
    
    
mf_1="yz_jules"
mf_2="kv_visit2"
plot_diff_2(mf_1, mf_2,delta_t_val=10,model_cols=model_cols)

# +
model_cols={
    "yz_jules": "blue",
    "kv_visit2": "orange",
}
from scipy.interpolate import interp1d, splprep

def plot_diff_3(mf_1, mf_2, delta_t_val, model_cols):
    
    part=30
    start_min_1,stop_max_1=gh.min_max_index(mf_1,delta_t_val,*gh.t_min_tmax_overlap([mf_1,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_1,stop_1 = start_min_1, int(start_min_1+(stop_max_1-start_min_1)/part)
    itr_1=tracebility_iterator(mf_1,delta_t_val)
    vals_1=itr_1[start_1:stop_1]
    times_1=gh.times_in_days_aD(mf_1,delta_t_val)[start_1:stop_1]/365
    
    start_min_2,stop_max_2=gh.min_max_index(mf_2,delta_t_val,*gh.t_min_tmax_overlap([mf_2,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_2,stop_2 = start_min_2, int(start_min_2+(stop_max_2-start_min_2)/part)
    itr_2=tracebility_iterator(mf_2,delta_t_val)
    vals_2=itr_2[start_2:stop_2]
    times_2=gh.times_in_days_aD(mf_2,delta_t_val)[start_2:stop_2]/365
    
    # Since the two models do not necessarily share the same point in time and not
    # even the same stepsize or number of steps we compute interpolating functions
    # to make them comparable
    # create interpolation functions for common timeline
    x_c_1=interp1d(times_1,vals_1.x_c)
    x_c_2=interp1d(times_2,vals_2.x_c)
    u_1=interp1d(times_1,vals_1.u)
    u_2=interp1d(times_2,vals_2.u)
    rt_1=interp1d(times_1,vals_1.rt)
    rt_2=interp1d(times_2,vals_2.rt)
    
    # common plot times
    start=max(times_1.min(),times_2.min())
    stop=min(times_1.max(),times_2.max())
    nstep=min(len(times_1),len(times_2))
    times=np.linspace(start,stop,nstep)
    
    # values for plots
    delta_u=u_1(times)-u_2(times)
    delta_rt=rt_1(times)-rt_2(times)
    delta_x_c=x_c_1(times)-x_c_2(times)
    
    
    fig=plt.figure(figsize=(10,10))
    axs=fig.subplots(1,1)
    ####################################################
    ## plot x_c_1, x_c_2,  RT_1*delta_u and delta_rt  
    ####################################################
    ## rt1*delta_u
    
    ax=axs#[2,1]
    ax.plot(times,x_c_1(times),label="$x_c$ {}".format(mf_1),color=model_cols[mf_1])
    ax.stackplot(
        times,
        [
            x_c_2(times),
            delta_rt*u_1(times),
            delta_u*rt_1(times),
        ],
        labels=[
            "$x_c$ {}".format(mf_2),
            "$\Delta rt * u_1$",
            "$\Delta u * rt_1$",
        ],
        colors=[
            model_cols[mf_2],
            "red",
            "black",
        ]
    )
    #ax.plot(
    #    times,
    #    delta_u*rt_1(times)+delta_rt*u_1(times),
    #    label="Approximate $\Delta  X_c= \Delta u * rt_1 + \Delta rt * u_1$",
    #    color="blue"
    #)
    #ax.plot(times,delta_x_c,label="Exact $\Delta X_c= X_1 - X_2$",color="green")
    ax.set_title("$Approximation$")
    ax.legend()
    
    
mf_1="yz_jules"
mf_2="kv_visit2"
plot_diff_3(mf_1, mf_2,delta_t_val=10,model_cols=model_cols)
# -


