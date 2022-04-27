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

# ## Model comparison using traceability analysis
# We use the infrastructure built so far to compare two ore more models.
# The workhorse will be an iterator (returned by a general function). That allows us to compute and easily access the desired timelines with python index notation it[2:5] will return the values from position 2 to 5 of the solution (and desired variables).
# The notebook also contains some functions to compute where the times of two models overlap and some plot functions. 

from IPython.display import Markdown, display
display(Markdown("TracebilityText.md"))

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
    epa_t=ta.epa_0
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
    
model_folders=['yz_jules','kv_visit2']#, 'Aneesh_SDGVM', 'kv_ft_dlem', 'jon_yib']#,'Aneesh_SDGVM','cable-pop','yz_jules']#,]
mf=model_folders[0]


# +
def times_in_days_aD(mf,delta_t_val):
    n_months=len(gh.test_args(mf).dvs[0])
    n_days=n_months*30
    n_iter=int(n_days/delta_t_val)
    days_after_sim_start=delta_t_val*np.arange(n_iter)
    return np.array(tuple(map(sim_day_2_day_aD_func(mf),days_after_sim_start)))

delta_t_val=30 # assuming the same step size for every model (could be made model specific by an additional testarg)


# +
def t_min_tmax(model_folders,delta_t_val):
    td={
        mf: times_in_days_aD(mf,delta_t_val)
        for mf in model_folders
    }
    t_min = max([t.min() for t in td.values()])
    t_max = min([t.max() for t in td.values()])
    return (t_min,t_max)

t_min_tmax(model_folders,delta_t_val)

# +
#find the timesteps corresponding to shared times
from functools import reduce
def min_max_index(mf,delta_t_val,t_min,t_max):
    ts=times_in_days_aD(mf,delta_t_val)
    def count(acc,i):
        min_i,max_i = acc
        t=ts[i]
        min_i = min_i+1 if t < t_min else min_i 
        max_i = max_i+1 if t < t_max else max_i 
        return (min_i,max_i)
    
    return reduce(count,range(len(ts)),(0,0)) 

s=slice(*min_max_index("yz_jules",delta_t_val,*t_min_tmax(model_folders,delta_t_val)))
s.step is None 

#ind_d={mf: min_max_index(mf) for mf in model_folders}
# -


# You can use python index notation on the iterator
# it works the same way as  with a list
l=[1,2,3,4,5,6,7,8,9]
l[3:8:2]

itr=tracebility_iterator(mf,delta_t_val)
res=itr[3:8:2]
# take a peek
type(res), res._fields, res.X_c.shape, res.X_dot

#mf='kv_visit2'
mf="yz_jules"
start,stop=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
start,stop

times=times_in_days_aD(mf,delta_t_val)[start:stop]/365
vals=itr[start:stop]
vals.X_c.shape,times.shape 
#vals

# +
import matplotlib.pyplot as plt
def plot_sums(model_folders,delta_t_val):
    n = len(model_folders)
    fig=plt.figure(figsize=((n+1)*10,20), dpi = 400)
    axs=fig.subplots(1,n)
    plt.rcParams['font.size'] = 20
    names=['X','X_c','X_p']
    cd={
            'X':"green",
            'X_c':"orange",
            'X_p':"blue"
    }
    for i,mf in enumerate(model_folders):
        itr=tracebility_iterator(mf,delta_t_val)
        start,stop=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
        #tups=values(itr,start,stop)
        #vals=values_2_TraceTuple(tups)
        vals=itr[start:stop]
        times=times_in_days_aD(mf,delta_t_val)[start:stop]/365
        ax=axs[i]
        for name in names:
            ax.plot(
                times,
                vals.__getattribute__(name).sum(axis=1),
                label=name+"_sum",
                color=cd[name]
            )
        ax.legend()
        ax.set_title(mf)
    #plt.close(fig)
    
#plot_sums(model_folders,delta_t_val)


# -

# ## temporal averages
# In some cases you might want to see a yearly or monthly average.
# You could compute all the values and then compute the averages of parts of the arrays.
# The iterator can also do it for you (without storing the unnecessary potentially huge intermediate fine grained arrays). The argument for the averaging will be a list of (start,stop) tuples describing the ranges over which to average. 

# +
# we can build this partition by a little function 
def partitions(start,stop,nr_acc=1):
    diff=stop-start
    step=nr_acc
    number_of_steps=int(diff/step)
    last_start=start+number_of_steps*step
    last_tup=(last_start,stop)
    return [
        (
            start + step * i,
            start + step *(i+1)
        )
        for i in range(number_of_steps)
    ]+[last_tup]

# an example we want to partition 
partitions(0,10,nr_acc=3)
#len(partitions(start,stop,12))
# -

itr=tracebility_iterator(mf,delta_t_val)
start,stop=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
tavg_vals=itr.averaged_values(
    partitions(start,stop,12)
)

tavg_vals.X.shape


# +
def averaged_times(times,partitions):
    return np.array(
        [
            times[p[0]:p[1]].sum()/(p[1]-p[0]) for p in partitions
        ]
    )
    
averaged_times(
    times_in_days_aD(mf,delta_t_val),
    partitions(start,stop,12)
).shape


# +
def plot_yearly_avg_sums(model_folders,delta_t_val):
    n = len(model_folders)
    fig=plt.figure(figsize=((n+1)*10,20))#, dpi = 400)
    axs=fig.subplots(1,n)
    plt.rcParams['font.size'] = 20
    names=['X','X_c','X_p']
    cd={
            'X':"green",
            'X_c':"orange",
            'X_p':"blue"
    }
    for i,mf in enumerate(model_folders):
        itr=tracebility_iterator(mf,delta_t_val)
        start,stop=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
        parts=partitions(start,stop,12)
        times=averaged_times(
            times_in_days_aD(mf,delta_t_val)/365,
            parts
        )
        vals=itr.averaged_values(
            parts
        )
        ax=axs[i]
        for name in names:
            ax.plot(
                times,
                vals.__getattribute__(name).sum(axis=1),
                label=name+"_avg_sum",
                color=cd[name]
            )
        ax.legend()
        ax.set_title(mf)
        
#plot_yearly_avg_sums(model_folders,delta_t_val)


# -
# ## Application: Numerical experiment concerning the  attraction of the solution to wards X_c
#
# We check the claim that the total Carbon $X$ as vector as well as the sum over all pools.
# To this end we plot $X_p$ and the derivative $\frac{d}{d t}X$ 

# +
def plot_sums(model_folders,delta_t_val):
    n = len(model_folders)
    fig=plt.figure(figsize=((n+1)*10,20))#, dpi = 400)
    axs=fig.subplots(2,n)
    plt.rcParams['font.size'] = 20
    mass_names=['X','X_c','X_p']
    names_2=['X_dot']
    cd={
            'X':"green",
            'X_c':"orange",
            'X_p':"blue",
            'X_dot':"black"
    }
    for i,mf in enumerate(model_folders):
        itr=tracebility_iterator(mf,delta_t_val)
        start_min,stop_max=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
        # we do not want the whole interval but look at a smaller part to observe the dynamics
        start,stop = start_min, int(start_min+(stop_max-start_min)/30)
        vals=itr[start:stop]
        times=times_in_days_aD(mf,delta_t_val)[start:stop]/365
        ax=axs[0,i]
        for name in mass_names:
            ax.plot(
                times,
                vals.__getattribute__(name).sum(axis=1),
                label=name+"_sum",
                color=cd[name]
            )
        ax.legend()
        ax.set_title(mf)
        ax=axs[1,i]
        for name in names_2:
            ax.plot(
                times,
                vals.__getattribute__(name).sum(axis=1),
                label=name+"_sum",
                color=cd[name]
            )
            ax.plot(
                times,
                np.zeros_like(times),
                #label=name+"_sum",
                color=cd[name]
            )
        ax.legend()
        ax.set_title(mf)
    #plt.close(fig)
    
#plot_sums(model_folders,delta_t_val)


# -


# ### We see that the sum of the derivative is clearly positive all the time even if X_c_sum crosses the X_sum lines

# +
def plot_components(mf,delta_t_val):
    names_1=['X','X_c']#,'X_p']
    names_2=['X_dot']
    cd={
            'X':"green",
            'X_c':"orange",
            'X_p':"blue",
            'X_dot':"black"
    }
    itr=tracebility_iterator(mf,delta_t_val)
    start_min,stop_max=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start,stop = start_min, int(start_min+(stop_max-start_min)/30)
    vals=itr[start:stop]
    n = vals.X.shape[1]
    fig=plt.figure(figsize=(20,(n+1)*10))#, dpi = 400)
    axs=fig.subplots(2*n,1)
    plt.rcParams['font.size'] = 20
    times=times_in_days_aD(mf,delta_t_val)[start:stop]/365
    i=0
    for j in range(n):
        ax=axs[2*j]#,i]
        for name in names_1:
            time_line= vals.__getattribute__(name)
            ax.plot(
                times,
                time_line[:,j,0],
                label=name+str(j),
                color=cd[name]
            )
            ax.legend()
            ax.set_title(mf)
        
        ax=axs[2*j+1]#,i]
        for name in names_2:
            time_line= vals.__getattribute__(name)
            ax.plot(
                times,
                time_line[:,j,0],
                label=name+str(j),
                color=cd[name]
            )
            ax.plot(
                times,
                np.zeros_like(times),
                color=cd[name]
            )
            ax.legend()
            ax.set_title(mf)
    
#plot_components(model_folders[0],delta_t_val)
plot_components(model_folders[1],delta_t_val)
# -

[ i for i in range(0,10,2)]



