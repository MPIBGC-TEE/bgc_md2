# ---
# jupyter:
#   jupytext:
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
# #%load_ext autoreload
# #%autoreload 2
import json
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
# define some functions that yield the result depending on the model folder  mf
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


# +
def values(itr,start,stop,increment=1):
    from copy import copy
    itr=copy(itr)
    # run the iterator for start iterations
    for i in range(start):
        itr.__next__()
    
    # now collect the desired
    return tuple(
        (
            itr.__next__() 
            for i in range(stop-start)
            if i%increment==0
        )
    )

#check with a list
l=[1,2,3,4]
start,stop=1,3
print(l[start:stop])
itr=l.__iter__()
values(itr,start,stop)

# -

# now with 
itr=tracebility_iterator(mf,delta_t_val)
#values(itr,3,4)

# +
def values_2_TraceTuple(tups):
    # instead of the collection of TraceTuples that the iterator returns
    # we want a Tracetuple of arrays whith the added time dimension
    return gh.TraceTuple(*(
        np.stack(
            tuple((tup.__getattribute__(name)  for tup in tups))
        )
        for name in gh.TraceTuple._fields
    ))
    
#number_of_iterations=100
mf='kv_visit2'
mf="yz_jules"
start,stop=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
start,stop

# +
tups=values(itr,start,stop)
#tt=tups[0]
#tt.__getattribute__("X")

vals=values_2_TraceTuple(tups)
times=times_in_days_aD(mf,delta_t_val)[start:stop]/365
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
        tups=values(itr,start,stop)
        vals=values_2_TraceTuple(tups)
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
    
plot_sums(model_folders,delta_t_val)


# -

# check how to sum a TraceTuple of arrays 
def tt_sum(tt):
    return gh.TraceTuple( *(
            tt.__getattribute__(name).sum(axis=0)
            for name in tt._fields 
        )
    )
#tt_sum(vals)


# +
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

#partitions(start,stop,12)


# -

#temporaly averaged values
def averaged_values(itr,partitions):
    start=partitions[0][0]
    def tt_avg(tups):
        l=len(tups)
        return gh.TraceTuple(*(
                np.stack(
                    [
                        tup.__getattribute__(name) 
                        for tup in tups
                    ],
                    axis=0
                ).sum(axis=0)/l
                for name in gh.TraceTuple._fields
            )
        )
                             
    # move to the start
    for i in range(start):
        itr.__next__()
        
    tts=[
        tt_avg(
            [ 
                itr.__next__()
                for i in range(stop_p-start_p)
            ]
        )
        for (start_p,stop_p) in partitions
    ]
    return values_2_TraceTuple(tts)



itr=tracebility_iterator(mf,delta_t_val)
start,stop=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
tavg_vals=averaged_values(
    itr,
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
def plot_avg_sums(model_folders,delta_t_val):
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
        parts=partitions(start,stop,12)
        times=averaged_times(
            times_in_days_aD(mf,delta_t_val)/365,
            parts
        )
        vals=averaged_values(
            itr,
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
        
plot_avg_sums(model_folders,delta_t_val)
# -




