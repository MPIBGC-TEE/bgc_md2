# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
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
    func_dict=gh.msh(mf).make_func_dict(dvs_t,cpa_t,epa_t)
    
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


# +
def times_in_days_aD(mf,delta_t_val):
    n_months=len(gh.test_args(mf).dvs[0])
    #n_days=n_months*30
    n_days=gh.month_2_day_index([n_months],gh.msh(mf).start_date())[0]
    n_iter=int(n_days/delta_t_val)
    #days_after_sim_start=delta_t_val*np.arange(n_iter)
    #print(n_months,n_days,days_after_sim_start)
    return np.array(
        tuple((
            gh.days_since_AD(i,delta_t_val,gh.msh(mf).start_date()) 
            for i in np.arange(n_iter)#days_after_sim_start
        ))
    )

delta_t_val=30 # assuming the same step size for every model (could be made model specific by an additional testarg)
times_in_days_aD(mf,delta_t_val)


# +
def t_min_tmax(model_folders,delta_t_val):
    td={
        mf: times_in_days_aD(mf,delta_t_val)
        for mf in model_folders
    }
    t_min = max([t.min() for t in td.values()])
    t_max = min([t.max() for t in td.values()])
    return (t_min,t_max)

#t_min_tmax(model_folders,delta_t_val)
t_min_tmax(model_folders[0:1],delta_t_val)
# -

#gh.month_2_day_index([13])


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

s=slice(
    *min_max_index(
            "yz_jules",
             delta_t_val,
            *t_min_tmax(model_folders[0:1],delta_t_val)
    )
)
s
#ind_d={mf: min_max_index(mf) for mf in model_folders}
# -


# You can use python index notation on the iterator
# it works the same way as  with a list
l=[1,2,3,4,5,6,7,8,9]
l[3:8:2]

itr=tracebility_iterator(mf,delta_t_val)
res=itr[3:8:2]
# take a peek
#type(res), res._fields, res.X_c.shape, res.X_dot

#mf='kv_visit2'
mf="yz_jules"
start,stop=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders[0:1],delta_t_val))
start,stop

#for better readability we convert the times to (approximate) years 
times=times_in_days_aD(mf,delta_t_val)[start:stop]/gh.days_per_year()
vals=itr[start:stop]
vals.X_c.shape,times.shape 
#vals

# +
import matplotlib.pyplot as plt
def plot_sums(model_folders,delta_t_val,cd):
    n = len(model_folders)
    fig=plt.figure(figsize=((n+1)*10,20))
    axs=fig.subplots(1,n)
    plt.rcParams['font.size'] = 20
    names=['X','X_c','X_p']
    for i,mf in enumerate(model_folders):
        itr=tracebility_iterator(mf,delta_t_val)
        start,stop=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
        #tups=values(itr,start,stop)
        #vals=values_2_TraceTuple(tups)
        vals=itr[start:stop]
        times=times_in_days_aD(mf,delta_t_val)[start:stop]/gh.days_per_year()
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
    
plot_sums(model_folders[0:2],delta_t_val, cd)
