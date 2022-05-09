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
# We use the infrastructure built so far to compare two or more models.
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
model_cols={
    "yz_jules": "blue",
    "kv_visit2": "orange",
}
from scipy.interpolate import interp1d, splprep

def plot_diff(mf_1, mf_2, delta_t_val, model_cols):
    n=2 #for 2 models
    part=30
    start_min_1,stop_max_1=min_max_index(mf_1,delta_t_val,*t_min_tmax([mf_1,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_1,stop_1 = start_min_1, int(start_min_1+(stop_max_1-start_min_1)/part)
    itr_1=tracebility_iterator(mf_1,delta_t_val)
    vals_1=itr_1[start_1:stop_1]
    times_1=times_in_days_aD(mf_1,delta_t_val)[start_1:stop_1]/365
    
    start_min_2,stop_max_2=min_max_index(mf_2,delta_t_val,*t_min_tmax([mf_2,mf_2],delta_t_val))
    # we do not want the whole interval but look at a smaller part to observe the dynamics
    start_2,stop_2 = start_min_2, int(start_min_2+(stop_max_2-start_min_2)/part)
    itr_2=tracebility_iterator(mf_2,delta_t_val)
    vals_2=itr_2[start_2:stop_2]
    times_2=times_in_days_aD(mf_2,delta_t_val)[start_2:stop_2]/365
    fig=plt.figure(figsize=((n)*10,20))
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
    # Since the two models do not necessarily share the same point in time and not
    # even the same stepsize or number of steps we compute interpolating functions
    # to make them comparable
    def diffp(ax,name):
        f1=interp1d(times_1,vals_1.__getattribute__(name))
        f2=interp1d(times_2,vals_2.__getattribute__(name))
        # chose the interval covered by both to avoid extrapolation
        start=max(times_1.min(),times_2.min())
        stop=min(times_1.max(),times_2.max())
        nstep=min(len(times_1),len(times_2))
        times=np.linspace(start,stop,nstep)
        
        diff=f1(times)-f2(times)
        ax.plot(times,diff,color="black")
        ax.set_title("{0}_{1}-{2}".format(name,mf_1,mf_2))
        
    diffp(axs[1,0],"x_c")
    diffp(axs[1,1],"u")
    diffp(axs[1,2],"rt")
    
mf_1="yz_jules"
mf_2="kv_visit2"
plot_diff(mf_1, mf_2,delta_t_val=10,model_cols=model_cols)
# -

sorted([[1,2],[5]],key=len)

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
    
plot_sums(model_folders,delta_t_val, cd)


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
#partitions(0,10,nr_acc=3)
partitions(start,stop,12)
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
def plot_yearly_avg_sums(model_folders,delta_t_val,cd):
    n = len(model_folders)
    fig=plt.figure(figsize=((n+1)*10,20))#, dpi = 400)
    axs=fig.subplots(1,n)
    plt.rcParams['font.size'] = 20
    #names=['X','X_c']#,'X_p']
    names=['x','x_c']#,'x_p']
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
                vals.__getattribute__(name),
                label=name+"_avg_sum",
                color=cd[name]
            )
        ax.legend()
        ax.set_title(mf)
        
plot_yearly_avg_sums(model_folders, delta_t_val, cd)


# -
# ## Application: Numerical experiment concerning the  attraction of the solution to wards $X_c$.
# In  Yiqis 2017 Biogeosciences paper: "Transient dynamics of terrestrial carbon storage: mathematical foundation and its applications" the carbon storage potential $\mathbf{X_c}(t)$ is called an "attractor" for the solution $\mathbf{X}(t)$ "at ecosystem scale".
# We will try to visualize this for our models.
# Before we can plot anything we have to specify what we mean concretely.
# We will plot: 
# 1. The components of $(\mathbf{X_c})_p \; p \in pools $ with the components of the stocks     
# $(\mathbf{X})_p \; p \in pools$ and the components of the derivative  $(\dot{\mathbf{X}})_p \; p \in pools$
# 2. The  sums of the components e.g. $\sum_{p \in pools}(\mathbf{X_c})_p$
# 3. The scalar variables $x_c,x,\dot{x}$ derived for a one pool surrogate system for the combined mass of all the pools and the combined inputs.
# $$
# \dot{x}=u(t)-m(t)x
# $$ 
# with m(t) specified as follows:
#
# We start with the special case of a linear but nonautonoumous System:
# $$
# \frac{d \mathbf{X}}{d t}= \mathbf{I}(t) - M(t) \mathbf{X} 
# $$
# Taking the sum over all pools yields.
# $$
# \sum_{p \in pools} \left( \frac{d \mathbf{X}}{d t} \right)_p
# =
# \left( \mathbf{I}(t) - M(t) \mathbf{X} \right)_p
# $$
#
# With: 
# $$
# u=\sum_{p \in pools} (\mathbf{I})_p, 
# $$
# $$
# x = \sum_{p \in pools} (\mathbf{X})_p
# \text{ and }
# $$ 
# $$
# \sum_{p \in pools} \left( \frac{d \mathbf{X}}{d t} \right)_p
# =\frac{d}{d t}\sum_{p \in pools} (\mathbf X )_p
# =\frac{d}{d t} x
# $$
# We can now try to costruct our new system for the combined mass $x$, in particular we want to find a function for the time dependent rate $m(t)$ such that. 
# $$
# \dot{x}
# =u(t)-m(t) x 
# =\sum_{p \in pools} \left( \mathbf{I}(t) - M(t) \mathbf{X} \right)_p
# =u(t)-\sum_{p \in pools} ( M(t) \mathbf{X} )_p
# $$
# This yields: 
# $$
# m(t) = \frac{
#     \sum_{p \in pools} ( M(t) \mathbf{X} )_p
#     }{
#     \sum_{p \in pools} (\mathbf{X})_p
#     }
# $$
# We can even extend this treatment to nonlinear systems:
# $$
# \frac{d \mathbf{X}}{d t}= \mathbf{I}(\mathbf{X},t) - M(\mathbf{X},t) \mathbf{X} 
# $$
# Assume that we first solve the system numerically and therefore have $\mathbf{X}(t)$ available.
# Substituting the solution we get:
# $$
#
# $$
#

# +
def plot_sums(model_folders,delta_t_val,cd):
    n = len(model_folders)
    fig=plt.figure(figsize=((n+1)*10,20))#, dpi = 400)
    axs=fig.subplots(3,n)
    plt.rcParams['font.size'] = 20
    names_0=['X','X_c','X_p']
    names_1=['x','x_c','x_p']
    names_2=['X_dot']
    for i,mf in enumerate(model_folders):
        itr=tracebility_iterator(mf,delta_t_val)
        start_min,stop_max=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
        # we do not want the whole interval but look at a smaller part to observe the dynamics
        start,stop = start_min, int(start_min+(stop_max-start_min)/30)
        vals=itr[start:stop]
        times=times_in_days_aD(mf,delta_t_val)[start:stop]/365
        ax=axs[0,i]
        for name in names_0:
            ax.plot(
                times,
                vals.__getattribute__(name).sum(axis=1),
                label=name+"_sum",
                color=cd[name]
            )
        ax.legend()
        ax.set_title(mf)
        
        ax=axs[1,i]
        for name in names_1:
            ax.plot(
                times,
                vals.__getattribute__(name),
                label=name+"1_pool_surrogate",
                color=cd[name]
            )
        ax.legend()
        ax.plot(
            times,
            np.zeros_like(times),
            #label=name+"_sum",
            color=cd[name]
        )
        ax.set_title(mf)
        
        ax=axs[2,i]
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
    
plot_sums(model_folders,delta_t_val, cd)


# -


# We see that the sum of the derivative is clearly positive all 
# the time even if $\sum_{p \in pools}(\mathbf{X_c})_p$ 
# crosses the $\sum_{p \in pools}( \mathbf{X})_p$ lines 
# which shows that  $\sum_{p \in pools}(\mathbf{X_c})_p$ is NOT ALWAYS 
# attractive in the sense that the summed derivative  $\sum_{p \in pools}
# (\mathbf{\dot{X}})_p$ points in the same direction as the difference
# $\sum_{p \in pools}(\mathbf{X_p})_p$ 
#

# +
def plot_components(mf,delta_t_val,cd):
    names_1=['X','X_c']#,'X_p']
    names_2=['X_dot']
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
plot_components(model_folders[1],delta_t_val,cd)
# -

[ i for i in range(0,10,2)]


# ### Numerical expiriment

# +
from scipy.integrate import solve_ivp

def u(t):
    return (np.cos(t)+1)
def M(t):
    #return (np.sin(t)+2)*0.1
    return 0.1

def X_dot(t, X): 
    return u(t) - M(t)*X

X_0=0.1
t_start=0
t_end=4*np.pi


sol = solve_ivp(X_dot, t_span=[t_start,t_end],t_eval=np.linspace(t_start,t_end,100),y0=np.array([X_0]))

times=sol.t
n=len(times)
Xs=sol.y.transpose().reshape(n,)

def inv_M(t):
    return 1/M(t)

def X_c(t):
    return inv_M(t)*u(t)
X_cs = np.array(list(map(X_c,times)))
us= np.array(list(map(u,times)))
X_ps = X_cs-Xs
X_dot_ts=np.array(list(map(lambda i:X_dot(times[i],Xs[i]),range(n))))    
Xs.shape,X_cs.shape,X_dot_ts.shape
d={
    "X": Xs, 
    "X_c": X_cs,
    "X_p": X_ps,
    "X_dot": X_dot_ts,
    "u": us 
}
names_1=['X','X_c','X_p']
names_2=['X_dot','u']
f=plt.figure(figsize=(20,20))
axs=f.subplots(2,1)
for name in names_1:
    axs[0].plot(times,d[name],color=cd[name],label=name)
axs[0].legend()

for name in names_2:
    axs[1].plot(times,d[name],color=cd[name],label=name)
    axs[1].plot(times,np.zeros_like(times),color='black')
axs[1].legend()

# +
# np.array([1,2]).reshape?
# -

np.array([1,2]).reshape

np.array([1,2]).reshape

# +
# np.reshape?
# -

np.array


