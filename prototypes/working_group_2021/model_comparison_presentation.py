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

# # Matrix Model Intercomparison Project: Traceability Analysis
# <p> We use the <em>bgc_md2</em> package infrastructure to compare two or more models using traceability analysis approach (Jia et al., Luo et al., Jiang et al.) modified for transient simulations of TRENDY model intercomparizon project. </p>
# <p>The biogeochemical models from TRENDY project have been reconstructed from literature using matrix approach (Luo et al.) and <em>sympy</em> python package (reference) in the <em>bgc_md2</em> database. Currently the models are simplified to be driven by NPP and thus omitting explicit simulation of leaf processes and autotrophic respiration. Data assimilation was used to optimize parameters of reconstructed models to fit TRENDY output. </p>
# <ol>
# <li>We start by comparing model outputs - <em><strong>C storage (X)</strong></em> over time. Then we compute and compare traceable components to investigate sources of discrepancy between model predictions of C storage, and thus sources of uncertainty in our understanding of global C dynamics. </li>
# <li>We compute <em><strong>C storage capacity (X<sub>C</sub>)</strong></em> for each model. <em><strong>X<sub>C</sub></strong></em> for each point in time shows what the system <em><strong>X</strong></em> would be if the the system was in the steady state. <em><strong>C storage potential (X<sub>P</sub>)</strong></em> is the difference between <em><strong>X<sub>C</sub></strong></em> and <em><strong>X</strong></em> at each point in time - it shows how far the current C storage of the system is from the steady state.</li> 
# <li><em><strong>X<sub>C</sub></strong></em> of each model depends on <em><strong>C input</strong></em> (in our case - <em><strong>NPP</strong></em>) and <em><strong>Residence Time</strong></em> (a characteristic of model structure). We compare <em><strong>NPP</strong></em> and <em><strong>RT</strong></em> for each model and attribute the discrepancy of <em><strong>X<sub>C</sub></strong></em> to <em><strong>NPP</strong></em> and <em><strong>RT</strong></em>.</li>
# <li><em><strong>RT</strong></em> is a dynamic characteristic of a model that depends on the model structure and on its sensitivity to environmental factors. We perform sensitivity analysis to determine how <em><strong>RT</strong></em> of each model depends on <em><strong>temperature</strong></em> and <em><strong>moisture</strong></em>.</li>
# </ol>
# <p>The analysis is currently performed for total global C storage, including vegetation and soil parts. <br> Next steps to expand the analysis may include the following: </p>
# <ul>
# <li>Explore temperature and moisture sensitivities of major fluxes (e.g. litter decomposition, soil respiration);</li> 
# <li>Separetely compare vegetation and soil components for each model;</li> 
# <li>Investigate differences between model predictions over different land covers / plant functional types;</li>
# <li>Expand models by explicitely including autotrophic processes;</li>
# <li>Include sensitivity to more environmental factors and anthropogenic disturbances</li> 
# </ul>
# <p>The short description of methodology for deriving traceable components is given below. </p>

from IPython.display import Markdown, display
display(Markdown("TracebilityText.md"))

# ### Loading required packages  and functions

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

# ### Selecting models to compare

# define some shortcut functions that yield the result depending on the model folder  mf
model_folders=['yz_jules', 
               'kv_visit2', 
               'yz_jules2'  # placeholder for a 3rd model - copy folder yz_jules and call it "yz_jules2" to run this
              ]#, 'Aneesh_SDGVM']#, 'kv_ft_dlem', 'jon_yib']#,'cable-pop']#,]
model_names={
    "yz_jules": "YULES",
    "kv_visit2": "VISIT",
    "yz_jules2": "YULES_2",
}
# selecting colors for plotting models
model_cols={
    "yz_jules": "blue",
    "kv_visit2": "orange",
    "yz_jules2": "green",
}

# ### Exploring model structures

# +
# bgc_md2 automatically generates a flow diagram, input allocation vector and compartmental matrix for each model
mf=model_folders[0]
import_module("{}.source".format(mf))
def mvs(mf):
    return import_module("{}.source".format(mf)).mvs

tups=[(mf,mvs(mf))
      for mf in model_folders
]

dh.table(
    tups=tups,
    types_compact = [],
    types_expanded= [InputTuple,CompartmentalMatrix]

)


# -

# ### Functions to standartize model outputs

# +
# functions to synchronize model outputs to the scale of days since AD
def sim_day_2_day_aD_func(mf): #->function
    return gh.msh(mf).make_sim_day_2_day_since_a_D(gh.confDict(mf))

def times_in_days_aD(mf,delta_t_val):
    n_months=len(gh.test_args(mf).dvs[0])
    n_days=n_months*30
    n_iter=int(n_days/delta_t_val)
    days_after_sim_start=delta_t_val*np.arange(n_iter)
    return np.array(tuple(map(sim_day_2_day_aD_func(mf),days_after_sim_start))) 

# function to determine overlapping time frames for models simulations 
def t_min_tmax(model_folders,delta_t_val):
    td={
        mf: times_in_days_aD(mf,delta_t_val)
        for mf in model_folders
    }
    t_min = max([t.min() for t in td.values()])
    t_max = min([t.max() for t in td.values()])
    return (t_min,t_max)

# function find the timesteps corresponding to shared times
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

def averaged_times(times,partitions):
    return np.array(
        [
            times[p[0]:p[1]].sum()/(p[1]-p[0]) for p in partitions
        ]
    )


# -

# ### The iterator for running model simulations and outputting traceable components. 

# +
# The iterator allows us to compute and easily access the desired timelines with python index notation [2:5] 
# it will return the values from position 2 to 5 of the solution (and desired variables).
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

# colors for plotting traceable components
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
# -
# ### Plots of traceable components

# +
# assuming same step size for each model
delta_t_val=30

from scipy.interpolate import interp1d, splprep

def plot_yearly_components(model_folders, delta_t_val, model_cols):
    variables=['x', 'x_p','u','rt']
    names={
        'x':'Total Carbon (X) and Carbon Storage Capacity (X_c)',
        'x_p':'Carbon Storage Potential (X_p)',
        'u':'Carbon Input (NPP)',
        'rt':'System Residense Time (RT)'
    }
    n=len(variables)
    yr=int(365/delta_t_val)
    
    fig=plt.figure(figsize=(30,n*10))
    axs=fig.subplots(n,1)
    
    # calculate output array length for the average base on the 1st model
    start,stop=min_max_index(model_folders[0],delta_t_val,*t_min_tmax(model_folders,delta_t_val))
    output_length=int((stop-start)/12)+1
    
    for i,name in enumerate(variables):
        sum=np.zeros(output_length)
        for mf in model_folders:
            itr=tracebility_iterator(mf,delta_t_val)
            start,stop=min_max_index(mf,delta_t_val,*t_min_tmax(model_folders,delta_t_val))
            parts=partitions(start,stop,yr)
            times=averaged_times(
                times_in_days_aD(mf,delta_t_val)/365,
                parts
            )
            vals=itr.averaged_values(
                parts
            )
            sum+=vals.__getattribute__(name) # for calculating avarage
            ax=axs[i]
            ax.plot(
                times,
                vals.__getattribute__(name),
                label=model_names[mf]+' - X',
                color=model_cols[mf],
            )
            if name=='x':  # we plot x together with x_c
                ax.plot(
                    times,
                    vals.__getattribute__('x_c'),
                    label=model_names[mf]+' - X_c',
                    color=model_cols[mf],
                    linestyle = 'dashed'
                )
            if name=='x_p':  # 0 line for X_p
                ax.plot(
                    times,
                    np.zeros_like(times),
                    color="black",
                    linestyle = 'dotted',
                    alpha=0.5
        )
        if name!='x_p':     # plot avarage for all models (doesn't make sense for X_p)   
            ax.plot(
                times,
                sum/len(model_names),
                label='Average (all models)',
                color="black",
                linestyle = 'dotted',                
                alpha=0.5
            )                
            ax.legend()
            ax.set_title(names[name])
                           
plot_yearly_components(model_folders=model_folders,delta_t_val=delta_t_val,model_cols=model_cols)
# -

# ## Below are plots in development

# +
####################################################
## plot delta_x_c delta_u and delta_rt for all possible pairs of models 
####################################################
# Since the two models do not necessarily share the same point in time and not
# even the same stepsize or number of steps we compute interpolating functions
# to make them comparable
model_folders=['yz_jules', 
               'kv_visit2', 
               'yz_jules2'  # placeholder for a 3rd model - copy folder yz_jules and call it "yz_jules2" to run this
              ]#, 'Aneesh_SDGVM']#, 'kv_ft_dlem', 'jon_yib']#,'cable-pop']#,]

from scipy.interpolate import interp1d, splprep

def diffp(ax,name,mf_1,times_1,vals_1,mf_2,times_2,vals_2):
    f1=interp1d(times_1,vals_1.__getattribute__(name))
    f2=interp1d(times_2,vals_2.__getattribute__(name))
    # chose the interval covered by both to avoid extrapolation
    start=max(times_1.min(),times_2.min())
    stop=min(times_1.max(),times_2.max())
    nstep=min(len(times_1),len(times_2))
    times=np.linspace(start,stop,nstep)
        
    diff=f1(times)-f2(times)
    #diff=vals_1.__getattribute__(name)-vals_2.__getattribute__(name) # if time step is same, this should work instead of interpolation
    ax.plot(times,diff,color="black")
    ax.set_title("{0}_{1}-{2}".format(name,mf_1,mf_2))    
    
def plot_yearly_diff(model_folders, delta_t_val):
    variables=['x_c', 'u','rt']
    names={
        'x_c':'Carbon Storage Potential (X_p)',
        'u':'Carbon Input (NPP)',
        'rt':'System Residense Time (RT)'
    }
    n=int(len(variables)*( len(model_folders) * (len(model_folders)-1) / 2))
    yr=int(365/delta_t_val)
    
    fig=plt.figure(figsize=(30,n*10))
    axs=fig.subplots(n,1)
    part=1 # for the whole interval
    plot_number=0
    for i,name in enumerate(variables):
        for j, mf_1 in enumerate(model_folders[0:len(model_folders)-1]):            
            for mf_2 in model_folders[j+1:len(model_folders)]:
                start_min_1,stop_max_1=min_max_index(mf_1,delta_t_val,*t_min_tmax([mf_1,mf_2],delta_t_val))
                # we do not want the whole interval but look at a smaller part to observe the dynamics
                start_1,stop_1 = start_min_1, int(start_min_1+(stop_max_1-start_min_1)/part)                               
                itr_1=tracebility_iterator(mf_1,delta_t_val)
                
                parts_1=partitions(start_1,stop_1,yr)
                times_1=averaged_times(
                    times_in_days_aD(mf_1,delta_t_val)/365,
                    parts_1
                )
                vals_1=itr_1.averaged_values(
                    parts_1
                )        
                
                start_min_2,stop_max_2=min_max_index(mf_2,delta_t_val,*t_min_tmax([mf_2,mf_2],delta_t_val))
                # we do not want the whole interval but look at a smaller part to observe the dynamics
                start_2,stop_2 = start_min_2, int(start_min_2+(stop_max_2-start_min_2)/part)
                itr_2=tracebility_iterator(mf_2,delta_t_val)                
                parts_2=partitions(start_2,stop_2,yr)
                times_2=averaged_times(
                    times_in_days_aD(mf_2,delta_t_val)/365,
                    parts_2
                )
                vals_2=itr_2.averaged_values(
                    parts_2
                )   
                print ("Plotting "+str(plot_number+1)+" out of "+str(n))
                diffp(axs[plot_number],name,mf_1,times_1,vals_1,mf_2,times_2,vals_2)
                plot_number+=1 
                
plot_yearly_diff(model_folders=model_folders,delta_t_val=delta_t_val)

# +
# this is to show full duration for YULES and VISIT.
# JULES starts in 1700 and VISIT - in 1860.
# X and NPP are very different for both models at the start of the simulation.
# Therefore, differences also can be traced to the spin-up process to NPP estimion by each model

# assuming same step size for each model
delta_t_val=30

from scipy.interpolate import interp1d, splprep

def plot_yearly_diff(mf_1, mf_2, delta_t_val, model_cols):
    n=len(model_folders)

    itr_1=tracebility_iterator(mf_1,delta_t_val)
    start_1,stop_1=min_max_index(mf_1,delta_t_val,*t_min_tmax([mf_1,mf_1],delta_t_val))
    parts_1=partitions(start_1,stop_1,12)
    times_1=averaged_times(
        times_in_days_aD(mf_1,delta_t_val)/365,
        parts_1
        )
    vals_1=itr_1.averaged_values(
        parts_1
        )        
    

    itr_2=tracebility_iterator(mf_2,delta_t_val)
    start_2,stop_2=min_max_index(mf_2,delta_t_val,*t_min_tmax([mf_1,mf_2],delta_t_val))    
    parts_2=partitions(start_2,stop_2,12)
    times_2=averaged_times(
        times_in_days_aD(mf_2,delta_t_val)/365,
        parts_2
        )
    vals_2=itr_2.averaged_values(
        parts_2
        )
    
    fig=plt.figure(figsize=((n)*10,40))
    axs=fig.subplots(5,1)
    
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

    subp(axs[0],"x")        
    subp(axs[1],"x_c")
    subp(axs[2],"x_p")
    subp(axs[3],"u")
    subp(axs[4],"rt")
    
mf_1="yz_jules"
mf_2="kv_visit2"
plot_yearly_diff(mf_1, mf_2,delta_t_val=delta_t_val,model_cols=model_cols)
# -


