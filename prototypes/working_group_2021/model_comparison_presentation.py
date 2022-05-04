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
# <p> We use the <a href="https://github.com/MPIBGC-TEE/bgc_md2">bgc_md2</a> package infrastructure to compare two or more models using traceability analysis approach (<a href="https://doi.org/10.1111/gcb.12172">Xia et al. 2013</a>, <a href="https://doi.org/10.5194/bg-14-145-2017">Luo et al. 2017</a>, <a href="https://doi.org/10.1002/2017MS001004">Jiang et al. 2017</a>) modified for transient simulations of TRENDY model intercomparison project. </p>
# <p>The biogeochemical models from TRENDY project have been reconstructed from literature using matrix approach (<a href="https://doi.org/10.1029/2002GB001923">Luo et al. 2003</a>, <a href=" https://doi.org/10.1111/gcb.12766">2015</a>, <a href="https://doi.org/10.5194/gmd-5-1045-2012">Sierra et al. 2012</a>) and <a href="https://www.sympy.org/en/index.html">sympy</a> python package (reference) in the <a href="https://github.com/MPIBGC-TEE/bgc_md2">bgc_md2</a> database. Currently the models are simplified to be driven by NPP and thus omitting explicit simulation of leaf processes and autotrophic respiration. Data assimilation (<a href=" https://doi.org/10.1890/09-1275.1">Luo et al. 2011</a>, <a href="https://doi.org/10.1002/2015GB005239">2015</a>) was used to optimize parameters of reconstructed models to fit TRENDY output. </p>
# <ol>
# <li>We start by comparing model outputs - <em><strong>C storage (X)</strong></em> over time. Then we compute and compare traceable components to investigate sources of discrepancy between model predictions of C storage, and thus sources of uncertainty in our understanding of global C dynamics. </li>
# <li>We compute <em><strong>C storage capacity (X<sub>C</sub>)</strong></em> for each model. <em><strong>X<sub>C</sub></strong></em> for each point in time shows what the system <em><strong>X</strong></em> would be if the the system was in the steady state. <em><strong>C storage potential (X<sub>P</sub>)</strong></em> is the difference between <em><strong>X<sub>C</sub></strong></em> and <em><strong>X</strong></em> at each point in time - it shows how far the current C storage of the system is from the steady state.</li> 
# <li><em><strong>X<sub>C</sub></strong></em> of each model depends on <em><strong>C input</strong></em> (in our case - <em><strong>NPP</strong></em>) and <em><strong>Residence Time</strong></em> (a characteristic of model structure). We compare <em><strong>NPP</strong></em> and <em><strong>RT</strong></em> for each model and attribute the discrepancy of <em><strong>X<sub>C</sub></strong></em> to <em><strong>NPP</strong></em> and <em><strong>RT</strong></em>.</li>
# <li><em><strong>RT</strong></em> is a dynamic characteristic of a model that depends on the model structure and on its sensitivity to environmental factors. We perform sensitivity analysis to determine how <em><strong>RT</strong></em> of each model depends on <em><strong>temperature</strong></em> and <em><strong>moisture</strong></em>. Then we run the models with same and constant temperature and moisture to determine differences in <em><strong>baseline RT</strong></em>.</li>
# </ol>
# <p>The analysis is currently performed for total global C storage, including vegetation and soil parts. <br> Next steps to expand the analysis may include the following: </p>
# <ul>
# <li>Explore temperature and moisture sensitivities of major fluxes (e.g. litter decomposition, soil respiration);</li> 
# <li>Separetely compare vegetation and soil components for each model;</li> 
# <li>Investigate differences between model predictions over different biomes;</li>
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

# +
# define models to compare as folder name : model name
model_names={ 
    "yz_jules": "JULES",
    "kv_visit2": "VISIT",
    "yz_jules2": "JULES_2", # placeholder for a 3rd model - copy folder yz_jules and call it "yz_jules2" to run this
}
model_folders=[(k) for k in model_names]

# selecting colors for plotting models
model_cols={
    "yz_jules": "blue",
    "kv_visit2": "orange",
    "yz_jules2": "green",
}
# -

# ### Exploring model structures

# bgc_md2 automatically generates a flow diagram, input allocation vector and compartmental matrix for each model
comparison_table=gh.model_table(model_names)
comparison_table

# ### Plots of traceable components

# +
# assuming same step size for each model
delta_t_val=30

# variables to plot
var_names={
        'x':'Total Carbon (X) and Carbon Storage Capacity (X_c)',
        'x_p':'Carbon Storage Potential (X_p)',
        'u':'Carbon Input (NPP)',
        'rt':'Residense Time (RT)'
    } 

gh.plot_yearly_components(model_names=model_names, var_names=var_names, delta_t_val=delta_t_val, model_cols=model_cols)
# -

# ## Below are plots in development

# +
####################################################
## plot delta_x_c delta_u and delta_rt for all possible pairs of models 
####################################################

# define models to compare as folder name : model name
model_names={ 
    "yz_jules": "JULES",
    "kv_visit2": "VISIT",
    "yz_jules2": "JULES_2", # placeholder for a 3rd model - copy folder yz_jules and call it "yz_jules2" to run this
}
model_folders=[(k) for k in model_names]

var_names={
        'x_c':'Carbon Storage Capacity (X_c)',
        'u':'Carbon Input (NPP)',
        'rt':'Residense Time (RT)'
    } 
                
gh.plot_yearly_diff(model_names=model_names, var_names=var_names, delta_t_val=delta_t_val)

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

    itr_1=gh.traceability_iterator_instance(mf_1,delta_t_val)
    start_1,stop_1=gh.min_max_index(mf_1,delta_t_val,*gh.t_min_tmax([mf_1,mf_1],delta_t_val))
    parts_1=gh.partitions(start_1,stop_1,12)
    times_1=gh.averaged_times(
        gh.times_in_days_aD(mf_1,delta_t_val)/365,
        parts_1
        )
    vals_1=itr_1.averaged_values(
        parts_1
        )        
    

    itr_2=gh.traceability_iterator_instance(mf_2,delta_t_val)
    start_2,stop_2=gh.min_max_index(mf_2,delta_t_val,*gh.t_min_tmax([mf_1,mf_2],delta_t_val))    
    parts_2=gh.partitions(start_2,stop_2,12)
    times_2=gh.averaged_times(
        gh.times_in_days_aD(mf_2,delta_t_val)/365,
        parts_2
        )
    vals_2=itr_2.averaged_values(
        parts_2
        )
    
    fig=plt.figure(figsize=(20,(n)*10))
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


