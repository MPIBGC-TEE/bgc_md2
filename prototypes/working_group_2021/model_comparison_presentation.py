# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Matrix Model Intercomparison Project: Traceability Analysis
# <p>We use the <a href="https://github.com/MPIBGC-TEE/bgc_md2">bgc_md2</a> package infrastructure to compare two or more models using traceability analysis approach (<a href="https://doi.org/10.1111/gcb.12172">Xia et al. 2013</a>, <a href="https://doi.org/10.5194/bg-14-145-2017">Luo et al. 2017</a>, <a href="https://doi.org/10.1002/2017MS001004">Jiang et al. 2017</a>) modified for transient simulations of the TRENDY model intercomparison project. </p>
# <p>The biogeochemical models from the TRENDY project have been reconstructed from literature using matrix approach (<a href="https://doi.org/10.1029/2002GB001923">Luo et al. 2003</a>, <a href=" https://doi.org/10.1111/gcb.12766">2015</a>, <a href="https://doi.org/10.5194/gmd-5-1045-2012">Sierra et al. 2012</a>) and <a href="https://www.sympy.org/en/index.html">sympy</a> python package in the <a href="https://github.com/MPIBGC-TEE/bgc_md2">bgc_md2</a> database. Currently the models are simplified to be driven by NPP and thus omitting explicit simulation of leaf processes and autotrophic respiration. Data assimilation (<a href=" https://doi.org/10.1890/09-1275.1">Luo et al. 2011</a>, <a href="https://doi.org/10.1002/2015GB005239">2015</a>) was used to optimize parameters of reconstructed models to fit TRENDY output.</p>
# <p>Currently our analysis includes the following steps:
# <ol>
# <li>We start by comparing model outputs - <em><strong>C storage</strong></em> ($X$) over time. Then we compute and compare traceable components to investigate sources of discrepancy between model predictions of C storage, and thus sources of uncertainty in our understanding of global C dynamics. </li>
# <li>We compute <em><strong>C storage capacity </strong></em> $X_{C}$ for each model. $X_{C}$ (as a function of time) represents the maximum C storage for the system under current conditions. 
# It is a theoretically predicted steady state of the autonomous system that we would obtain if we froze our original non-autonomous system at time point $t_{freeze}$ keeping the conditions (C input and environmental factors) constant after $t_{freeze}$. While $X_{C}$ is not a correct prediction for the solution $X$ of the original time dependent system, it is rather a momentary characteristic of the system that depends on conditions at each point in time.  
# $X_{C}$ is an attractor of the frozen autonomous system, although it is not attracting the real solution  $X$. In a non-autonomous system with variable C input and environmental factors (such as temperature and moisture)  <em><strong>ecosystem</strong></em> $X$ (if we consider ecosystem as a whole - as a <em>"surrogate" 1-pool system</em>, not pool-wise) is constantly chasing $X_{C}$, but cannot consistently reach it due to variable conditions: when $X_{C}$ &gt $X$, $X$ is increasing, and when $X_{C}$ &lt $X$, $X$ is decreasing.</li>
# <li>The difference between $X_{C}$ and $X$ at each point in time is called  <em><strong>C storage potential </strong></em>($X_{P}$): it shows how far the current C storage of the system is from the theoretical steady state (in a positive or negative direction).</li> 
# <li>$X_{C}$ of each model depends on <em><strong>C input</strong></em> (in our case - $NPP$) and <em><strong>Equilibrium Residence Time </strong></em>($RT$). We compare $NPP$ and $RT$ for each model and attribute the contribution of $NPP$ and $RT$ to the discrepancy of $X_{C}$.</li>
# <li>$RT$ represents the time which it would take on average for C particle to exit the system after it entered it, if the system was at the equilibrium. In a transient simulation, where the system is not at the steady state, $RT$ becomes a dynamic model characteristic that can be interpreted as a measure of the inverse C turnover rate of the system at each point in time. $RT$ depends on the model structure and on environmental factors. To single out environmental effects we determine <em><strong>temperature</strong></em> and <em><strong>moisture sensitivity</strong></em> of $RT$, and compare it between models. We also compare $RT$ at fixed temperatures and moisture conditions (gauge conditions) including minimum, maximum and mean conditions across the period of the simulation.</li>
# <li> Based on the traceable components: $NPP$, $RT$, <em><strong>temperature</strong></em> and <em><strong>moisture sensitivity</strong></em> of $RT$, we can make conclusions regarding the inherent similarity/dissimilarity between the models and single out mechanisms that contribute most to the discrepancy of C predictions. 
# </ol>
# <p>The analysis is currently performed for total global C storage, including vegetation and soil parts. <br> Next steps to expand the analysis may include the following: </p>
# <ul>
# <li>Explore temperature and moisture sensitivities of major fluxes (e.g. litter decomposition, soil respiration);</li> 
# <li>Separetely compare vegetation and soil components for each model;</li> 
# <li>Investigate differences over different biomes;</li>
# <li>Include additional diagnostic variables (e.g. transit time, carbon use efficiency, etc.);</li>
# <li>Expand models by explicitly including autotrophic processes (start from GPP as C input);</li>
# <li>Include sensitivity to more environmental factors and anthropogenic disturbances</li>
# </ul>
# <p>The short description of methodology for deriving traceable components is given below. </p>

from IPython.display import Markdown, display
display(Markdown("TracebilityText.md"))

# ### Loading required packages  and functions

# %load_ext autoreload
# %autoreload 2
import matplotlib.pyplot as plt
import numpy as np
from functools import lru_cache
import general_helpers as gh
import MIP_output_helpers as moh
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
    #"Aneesh_SDGVM":"SDGVM",
    "cj_isam": "ISAM",
    ##"bian_ibis2":"IBIS",
    ##"ORCHIDEE-V2":"OCN",
}

# ### Exploring model structures

# bgc_md2 automatically generates a flow diagram, input allocation vector and compartmental matrix for each model
comparison_table=gh.model_table(model_names)
comparison_table

# ### Loading TRENDY data and model parameters

# define same step size for each model (in days)
delta_t_val=30

# load data and parameters
model_folders=[(k) for k in model_names]
#test_arg_list=gh.get_test_arg_list(model_folders)

# fixme mm 8-12: 
# it would make sense to create a dictionary indexed by the model name 
# so it can be used for model comparisons by name like this one
test_arg_dict = {mf: gh.test_args(mf) for mf in model_folders}
#test_arg_dict

# ### Checking data assimilation quality for all models

plt.rcParams.update({'font.size': 32})

# +
import matplotlib.lines as mlines

#for j, mf in enumerate(model_folders):
for j, mf in enumerate(["cj_isam"]):
    print ('\033[1m'+"Matrix version vs original Trendy output: "+ mf)
    print ('\033[0m')
    fig = plt.figure(figsize=(50,10))
    axs=fig.subplots(1, len(gh.msh(mf).Observables._fields))
  
    mvs=test_arg_dict[mf].mvs
    dvs=test_arg_dict[mf].dvs
    cpa=test_arg_dict[mf].cpa
    epa_opt=test_arg_dict[mf].epa_opt
  
    param2res_sym = gh.msh(mf).make_param2res_sym(mvs,cpa,dvs)
    out_simu=param2res_sym(epa_opt)._asdict()
    
    obs=test_arg_dict[mf].svs._asdict()

    
    print ("Amount of variance explined: ")
    for i,f in enumerate(gh.msh(mf).Observables._fields):
        resid=out_simu[f]-obs[f]
        mean_obs = obs[f].mean()
        mean_centered_obs = obs[f] - mean_obs
        AVE=1 - np.sum( resid**2) / np.sum( mean_centered_obs**2 )
        print(f+ " : " + str(round(AVE,3)) )
    for i,f in enumerate(gh.msh(mf).Observables._fields):
        axs[i].scatter(out_simu[f], obs[f], c='black')
        line = mlines.Line2D([0, 1], [0, 1], color='red')
        transform = axs[i].transAxes
        line.set_transform(transform)
        axs[i].add_line(line)
        axs[i].set_title(f)
        axs[i].set_xlabel('Matrix output')
        axs[i].set_ylabel('Original Trendy output')
        axs[i].grid()
    plt.show()
# -


for k in obs.keys():
    print(k, obs[k].shape, out_simu[k].shape)


# ### Plots of traceable components

plt.rcParams.update({'font.size': 16})

model_cols={
    "JULES": "blue",
    "VISIT": "orange",
    "YIBs": "green",
    "DLEM": "red",
    "SDGVM":"yellow",
    "ISAM": "purple",
    "IBIS":"magenta",
    "OCN":"teal",
}

all_comp_dict= moh.get_traceable_components(model_names=model_names,
             test_arg_list=[test_arg_dict[k] for k in model_names],
             delta_t_val=delta_t_val, 
             model_cols=model_cols,
             part=1,
             averaging=12*30//delta_t_val, # yearly averaging
             #averaging=30//delta_t_val, # monthly averaging
             overlap=True
             )

x_x_c,sigma_x_x_c=moh.plot_traceable_component(
    all_comp_dict,
    "x_x_c",
    model_cols,
    #delta=True,
)

x1,sigma_x1=moh.plot_traceable_component(
    all_comp_dict,
    "x",
    model_cols,
    #delta=True,
)
x2,sigma_x2=moh.plot_traceable_component(
    all_comp_dict,
    "x",
    model_cols,
    delta=True,
)

# +
times=all_comp_dict["Times"]
var=sigma_x1
moh.plot_single_trend(var,times,3,"Standard deviation of X over time")

times=all_comp_dict["Times"]
var=sigma_x2
moh.plot_single_trend(var,times,3, "Standard deviation of $\Delta$X over time")
# -

x_c1,sigma_x_c1=moh.plot_traceable_component(
    all_comp_dict,
    "x_c",
    model_cols,
    #delta=True,
)
x_c2,sigma_x_c2=moh.plot_traceable_component(
    all_comp_dict,
    "x_c",
    model_cols,
    delta=True,
)

# +
times=all_comp_dict["Times"]
var=sigma_x_c1
moh.plot_single_trend(var,times,3,"Standard deviation of X_c over time")

times=all_comp_dict["Times"]
var=sigma_x_c2
moh.plot_single_trend(var,times,3, "Standard deviation of $\Delta$X_c over time")
# -

x_p1,sigma_x_p1=moh.plot_traceable_component(
    all_comp_dict,
    "x_p",
    model_cols,
    #delta=True,
)
# x_p2,sigma_x_p2=moh.plot_traceable_component(
#     all_comp_dict,
#     "x_p",
#     model_cols,
#     delta=True,
# )

# +
# times=all_comp_dict["Times"]
# var=sigma_x_p1
# moh.plot_single_trend(var,times,3,"Uncertainty of X_p over time")

# times=all_comp_dict["Times"]
# var=sigma_x_p2
# moh.plot_single_trend(var,times,3, "Uncertainty of $\Delta$X_p over time")
# -

times=all_comp_dict["Times"]
var=x_p1
moh.plot_single_trend(var,times,3,"Mean X_p over time")

u1, sigma_u1=moh.plot_traceable_component(
    all_comp_dict,
    "u",
    model_cols,
    #delta=True,
)
u2, sigma_u2=moh.plot_traceable_component(
    all_comp_dict,
    "u",
    model_cols,
    delta=True,
)

# +
times=all_comp_dict["Times"]
var=sigma_u1
moh.plot_single_trend(var,times,3,"Standard deviation of NPP over time")

times=all_comp_dict["Times"]
var=sigma_u2
moh.plot_single_trend(var,times,3, "Standard deviation of $\Delta$NPP over time")
# -

rt1,sigma_rt1=moh.plot_traceable_component(
    all_comp_dict,
    "rt",
    model_cols,
    #delta=True,
)
rt2,sigma_rt2=moh.plot_traceable_component(
    all_comp_dict,
    "rt",
    model_cols,
    delta=True,
)

# +
times=all_comp_dict["Times"]
var=sigma_rt1
moh.plot_single_trend(var,times,3,"Standard deviation of RT over time")

times=all_comp_dict["Times"]
var=sigma_rt2
moh.plot_single_trend(var,times,3, "Standard deviation of $\Delta$RT over time")
# -

# ## Contribution of traceable components to the uncertainty of C storage (X)

plt.rcParams.update({'font.size': 12})

all_comp_dict['VISIT'].keys()

"u" in all_comp_dict['VISIT'].keys()

moh.plot_attribution_sum (
    all_comp_dict=all_comp_dict,
    percent=True,
    part=1,
)

# ### Per model

moh.plot_attribution_per_model(
    all_comp_dict=all_comp_dict,
    #percent=True,    
)


