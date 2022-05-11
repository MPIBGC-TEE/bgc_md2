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

# # SDGVM model

# +
from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
# %load_ext autoreload
# %autoreload 2


import numpy as np
from pathlib import Path
import json
import sys
sys.path.insert(0,'..') # necessary to import general_helpers
import general_helpers as gh

# +
#model description sourced out
from source import mvs

import model_specific_helpers_2 as msh
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h
# -

h.compartmental_graph(mvs)

dh.mass_balance_equation(mvs)

# Nothing has changed in the model description but be have some more symbols to work with.
# We could type `C_leaf` somewhere without an error since it is now known as a variable.
# In the next step we replace all occurences of `vl` by `C_leaf` 

# +
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

msh.download_my_TRENDY_output(conf_dict)
# -

# Before we build a function to load the data lets look at it to get an idea.
#

svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))

# +

svs_0=msh.Observables(*map(lambda v: v[0],svs))
cpa=msh.Constants(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 cRoot_0 = svs_0.cRoot,
 npp_0=dvs.npp[0],
 rh_0=svs_0.rh,
 number_of_months=len(svs.rh)
)
# -

cpa

# +
#r_C_leaf2abvstrlit=k_C_leaf *f_leaf2abvstrlit
#r_leaf2abvmetlit = k_C_leaf *f_leaf2abvmetlit
#r_wood2abvstrlit = k_C_wood *f_wood2abvstrlit
#r_wood2abvmetlit = k_C_wood *f_wood2abvmetlit
#r_root2belowstrlit = k_C_root * f_root2belowstrlit
#r_root2belowmetlit = k_C_root * f_root2belowmetlit
#r_abvstrlit2surface_microbe = k_C_abvstrlit *f_abvstrlit2surface_microbe
#r_abvmetlit2surface_microbe = k_C_abvmetlit *f_abvmetlit2surface_microbe
#r_abvstrlit2slowsom = k_C_abvstrlit*f_abvstrlit2slowsom
#r_belowstrlit2soil_microbe = k_C_belowstrlit * f_belowstrlit2soil_microbe
#r_belowmetlit2soil_microbe = k_C_belowmetlit  * f_belowmetlit2soil_microbe
#r_belowstrlit2slowsom = k_C_belowstrlit *f_belowstrlit2slowsom
#r_surface_microbe2slowsom = k_C_surface_microbe*f_surface_microbe2slowsom
#r_soil_microbe2slowsom = k_C_soil_microbe *f_soil_microbe2slowsom
#r_slowsom2soil_microbe = k_C_slowsom *f_slowsom2soil_microbe
#r_soil_microbe2passsom = k_C_soil_microbe*f_soil_microbe2passsom
#r_slowsom2passsom = k_C_slowsom*f_slowsom2passsom
#r_passsom2soil_microbe =k_C_passsom * f_passsom2soil_microbe
#r_leached = leached = (leachedwtr30)/18 * (0.01 + 0.04* (1- silt_clay))

# r_C_leaf2abvstrlit=1/360*60*0.2 #k_C_leaf *f_leaf2abvstrlit
# r_leaf2abvmetlit = 1/360*60*0.8 #k_C_leaf *f_leaf2abvmetlit
# r_wood2abvstrlit = 1/360*30*0.2 #k_C_wood *f_wood2abvstrlit
# r_wood2abvmetlit = 1/(360*30)*0.8 #k_C_wood *f_wood2abvmetlit
# r_root2belowstrlit = 1/(360*22)*0.2#k_C_root * f_root2belowstrlit
# r_root2belowmetlit = 1/(360*60)*0.8 #k_C_root * f_root2belowmetlit
# r_abvstrlit2surface_microbe = #k_C_abvstrlit *f_abvstrlit2surface_microbe
# r_abvmetlit2surface_microbe = k_C_abvmetlit *f_abvmetlit2surface_microbe
# r_abvstrlit2slowsom = k_C_abvstrlit*f_abvstrlit2slowsom
# r_belowstrlit2soil_microbe = k_C_belowstrlit * f_belowstrlit2soil_microbe
# r_belowmetlit2soil_microbe = k_C_belowmetlit  * f_belowmetlit2soil_microbe
# r_belowstrlit2slowsom = k_C_belowstrlit *f_belowstrlit2slowsom
# r_surface_microbe2slowsom = k_C_surface_microbe*f_surface_microbe2slowsom
# r_soil_microbe2slowsom = k_C_soil_microbe *f_soil_microbe2slowsom
# r_slowsom2soil_microbe = k_C_slowsom *f_slowsom2soil_microbe
# r_soil_microbe2passsom = k_C_soil_microbe*f_soil_microbe2passsom
# r_slowsom2passsom = k_C_slowsom*f_slowsom2passsom
# r_passsom2soil_microbe =k_C_passsom * f_passsom2soil_microbe
# r_leached = leached = (leachedwtr30)/18 * (0.01 + 0.04* (1- silt_clay))
# -

epa_0=msh.EstimatedParameters(
    
     beta_leaf=0.5, 
     beta_wood=0.4, 
     r_C_leaf2abvstrlit= 0.0003*0.85,
     r_C_leaf2abvmetlit=0.0003*0.85,
     r_C_wood2abvmetlit=0.00006*1.9,
     r_C_wood2abvstrlit=0.00006*1.9,
     r_C_root2belowmetlit=0.000034*2.105,
     r_C_root2belowstrlit=0.000032*2.105,
    
     r_C_abvstrlit2slowsom=0.00003*1.8,
     r_C_abvstrlit2surface_microbe=0.00005*2,
     r_C_abvmetlit2surface_microbe=0.00005,
     r_C_belowmetlit2soil_microbe=0.00002*50,
     r_C_belowstrlit2slowsom=0.00006*60,
     r_C_belowstrlit2soil_microbe=0.0003*60,
    
     #r_C_leached=0.0001,
    
     r_C_passsom2soil_microbe=0.00005*100,
     r_C_slowsom2passsom=0.0005*100,
     r_C_slowsom2soil_microbe=0.0005*100,
     r_C_soil_microbe2passsom=0.006*100,
     r_C_soil_microbe2slowsom=0.006*105,
     r_C_surface_microbe2slowsom=0.006*100,
    
     r_C_abvstrlit_rh=0.00975/11,
     r_C_abvmetlit_rh=0.024667/10,
     r_C_belowstrlit_rh=0.011333/10,
     r_C_belowmetlit_rh=0.028264/10,

     r_C_soil_microbe_rh=0.0003*2.2,
     r_C_slowsom_rh=0.00004*2,
     r_C_passsom_rh=0.000006875*2,
     r_C_surface_microbe_rh=0.000003*2.2,
    
     C_leaf_0=svs_0.cVeg/3,
     C_abvstrlit_0= 0.59,    #svs_0.cLitter/4,
     C_abvmetlit_0= 0.24,    #svs_0.cLitter/4,
     C_blwstrlit_0= 0.00287261,    #svs_0.cLitter/4,
     C_surfacemic_0= 9.4e-5,    #svs_0.cSoil/4,
     C_soilmic_0= 0.049,      #svs_0.cSoil/4,
     C_slow_0= 0.29828058     #svs_0.cSoil/4
)

# +
# now test it 
import matplotlib.pyplot as plt
from general_helpers import plot_solutions
const_params = cpa

param2res_sym = msh.make_param2res_sym(mvs,const_params,dvs)
obs_0 = param2res_sym(epa_0)

# +
#day_indices=month_2_day_index(range(cpa.number_of_months)),

out_simu_d=obs_0._asdict()
obs_d=svs._asdict()

print("Forward run with initial parameters (blue) vs TRENDY output (orange)")
obs_d['cVeg'],svs.cVeg

fig = plt.figure(figsize=(10,50))
axs=fig.subplots(5,1)


for ind,f in enumerate(msh.Observables._fields):
    val_sim=out_simu_d[f]
    val_obs=obs_d[f]
    axs[ind].plot(range(len(val_sim)),val_sim,label=f+"_sim")
    axs[ind].plot(range(len(val_obs)),val_obs,label=f+"_obs")
    axs[ind].legend()
    
fig.savefig('solutions_SDGVM.pdf')

# +
epa_min=np.array(
    msh.EstimatedParameters(
     beta_leaf=0, 
     beta_wood=0, 
     r_C_leaf2abvstrlit= epa_0.r_C_leaf2abvmetlit/100,
     r_C_abvmetlit2surface_microbe= epa_0.r_C_abvmetlit2surface_microbe/100,
     r_C_abvstrlit2slowsom=epa_0.r_C_abvstrlit2slowsom/100,
     r_C_abvstrlit2surface_microbe=epa_0.r_C_abvstrlit2surface_microbe/100,
     r_C_belowmetlit2soil_microbe=epa_0.r_C_belowmetlit2soil_microbe/100,
     r_C_belowstrlit2slowsom=epa_0.r_C_belowstrlit2slowsom/100,
     r_C_belowstrlit2soil_microbe=epa_0.r_C_belowstrlit2soil_microbe/100,
     #r_C_leached=epa_0.r_C_leached/100,
     r_C_leaf2abvmetlit=epa_0.r_C_leaf2abvmetlit/100,
     r_C_passsom2soil_microbe=epa_0.r_C_passsom2soil_microbe/100,
     r_C_root2belowmetlit=epa_0.r_C_root2belowmetlit/100,
     r_C_root2belowstrlit=epa_0.r_C_root2belowstrlit/100,
     r_C_slowsom2passsom=epa_0.r_C_slowsom2passsom/100,
     r_C_slowsom2soil_microbe=epa_0.r_C_slowsom2soil_microbe/100,
     r_C_soil_microbe2passsom=epa_0.r_C_soil_microbe2passsom/100,
     r_C_soil_microbe2slowsom=epa_0.r_C_soil_microbe2slowsom/100,
     r_C_surface_microbe2slowsom=epa_0.r_C_surface_microbe2slowsom/100,
     r_C_wood2abvmetlit=epa_0.r_C_wood2abvmetlit/100,
     r_C_wood2abvstrlit=epa_0.r_C_wood2abvstrlit/100,
     r_C_abvstrlit_rh=epa_0.r_C_abvstrlit_rh/100,
     r_C_abvmetlit_rh=epa_0.r_C_abvmetlit_rh/100,
     r_C_belowstrlit_rh=epa_0.r_C_belowstrlit_rh/100,
     r_C_belowmetlit_rh=epa_0.r_C_belowmetlit_rh/100,
     r_C_soil_microbe_rh=epa_0.r_C_soil_microbe_rh/100,
     r_C_surface_microbe_rh=epa_0.r_C_surface_microbe_rh/100,
     r_C_slowsom_rh=epa_0.r_C_slowsom_rh/100,
     r_C_passsom_rh=epa_0.r_C_passsom_rh/100,
     C_leaf_0=0,
     #C_root_0=svs_0.cVeg/3,
     C_abvstrlit_0=0,
     C_abvmetlit_0=0,
     C_blwstrlit_0=0,
     C_surfacemic_0=0,
     C_soilmic_0=0,
     C_slow_0=0
    )
)

epa_max=np.array(
    msh.EstimatedParameters(
     beta_leaf=1, 
     beta_wood=1, 
     r_C_leaf2abvstrlit= epa_0.r_C_leaf2abvmetlit*100,
     r_C_abvmetlit2surface_microbe= epa_0.r_C_abvmetlit2surface_microbe*100,
     r_C_abvstrlit2slowsom=epa_0.r_C_abvstrlit2slowsom*100,
     r_C_abvstrlit2surface_microbe=epa_0.r_C_abvstrlit2surface_microbe*100,
     r_C_belowmetlit2soil_microbe=epa_0.r_C_belowmetlit2soil_microbe*100,
     r_C_belowstrlit2slowsom=epa_0.r_C_belowstrlit2slowsom*100,
     r_C_belowstrlit2soil_microbe=epa_0.r_C_belowstrlit2soil_microbe*100,
     #r_C_leached=epa_0.r_C_leached*100,
     r_C_leaf2abvmetlit=epa_0.r_C_leaf2abvmetlit*100,
     r_C_passsom2soil_microbe=epa_0.r_C_passsom2soil_microbe*100,
     r_C_root2belowmetlit=epa_0.r_C_root2belowmetlit*100,
     r_C_root2belowstrlit=epa_0.r_C_root2belowstrlit*100,
     r_C_slowsom2passsom=epa_0.r_C_slowsom2passsom*100,
     r_C_slowsom2soil_microbe=epa_0.r_C_slowsom2soil_microbe*100,
     r_C_soil_microbe2passsom=epa_0.r_C_soil_microbe2passsom*100,
     r_C_soil_microbe2slowsom=epa_0.r_C_soil_microbe2slowsom*100,
     r_C_surface_microbe2slowsom=epa_0.r_C_surface_microbe2slowsom*100,
     r_C_wood2abvmetlit=epa_0.r_C_wood2abvmetlit*100,
     r_C_wood2abvstrlit=epa_0.r_C_wood2abvstrlit*100,
     r_C_abvstrlit_rh=epa_0.r_C_abvstrlit_rh*100,
     r_C_abvmetlit_rh=epa_0.r_C_abvmetlit_rh*100,
     r_C_belowstrlit_rh=epa_0.r_C_belowstrlit_rh*100,
     r_C_belowmetlit_rh=epa_0.r_C_belowmetlit_rh*100,
     r_C_soil_microbe_rh=epa_0.r_C_soil_microbe_rh*100,
     r_C_surface_microbe_rh=epa_0.r_C_surface_microbe_rh*100,
     r_C_slowsom_rh=epa_0.r_C_slowsom_rh*100,
     r_C_passsom_rh=epa_0.r_C_passsom_rh*100,
     C_leaf_0=svs_0.cVeg,
     #C_root_0=svs_0_0.cVeg/3,
     C_abvstrlit_0=svs_0.cLitter,
     C_abvmetlit_0=svs_0.cLitter,
     C_blwstrlit_0=svs_0.cLitter,
     C_surfacemic_0=svs_0.cSoil,
     C_soilmic_0=svs_0.cSoil,
     C_slow_0=svs_0.cSoil
    )
)

# +
from general_helpers import autostep_mcmc, make_param_filter_func

isQualified = msh.make_param_filter_func(epa_max, epa_min)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa_0,
    filter_func=isQualified,
    param2res=param2res,
    costfunction=msh.make_weighted_cost_func(svs),
    nsimu=2000, # for testing and tuning mcmc
    #nsimu=20000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=15,   # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=50,   # default value | increase value to reduce initial step size
    K=1 # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

# +
# optimized parameter set (lowest cost function)
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields),1),axis=1)
epa_opt=msh.EstimatedParameters(*par_opt)
mod_opt = param2res(epa_opt)  
DA = mod_opt._asdict()
obs_d=svs._asdict()
print("Forward run with optimized parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(10, 50), dpi=80)
axs_DA = fig.subplots(5,1)

for a, b in enumerate(msh.Observables._fields):
    val_DA=DA[b]
    val_obs=obs_d[b]
    axs_DA[a].plot(range(len(val_DA)), val_DA, label=b+"_DA")
    axs_DA[a].plot(range(len(val_obs)),val_obs, label=b+"_obs")
    axs_DA[a].legend()

fig.savefig('solution_DA.pdf')

fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

import pandas as pd
pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('SDGVM_da_aa.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('SDGVM_da_j_aa.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('SDGVM_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('SDGVM_optimized_solutions.csv'), sep=',')
# +
it_sym_trace = msh.make_traceability_iterator(mvs,dvs,cpa,epa_opt)
ns=1500
StartVectorTrace=gh.make_StartVectorTrace(mvs)
nv=len(StartVectorTrace._fields)
res_trace= np.zeros((ns,nv))
for i in range(ns):
    res_trace[i,:]=it_sym_trace.__next__().reshape(nv)
#res_trace

import matplotlib.pyplot as plt
n=len(mvs.get_StateVariableTuple())
fig=plt.figure(figsize=(20,(n+1)*10), dpi=80)
axs=fig.subplots(n+1,2)
days=list(range(ns))


for i in range(n):
    
    ax = axs[i,0]
    #  the solution
    pos=i
    ax.plot(
        days,
        res_trace[:,i],
        label=StartVectorTrace._fields[pos],
        color='blue'
    )
    # X_p
    pos=i+n
    ax.plot(
        days,
        res_trace[:,pos],
        label=StartVectorTrace._fields[pos],
        color='red'
    )
    # X_c
    pos=i+2*n
    ax.plot(
        days,
        res_trace[:,pos],
        label=StartVectorTrace._fields[pos],
        color='yellow'
    )
    ax.legend()
    
    ax = axs[i,1]
    # RT
    pos=i+3*n
    ax.plot(
        days,
        res_trace[:,pos],
        label=StartVectorTrace._fields[pos],
        color='black'
    )
    ax.legend()
    
axs[n,0].plot(
    days,
    [msh.make_npp_func(dvs)(d) for d in days],
    label='NPP',
    color='green'
)
axs[n,0].legend()

fig.savefig('solution_Traceability.pdf')

# -


