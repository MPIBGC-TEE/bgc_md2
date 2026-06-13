# %load_ext autoreload
# %autoreload 2
from pathlib import Path


import netCDF4 as nc
import numpy as np
import dask.array as da
from tqdm import tqdm
from functools import reduce
import sys
import json
import scipy
from CompartmentalSystems.smooth_model_run import SmoothModelRun
from CompartmentalSystems.start_distributions import (
    start_age_moments_from_empty_spinup,
    start_age_moments_from_steady_state,
    start_age_moments_from_zero_initial_content,
    compute_fixedpoint_numerically,
    start_age_distributions_from_steady_state,
    start_age_distributions_from_empty_spinup,
    start_age_distributions_from_zero_initial_content,
)
import CompartmentalSystems.helpers_reservoir as hr
import unittest
import pathlib
import inspect
import shutil
import matplotlib.pyplot as plt
import numpy as np
from unittest.case import TestCase, skip
from importlib import import_module
from collections import OrderedDict, namedtuple
from sympy import (
    Symbol,
    Function,
    sympify,
    simplify,
    lambdify,
    Function,
    Symbol,
    diff,
    exp,
    diag,
)
import shutil
import matplotlib.pyplot as plt
from plotly.offline import plot
from importlib.resources import files as mod_files


from ComputabilityGraphs.CMTVS import CMTVS



from bgc_md2.resolve.mvars import (
    NumericParameterization,
    NumericSimulationTimes,
    NumericStartValueArray,
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.helper as h
import bgc_md2.resolve.computers as bgc_c

from testinfrastructure.InDirTest import InDirTest
from trendy9helpers import general_helpers as gh
#import MIP_output_helpers as moh
model_folders = [
    "kv_visit2",
    "jon_yib",
    #"Aneesh_SDGVM", # very fast drop in soil and system tot
    ##"cable-pop", # has not EstimatedParameters
    ##"cj_isam", # msh.numericX0 also yields a negative pool value for the last pool
    "yz_jules",
    ##"kv_ft_dlem",
    ##"bian_ibis2",
]
model_names = {
    "yz_jules": "JULES",
    "cable-pop": "CABLE",
    "kv_visit2": "VISIT",
    "jon_yib": "YIBs",
    "kv_ft_dlem": "DLEM",
    "Aneesh_SDGVM": "SDGVM",
    "cj_isam": "ISAM",  
    "bian_ibis2": "IBIS",
    "ORCHIDEE-V2": "OCN",
}
# here we assume that all models started from equilibrium at
# t_min (which is not correct)
# according to the S2 experiment they started from equilibrium
# but at different times
# we want the equilibrium to be computed from Temperatur
# and Moisture data in the spring. Since all models start
# in January we start 120 days after the model start of all
# models in days since THIS MODELS start
start_shift = 120
#start_sAD, stop_sAD = gh.t_min_tmax_overlap_gm(model_folders, delta_t_val, start_shift=start_shift)
delta_t_val = 15 #SDGVM needs small timestep since the rates are so high

# ### @Kostia:
# #### To run any model we only need a parameterization which consists of 
# - a parameter_dict (stored as a json file)
# - a func_dict
#     - build from a dvs cached as netcdf file and  
#     - for some models (jules) and additional func_dict_param_dict to build the func_dict) 
# - a start_value_dict (stored as a json file)
# - a model specific class called CachedParameterization which handles reading, and writing 
#   the files and building the more abstract ingredients:
#   - parameter_dict
#   - func_dict and 
#   - start_value_dict 
#   
#   from these files 
#
# The former distinction in epa and cpa was specific to a certain kind of data assimilation, which is irrelevant for running the model, and had to be abandoned since we now have the
# possibility of different da schemes per model (and therefore different sets of parameters
# we consider constant or to be estimated)
#
# #### How to get such a parameterization for the models?
# There are three approaches:
# 1. from data assimilation 
#    We specify a da module (implemented in the trendy9helpers package) which implies a 
#    certain definition of estimated parameters and constants and assumes the existence of  
#    directories and files containing the such parameters including hyper parameters for the 
#    mcmc runs. 
#    (detailes explained below)
#    Although the results of different data assimilation schemes imply different optimized parameters 
#    (not only values) there is always a way to combine the constants and estimated 
#    parameters to build a parameterization. The da functions now return not only a specific epa_opt 
#    but also a parameterization.
#
# 1. from (possibly hand tuned) model specific parameter files (which assumes the existence of directories containing the json files for the parameter_dicts and X_0_dicts, and constructing the CachedParameterization instance from the global_means which are part of the trendy9helpers package.
#
# 1. from (model specific) complete CachedParameterization instances (written to disk in a named directory). This would save an extra netcdf file for the drivers and thus be independent from the trendy9helpers package, but wasting space by saving the (same) drivers multiple times. This mehtod is used for examples for the core framework where the dependency of trendy9helpers is to be avoided...
#   
#
# Approaches 1 and 2 (combined with the fiddeling of either the parameters for the da (1) or directly with the model parameters (2) are intended to find a reasonable parameterization,
# whereas (3) is intended to store the result of such fiddeling independently from the data
# used to get it. 
#
# The following code demonstrates all approaches and also contains some 
# diagnostic plots to check that the parametrization in question is sane.
# So the first part is similar to the old "inspect_model.py".
# You can easily remove it from here, put it in an extra notebook or script(s), and keep only the essential second part as test3.py

# +
# variant 1.) Model parameters from data assimilation.
# chose da (d"ata "a"ssimilation") module of the trendy9helpers package to be used.
# da_1 exists for every model as a submodule in the trendy9helpers package under the model folder
# da_2 exists only for kv_visit2. 
# the parset name designates a subfolder under 
# trendy9helpers/model_name/da_scheme/
# and contains two directories "in" and "out" 
# (with "out" appearing after running da)
# e.g. for visit trendy9helpers/src/trendy9helpers/kv_visit/da_2/par_1/in or
# e.g. for visit trendy9helpers/src/trendy9helpers/kv_visit/da_2/par_2 /in
# it represents different  parameters (start, min, max, hyper-parameters ) 
# for an mcmc run
    
# set a directory where to create all the foldesr and store all the files  (parameters, plots..)
p = Path(".")
print(p.absolute())
mf="kv_visit2"
da_scheme="da_2"
par_dir="par_1"
# this implies the existence of a directory p/da_1/par_1/in/
# containing epa0.json, epa_min.json, epa_max.json, cpa.json, hyper.json
Cs, Js, epa_opt, test_cp1 = gh.gm_da_from_folder_names(p, mf, da_scheme, par_dir)
#lets also store the parameterization for later in a dictionary, that we will expand later
all_cps={mf: test_cp1}
# -


# some diagnostic plots for the data assimilation 
fig=plt.figure(figsize=(15,15))
axs=fig.subplots(1,2)
ax=axs[0]
ax.plot(Js[0,:])
ax.set_title("#iterations over #accepted")
ax=axs[1]
ax.plot(Js[1,:])
ax.set_title("costfunction values over #accepted")
n_par=Cs.shape[0]
n_par
fig=plt.figure(figsize=(15,n_par*15))
axs=fig.subplots(n_par)
for i in range(n_par):
    ax=axs[i]
    ax.plot(Cs[i,:])
    ax.set_title(f"accepted {epa_opt._fields[i]}")

mf,da_scheme

# +
# look at the output of the param2res function (for da_2) for epa0 and epa_opt
# and compare it to the observations
cpa = gh.da_mod(mf,da_scheme).FreeConstants(
    **h.load_dict_from_json_path(gh.da_param_path(p,mf,da_scheme,par_dir).joinpath("cpa.json"))
)
epa_0 = gh.da_mod(mf,da_scheme).EstimatedParameters(
    **h.load_dict_from_json_path(gh.da_param_path(p,mf,da_scheme,par_dir).joinpath("epa_0.json"))
)
svs, dvs = gh.msh(mf).get_global_mean_vars(gh.data_path(mf), gh.target_path(p,mf), flash_cache=False)
param2res= gh.da_mod(mf,da_scheme).make_param2res_sym(gh.mvs(mf),cpa,dvs,svs)
epa_opt = gh.da_mod(mf,da_scheme).EstimatedParameters(
        **h.load_dict_from_json_path(
            gh.output_cache_path(p,mf,da_scheme,par_dir).joinpath("epa_opt.json")))   
sim_0 = param2res(epa_0)
sim_opt =  param2res(epa_opt)

fig = plt.figure(figsize=(10,50))
axs=fig.subplots(len(svs._fields),1)


for ind,f in enumerate(svs._fields):
    val_sim_0=sim_0.__getattribute__(f)
    val_sim_opt=sim_opt.__getattribute__(f)
    val_obs=svs.__getattribute__(f)
    axs[ind].plot(range(len(val_sim_0)),val_sim_0,label=f+"_sim_0")
    axs[ind].plot(range(len(val_sim_opt)),val_sim_opt,label=f+"_sim_opt")
    axs[ind].plot(range(len(val_obs)),val_obs,label=f+"_obs")
    axs[ind].legend()
    
fig.savefig(p.joinpath(mf,'param2res.pdf'))


# +
# read hand-tuned Parameter files from named directories.
# for mf in model_folders:
def cp_from_parameter_dir(p,mf,sub_dir_path):
    CP=import_module(f"bgc_md2.models.{mf}.CachedParameterization").CachedParameterization
    htpi=p.joinpath(mf,sub_dir_path)
    svs,dvs=gh.msh(mf).get_global_mean_vars(gh.data_path(mf))
    cp=CP(
        parameter_dict=CP.parameter_dict_from_path(htpi),
        drivers=dvs,
        X_0_dict=CP.X_0_dict_from_path(htpi),
        func_dict_param_dict=CP.func_dict_param_dict_from_path(htpi)
    )
    return cp

test_cp2A=cp_from_parameter_dir(p,"kv_visit2",Path("hand_tuned_1").joinpath("in"))
test_cp2A.X_0_dict,test_cp2A.parameter_dict
# +
def cp_from_model_dir(mf):
    # alternatively method 2B reads a complete Parameterization which 
    # includes a Drivers.nc file.
    # (here from the folder in which the source.py file resides)
    # this is usefull for ONE example complete parameterization
    # and does not need the trendy9helper package.

    cpp=mod_files("bgc_md2.models").joinpath(mf,"parameterization_from_test_args")
    #cp.write(cpp) #uncomment if the parameterization should be available without the original driver data (which would be duplicated)
    CP=import_module(f"bgc_md2.models.{mf}.CachedParameterization").CachedParameterization
    return CP.from_path(cpp)

test_cp2B=cp_from_model_dir("kv_visit2")
# -
test_cp2B



# +
# sanity checke (reproducing observables)
# make some basic test plots ensuring that we don't have negative pool values 
# and so on.
#cp=test_cp1
#cp=test_cp2A
cp=test_cp2A

mvs=import_module(f"bgc_md2.models.{mf}.source").mvs
synth_obs=gh.msh(mf).synthetic_observables(
   mvs,
   np.array([cp.X_0_dict[v] for v in mvs.get_StateVariableTuple()]),
   cp.parameter_dict,
   cp.func_dict,
   dvs
)


fig = plt.figure(figsize=(10,50))
axs=fig.subplots(len(svs._fields),1)
for ind,f in enumerate(svs._fields):
   val_sim=synth_obs.__getattribute__(f)
   val_obs=svs.__getattribute__(f)
   axs[ind].plot(range(len(val_sim)),val_sim,label=f+"_sim")
   axs[ind].plot(range(len(val_obs)),val_obs,label=f+"_obs")
   axs[ind].legend()
   
fig.savefig(p.joinpath(mf,'sythetic_observables.pdf'))
# -

# sanety check solutions 
import bgc_md2.resolve.mvars as mvars
mvs.provided_mvar_types
mvs=mvs.remove(
   [
       mvars.NumericParameterization,
       mvars.StartConditionMaker,
       mvars.NumericSimulationTimes
   ]
)
dpy = h.date.days_per_year
dpm = h.date.days_per_month
times=np.arange(0,dpm*len(cp.drivers[0]),dpm/2)
times

td_AD = h.date.days_since_AD(gh.msh(mf).start_dt())
ad_days= td_AD + times
ad_times= ad_days / dpy 
mvs=mvs.update(
   {
       mvars.NumericParameterization(
           par_dict=cp.parameter_dict,
           func_dict=cp.func_dict
       ),
       mvars.NumericStartValueDict(cp.X_0_dict),
       mvars.NumericSimulationTimes(times)
   }    
)
sv = mvs.get_StateVariableTuple()
n_pools = len(sv)
fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
axs = fig1.subplots(n_pools, 1)
sol_arr2 = mvs.get_NumericSolutionArray()
#vcsv = mvs.get_VegetationCarbonStateVariableTuple()
#veg_m_a_arr2 = mvs.get_NumericVegetationCarbonMeanAgeSolutionArray()
#start_ind=200
start_ind=0
for i in range(n_pools):
   ax = axs[i,]
   ax.plot(ad_times, sol_arr2[:, i], label="sol")
   ax.legend()
   ax.set_title(f"{sv[i]} solution")






def timelines(mf,cp):
   stride = 1  # this does not affect the precision of the iterator but of the averages
   ## but makes it more effiecient (the values in between the strides
#
   ## Every model has it's own timeline counted from start_shift in days 
   ## till 
   n_days = 30 * gh.msh(mf).n_months()
#
   
   vals = gh.all_timelines_starting_at_steady_state(
       gh.mvs(mf),
       cp.func_dict,
       cp.parameter_dict,
       t_min=start_shift,
       index_slice=slice(0, int((n_days-start_shift)/delta_t_val), stride),
       delta_t_val=delta_t_val,
   )
   return vals


all_values2 = {
   mf : timelines(mf,cp_from_parameter_dir(p,mf,Path("hand_tuned_1").joinpath("in"))) 
   for mf in model_folders
}

# +
# you could also handarrange which kind of parameterisation is to be used for which model

#all_values3 = {mf : timelines(mf,gh.gm_cp_from_folder_names(p,mf,da_schemes[mf],parset_names[mf])) for mf,tup in da_res.items()}
# -

all_values2[mf]['t']

# +
# just plot the values for checking just don't run the cell or comment it
# if you dont want to 
# we see that Aneeshs model is not in good shape
mf='Aneesh_SDGVM'
Xs=all_values2[mf]['X']
mvs=gh.mvs(mf)
svt=mvs.get_StateVariableTuple()
n=len(svt)
fig=plt.figure(figsize=(10,n/2*10))
axs=fig.subplots(n)
td_AD = h.date.days_since_AD(gh.msh(mf).start_dt())
c_times = [(td_AD + mt) / 360 for mt in all_values2[mf]['t']]
sl=slice(0,None,None)
for i in range(n):
   ax=axs[i]
   ax.plot(c_times[sl], Xs[sl,i])
   ax.set_title(str(svt[i]))
Xs.shape,len(c_times)    



# -

obs_0._fields


def yearly_averages(vals):
   n_days = vals.t.shape[0]
   step = int(365.25 / delta_t_val)
   parts = hr.partitions(0, n_days, step)
   return vals.averaged_values(parts)

#- 
# If you can afford the memory you can cache all the averages
# (usully this is not necessary)
all_averaged_values2 = {mf : yearly_averages(vals) for mf,vals in all_values2.items()}

#all_values2['kv_visit2'].t.shape
#
#all_averaged_values2['kv_visit2'].system_tot
#
#
from bgc_md2 import helper as h
def plot_time_lines_one_plot_per_model(
   value_dict,
   title_dict,
   desired_keys,
   style_dict,
   fig,
   limit=None
):
   title_keys = title_dict.keys()
   axs = fig.subplots(len(value_dict.keys()), 1)  # Different timelines no sharing, sharex=True)
   for i, mf in enumerate(value_dict.keys()):
       vals=value_dict[mf]
       try:
           ax=axs[i]
       except:
           IndexError
           ax=axs

       ax.set_title(title_dict[mf] if mf in title_keys else mf)
       # from IPython import embed; embed()
       # transform the times of the individual iterators back to
       # the common format (days since aD and then to years)
       td_AD = h.date.days_since_AD(gh.msh(mf).start_dt())
       c_times = [(td_AD + mt) / 360 for mt in vals.t]
       for i,key in enumerate(desired_keys):
           y=vals[key]
           if limit is not None:
               y=y[0:limit]
               c_times=c_times[0:limit]
           ax.plot(
               c_times,
               y,
               label=key,
               **style_dict[key]
           )
   
   handles, labels = ax.get_legend_handles_labels()
   fig.legend(
       handles,
       labels,
       #loc='upper center'
       loc='upper left'
   )

sub_system_cols= {
  "veg" : "green",
  "soil" : "brown",
  "system" : "black",
}

marker_dict = {
   "RT_sum": "*",
   "continuous_mean_btt": "+", 
   "tot": "o"
}
style_dict = {
   f"{k1}_{k2}": {'color': v1, 'marker': v2} 
   for k1,v1 in sub_system_cols.items()
   for k2,v2 in marker_dict.items()
}
desired_keys = [
   "system_continuous_mean_btt",
   #"veg_continuous_mean_btt",
   #"soil_continuous_mean_btt",
   "system_RT_sum",
   "system_tot", 
   #"veg_tot",   
   #"soil_tot",   
]
fontsize=16
fsx=15
fsy=25
fig = plt.figure(figsize=(fsx,fsy))
fig.suptitle("Daily timelines of Transit Times and Approximations", fontsize=fontsize)
plot_time_lines_one_plot_per_model(
   value_dict=all_values2, 
   title_dict=model_names,    
   desired_keys=desired_keys,
   style_dict=style_dict,
   fig=fig,
   limit=int(5*360/delta_t_val) # 5 years
)
fig.subplots_adjust(
   left=0.1,
   bottom=0.1,
   right=0.9,
   top=0.95,
   wspace=0.4,
   hspace=0.3
)
fig.savefig(p.joinpath("fine.pdf"))
style_dict

# +
from copy import copy
all_averaged_values3=all_averaged_values2.copy()

all_averaged_values3['kv_visit2'].system_tot[0:5]
all_averaged_values4 = {}
#print({mf : vals for mf,vals in all_averaged_values3.items()})
for mf,vals in all_averaged_values3.items():
   if mf=='kv_visit2': 
       start=60
       end=159
   else: 
       start=220
       end=319        
   vals.system_tot=vals.system_tot[start:end]
   vals2={'system_tot':vals.system_tot[start:end]/365, # in yr
          'system_RT_sum':vals.system_RT_sum[start:end]/365, # in yr
          'system_continuous_mean_btt':vals.system_continuous_mean_btt[start:end]/365, # in yr
          't':vals.t[start:end],
          'x':vals.x[start:end]*148.94, # Pg C global
          'x_c':vals.x_c[start:end]*148.94, # Pg C global
          'x_approx':vals.u[start:end]*vals.system_tot[start:end]*148.94, # Pg C global
         }
   print(type(vals))
   print(type(vals2))
   all_averaged_values4[mf] = vals2


# -

def plot_time_lines_one_plot_per_model2(
   value_dict,
   title_dict,
   desired_keys,
   style_dict,
   fig,
   limit=None
):
   title_keys = title_dict.keys()
   axs = fig.subplots(len(value_dict.keys()), 1)  # Different timelines no sharing, sharex=True)
   for i, mf in enumerate(value_dict.keys()):
       vals=value_dict[mf]
       try:
           ax=axs[i]
       except:
           IndexError
           ax=axs
       ax.grid()
       ax.set_title(title_dict[mf] if mf in title_keys else mf)
       # from IPython import embed; embed()
       # transform the times of the individual iterators back to
       # the common format (days since aD and then to years)
       td_AD = h.date.days_since_AD(gh.msh(mf).start_dt())
       c_times = [(td_AD + mt) / 360 for mt in vals['t']]
       for i,key in enumerate(desired_keys):
           y=vals[key]
           if limit is not None:
               y=y[0:limit]
               c_times=c_times[0:limit]
           if key in ['system_RT_sum','system_tot','system_continuous_mean_btt']:
               ax.set_ylabel('$yr$')
           elif key in ['x', 'x_c', 'x_approx']:
               ax.set_ylabel('$Pg$ $C$')
           ax.plot(
               c_times,
               y,
               label=key,
               **style_dict[key]
           )
   
   handles, labels = ax.get_legend_handles_labels()
   fig.legend(
       handles,
       labels,
       #loc='upper center'
       loc='upper left'
   )


for mf in ['kv_visit2']:
   print(mf)
   slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(all_averaged_values4[mf]['system_RT_sum'], all_averaged_values4[mf]['system_tot'])
   print("system_RT_sum vs system_tot")
   print("r2 = "+str(r_value*r_value))
   print("std_err = "+str(std_err))
   print("slope = "+str(slope))
   print("intercept = "+str(intercept))
   slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(all_averaged_values4[mf]['system_RT_sum'], all_averaged_values4[mf]['system_continuous_mean_btt'])
   print("system_RT_sum vs system_continuous_mean_btt")
   print("r2 = "+str(r_value*r_value))
   print("std_err = "+str(std_err))
   print("slope = "+str(slope))
   print("intercept = "+str(intercept))
   slope3, intercept3, r_value3, p_value3, std_err3 = scipy.stats.linregress(all_averaged_values4[mf]['system_tot'], all_averaged_values4[mf]['system_continuous_mean_btt'])
   print("system_tot vs system_continuous_mean_btt")
   print("r2 = "+str(r_value3*r_value3))
   print("std_err = "+str(std_err3))
   print("slope = "+str(slope3))
   print("intercept = "+str(intercept3))

fig = plt.figure(figsize=(15, 5))
axs = fig.subplots(1, 3)
#create basic scatterplot
ax0=axs[2]
x=all_averaged_values4[mf]['system_RT_sum']
y=all_averaged_values4[mf]['system_tot']
ax0.plot(x, y, 'o')
#obtain m (slope) and b(intercept) of linear regression line
m, b = np.polyfit(x, y, 1)
#add linear regression line to scatterplot 
ax0.plot(x, m*x+b)
ax0.grid()
ax0.set_title('Residence time vs Turnover time')
ax0.set_xlabel('$RT,$ $yr$')
ax0.set_ylabel('$τ,$ $yr$')
#create basic scatterplot
ax1=axs[1]
x=all_averaged_values4[mf]['system_RT_sum']
y=all_averaged_values4[mf]['system_continuous_mean_btt']
ax1.plot(x, y, 'o')
#obtain m (slope) and b(intercept) of linear regression line
m, b = np.polyfit(x, y, 1)
#add linear regression line to scatterplot 
ax1.plot(x, m*x+b)
ax1.grid()
ax1.set_title('Residence time vs Transit time')
ax1.set_xlabel('$RT,$ $yr$')
ax1.set_ylabel('$TR,$ $yr$')
#create basic scatterplot
ax2=axs[0]
x=all_averaged_values4[mf]['system_tot']
y=all_averaged_values4[mf]['system_continuous_mean_btt']
ax2.plot(x, y, 'o')
#obtain m (slope) and b(intercept) of linear regression line
m, b = np.polyfit(x, y, 1)
#add linear regression line to scatterplot 
ax2.plot(x, m*x+b)
ax2.grid()
ax2.set_title('Turnover time vs Transit time')
ax2.set_xlabel('$τ,$ $yr$')
ax2.set_ylabel('$TR,$ $yr$')

desired_keys = [
   "system_continuous_mean_btt",
   "system_RT_sum",
   "system_tot",      
]
style_dict = {'veg_RT_sum': {'color': 'green', 'marker': '*'},
'veg_continuous_mean_btt': {'color': 'green', 'marker': '+'},
'veg_tot': {'color': 'green', 'marker': 'o'},
'soil_RT_sum': {'color': 'brown', 'marker': '*'},
'soil_continuous_mean_btt': {'color': 'brown', 'marker': '+'},
'soil_tot': {'color': 'brown', 'marker': 'o'},
'system_RT_sum': {'color': 'orange', 'marker': '.'},
'system_continuous_mean_btt': {'color': 'green', 'marker': '.'},
'system_tot': {'color': 'blue', 'marker': '.'}}
fig = plt.figure(figsize=(fsx,fsy))
fig.suptitle("Yearly averages of Transit Times and Approximations", fontsize=fontsize)
plot_time_lines_one_plot_per_model2(
   value_dict=all_averaged_values4, #all_averaged_values2,
   title_dict=model_names,
   desired_keys=desired_keys,
   style_dict=style_dict,
   fig=fig,
)
fig.subplots_adjust(
   left=0.1,
   bottom=0.1,
   right=0.9,
   top=0.95,
   wspace=0.4,
   hspace=0.3
)
fig.savefig(p.joinpath("yearly.pdf"))

# +
for mf in ['kv_visit2']:
   print(mf)
   slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(all_averaged_values4[mf]['x_c'], all_averaged_values4[mf]['x_approx'])
   print("x_c vs x_approx")
   print("r2 = "+str(r_value*r_value))
   print("std_err = "+str(std_err))
   print("slope = "+str(slope))
   print("intercept = "+str(intercept))
   slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(all_averaged_values4[mf]['x_c'], all_averaged_values4[mf]['x'])
   print("x_c vs x")
   print("r2 = "+str(r_value*r_value))
   print("std_err = "+str(std_err))
   print("slope = "+str(slope))
   print("intercept = "+str(intercept))
   slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(all_averaged_values4[mf]['x_approx'], all_averaged_values4[mf]['x'])
   print("x_approx vs x")
   print("r2 = "+str(r_value*r_value))
   print("std_err = "+str(std_err))
   print("slope = "+str(slope))
   print("intercept = "+str(intercept))

# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(all_averaged_values4['jon_yib']['x_c'], all_averaged_values4['jon_yib']['x_approx'])
# print(r_value)
# print(std_err)
# print(slope)
# print(intercept)
# -

fig = plt.figure(figsize=(15, 5))
axs = fig.subplots(1, 3)
#create basic scatterplot
ax0=axs[2]
x=all_averaged_values4[mf]['x_c']
y=all_averaged_values4[mf]['x_approx']
ax0.plot(x, y, 'o')
#obtain m (slope) and b(intercept) of linear regression line
m, b = np.polyfit(x, y, 1)
#add linear regression line to scatterplot 
ax0.plot(x, m*x+b)
ax0.grid()
ax0.set_title('C Storage Capacity vs C Approximation')
ax0.set_xlabel('$X_{approx}$, $Pg$ $C$')
ax0.set_ylabel('$X_{ecosystem},$ $Pg$ $C$')
#create basic scatterplot
ax1=axs[1]
x=all_averaged_values4[mf]['x_c']
y=all_averaged_values4[mf]['x']
ax1.plot(x, y, 'o')
#obtain m (slope) and b(intercept) of linear regression line
m, b = np.polyfit(x, y, 1)
#add linear regression line to scatterplot 
ax1.plot(x, m*x+b)
ax1.grid()
ax1.set_title('C Storage Capacity vs Ecosystem C')
ax1.set_xlabel('$X_{approx}$, $Pg$ $C$')
ax1.set_ylabel('$X_{ecosystem},$ $Pg$ $C$')
#create basic scatterplot
ax2=axs[0]
x=all_averaged_values4[mf]['x_approx']
y=all_averaged_values4[mf]['x']
ax2.plot(x, y, 'o')
#obtain m (slope) and b(intercept) of linear regression line
m, b = np.polyfit(x, y, 1)
#add linear regression line to scatterplot 
ax2.plot(x, m*x+b)
ax2.grid()
ax2.set_title('C Approximation vs Ecosystem C')
ax2.set_xlabel('$X_{approx}$, $Pg$ $C$')
ax2.set_ylabel('$X_{ecosystem},$ $Pg$ $C$')

desired_keys = [
   "x",
   "x_c",
   "x_approx"
]
style_dict = {'veg_RT_sum': {'color': 'green', 'marker': '*'},
'veg_continuous_mean_btt': {'color': 'green', 'marker': '+'},
'veg_tot': {'color': 'green', 'marker': 'o'},
'soil_RT_sum': {'color': 'brown', 'marker': '*'},
'soil_continuous_mean_btt': {'color': 'brown', 'marker': '+'},
'soil_tot': {'color': 'brown', 'marker': 'o'},
'x_c': {'color': 'orange', 'marker': '.'},
'x': {'color': 'green', 'marker': '.'},
'x_approx': {'color': 'blue', 'marker': '.'}}
fig = plt.figure(figsize=(fsx,fsy))
fig.suptitle("Yearly averages of Transit Times and Approximations", fontsize=fontsize)
plot_time_lines_one_plot_per_model2(
   value_dict=all_averaged_values4, #all_averaged_values2,
   title_dict=model_names,
   desired_keys=desired_keys,
   style_dict=style_dict,
   fig=fig,
)
fig.subplots_adjust(
   left=0.1,
   bottom=0.1,
   right=0.9,
   top=0.95,
   wspace=0.4,
   hspace=0.3
)
fig.savefig(p.joinpath("yearly2.pdf"))

1* 148940000* 1000000 * 0.000000000001




