# %load_ext autoreload
# %autoreload 2
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
from pathlib import Path
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
import matplotlib.pyplot as plt
from plotly.offline import plot

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
import general_helpers as gh
import MIP_output_helpers as moh
model_folders = [
    "kv_visit2",
    "jon_yib",
    #"Aneesh_SDGVM", # very fast drop in soil and system tot
    #"cable-pop", # has not EstimatedParameters
    #"cj_isam", # msh.numericX0 also yields a negative pool value for the last pool
    "yz_jules",
    #"kv_ft_dlem",
    #"bian_ibis2",
]
delta_t_val = 1
# test_arg_dict = gh.get_test_arg_dict(model_folders)
# t_min, t_max=gh.t_min_tmax_overlap_2(test_arg_dict, delta_t_val)
# here we assume that all models started from equilibrium at
# t_min (which is not correct)
# according to the S2 experiment they started from equilibrium
# but at different times
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
# compute the overlap of the models
td = gh.get_test_arg_dict(model_folders)
# we want the equilibrium to be computed from Temperatur
# and Moisture data in the spring. Since all models start
# in January we start 120 days after the model start of all
# models in days since THIS MODELS start
start_shift = 120
start_sAD, stop_sAD = gh.t_min_tmax_overlap_2(td, delta_t_val, start_shift=start_shift)

def timelines_from_model_folder(mf):
    test_args = gh.test_args_2(mf)
    msh = gh.msh(mf)
    mvs = gh.mvs_2(mf)
    state_vector = mvs.get_StateVariableTuple()
    # load global mean vars
    target_path = Path(mf).joinpath("global_means")
    data_path=Path(gh.confDict(mf)["dataPath"])
    svs, dvs = msh.get_global_mean_vars(data_path, target_path, flash_cache=False)
    #cpa = test_args.cpa
    #epa = test_args.epa_opt
    #dvs = test_args.dvs
    
    ##read parameters for data-assimilation
    #tr_path = Path(mf).joinpath(
    #    "data_assimilation_parameters_from_test_args"
    #)
    #cpa = msh.Constants(
    #    **h.load_dict_from_json_path(tr_path.joinpath("cpa.json"))
    #)
##
    #epa_min, epa_max, epa_0 = tuple(
    #    map(
    #        lambda p: msh.EstimatedParameters(
    #            **h.load_dict_from_json_path(p)
    #        ),
    #        [
    #            tr_path.joinpath(f"{s}.json")
    #            for s in ["epa_min", "epa_max", "epa_0"]
    #        ],
    #    )
    #)
    #perform da (or read results from cache)
    dir_path = Path(mf).joinpath("output")
    cp, Cs, Js, epa_opt = msh.gm_da_res_1(output_cache_path=dir_path)
    #read the dataassimilation results from a cache directory
    #CP = import_module( f"{msh.model_mod}.CachedParameterization" ).CachedParameterization
    CP = msh.CachedParameterization
    cp = CP.from_path(dir_path)
    X_0 = np.array(
        [cp.X_0_dict[str(s)] for s in mvs.get_StateVariableTuple()]
    )
    func_dict = cp.func_dict
    par_dict = cp.parameter_dict
    stride = 1  # this does not affect the precision of the iterator but of the averages
    # but makes it more effiecient (the values in between the strides
    func_dict = msh.make_func_dict(dvs, cpa=cpa, epa=epa)

    # Every model has it's own timeline counted from start_shift in days 
    # till 
    n_days=len(dvs[0])*30

    
    vals = gh.all_timelines_starting_at_steady_state(
        mvs,
        func_dict,
        par_dict,
        #dvs,
        #cpa,
        #epa,
        t_min=start_shift,
        index_slice=slice(0, int((n_days-start_shift)/delta_t_val), stride),
        delta_t_val=delta_t_val,
    )
    return vals


def yearly_averages(vals):
    n_days = vals.t.shape[0]
    step = int(360 / delta_t_val)
    parts = hr.partitions(0, n_days, step)
    return vals.averaged_values(parts)


all_values2 = {mf : timelines_from_model_folder(mf) for mf in model_folders}
all_averaged_values2 = {mf : yearly_averages(vals) for mf,vals in all_values2.items()}

all_values2['kv_visit2'].t.shape

all_averaged_values2['kv_visit2'].system_tot


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
        td_AD = gh.td_AD(gh.test_args_2(mf).start_date)
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
fig.savefig("test2_fine.pdf")
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
        td_AD = gh.td_AD(gh.test_args_2(mf).start_date)
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
fig.savefig("test2_yearly.pdf")

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
fig.savefig("test2_yearly.pdf")

1* 148940000* 1000000 * 0.000000000001



