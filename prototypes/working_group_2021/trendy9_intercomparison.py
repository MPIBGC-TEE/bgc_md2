# %matplotlib inline
# %load_ext autoreload
# %autoreload 2

""
import netCDF4 as nc
import numpy as np
import dask.array as da
from tqdm import tqdm
from functools import reduce
import sys
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
import inspect
import shutil
import matplotlib.pyplot as plt
import numpy as np
from unittest.case import TestCase, skip
from pathlib import Path
from importlib import import_module
from importlib.resources import files as mod_files
from collections import OrderedDict, namedtuple
import json
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
from trendy9helpers import general_helpers as gh
#import MIP_output_helpers as moh

# model_folders = self.model_folders
output_path=Path("trendy9_intercomparison")
output_path.mkdir(exist_ok=True)

model_folders = [
    "kv_visit2",
    "jon_yib",
    "yz_jules",
    #"Aneesh_SDGVM", # very fast drop in soil and system tot
    #"cj_isam", # msh.numericX0 also yields a negative pool value for the last pool
    #"kv_ft_dlem", # not yet converted to the new format
    #"bian_ibis2", # not yet converted to the new format
    #"cable-pop", # has not EstimatedParameters
]
delta_t_val = 15.22668056279312
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
model_cols_by_mf = {
     "yz_jules": "blue",
     "cable-pop": "violet",
     "kv_visit2": "orange",
     "jon_yib": "green",
     "kv_ft_dlem": "red",
     "Aneesh_SDGVM": "yellow",
     "cj_isam": "purple",
     "bian_ibis2": "magenta",
     "ORCHIDEE-V2": "teal",
}


# compute the overlap of the models

# we want the equilibrium to be computed from Temperatur
# and Moisture data in the spring. Since all models start
# in January we start 120 days after.
# Note:
# Since the S2 experiments claims to start in equilibrium a correct replication
# suggests 0 days
# For different models these might be different years
start_shift = 120

# now we compute the overlaping time for all maodels 
# this will be the time for which we PLOT the results TOGETHER
# so this times relative to a model independent start date
# here AD = Anno Domini = year 0
# we first have to SOLVE all IVPs from THEIR respective start times (which differ by many years)
start_sAD, stop_sAD = gh.t_min_max_overlap_gm(
    model_folders,
    delta_t_val,
    start_shift=start_shift
)

""
start_sAD,stop_sAD


""
def plot_syncronized_timeline(
        mf, 
        func, 
        start_dsAD, # start in days after AD
        stop_dsAD, # start in days after AD
        ax,
        **plt_kwargs
    ):
    msh = gh.msh(mf)
    mvs = gh.mvs(mf)
    nupa=mvs.get_NumericParameterization()
    state_vector = mvs.get_StateVariableTuple()
    
    # this does not affect the precision of the iterator but of the averages
    # but makes it more effiecient (the values in between the strides
    #stride = 15  
    stride = 5  

    # Every model has it's own timeline counted from 0 in days starting
    # from the first day of where IT'S data is available

    model_start_dsAD=h.date.days_since_AD(gh.msh(mf).start_dt())

    overlap_times_since_model_start=np.arange(start_dsAD,stop_dsAD,stride*delta_t_val) - model_start_dsAD + start_shift 
    # We add the start time which is the same for every model to the array if it is not the first value anyway.
    # since the models are supposed to start from equilibrium at their model specific start date
    first = overlap_times_since_model_start[0] 
    cond = (first == start_shift)
    
    compute_times = overlap_times_since_model_start if cond else np.append(
        np.array([start_shift]),
        overlap_times_since_model_start
    ) 
    mvs=mvs.update([NumericSimulationTimes(compute_times)])
    vals=func(mvs)
    overlap_vals= vals if cond else vals[1:]

    dpy=h.date.days_per_year
    c_times = [(model_start_dsAD + mt) / dpy for mt in overlap_times_since_model_start]
    
    ax.plot(
        c_times,
        overlap_vals,
        **plt_kwargs
    )



""
diff_d = stop_sAD - start_sAD
diff_d_plot = diff_d / 4
fig0 = plt.figure(figsize=(15,15))
ax=fig0.subplots(1,1)
for mf in model_folders:
    plot_syncronized_timeline(
        mf,
        lambda mvs: mvs.get_NumericVegetationCarbonMeanBackwardTransitTimeSolution(),
        start_sAD,
        start_sAD + diff_d_plot,
        ax,
        label=f"{model_names[mf]}_veg",
        color=model_cols_by_mf[mf],
        linestyle="dotted"
    )
    plot_syncronized_timeline(
        mf,
        lambda mvs: mvs.get_NumericSoilCarbonMeanBackwardTransitTimeSolution(),
        start_sAD,
        start_sAD + diff_d_plot,
        ax,
        label=f"{model_names[mf]}_soil",
        color=model_cols_by_mf[mf],
        linestyle="dashed"
    )
    plot_syncronized_timeline(
        mf,
        lambda mvs: mvs.get_NumericMeanBackwardTransitTimeSolution(),
        start_sAD,
        start_sAD + diff_d_plot,
        ax,
        label=f"{model_names[mf]}_system",
        color=model_cols_by_mf[mf],
        linestyle="solid"
    )
ax.legend()
ax.set_title("Carbon Mean Transit Times for Vegetation, Soil and the System")

fig0.savefig(output_path.joinpath("continuous_mean_btts.pdf"))
#from IPython import embed; embed()



""
#def syncronized_timelines_from_model_folder(mf):
#    msh = gh.msh(mf)
#    mvs = gh.mvs(mf)
#    nupa=mvs.get_NumericParameterization()
#    state_vector = mvs.get_StateVariableTuple()
#    stride = 15  # this does not affect the precision of the iterator but of the averages
#    # but makes it more effiecient (the values in between the strides
#
#    # Every model has it's own timeline counted from 0 in days starting
#    # from the first day of where IT'S data is available
#
#    # To compare the output of one model to the simultaneous output of
#    # another model to compute the indices of the iterator timesteps
#    # for the appropriate common time tc
#    start, stop = gh.min_max_index_2(
#        mf,
#	delta_t_val,
#	start_sAD,
#	stop_sAD,
#	start_shift
#    )
#    vals = gh.all_timelines_starting_at_steady_state(
#        mvs,
#        nupa.func_dict,
#        nupa.par_dict,
#        # we want the equilibrium to be computed from Temperatur
#        # and Moisture data in the spring. Since all models start
#        # in January we start 120 days after the model start of all
#        # models in days since THIS MODELS start
#        t_min=start_shift,
#        # - t_max in THIS MODELS' time  is (len(dvs[0]))* 30 (in days),
#        # - Step 0 of the iterator refers to t_min
#        # The maximum index of the iterator < t_max is therefore (t_max-t_min)/delta_t_val
#        index_slice=slice(start, stop, stride),
#        delta_t_val=delta_t_val,
#    )
#    return vals
#
#""
#
#
#""
## for caching all values (dont do this if you have  RAM shortage}
## takes some minutes to compute
#all_values = {mf : syncronized_timelines_from_model_folder(mf) for mf in model_folders}
#
#""
#mf="kv_visit2"
#stride = 1  # this does not affect the precision of the iterator but of the averages
## but makes it more effiecient (the values in between the strides
#
## Every model has it's own timeline counted from 0 in days starting
## from the first day of where IT'S data is available
#
## To compare the output of one model to the simultaneous output of
## another model to compute the indices of the iterator timesteps
## for the appropriate common time tc
#start, stop = gh.min_max_index_2(
#    mf, delta_t_val, start_sAD, stop_sAD, start_shift
#)
#start, stop
#
#""
#vals=all_values[mf]
#model_mod=f'bgc_md2.models.{mf}'
#mvs = import_module(f"{model_mod}.source").mvs
#sv = mvs.get_StateVariableTuple()
#n_pools = len(sv)
#vals.t
#
#""
#fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
#axs = fig1.subplots(n_pools, 2)
#dpy=h.date.days_per_year
#td_AD = h.date.days_since_AD(gh.msh(mf).start_dt())
#c_times = [(td_AD + mt) / dpy for mt in vals['t']]
#for i in range(n_pools):
#    ax = axs[i, 0]
#    ax.plot(c_times, vals.X[:,i], label="bit")
#    ax.legend()
#    ax.set_title(f"{sv[i]} solution")
#    
#    #ax = axs[i, 1]
#    #ax.plot(cty, m_a_arr[:n_steps, i]/dpy)
#    #ax.set_title(f"{sv[i]} mean_age")
#
#fig1.savefig(Path(mf).joinpath("poolwise.pdf"))
#
#
#""
#def yearly_averages(vals):
#    n_days = vals.t.shape[0]
#    step = int(360 / delta_t_val)
#    parts = hr.partitions(0, n_days, step)
#    return vals.averaged_values(parts)
#
## for caching all values (dont do this if you have  RAM shortage}
## takes some minutes to compute
#all_averaged_values = {mf : yearly_averages(vals) for mf,vals in all_values.items()}
#value_dicts = [all_values, all_averaged_values]
#
#
#""
#def plot_time_lines_one_plot_per_property(
#    value_dict,
#    title_dict,
#    y_label_dict,
#    desired_keys,
#    style_dict,
#    fig,
#    limit=None,
#    ensemble_mean=False # only works for arrays that have the same dimensions
#    # for all models NOT for vectors or matrices
#):
#    title_keys = title_dict.keys()
#    axs = fig.subplots(len(desired_keys), 1)#, sharex=True)
#    for mf, vals in value_dict.items():
#        # from IPython import embed; embed()
#        # transform the times of the individual iterators back to
#        # the common format (days since aD and then to years)
#        td_AD = h.date.days_since_AD(gh.msh(mf).start_dt())
#        c_times = [(td_AD + mt) / 365.25 for mt in vals['t']]
#        for i,key in enumerate(desired_keys):
#            ax = axs[i]
#            ax.set_title(title_dict[key] if key in title_keys else key)
#            y=vals[key]
#            if limit is not None:
#                y=y[0:limit]
#                c_times=c_times[0:limit]
#            ax.plot(
#                c_times,
#                y,
#                **style_dict[mf],
#                label=model_names[mf],
#                # marker=marker_dict["system"],
#            )
#            if key in y_label_dict.keys():
#                ax.set_ylabel(y_label_dict[key])
#            #ax.legend()
#
#    if ensemble_mean:
#        style_dict.update({"ensemble_mean": {'color': "gray"}} )
#        def array_maen(array_list):
#            return sum(array_list)/len(array_list) 
#        ensemble_means = {
#            key: array_maen([ value_dict[mf][key]  for mf in model_folders]) 
#            for key in desired_keys
#        }
#        def array_variance(array_list,mean_arr):
#            return sum(
#                [(arr-mean_arr)**2 for arr in array_list]
#            )/(len(array_list)-1) 
#
#        ensemble_sigma = {
#            key: np.sqrt(array_variance(
#                [ value_dict[mf][key]  for mf in model_folders],
#                ensemble_means[key]
#            )) 
#            for key in desired_keys
#        }
#        for i, key in enumerate(desired_keys):
#            ax = axs[i]
#            y =ensemble_means[key]
#            if limit is not None:
#                y = y[0:limit]
#                c_times = c_times[0:limit]
#            ax.plot(
#                c_times,
#                y,
#                **style_dict["ensemble_mean"],
#                label="ensemble mean"
#                # marker=marker_dict["system"],
#            )
#            y1 =ensemble_means[key]-ensemble_sigma[key]
#            y2 =ensemble_means[key]+ensemble_sigma[key]
#            if limit is not None:
#                y1 = y1[0:limit]
#                y2 = y2[0:limit]
#                c_times = c_times[0:limit]
#            ax.fill_between(
#                c_times,
#                y1,
#                y2,
#                **style_dict["ensemble_mean"],
#                alpha=0.2,
#                label="$\sigma$ "
#
#                # marker=marker_dict["system"],
#            )
#
#    handles, labels = ax.get_legend_handles_labels()
#    fig.legend(
#        handles,
#        labels,
#        #loc='upper center'
#        loc='upper left'
#    )
#
#
#
#title_dict={
#    "system_continuous_mean_btt": "System backward transit time",
#    "veg_continuous_mean_btt": "Vegetation subsystem backward transit time",
#    "soil_continuous_mean_btt": "Soil subsystem backward transit time",
#    "system_RT_sum": "System $\sum_i (RT)_i$, (Luo Equilibrium) Residence Time",
#    "system_tot":"System turn over time", 
#}
#
#y_label_dict={
#    "system_continuous_mean_btt": "days",
#    "veg_continuous_mean_btt": "days",
#    "soil_continuous_mean_btt": "days",
#    "system RT_sum": "days",
#}
#
#desired_keys=[
#    "x",
#    "x_dot",
#    "veg_x",
#    "soil_x",
#    "system_continuous_mean_btt",
#    "veg_continuous_mean_btt",
#    "soil_continuous_mean_btt",
#    'out_2_veg',
#    'veg_2_soil',
#    'veg_2_out',
#    'soil_2_out',
#    "system_RT_sum",
#    "system_tot", 
#    "veg_tot",
#    "soil_tot",
#]
#fontsize=16
#fsx=15
#fsy=25
#style_dict = {mf: {'color': model_cols[model_names[mf]]} for mf in model_folders}
#""
#fig = plt.figure(figsize=(fsx,10))
#plot_time_lines_one_plot_per_property(
#   value_dict=all_values, 
#   title_dict=title_dict,
#   y_label_dict=y_label_dict,
#   desired_keys=[
#       "system_continuous_mean_btt",
#       "system_RT_sum",
#       "system_tot", 
#       # "veg_continuous_mean_btt",
#       # "veg_tot",   
#       # "soil_continuous_mean_btt",
#       # "soil_tot",
#   ],
#   style_dict=style_dict,
#   fig=fig,
#   limit=int(5*360/delta_t_val), # 5 years
#   #ensemble_mean=True
#)
#fig.suptitle("""Daily timelines of different measures of carbon transit""", fontsize=fontsize)
##fig.tight_layout()
#fig.subplots_adjust(
#   left=0.1,
#   bottom=0.1,
#   right=0.9,
#   top=0.90,
#   wspace=0.4,
#   hspace=0.7
#)
#fig.savefig("test.pdf")
#
#""
#""
#fig = plt.figure(figsize=(fsx,7))
#plot_time_lines_one_plot_per_property(
#    value_dict=all_values, 
#    title_dict=title_dict,
#    y_label_dict=y_label_dict,
#    desired_keys=[
#         "veg_continuous_mean_btt",
#         #"veg_tot",   
#         "soil_continuous_mean_btt",
#         #"soil_tot",
#    ],
#    style_dict=style_dict,
#    fig=fig,
#    limit=int(5*360/delta_t_val), # 5 years
#    #ensemble_mean=True
#)
#fig.suptitle("""Daily timelines carbon transit times for vegetation and soil sub systems""", fontsize=fontsize)
##fig.tight_layout()
#fig.subplots_adjust(
#    left=0.1,
#    bottom=0.1,
#    right=0.9,
#    top=0.90,
#    wspace=0.4,
#    hspace=0.7
#)
#fig.savefig("test_veg_soil.pdf")
#
#############################################################################
## ###########################################################################
#desired_keys=[
#    "x",
#    "veg_x",
#    "soil_x",
#]
#
#def array_diff0(arr):
#    return arr-arr[0]
#
#all_averaged_values_zero = {
#    mf:{key: array_diff0(vals[key]) if key !="t" else vals[key]
#    for key in vals.keys() } 
#    for mf,vals in all_averaged_values.items()
#}
#fig = plt.figure(figsize=(fsx,9))
#plot_time_lines_one_plot_per_property(
#    value_dict=all_averaged_values_zero, 
#    title_dict=title_dict,
#    y_label_dict=y_label_dict,
#    desired_keys=desired_keys,
#    style_dict=style_dict,
#    fig=fig,
#    #limit=int(5*360/delta_t_val), # 5 years
#    ensemble_mean=True
#)
#fig.suptitle("Yearly averaged timelines of stocks spread ", fontsize=fontsize)
##fig.tight_layout()
#fig.subplots_adjust(
#    left=0.1,
#    bottom=0.1,
#    right=0.9,
#    top=0.9,
#    wspace=0.4,
#    hspace=0.7
#)
#fig.savefig("test_stock_mean.pdf")
#
#""
#fig = plt.figure(figsize=(fsx,fsy))
#plot_time_lines_one_plot_per_property(
#   value_dict=all_averaged_values, 
#   title_dict=title_dict,
#   y_label_dict=y_label_dict,
#   desired_keys=desired_keys,
#   style_dict=style_dict,
#   fig=fig
#)
#fig.subplots_adjust(
#   left=0.1,
#   bottom=0.1,
#   right=0.9,
#   top=0.95,
#   wspace=0.4,
#   hspace=0.7
#)
#fig.suptitle("Yearly averaged timelines of diagnostic model properties", fontsize=fontsize)
##fig.tight_layout()
#fig.savefig("test_yearly.pdf")
#""
#
#""
#
