%load_ext autoreload
%autoreload 2
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
import pathlib
import inspect
import shutil
import matplotlib.pyplot as plt
import numpy as np
from unittest.case import TestCase, skip
from pathlib import Path
from importlib import import_module
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
import general_helpers as gh
import MIP_output_helpers as moh

# model_folders = self.model_folders

model_folders = [
    "kv_visit2",
    "jon_yib",
    #"Aneesh_SDGVM", # very fast drop in soil and system tot
    #"cable-pop", # has not EstimatedParameters
    #"cj_isam", # msh.numericX0 also yields a negative pool value for the last pool
    "yz_jules",
    "kv_ft_dlem",
    "bian_ibis2",
]
delta_t_val = 15
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


def syncronized_timelines_from_model_folder(mf):
    test_args = gh.test_args(mf)
    msh = gh.msh(mf)
    mvs = gh.mvs(mf)
    state_vector = mvs.get_StateVariableTuple()
    cpa = test_args.cpa
    epa = test_args.epa_opt
    dvs = test_args.dvs
    stride = 1  # this does not affect the precision of the iterator but of the averages
    # but makes it more effiecient (the values in between the strides
    func_dict = msh.make_func_dict(mvs, dvs, cpa, epa)

    # Every model has it's own timeline counted from 0 in days starting
    # from the first day of where IT'S data is available

    # To compare the output of one model to the simultaneous output of
    # another model to compute the indices of the iterator timesteps
    # for the appropriate common time tc
    start, stop = gh.min_max_index(
        test_args, delta_t_val, start_sAD, stop_sAD, start_shift
    )
    # from IPython import embed; embed()
    vals = gh.all_timelines_starting_at_steady_state(
        mvs,
        func_dict,
        dvs,
        cpa,
        epa,
        # we want the equilibrium to be computed from Temperatur
        # and Moisture data in the spring. Since all models start
        # in January we start 120 days after the model start of all
        # models in days since THIS MODELS start
        t_min=start_shift,
        # - t_max in THIS MODELS' time  is (len(dvs[0]))* 30 (in days),
        # - Step 0 of the iterator refers to t_min
        # The maximum index of the iterator < t_max is therefore (t_max-t_min)/delta_t_val
        index_slice=slice(start, stop, stride),
        delta_t_val=delta_t_val,
    )
    return vals


def yearly_averages(vals):
    n_days = vals.t.shape[0]
    step = int(360 / delta_t_val)
    parts = hr.partitions(0, n_days, step)
    return vals.averaged_values(parts)

all_values = {mf : syncronized_timelines_from_model_folder(mf) for mf in model_folders}
all_averaged_values = {mf : yearly_averages(vals) for mf,vals in all_values.items()}
value_dicts = [all_values, all_averaged_values]

###########################################################################################
def plot_time_lines_one_plot_per_property(
    value_dict,
	 title_dict,
	 y_label_dict,
	 desired_keys,
	 style_dict,
	 fig,
     limit=None
    ):
    title_keys = title_dict.keys()
    axs = fig.subplots(len(desired_keys), 1)#, sharex=True)
    for mf, vals in value_dict.items():
        # from IPython import embed; embed()
        # transform the times of the individual iterators back to
        # the common format (days since aD and then to years)
        td_AD = gh.td_AD(gh.test_args(mf).start_date)
        c_times = [(td_AD + mt) / 360 for mt in vals.t]
        for i,key in enumerate(desired_keys):
            ax = axs[i]
            ax.set_title(title_dict[key] if key in title_keys else key)
            y=vals[key]
            if limit is not None:
                y=y[0:limit]
                c_times=c_times[0:limit]
            ax.plot(
                c_times,
                y,
                **style_dict[mf],
                label=model_names[mf],
                # marker=marker_dict["system"],
            )
            #ax.legend()



    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        #loc='upper center'
        loc='upper left'
    )


# for caching all values (dont do this if you have  RAM shortage}
# take some minutes to compute

model_cols = {
    "JULES": "blue",
    "VISIT": "orange",
    "YIBs": "green",
    "DLEM": "red",
    "SDGVM": "yellow",
    "ISAM": "purple",
    "IBIS": "magenta",
    "OCN": "teal",
    "CABLE": "violet",
}
title_dict={
    "system_continuous_mean_btt": "System backward transit time",
    "veg_continuous_mean_btt": "Vegetation subsystem backward transit time",
    "soil_continuous_mean_btt": "Soil subsystem backward transit time",
    "system RT_sum": "System $\sum_i (RT)_i$",
}

y_label_dict={
    "system_continuous_mean_btt": "days",
    "veg_continuous_mean_btt": "days",
    "soil_continuous_mean_btt": "days",
    "system RT_sum": "days",
}

desired_keys=[
    "x",
    "x_dot",
    "veg_x",
    "soil_x",
    "veg_continuous_mean_btt",
    "veg_continuous_mean_btt",
    "soil_continuous_mean_btt",
    'out_2_veg',
    'veg_2_soil',
    'veg_2_out',
    'soil_2_out',
    "system_RT_sum",
    "system_tot", 
    "veg_tot",   
    "soil_tot",

]
fontsize=16
fsx=15
fsy=25
style_dict = {mf: {'color': model_cols[model_names[mf]]} for mf in model_folders} 
fig = plt.figure(figsize=(fsx,fsy))
plot_time_lines_one_plot_per_property(
    value_dict=all_values, 
    title_dict=title_dict,
    y_label_dict=y_label_dict,
    desired_keys=desired_keys,
    style_dict=style_dict,
    fig=fig,
    limit=int(5*360/delta_t_val) # 5 years
)
fig.suptitle("Daily timelines of diagnostic model properties", fontsize=fontsize)
#fig.tight_layout()
fig.subplots_adjust(
    left=0.1,
    bottom=0.1,
    right=0.9,
    top=0.95,
    wspace=0.4,
    hspace=0.7
)
fig.savefig("test.pdf")

fig = plt.figure(figsize=(fsx,fsy))
plot_time_lines_one_plot_per_property(
    value_dict=all_averaged_values, 
    title_dict=title_dict,
    y_label_dict=y_label_dict,
    desired_keys=desired_keys,
    style_dict=style_dict,
    fig=fig
)
fig.subplots_adjust(
    left=0.1,
    bottom=0.1,
    right=0.9,
    top=0.95,
    wspace=0.4,
    hspace=0.7
)
fig.suptitle("Yearly averaged timelines of diagnostic model properties", fontsize=fontsize)
#fig.tight_layout()
fig.savefig("test_yearly.pdf")
###########################################################################################