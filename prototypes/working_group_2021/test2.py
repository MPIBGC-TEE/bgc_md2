# %load_ext autoreload
# %autoreload 2
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
from trendy9helpers import general_helpers as gh
#import MIP_output_helpers as moh
model_folders = [
    "kv_visit2",
    #"jon_yib",
    #"Aneesh_SDGVM", # very fast drop in soil and system tot
    #"cable-pop", # has not EstimatedParameters
    #"cj_isam", # msh.numericX0 also yields a negative pool value for the last pool
    #"yz_jules",
    #"kv_ft_dlem",
    #"bian_ibis2",
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
# we want the equilibrium to be computed from Temperatur
# and Moisture data in the spring. Since all models start
# in January we start 120 days after the model start of all
# models in days since THIS MODELS start
start_shift = 120
start_sAD, stop_sAD = gh.t_min_max_overlap_gm(model_folders, delta_t_val, start_shift=start_shift)

start_sAD,stop_sAD

import bgc_md2.helper as h
h.date.days_per_month


# +
def timelines_from_model_folder(mf):
    msh = gh.msh(mf)
    mvs = gh.mvs(mf)
    nupa=mvs.get_NumericParameterization()
    # The iterator makes timesteps with length delta_t_val.
    # iteration 0 refers to t_0 
    # iteration 1 refers to t_0 + 1
    # The iterator result supports index slicing in the same way as a list [start,stop,stride]
    # All three numbers are integers 
    # This means the first reported value will be the iterator result after start iterations 
    # with timestep delta_t_val
    n_days=h.date.days_per_month* msh.n_months()
    start_index=0
    stride = 1   
    index_slice=slice(start_index, int((n_days-start_shift)/delta_t_val), stride)


    vals = gh.all_timelines_starting_at_steady_state(
        mvs,
        nupa.func_dict,
        nupa.par_dict,
        t_min=start_shift,
        index_slice=index_slice,
        delta_t_val=delta_t_val,
    )
    return vals

all_values2 = {mf : timelines_from_model_folder(mf) for mf in model_folders}

# +
## test plot for visit
#mf="kv_visit2"
mf="jon_yib"
model_mod=f'bgc_md2.models.{mf}'
#
mvs = import_module(f"{model_mod}.source").mvs
fig=plt.figure()
ax=fig.subplots(1,1)
g=mvs.closure_graph_plot(ax)

#smr = mvs.get_SmoothModelRun()
#smr.initialize_state_transition_operator_cache(lru_maxsize=None)
#start_mean_age_vec = mvs.get_NumericStartMeanAgeVector()
#sv = mvs.get_StateVariableTuple()
#n_pools = len(sv)
#vals=values2["kv_visit2"]
#fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
#axs = fig1.subplots(n_pools, 2)
#dpy=365.25
#cty=cut_times/dpy
#for i in range(n_pools):
#    ax = axs[i, 0]
#    ax.plot(cty, sol_arr[:n_steps, i], label="sol")
#    ax.plot(vals.t/dpy, vals.X[:n_steps, i], label="bit")
#    ax.legend()
#    ax.set_title(f"{sv[i]} solution")
#    
#    ax = axs[i, 1]
#    ax.plot(cty, m_a_arr[:n_steps, i]/dpy)
#    ax.set_title(f"{sv[i]} mean_age")
#
#fig1.savefig(Path(mf).joinpath("poolwise.pdf"))
# -

def yearly_averages(vals):
    n_days = vals.t.shape[0]
    step = int(360 / delta_t_val)
    parts = hr.partitions(0, n_days, step)
    #print(parts)
    return vals.averaged_values(parts)
all_averaged_values2 = {mf : yearly_averages(vals) for mf,vals in all_values2.items()}



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
        c_times = [(td_AD + mt) / h.date.days_per_year for mt in vals.t]
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
    "veg_continuous_mean_btt",
    "soil_continuous_mean_btt",
    "system_RT_sum",
    "system_tot", 
    "veg_tot",   
    "soil_tot",   
]
fontsize=16
fsx=12
fsy=12
fig = plt.figure(figsize=(fsx,fsy))
fig.suptitle("Daily timelines of Transit Times and Approximations", fontsize=fontsize)
plot_time_lines_one_plot_per_model(
    value_dict=all_values2, 
    title_dict=model_names,    
    desired_keys=desired_keys,
    style_dict=style_dict,
    fig=fig,
    limit=int(5*h.date.days_per_year/delta_t_val) # 5 years
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

fig = plt.figure(figsize=(fsx,fsy))
fig.suptitle("Yearly averages of Transit Times and Approximations", fontsize=fontsize)
plot_time_lines_one_plot_per_model(
    value_dict=all_averaged_values2,
    title_dict=model_names,
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
    hspace=0.3
)
fig.savefig("test2_yearly.pdf")


