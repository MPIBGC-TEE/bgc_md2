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

# +
import shutil
import re
from unittest.case import TestCase, skip
from pathlib import Path
from importlib import import_module
from importlib.resources import files as mod_files
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
from plotly.offline import plot
import matplotlib.pyplot as plt
import numpy as np

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
from CompartmentalSystems.ArrayDictResult import ArrayDictResult
import CompartmentalSystems.helpers_reservoir as hr

from ComputabilityGraphs.CMTVS import CMTVS
from testinfrastructure.InDirTest import InDirTest

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


# +
# get everything from mvs
# This is the blueprint for a notebook
# that does not depend on any data exept those
# provided in the model folder
#model_folders = [
#            "kv_visit2",
#            "jon_yib",
#            "yz_jules",
#            "Aneesh_SDGVM"
#]
#for mf in set(model_folders):
mf="kv_visit2"
#mf="jon_yib"
model_mod=f'bgc_md2.models.{mf}'

mvs = import_module(f"{model_mod}.source").mvs

smr = mvs.get_SmoothModelRun()
smr.initialize_state_transition_operator_cache(lru_maxsize=None)
start_mean_age_vec = mvs.get_NumericStartMeanAgeVector()
sv = mvs.get_StateVariableTuple()
n_pools = len(sv)
order = 1
arr, func = smr._solve_age_moment_system( order, start_mean_age_vec.reshape(1, n_pools)) 
# -

times = mvs.get_NumericSimulationTimes()
t0 = times[0]
n_steps = 2800 # for testing
#n_steps = len(times)
cut_times=times[:n_steps]
# compute the indeces for the iterator to test the resulst
stride=2
#assuming an equidistant times array
delta_t_val=(times[1]-times[0])/stride
sol_arr = arr[:n_steps,0:n_pools]
m_a_arr = arr[:n_steps,n_pools:2*n_pools]
# the colums n+1..2n are the first age_moments (mean)
## plot the continuous solution (by ODE) solver against the  iterator
## generated one.

mvs.get_NumericStartValueArray()

start_mean_age_vec

# +
from trendy9helpers import general_helpers as gh
par_dict = mvs.get_NumericParameterization().par_dict
func_dict = mvs.get_NumericParameterization().func_dict
X_0=mvs.get_NumericStartValueArray()

bit = ArrayDictResult(
    gh.traceability_iterator_internal(
        mvs, X_0, par_dict, func_dict, delta_t_val=delta_t_val, t_0=t0
    )
)
max_iter=int((cut_times[-1]-t0)/delta_t_val)
max_iter
vals=bit[0:max_iter:stride]
vals.t,cut_times

fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
axs = fig1.subplots(n_pools, 2)
dpy=364.25
cty=cut_times/dpy
for i in range(n_pools):
    ax = axs[i, 0]
    ax.plot(cty, sol_arr[:n_steps, i], label="sol")
    ax.plot(vals.t/dpy, vals.X[:n_steps, i], label="bit")
    ax.legend()
    ax.set_title(f"{sv[i]} solution")
    
    ax = axs[i, 1]
    ax.plot(cty, m_a_arr[:n_steps, i]/dpy)
    ax.set_title(f"{sv[i]} mean_age")

fig1.savefig(Path(mf).joinpath("poolwise.pdf"))

# +
sm_a_arr = smr.system_age_moment(order,start_mean_age_vec.reshape(1, n_pools))
fig1 = plt.figure(figsize=(2 * 10,  10))
axs = fig1.subplots(1, 2)
dpy=364.25
ax = axs[0]
#ax.plot(times, vals.X[:, i], label="bit")
ax.plot(cut_times/dpy, np.sum(sol_arr,axis=1), label="sol")
ax.plot(vals.t/dpy, vals.x, label="sol")
ax.legend()
ax.set_title("cumulativ solution")

ax = axs[ 1]
ax.plot(times/dpy, sm_a_arr/dpy)
ax.set_title(f"system mean_age")

fig1.savefig(Path(mf).joinpath("system_mean_age.pdf"))
# -

fig2 = plt.figure(figsize=(10, 10))
axs2 = fig2.subplots(1, 1)
mean_btts = smr.backward_transit_time_moment(
    order=1, start_age_moments=start_mean_age_vec.reshape(1, n_pools)
)
ax = axs2
ax.plot(times/dpy, mean_btts, label="mean backward transit time")
#ax.plot(
#    times, vals.system_RT_sum, label="$\sum_i (RT)_i$"
#)  # steady state transit times
#ax.plot(
#    times, vals.rt, label="rt of surrogate one pool system"
#)  # steady state transit times
ax.legend()
fig2.savefig(Path(mf).joinpath("system.pdf"))

# +
# construct a function p that takes an age array "ages" as argument
# and gives back a three-dimensional ndarray (ages x times x pools)
# from the a array-valued function of a single age a_dens_function
srm = mvs.get_SmoothReservoirModel()
a_dens_function, X_fix = start_age_distributions_from_steady_state(
    srm, t0=t0, parameter_dict=par_dict, func_set=func_dict, x0=X_0
)
p = smr.pool_age_densities_func(a_dens_function)
ages = np.linspace(
    0,
    (np.array(start_mean_age_vec, dtype=float).reshape(-1)).max() * 2,
    2,
)  # 21)
age_densities = p(ages)

for n in range(srm.nr_pools):
    max_ind = np.argmin(ages < start_mean_age_vec[n] * 2)
    fig = smr.plot_3d_density_plotly(
        "age distribution pool {0}".format(sv[n]),
        age_densities[0:max_ind, :, n],
        ages[0:max_ind],
    )
    # plot the computed start age density for t0 on top
    fig.add_scatter3d(
        x=np.array([-t0 for a in ages]),
        y=np.array([a for a in ages]),
        z=np.array([a_dens_function(a)[n] for a in ages]),
        mode="lines",
        line=dict(color="#FF0000", width=15),
    )
    smr.add_line_to_density_plot_plotly(
        fig,
        data=m_a_arr[:, n],
        color="#FF0000",
        name="mean age",
        time_stride=1,
        on_surface=True,
        bottom=True,
        legend_on_surface=True,
        legend_bottom=False,
    )

    plot(
        fig,
        filename=str(
            "age_distribution_{0}.html".format(sv[n])
        ),
        auto_open=False,
    )

btt_dens = smr.backward_transit_time_density(age_densities)
fig_btt = smr.plot_3d_density_plotly(
    "backward_transit_time_density_steady_state",
    btt_dens,
    ages,
    y_label="transit time",
)
smr.add_line_to_density_plot_plotly(
    fig_btt,
    data=mean_btts,
    color="#FF0000",
    name="mean age",
    time_stride=1,
    on_surface=True,
    bottom=True,
    legend_on_surface=True,
    legend_bottom=False,
)
plot(
    fig_btt,
    filename="btt_distribution.html",
    auto_open=False,
)
# -




