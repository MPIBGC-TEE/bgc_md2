# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Compute 14C for a single site in northern Sweden based on CARDAMOM output data (prob = 0)

# +
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from sympy import symbols

from CompartmentalSystems.helpers_reservoir import DECAY_RATE_14C_YEARLY, ALPHA_14C, DECAY_RATE_14C_DAILY
from CompartmentalSystems.pwc_model_run_fd import PWCModelRunFD

from bgc_md2.models.CARDAMOM import CARDAMOMlib

# +
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
netCDF_filestem = "sol_acc_age_btt"

CARDAMOM_path = Path("/home/data/CARDAMOM/")
intcal20_path = CARDAMOM_path.joinpath("IntCal20_Year_Delta14C.csv")

data_combinations = [
    ("monthly", None, "discrete"),
    ("monthly", None, "continuous"), # only first prob computed so far
    ("yearly", 0, "continuous"),
    ("yearly", 6, "continuous")
]

datasets = dict()
for dc in data_combinations:
    time_resolution, delay_in_months, model_type = dc
    params = CARDAMOMlib.load_params(time_resolution, delay_in_months)
    output_path = data_path.joinpath(data_path.joinpath(params["output_folder"]))
    project_path = output_path.joinpath(model_type)

    ds_path = project_path.joinpath(netCDF_filestem)
    print(dc, ds_path)
    datasets[dc] = xr.open_mfdataset(str(ds_path) + "*.nc")
# -

#ds = datasets[("monthly", None, "discrete")]
ds = datasets[("monthly", None, "continuous")]
ds

# choose a site in northern Sweden, ensemble member prob=0
(lat, lon, prob) = (28, 52, 0)
sub_ds = ds.isel(lat=lat, lon=lon, prob=prob)
sub_ds

# +
# load the corresponding model run

# convert monthly time steps to daily time steps, each month comprising 31 days
data_times = np.arange(len(sub_ds["time"])) * 31

start_values = sub_ds["start_values"]
us = sub_ds["us"]
Bs = sub_ds["Bs"]

pwc_mr = PWCModelRunFD.from_Bs_and_us(
    symbols("t"),
    data_times,
    start_values.values,
    Bs[:-1],
    us[:-1],
    state_variables=sub_ds.pool.values
)
# -

fig, ax = plt.subplots(figsize=(8, 8))
pwc_mr.model.plot_pools_and_fluxes(ax, legend=False)
_ = ax.set_title("CARDAMOM")

# compute the solution
#soln = pwc_mr.solve()
soln = sub_ds["solution"].values

# +
fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(18, 12))

for nr, (pool_name, ax) in enumerate(zip(sub_ds.pool.values, axes.flatten())):
    ax.plot(sub_ds.time, soln[:, nr])
    ax.set_title(pool_name)
    ax.set_xlabel("year")
    ax.set_xlim([sub_ds.time[0], sub_ds.time[-1]])
    ax.set_ylabel(r"stock ($%s$)" % start_values.attrs["units"])
    
plt.suptitle("CARDAMOM model reconstruction (solution)")
plt.tight_layout()
# -

# Now we fake an equilibrium age distribution of the pools in the year 1920 by averaging over us and Bs for the first ten years and plugging the resulting averages in the equilibrium formulas.

# +
nr_time_steps = 10 * 12
p0 = pwc_mr.fake_start_age_densities(nr_time_steps)
age_moments = pwc_mr.fake_start_age_moments(nr_time_steps, 2)
mean_age_vector = age_moments[0, :] / 31 / 12 # convert from days to years
age_sd_vector = np.sqrt(age_moments[1, :]) / 31 / 12

fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(18, 12))
for nr, (pool_name, ax) in enumerate(zip(sub_ds.pool.values, axes.flatten())):
    ages = np.linspace(0, age_sd_vector[nr]*2, 100)
    pool_age_density = np.vectorize(lambda a: p0(a*31*12)[nr])
    ax.plot(ages, pool_age_density(ages))

    ax.set_title(pool_name)
    ax.set_xlabel("age (yr)")
    ax.set_xlim([ages[0], ages[-1]])
    ax.set_ylabel(r"density ($%s/yr$)" % start_values.attrs["units"])
    ax.set_ylim([0, ax.get_ylim()[-1]])

    ax.axvline(mean_age_vector[nr], ymax=pool_age_density(mean_age_vector[nr]) / ax.get_ylim()[-1], ls="--")
    
plt.suptitle("(Fake) equilibrium age distributions")
plt.tight_layout()
# -

# ## As an example show how to compute initial 14C stocks in 1920 for the Soil pool

# +
intcal20 = np.loadtxt(intcal20_path, delimiter=" ", skiprows=1, usecols=(1, 2)).transpose()

left_val = intcal20[1][np.argmin(intcal20[0])]
F_atm = interp1d(
    intcal20[0],
    intcal20[1],
#    kind="cubic",
    bounds_error=False,
    fill_value=(left_val, -1000.0) # no 14C oustide of dataset
#    fill_value="extrapolate"
)

# +
# at t0 = 1920 begins the CARDAMOM data
t0 = 1920

pool_name = "Soil"
pool = sub_ds.pool.values.tolist().index(pool_name)

# vectorize p0 in order to shorten the syntax in the plotting part
f = np.vectorize(lambda a: p0(a*31*12)[pool]) # convert a from years to days

# 14C inputs based on the age of the steady state pool value and the atmospheric Delta 14C at time t0-a
p0_14C = lambda a: (F_atm(t0-a)/1000.0 + 1) * f(a) * ALPHA_14C

# what remains from p0_14C until t0 (removing 14C particles that decayed in time interval [t0-a, t0])
p0_14C_with_decay = lambda a: p0_14C(a) * np.exp(-DECAY_RATE_14C_YEARLY*a)

# +
times = np.linspace(t0-age_sd_vector[pool]*2, t0, 200)

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(18, 12), sharex=True)
ax1.set_title("Intcal20", fontsize=24)
ax1.plot(times, F_atm(times))
ax1.set_ylabel(r"$\Delta^{14}$C (‰)", fontsize=24)
ax1.tick_params(axis='y', labelsize=20)
ax1.axhline(0, color="black", ls="--")
#ax1.axvline(np.min(intcal20[0]), c="black")

ax2.set_title(r"${}^{12}$C inputs ending up in $%s$ at time $t_0=%d$" % (pool_name, t0), fontsize=24)
ax2.plot(times, f(t0-times))
ax2.set_ylabel(r"mass (${}^{12}$C)", fontsize=24)
ax2.tick_params(axis='y', labelsize=20)
ax2.axvline(np.min(intcal20[0]), c="black")
ax2.set_ylim([0, ax2.get_ylim()[1]])
ax2.axvline(t0-mean_age_vector[pool], ymax=f(mean_age_vector[pool])/ax2.get_ylim()[1], ls="--")

ax3.set_title(r"Corresponding ${}^{14}$C inputs", fontsize=24)
ax3.plot(times, p0_14C(t0-times), label="without decay", color="purple")
ax3.plot(times, p0_14C_with_decay(t0-times), label="with decay", color="orange")
ax3.set_xlabel("year (AD)", fontsize=24)
ax3.set_ylabel(r"mass (${}^{14}$C)", fontsize=24)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim([times[0], times[-1]])
#ax3.axvline(np.min(intcal20[0]), c="black")
ax3.legend(fontsize=20)
ax3.set_ylim([0, ax3.get_ylim()[1]])

xs = np.linspace(times[0], times[-1], 50)
for x in xs:
    y = p0_14C_with_decay(t0-x) / ax3.get_ylim()[1]
    ax3.axvline(x, ymin=0, ymax=y, color="orange", lw=2)
    y1 = p0_14C(t0-x) / ax3.get_ylim()[1]
    ax3.axvline(x, ymin=y, ymax=y1, lw=2, color="purple")

fig.tight_layout()
# -

# The area under the blue curve is the initial 12C Soil content at time t0=1920, the purple area is 14C lost by radioactive decay in the system until t0, and the orange area is the 14C Soil content at time t0.

# ### Integration

# amount of material in the Soil in 1920
print("%2.2f PgC" % start_values[pool].values)

# +
# integration by solving an initial value problem is so much faster and less picky than using scipy.integrate.quad

soil_by_integration = solve_ivp(
    lambda a, _: f(a)*31*12,
    (0, 4500),#*31*12),
    np.array([0]),
)
soil_by_integration = soil_by_integration.y.reshape(-1)[-1]
soil_by_integration
# -

soil_14C_by_integration = solve_ivp(
    # if I divide by ALPHA_14C or not makes a stupid numerical difference for the integration 
    # which result is right?
    lambda a, _: p0_14C_with_decay(a)*31*12 / ALPHA_14C, 
    (0, 45000000),
    np.array([0])
)
soil_14C_by_integration = soil_14C_by_integration.y.reshape(-1)[-1] * ALPHA_14C
soil_14C_by_integration

# +
# Delta 14C in Soil in t0 = 1920 AD

print("Soil Delta 14C at t0 = 1920")
print("%2.2f ‰" % ((soil_14C_by_integration / soil_by_integration / ALPHA_14C - 1) * 1000))
# -

# ## Or use CARDAMOMlib

# +
start_values_14C = CARDAMOMlib.compute_start_values_14C(pwc_mr, nr_time_steps)
start_values_Delta_14C = (start_values_14C / start_values.values / ALPHA_14C - 1) * 1000

system_Delta_14C = (start_values_14C.sum() / start_values.values.sum() / ALPHA_14C - 1 ) * 1000

with np.printoptions(precision=2, suppress=True):
    print("Initial Delta14C (t0 = 1920)\n")
    print(pwc_mr.model.state_variables)
    print(start_values_Delta_14C, "‰")
    
    print("\nInitial Delta14C total system (t0 = 1920)\n")
    print("%2.2f ‰" % system_Delta_14C)
# -


