# -*- coding: utf-8 -*-
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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Compute 14C for a single site in northern Sweden based on CARDAMOM output data

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
from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR

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

ds_cont = datasets[("monthly", None, "continuous")]
ds_disc = datasets[("monthly", None, "discrete")]

# choose a site in northern Sweden, ensemble member prob=0
(lat, lon, prob) = (28, 52, 0)
sub_ds_cont = ds_cont.isel(lat=lat, lon=lon, prob=prob)
sub_ds_disc = ds_disc.isel(lat=lat, lon=lon, prob=prob)

# +
# load the corresponding model runs

# convert monthly time steps to daily time steps, each month comprising 31 days
data_times = np.arange(len(sub_ds_cont["time"])) * 31


# continuous mr
start_values_cont = sub_ds_cont["start_values"]
us_cont = sub_ds_cont["us"]
Bs_cont = sub_ds_cont["Bs"]

pwc_mr = PWCModelRunFD.from_Bs_and_us(
    symbols("t"),
    data_times,
    start_values_cont.values,
    Bs_cont[:-1],
    us_cont[:-1],
    state_variables=sub_ds_cont.pool.values
)

# discrete mr
start_values_disc = sub_ds_disc["start_values"]
Us_disc = sub_ds_disc["Us"]
Bs_disc = sub_ds_disc["Bs"]

dmr = DMR.from_Bs_and_net_Us(
    start_values_disc.values,
    data_times,
    Bs_disc[:-1].values,
    Us_disc[:-1].values
)

# -

fig, ax = plt.subplots(figsize=(8, 8))
pwc_mr.model.plot_pools_and_fluxes(ax, legend=False)
_ = ax.set_title("CARDAMOM")

# compute the solution
#soln = pwc_mr.solve()
soln_cont = sub_ds_cont["solution"].values
soln_disc = sub_ds_disc["solution"].values

# +
fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(18, 12))

for nr, (pool_name, ax) in enumerate(zip(sub_ds_cont.pool.values, axes.flatten())):
    ax.plot(sub_ds_cont.time, (soln_cont-soln_disc)[:, nr])
    ax.set_title(pool_name)
    ax.set_xlabel("year")
    ax.set_xlim([sub_ds_cont.time[0], sub_ds_cont.time[-1]])
    ax.set_ylabel(r"stock ($%s$)" % start_values_cont.attrs["units"])
    
plt.suptitle("CARDAMOM model reconstruction (solution cont-disc)")
plt.tight_layout()
# -

# Now we fake an equilibrium age distribution of the pools in the year 1920 by averaging over us/Us and Bs for the first ten years and plugging the resulting averages in the equilibrium formulas.

# +
nr_time_steps = 10 * 12

p0_cont = pwc_mr.fake_start_age_densities(nr_time_steps)
age_moments_cont = pwc_mr.fake_start_age_moments(nr_time_steps, 2)
mean_age_vector_cont = age_moments_cont[0, :] / 31 / 12 # convert from days to years
age_sd_vector_cont = np.sqrt(age_moments_cont[1, :]) / 31 / 12

p0_disc_ai = dmr.fake_start_age_masses(nr_time_steps) 
p0_disc = lambda a: p0_disc_ai(int(a/dmr.dt)) / dmr.dt # make a density from masses
age_moments_disc = dmr.fake_start_age_moments(nr_time_steps, 2)
mean_age_vector_disc = age_moments_disc[0, :] / 31 / 12 # convert from days to years
age_sd_vector_disc = np.sqrt(age_moments_disc[1, :]) / 31 / 12

# +
fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(18, 12))
for nr, (pool_name, ax) in enumerate(zip(sub_ds_cont.pool.values, axes.flatten())):
    ages = np.linspace(0, age_sd_vector_cont[nr]*2, 100)

    pool_age_density_cont = np.vectorize(lambda a: p0_cont(a*31*12)[nr])
    ax.plot(ages, pool_age_density_cont(ages), c="blue", label="cont")

    pool_age_density_disc = np.vectorize(lambda a: p0_disc(a*31*12)[nr])
    ax.plot(ages, pool_age_density_disc(ages), c="orange", label="disc")

    ax.set_title(pool_name)
    ax.set_xlabel("age (yr)")
    ax.set_xlim([ages[0], ages[-1]])
    ax.set_ylabel(r"density ($%s/yr$)" % start_values_cont.attrs["units"])
    ax.set_ylim([0, ax.get_ylim()[-1]])

    ax.axvline(mean_age_vector_cont[nr], ymax=pool_age_density_cont(mean_age_vector_cont[nr]) / ax.get_ylim()[-1], ls="--", c="blue")
    ax.axvline(mean_age_vector_disc[nr], ymax=pool_age_density_disc(mean_age_vector_disc[nr]) / ax.get_ylim()[-1], ls="--", c="orange")
    
    ax.legend()
    
plt.suptitle("(Fake) equilibrium age distributions")
plt.tight_layout()
# -

# ## As an example show how to compute initial 14C stocks in 1920 for the Soil pool

# +
intcal20 = np.loadtxt(intcal20_path, delimiter=" ", skiprows=1, usecols=(1, 2)).transpose()

left_val = intcal20[1][np.argmin(intcal20[0])]
F_atm_raw = interp1d(
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
pool = sub_ds_cont.pool.values.tolist().index(pool_name)

# vectorize p0 in order to shorten the syntax in the plotting part
f_cont = np.vectorize(lambda a: p0_cont(a*31*12)[pool]) # convert a from years to days
f_disc = np.vectorize(lambda a: p0_disc(a*31*12)[pool]) # convert a from years to days

# 14C inputs based on the age of the steady state pool value and the atmospheric Delta 14C at time t0-a
p0_14C_cont = lambda a: (F_atm_raw(t0-a)/1000.0 + 1) * f_cont(a) * ALPHA_14C
p0_14C_disc = lambda a: (F_atm_raw(t0-a)/1000.0 + 1) * f_disc(a) * ALPHA_14C

# what remains from p0_14C until t0 (removing 14C particles that decayed in time interval [t0-a, t0])
p0_14C_with_decay_cont = lambda a: p0_14C_cont(a) * np.exp(-DECAY_RATE_14C_YEARLY*a)
p0_14C_with_decay_disc = lambda a: p0_14C_disc(a) * np.exp(-DECAY_RATE_14C_YEARLY*a)

# +
times = np.linspace(t0-age_sd_vector_cont[pool]*2, t0, 100)

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(18, 12), sharex=True)
ax1.set_title("Intcal20", fontsize=24)
ax1.plot(times, F_atm_raw(times), c="black")
ax1.set_ylabel(r"$\Delta^{14}$C (‰)", fontsize=24)
ax1.tick_params(axis='y', labelsize=20)
ax1.axhline(0, color="black", ls="--")
#ax1.axvline(np.min(intcal20[0]), c="black")

ax2.set_title(r"${}^{12}$C inputs ending up in $%s$ at time $t_0=%d$" % (pool_name, t0), fontsize=24)
ax2.plot(times, f_cont(t0-times), c="blue")
ax2.plot(times, f_disc(t0-times), c="orange", ls="--")

ax2.set_ylabel(r"density (${}^{12}$C/yr)", fontsize=24)
ax2.tick_params(axis='y', labelsize=20)
ax2.axvline(np.min(intcal20[0]), c="black")
ax2.set_ylim([0, ax2.get_ylim()[1]])
ax2.axvline(t0-mean_age_vector_cont[pool], ymax=f_cont(mean_age_vector_cont[pool])/ax2.get_ylim()[1], c="blue")
ax2.axvline(t0-mean_age_vector_disc[pool], ymax=f_disc(mean_age_vector_disc[pool])/ax2.get_ylim()[1], c="orange", ls="--")

ax3.set_title(r"Corresponding ${}^{14}$C inputs", fontsize=24)
ax3.plot(times, p0_14C_cont(t0-times), label="without decay (cont)", color="blue")
ax3.plot(times, p0_14C_with_decay_cont(t0-times), label="with decay (cont)", color="blue", ls="--")

ax3.plot(times, p0_14C_disc(t0-times), label="without decay (disc)", color="orange")
ax3.plot(times, p0_14C_with_decay_disc(t0-times), label="with decay (disc)", color="orange", ls="--")

ax3.set_xlabel("year (AD)", fontsize=24)
ax3.set_ylabel(r"mass (${}^{14}$C)", fontsize=24)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim([times[0], times[-1]])
#ax3.axvline(np.min(intcal20[0]), c="black")
ax3.legend(fontsize=20)
ax3.set_ylim([0, ax3.get_ylim()[1]])

#xs = np.linspace(times[0], times[-1], 50)
#for x in xs:
#    y = p0_14C_with_decay(t0-x) / ax3.get_ylim()[1]
#    ax3.axvline(x, ymin=0, ymax=y, color="orange", lw=2)
#    y1 = p0_14C(t0-x) / ax3.get_ylim()[1]
#    ax3.axvline(x, ymin=y, ymax=y1, lw=2, color="purple")

fig.tight_layout()
# -

# The area under the blue curve is the initial 12C Soil content at time t0=1920, the purple area is 14C lost by radioactive decay in the system until t0, and the orange area is the 14C Soil content at time t0.

# ### Integration

# amount of material in the Soil in 1920
print("%2.2f PgC" % start_values_cont[pool].values)

# +
# integration by solving an initial value problem is so much faster and less picky than using scipy.integrate.quad

soil_by_integration_cont = solve_ivp(
    lambda a, _: f_cont(a)*31*12,
    (0, 4500),#*31*12),
    np.array([0]),
)
soil_by_integration_cont = soil_by_integration_cont.y.reshape(-1)[-1]

soil_by_integration_disc = solve_ivp(
    lambda a, _: f_disc(a)*31*12,
    (0, 4500),#*31*12),
    np.array([0]),
)
soil_by_integration_disc = soil_by_integration_disc.y.reshape(-1)[-1]

soil_by_integration_cont, soil_by_integration_disc

# +
soil_14C_by_integration_cont = solve_ivp(
    # if I divide by ALPHA_14C or not makes a stupid numerical difference for the integration 
    # which result is right?
    lambda a, _: p0_14C_with_decay_cont(a)*31*12 / ALPHA_14C, 
    (0, 45000000),
    np.array([0])
)
soil_14C_by_integration_cont = soil_14C_by_integration_cont.y.reshape(-1)[-1] * ALPHA_14C

soil_14C_by_integration_disc = solve_ivp(
    # if I divide by ALPHA_14C or not makes a stupid numerical difference for the integration 
    # which result is right?
    lambda a, _: p0_14C_with_decay_disc(a)*31*12 / ALPHA_14C, 
    (0, 45000000),
    np.array([0])
)
soil_14C_by_integration_disc = soil_14C_by_integration_disc.y.reshape(-1)[-1] * ALPHA_14C

soil_14C_by_integration_cont, soil_14C_by_integration_disc

# +
# Delta 14C in Soil in t0 = 1920 AD

print("Soil Delta 14C at t0 = 1920")
print("%2.2f ‰" % ((soil_14C_by_integration_disc / soil_by_integration_disc / ALPHA_14C - 1) * 1000))
# -

# ## Or use CARDAMOMlib

# +
start_values_14C_disc = CARDAMOMlib.compute_start_values_14C(dmr, nr_time_steps)

start_values_Delta_14C_disc = (start_values_14C_disc / start_values_disc.values / ALPHA_14C - 1) * 1000
with np.printoptions(precision=2, suppress=True):
    print("Initial Delta14C (t0 = 1920)\n")
    print(pwc_mr.model.state_variables)
    print(start_values_Delta_14C_disc, "‰")
# -


