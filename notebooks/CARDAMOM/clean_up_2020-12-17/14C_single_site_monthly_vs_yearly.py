# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Compute 14C for a single site in northern Sweden based on CARDAMOM output data
#
# We use monthly data, yearly data starting in January, and yearly data starting in July.

# +
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from sympy import symbols

from CompartmentalSystems.helpers_reservoir import DECAY_RATE_14C_YEARLY, ALPHA_14C, DECAY_RATE_14C_DAILY, F_Delta_14C

from CompartmentalSystems.pwc_model_run_fd import PWCModelRunFD
from CompartmentalSystems.pwc_model_run_14C import PWCModelRun_14C
from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR

from bgc_md2.models.CARDAMOM import CARDAMOMlib

# +
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
netCDF_filestem = "sol_acc_age_btt"

CARDAMOM_path = Path("/home/data/CARDAMOM/")
intcal20_path = CARDAMOM_path.joinpath("IntCal20_Year_Delta14C.csv")
Delta14C_atm_path = CARDAMOM_path.joinpath("C14Atm_NH.csv")

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

ds_m = datasets[("monthly", None, "continuous")]
ds_y0 = datasets[("yearly", 0, "continuous")]
ds_y6 = datasets[("yearly", 6, "continuous")]

# choose a site in northern Sweden, ensemble member prob=0
(lat, lon, prob) = (28, 52, 0)
sub_ds_m = ds_m.isel(lat=lat, lon=lon, prob=prob)
sub_ds_y0 = ds_y0.isel(lat=lat, lon=lon, prob=prob)
sub_ds_y6 = ds_y6.isel(lat=lat, lon=lon, prob=prob)

# +
# load the corresponding model runs

# convert monthly time steps to daily time steps, each month comprising 31 days
data_times_m = np.arange(len(sub_ds_m["time"])) * 31
data_times_y0 = np.arange(len(sub_ds_y0["time"])) * 31 * 12
data_times_y6 = np.arange(len(sub_ds_y6["time"])) * 31 * 12

# monthly
start_values_m = sub_ds_m["start_values"]
us_m = sub_ds_m["us"]
Bs_m = sub_ds_m["Bs"]

pwc_mr_m = PWCModelRunFD.from_Bs_and_us(
    symbols("t"),
    data_times_m,
    start_values_m.values,
    Bs_m[:-1],
    us_m[:-1],
    state_variables=sub_ds_m.pool.values
)

# yearly (Jan)
start_values_y0 = sub_ds_y0["start_values"]
us_y0 = sub_ds_y0["us"]
Bs_y0 = sub_ds_y0["Bs"]

pwc_mr_y0 = PWCModelRunFD.from_Bs_and_us(
    symbols("t"),
    data_times_y0,
    start_values_y0.values,
    Bs_y0[:-1],
    us_y0[:-1],
    state_variables=sub_ds_y0.pool.values
)

# yearly (Jul)
start_values_y6 = sub_ds_y6["start_values"]
us_y6 = sub_ds_y6["us"]
Bs_y6 = sub_ds_y6["Bs"]

pwc_mr_y6 = PWCModelRunFD.from_Bs_and_us(
    symbols("t"),
    data_times_y6,
    start_values_y6.values,
    Bs_y6[:-1],
    us_y6[:-1],
    state_variables=sub_ds_y6.pool.values
)
# -

fig, ax = plt.subplots(figsize=(8, 8))
pwc_mr_m.model.plot_pools_and_fluxes(ax, legend=False)
_ = ax.set_title("CARDAMOM")

# compute the solution
#soln = pwc_mr.solve()
soln_m = sub_ds_m["solution"].values
soln_y0 = sub_ds_y0["solution"].values
soln_y6 = sub_ds_y6["solution"].values

# +
fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(18, 12))

for nr, (pool_name, ax) in enumerate(zip(sub_ds_m.pool.values, axes.flatten())):
    ax.plot(sub_ds_m.time, soln_m[:, nr], label="monthly")
    ax.plot(sub_ds_y0.time, soln_y0[:, nr], label="yearly (Jan)")
    ax.plot(sub_ds_y6.time, soln_y6[:, nr], label="yearly (Jul)")

    ax.set_title(pool_name)
    ax.set_xlabel("year")
    ax.set_xlim([sub_ds_m.time[0], sub_ds_m.time[-1]])
    ax.set_ylabel(r"stock ($%s$)" % start_values_m.attrs["units"])
    ax.legend()
    
plt.suptitle("CARDAMOM model reconstruction (solution)")
plt.tight_layout()
# -

# Now we fake an equilibrium age distribution of the pools in the year 1920 by averaging over us and Bs for the first ten years and plugging the resulting averages in the equilibrium formulas.

# +
nr_time_steps_m = 10 * 12
nr_time_steps_y = 10

p0_m = pwc_mr_m.fake_start_age_densities(nr_time_steps_m)
age_moments_m = pwc_mr_m.fake_start_age_moments(nr_time_steps_m, 2)
mean_age_vector_m = age_moments_m[0, :] / 31 / 12 # convert from days to years
age_sd_vector_m = np.sqrt(age_moments_m[1, :]) / 31 / 12

p0_y0 = pwc_mr_y0.fake_start_age_densities(nr_time_steps_y)
age_moments_y0 = pwc_mr_y0.fake_start_age_moments(nr_time_steps_y, 2)
mean_age_vector_y0 = age_moments_y0[0, :] / 31 / 12 # convert from days to years
age_sd_vector_y0 = np.sqrt(age_moments_y0[1, :]) / 31 / 12

p0_y6 = pwc_mr_y6.fake_start_age_densities(nr_time_steps_y)
age_moments_y6 = pwc_mr_y6.fake_start_age_moments(nr_time_steps_y, 2)
mean_age_vector_y6 = age_moments_y6[0, :] / 31 / 12 # convert from days to years
age_sd_vector_y6 = np.sqrt(age_moments_y6[1, :]) / 31 / 12

# +
fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(18, 12))
for nr, (pool_name, ax) in enumerate(zip(sub_ds_m.pool.values, axes.flatten())):
    ages = np.linspace(0, age_sd_vector_m[nr]*2, 100)

    pool_age_density_m = np.vectorize(lambda a: p0_m(a*31*12)[nr])
    ax.plot(ages, pool_age_density_m(ages), c="blue", label="monthly")

    pool_age_density_y0 = np.vectorize(lambda a: p0_y0(a*31*12)[nr])
    ax.plot(ages, pool_age_density_y0(ages), c="orange", label="yearly (Jan)")

    pool_age_density_y6 = np.vectorize(lambda a: p0_y6(a*31*12)[nr])
    ax.plot(ages, pool_age_density_y6(ages), c="green", label="yearly (Jul)")

    ax.set_title(pool_name)
    ax.set_xlabel("age (yr)")
    ax.set_xlim([ages[0], ages[-1]])
    ax.set_ylabel(r"density ($%s/yr$)" % start_values_m.attrs["units"])
    ax.set_ylim([0, ax.get_ylim()[-1]])

    ax.axvline(mean_age_vector_m[nr], ymax=pool_age_density_m(mean_age_vector_m[nr]) / ax.get_ylim()[-1], ls="--", c="blue")
    ax.axvline(mean_age_vector_y0[nr], ymax=pool_age_density_y0(mean_age_vector_y0[nr]) / ax.get_ylim()[-1], ls="--", c="orange")
    ax.axvline(mean_age_vector_y6[nr], ymax=pool_age_density_y6(mean_age_vector_y6[nr]) / ax.get_ylim()[-1], ls="--", c="green")
    
    ax.legend()
    
plt.suptitle("(Fake) equilibrium age distributions")
plt.tight_layout()
# -

# ## Initial Delta 14C from (fake) initial age distributions

# +
start_values_14C_m = CARDAMOMlib.compute_start_values_14C(pwc_mr_m, nr_time_steps_m)
start_values_Delta_14C_m = (start_values_14C_m / start_values_m.values / ALPHA_14C - 1) * 1000
system_Delta_14C = (start_values_14C_m.sum() / start_values_m.values.sum() / ALPHA_14C - 1 ) * 1000
with np.printoptions(precision=2, suppress=True):
    print("Monthly initial Delta14C (t0 = 1920)\n")
    print(pwc_mr_m.model.state_variables)
    print(start_values_Delta_14C_m, "‰")
    print("System (1920): %2.2f ‰" % system_Delta_14C)

print("\n")
start_values_14C_y0 = CARDAMOMlib.compute_start_values_14C(pwc_mr_y0, nr_time_steps_y)
start_values_Delta_14C_y0 = (start_values_14C_y0 / start_values_y0.values / ALPHA_14C - 1) * 1000
system_Delta_14C_y0 = (start_values_14C_y0.sum() / start_values_y0.values.sum() / ALPHA_14C - 1 ) * 1000
with np.printoptions(precision=2, suppress=True):
    print("Yearly (Jan) initial Delta14C (t0 = 1920)\n")
    print(pwc_mr_y0.model.state_variables)
    print(start_values_Delta_14C_y0, "‰")
    print("System (1920): %2.2f ‰" % system_Delta_14C_y0)

print("\n")
start_values_14C_y6 = CARDAMOMlib.compute_start_values_14C(pwc_mr_y6, nr_time_steps_y)
start_values_Delta_14C_y6 = (start_values_14C_y6 / start_values_y6.values / ALPHA_14C - 1) * 1000
system_Delta_14C_y6 = (start_values_14C_y6.sum() / start_values_y6.values.sum() / ALPHA_14C - 1 ) * 1000
with np.printoptions(precision=2, suppress=True):
    print("Yearly (Jul) initial Delta14C (t0 = 1920)\n")
    print(pwc_mr_y6.model.state_variables)
    print(start_values_Delta_14C_y6, "‰")
    print("System (1920): %2.2f ‰" % system_Delta_14C_y6)
# -

# ## Construct a 14C model run and make a transient run from 1920 to 2009

# +
# load NH Delta14C dataset
# NOTE: apparently not completely consistent with intcal20 where they overlap
NH_Delta_14C = np.loadtxt(Delta14C_atm_path, delimiter=",", skiprows=1).transpose()

F_atm_NH = interp1d(
    NH_Delta_14C[0],
    NH_Delta_14C[1],
#    kind="cubic",
#    bounds_error=False,
#    fill_value=(left_val, -1000.0) # no 14C oustide of dataset
#    fill_value="extrapolate"
)

times = np.arange(1920, np.max(NH_Delta_14C[0]), 1)

fig, ax = plt.subplots(figsize=(18, 6))
ax.plot(times, F_atm_NH(times))
ax.set_xlim([times[0], times[-1]])
ax.axhline(0, c="black")
ax.set_ylabel(r"$\Delta^{14}$C (‰)")
_ = ax.set_title(r"${\Delta}^{14}$C")

# +
F_frac = lambda t: (F_atm_NH(t)/1000+1)*ALPHA_14C

fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(18, 8))
ax1.plot(times, F_frac(times))
ax1.set_xlim([times[0], times[-1]])
#ax.set_ylim([0, ax.get_ylim()[1]])
ax1.set_title(r"Fraction of ${}^{14}$C in the atmosphere (multiply external inputs to CARDAMOM with it)")
ax1.set_ylabel(r"${}^{14}$C fraction in ${}^{12}$C")

# rescale F_frac to model run time
t0 = 1920
F_frac_model = lambda t: F_frac(t0 + t/(31*12))
x = pwc_mr_m.times[pwc_mr_m.times / (31*12) < 2009 - 1920]
y = pwc_mr_y0.times[pwc_mr_y0.times / (31*12) < 2009 -1920]
ax2.plot(x, F_frac_model(x))
ax2.plot(y, F_frac_model(y), ls="--")
ax2.set_xlim([x[0], x[-1]])
#ax.set_ylim([0, ax.get_ylim()[1]])
ax2.set_ylabel(r"${}^{14}$C fraction in ${}^{12}$C")
_ = ax2.set_title(r"Fraction of ${}^{14}$C in the atmosphere (model run time: days after 1920-01-01)")

# +
# %%time

# construct a 14C model run from the 12C model run
pwc_mr_m.pwc_mr.times = x # cut off in 2009 because no 14C data later
pwc_mr_m_14C = PWCModelRun_14C(
    pwc_mr_m.pwc_mr,
    start_values_14C_m,
    F_frac_model,
    DECAY_RATE_14C_DAILY
)


# construct a 14C model run from the 12C model run
pwc_mr_y0.pwc_mr.times = y # cut off in 2009 because no 14C data later
pwc_mr_y0_14C = PWCModelRun_14C(
    pwc_mr_y0.pwc_mr,
    start_values_14C_y0,
    F_frac_model,
    DECAY_RATE_14C_DAILY
)


# construct a 14C model run from the 12C model run
pwc_mr_y6.pwc_mr.times = y # cut off in 2009 because no 14C data later
pwc_mr_y6_14C = PWCModelRun_14C(
    pwc_mr_y6.pwc_mr,
    start_values_14C_y6,
    F_frac_model,
    DECAY_RATE_14C_DAILY
)

# +
# %%time

soln_m_14C = pwc_mr_m_14C.solve()
soln_y0_14C = pwc_mr_y0_14C.solve()
soln_y6_14C = pwc_mr_y6_14C.solve()

# +
fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(18, 12))

cut_index_x = len(pwc_mr_m_14C.times)
cut_index_y = len(pwc_mr_y0_14C.times)
for nr, (pool_name, ax) in enumerate(zip(sub_ds_m.pool.values, axes.flatten())):
    ax.plot(sub_ds_m.time[:cut_index_x], F_Delta_14C(soln_m[:cut_index_x, nr], soln_m_14C[:, nr]), label="monthly")
    ax.plot(sub_ds_y0.time[:cut_index_y], F_Delta_14C(soln_y0[:cut_index_y, nr], soln_y0_14C[:, nr]), label="yearly (Jan)")
    ax.plot(sub_ds_y6.time[:cut_index_y], F_Delta_14C(soln_y6[:cut_index_y, nr], soln_y6_14C[:, nr]), label="yearly (Jul)")

    ax.set_title(pool_name)
    ax.set_xlabel("year")
    ax.set_xlim([sub_ds_m.time[0], sub_ds_m.time[-1]])
    ax.set_ylabel(r"$\Delta^{14}$C (‰)")
    ax.legend()
    
plt.suptitle(r"CARDAMOM $\Delta^{14}$C")
plt.tight_layout()
# -

# ## Now run 20 ensemble members 
#
# in new notebook, show also semi-transparent plots for age densities.


