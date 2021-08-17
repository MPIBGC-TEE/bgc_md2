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

# # Construct initial 14C stocks from initial 12C age distribution

# +
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
from sympy import Matrix, init_printing, symbols

from CompartmentalSystems.helpers_reservoir import DECAY_RATE_14C_YEARLY, ALPHA_14C
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from LAPM.linear_autonomous_pool_model import LinearAutonomousPoolModel as LAPM

init_printing()
# -

# Load the Intcal20 data.

CARDAMOM_path = Path("/home/data/CARDAMOM/")
intcal20_path = CARDAMOM_path.joinpath("IntCal20_Year_Delta14C.csv")

# Make a function from the data with constant continuation to the left (left_val) and no 14C continuation to the right (-1000).

# +
intcal20 = np.loadtxt(intcal20_path, delimiter=" ", skiprows=1, usecols=(1, 2)).transpose()

left_val = intcal20[1][np.argmin(intcal20[0])]
F_atm = interp1d(
    intcal20[0],
    intcal20[1],
#        kind="cubic",
    bounds_error=False,
    fill_value=(left_val, -1000.0) # no 14C oustide of dataset
)
# -

# Construct a simple linear autonomous pool model in equilibrium, and fix pool 0 as the pool of interest

# +
u = Matrix(2, 1, [1, 2])
B = Matrix(
    [[-0.00005,  0.0001],
     [ 0.00001, -0.0002]]
)
eq_model = LAPM(u, B, force_numerical=True)

t, x, y = symbols("t, x, y")
sv = Matrix(2, 1, [x, y])
srm = SmoothReservoirModel.from_B_u(sv, t, B, u)
fig, ax = plt.subplots(figsize=(6, 6))
srm.plot_pools_and_fluxes(ax, legend=False)

# +
# at t0 = 1920 begins the CARDAMOM data, for example
t0 = 1920

# steady state vector
xss = eq_model.xss

# initial/equilibrium age density in terms of mass --> intergrates to steady state value
def p0(a, pool):
    res = eq_model.a_density(a)[pool] * xss[pool] if a >= 0 else np.nan
    return np.float(res)

pool = 0

# vectorize p0 in order to shorten the syntax in the plotting part
f = np.vectorize(lambda a: p0(a, pool))

# 14C inputs based on the age of the steady state pool value and the atmospheric Delta 14C at time t0-a
p0_14C = lambda a: (F_atm(t0-a)/1000.0 + 1) * f(a) * ALPHA_14C

# what remains from p0_14C until t0 (removing 14C particles that decayed in time interval [t0-a, t0])
p0_14C_with_decay = lambda a: p0_14C(a) * np.exp(-DECAY_RATE_14C_YEARLY*a)

# +
times = np.linspace(-60_000, t0, 1000)

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(18, 12), sharex=True)
ax1.set_title("Intcal20", fontsize=24)
ax1.plot(times, F_atm(times))
ax1.set_ylabel(r"$\Delta^{14}$C (‰)", fontsize=24)
ax1.tick_params(axis='y', labelsize=20)
ax1.axhline(0, color="black", ls="--")
ax1.axvline(np.min(intcal20[0]), c="black")

ax2.set_title(r"${}^{12}$C inputs ending up in $x$ at time $t_0=%d$" % t0, fontsize=24)
ax2.plot(times, f(t0-times))
ax2.set_ylabel(r"mass (${}^{12}$C)", fontsize=24)
ax2.tick_params(axis='y', labelsize=20)
#ax2.axvline(t0-eq_model.a_expected_value[pool])#, ymin=0, ymax=f(eq_model.a_expected_value[pool]), ls="--")
ax2.axvline(np.min(intcal20[0]), c="black")
ax2.set_ylim([0, ax2.get_ylim()[1]])

xs = np.linspace(-60_000, t0, 200)
for x in xs:
    y = f(t0-x) / ax2.get_ylim()[1]
    ax2.axvline(x, ymin=0, ymax=y, lw=2)

ax3.set_title(r"Corresponding ${}^{14}$C inputs", fontsize=24)
ax3.plot(times, p0_14C(t0-times), label="without decay", color="purple")
ax3.plot(times, p0_14C_with_decay(t0-times), label="with decay", color="orange")
ax3.set_xlabel("year (AD)", fontsize=24)
ax3.set_ylabel(r"mass (${}^{14}$C)", fontsize=24)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim([times[0], times[-1]])
ax3.axvline(np.min(intcal20[0]), c="black")
ax3.legend(fontsize=20)
ax3.set_ylim([0, ax3.get_ylim()[1]])

xs = np.linspace(-60_000, t0, 200)
for x in xs:
    y = p0_14C_with_decay(t0-x) / ax3.get_ylim()[1]
    ax3.axvline(x, ymin=0, ymax=y, color="orange", lw=2)
    y1 = p0_14C(t0-x) / ax3.get_ylim()[1]
    ax3.axvline(x, ymin=y, ymax=y1, lw=2, color="purple")

fig.tight_layout()
# -

#
# The blue area is the initial 12C pool content x at time t0, the purple area is 14C lost by radioactive decay in the system until t0, and the orange area is the 14C pool content of x at time t0.


