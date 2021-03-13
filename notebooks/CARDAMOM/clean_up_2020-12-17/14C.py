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

# +
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy.integrate import quad
from scipy.interpolate import interp1d
from sympy import Matrix, init_printing

from CompartmentalSystems.helpers_reservoir import DECAY_RATE_14C_YEARLY, ALPHA_14C
from LAPM.linear_autonomous_pool_model import LinearAutonomousPoolModel as LAPM

init_printing()
# -

CARDAMOM_path = Path("/home/data/CARDAMOM/")
intcal20_path = CARDAMOM_path.joinpath("IntCal20_Year_Delta14C.csv")

intcal20 = np.loadtxt(intcal20_path, delimiter=" ", skiprows=1, usecols=(1, 2)).transpose()
#F_atm = interp1d(intcal20[0], intcal20[1], kind="linear")
left_val = intcal20[1][np.argmin(intcal20[0])]
F_atm = interp1d(
    intcal20[0],
    intcal20[1],
#        kind="cubic",
    bounds_error=False,
    fill_value=(left_val, -1000.0) # no 14C oustide of dataset
)


# +
u = Matrix(1, 1, [1])
B = Matrix([[-0.00005]])

eq_model = LAPM(u, B, force_numerical=True)
t0 = 1920
xss = eq_model.xss

def g(a, pool):
    res = eq_model.a_density(a)[pool] * xss[pool] if a >= 0 else np.nan
    return np.float(res)

pool = 0
f = np.vectorize(lambda a: g(a, pool))

p0_14C = lambda a: (F_atm(t0-a)/1000.0 + 1) * f(a) * ALPHA_14C
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

ax2.set_title(r"${}^{12}$C inputs remaining until $t_0$", fontsize=24)
ax2.plot(times, f(t0-times))
ax2.set_ylabel(r"mass (${}^{12}$C)", fontsize=24)
ax2.tick_params(axis='y', labelsize=20)
ax2.axvline(t0-eq_model.a_expected_value[pool], ymin=0, ymax=f(eq_model.a_expected_value[pool]), ls="--")
ax2.axvline(np.min(intcal20[0]), c="black")

ax3.set_title(r"Corresponding ${}^{14}$C inputs", fontsize=24)
ax3.plot(times, p0_14C(t0-times), label="without decay")
ax3.plot(times, p0_14C_with_decay(t0-times), label="with decay")
ax3.set_xlabel("year (AD)", fontsize=24)
ax3.set_ylabel(r"mass (${}^{14}$C)", fontsize=24)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim([times[0], times[-1]])
ax3.axvline(np.min(intcal20[0]), c="black")
ax3.legend(fontsize=20)

fig.tight_layout()
