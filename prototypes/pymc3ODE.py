# https://docs.pymc.io/notebooks/ODE_API_introduction.html
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pymc3 as pm
import theano

from pymc3.ode import DifferentialEquation
from scipy.integrate import odeint

plt.style.use("seaborn-darkgrid")
# For reproducibility
np.random.seed(20394)


def freefall(y, t, p):
    return 2.0 * p[1] - p[0] * y[0]


# Times for observation
times = np.arange(0, 10, 0.5)
gamma, g, y0, sigma = 0.4, 9.8, -2, 2
y = odeint(freefall, t=times, y0=y0, args=tuple([[gamma, g]]))
yobs = np.random.normal(y, 2)

fig, ax = plt.subplots()
ax.plot(times, yobs, label="observed speed", linestyle="dashed", marker="o", color="red")
ax.plot(times, y, label="True speed", color="k", alpha=0.5)
ax.legend()
ax.set_xlabel("Time (Seconds)")
ax.set_ylabel(r"$y(t)$")
fig.savefig('pymc3ODE.pdf')
ode_model = DifferentialEquation(func=freefall, times=times, n_states=1, n_theta=2, t0=0)

with pm.Model() as model:
    # Specify prior distributions for soem of our model parameters
    sigma = pm.HalfCauchy("sigma", 1)
    gamma = pm.Lognormal("gamma", 0, 1)

    # If we know one of the parameter values, we can simply pass the value.
    ode_solution = ode_model(y0=[0], theta=[gamma, 9.8])
    # The ode_solution has a shape of (n_times, n_states)

    Y = pm.Normal("Y", mu=ode_solution, sigma=sigma, observed=yobs)

    prior = pm.sample_prior_predictive()
    trace = pm.sample(2000, tune=1000, cores=1)
    posterior_predictive = pm.sample_posterior_predictive(trace)

    data = az.from_pymc3(trace=trace, prior=prior, posterior_predictive=posterior_predictive)


