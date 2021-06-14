#!/usr/bin/env python
# coding: utf-8

# # Transient age simulations, Emanuel model

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotly.offline import init_notebook_mode, iplot
from sympy import Matrix, symbols, Symbol, Function, latex, exp, log
from scipy.interpolate import interp1d
from LAPM.linear_autonomous_pool_model import LinearAutonomousPoolModel
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun
from CompartmentalSystems.helpers_reservoir import (
    numerical_function_from_expression
)


########## symbol definitions ##########

# time symbol
time_symbol = symbols("t")

# Terrestrial carbon pools
x_1, x_2, x_3, x_4, x_5 = symbols("x_1 x_2 x_3 x_4 x_5")

# GPP inputs
G = Function("G")(time_symbol)

# Time-dependet modifier for matrix B
xi = Function("xi")(time_symbol)

# Parameters
G_eq = symbols("G_eq")

########## model structure: equilibrium values and fluxes ##########
pool_names = [
    "Non-woody tree parts",
    "Woody tree parts",
    "Ground vegetation",
    "Detritus/Decomposers",
    "Soil",
]

# equilibrium values
x_1e, x_2e, x_3e, x_4e, x_5e = (37.0, 452.0, 69.0, 81.0, 1121.0)

state_vector = Matrix([x_1, x_2, x_3, x_4, x_5])

B_eq = Matrix(
    [
        [-(25 + 21 + 31) / 37, 0, 0, 0, 0],
        [31 / 37, -(14 + 15 + 2) / 452, 0, 0, 0],
        [0, 0, -(18 + 12 + 6) / 69, 0, 0],
        [21 / 37, 15 / 452, 12 / 69, -(45 + 3) / 81, 0],
        [0, 2 / 452, 6 / 69, 3 / 81, -11 / 1121],
    ]
)

B = xi * B_eq

u = G * Matrix(5, 1, [77 / G_eq, 0, 36 / G_eq, 0, 0])


srm = SmoothReservoirModel.from_B_u(state_vector, time_symbol, B, u)

# define the time and age windows of interest
start_year = 1851
end_year = 2014
max_age = 250

times = np.arange(start_year, end_year + 1, 1)
ages = np.arange(0, max_age + 1, 1)


# Constant GPP
G_emanuel = 113.0

# Time dependent GPP from Rasmussen et al. (2016, JMB). We assume that proportional increase in GPP is equal for both pools
T_s0 = 15
s_0 = 1
sigma = 4.5
alpha = 1
rho = 0.65
f_i = 1
x_a = start_year*exp(0.0305*time_symbol)/(start_year+exp(0.0305*time_symbol)-1)+284 
T_s = T_s0 + sigma/log(2)*log(x_a/285)
Gamma = 42.7 + 1.68*(T_s-25) + 0.012*(T_s-25)**2
beta = 3*rho*x_a*Gamma/((rho*x_a-Gamma)*(rho*x_a+2*Gamma))
s_i = f_i*alpha*s_0*(1+2.5*beta*log(x_a/285)) 

gpp_industrial = G_emanuel * s_i

u_func_numerical = numerical_function_from_expression(
    gpp_industrial, # sympy expression
    (time_symbol,), # tuple of symbols to be replaced by numerical values
    {}, # parameter dict
    {} # func dict
)

xi_b = 2
xi_t = xi_b**(0.1*T_s-1.5)

xi_func_numerical = numerical_function_from_expression(
    xi_t, # sympy expression
    (time_symbol,), # tuple of symbols to be replaced by numerical values
    {}, # parameter dict
    {} # func dict
)

# define a dictionary to connect the symbols with the according functions
func_set = {xi: xi_func_numerical,
             G: u_func_numerical}

# the system starts in equilibrium
start_values = np.array([x_1e, x_2e, x_3e, x_4e, x_5e])
# start_values = np.zeros(5)

# parameter dictionary
par_dict = {G_eq: G_emanuel}

# create the nonlinear model run and the according cache
smrs = []
smr = SmoothModelRun(srm, par_dict, start_values, times, func_set)
smr.initialize_state_transition_operator_cache(lru_maxsize=5000, size=100)
smrs.append(smr)

soln = smr.solve()
cstocks = pd.DataFrame(soln, columns=pool_names)
stocks = cstocks.join(pd.DataFrame({"Time": times}))
stocks.to_csv("stocks.csv", index=False)


##### load linear autonomous pool model in steady state #####

# no fossil fuel inputs
u_eq = Matrix(5, 1, [77, 0, 36, 0, 0])

# force purely numerical treatment of the LAPM
# symbolic treatment would be too slow here
LM = LinearAutonomousPoolModel(u_eq, B_eq, force_numerical=True)

## load equilibrium age densities ##

# the start age densities are given as a function of age that returns
# a vector of mass with that age
def start_age_densities(a):
    # we need to convert from sympy data types to numpy data types
    res = np.array(LM.a_density(a)).astype(np.float64).reshape((5,)) * start_values
    return res


# get the start mean ages
start_mean_ages = np.array(LM.a_expected_value).astype(np.float64).reshape((5,))
start_age_moments = start_mean_ages.reshape((1, 5))

# Equilibrium age distribution
y = np.array([start_age_densities(a) for a in ages])

start_age_dens = pd.DataFrame(y, columns=pool_names)
start_age_dens.to_csv("start_age_dens.csv", index=False)


age_densitiess = []
pool_age_densitiess = []
system_age_densities = []

p = smr.pool_age_densities_func(start_age_densities)
pool_age_densities = p(ages)

system_age_density = smr.system_age_density(pool_age_densities)

print("Saving age densities", flush=True)
smr.save_pools_and_system_density_csv(
    "Emmanuel_age_dens.csv", pool_age_densities, system_age_density, ages
)

# combine pool and system age densities to one numpy array
# age_densities = smr.age_densities(pool_age_densities, system_age_density)
# age_densitiess.append(age_densities)
#
# print('done', flush = True)


##### mean ages #####

pool_age_means = []
system_age_means = []

pool_age_mean = smr.age_moment_vector(1, start_age_moments)
system_age_mean = smr.system_age_moment(1, start_age_moments)

mean_pool_ages = pd.DataFrame(pool_age_mean, columns=pool_names)
mean_ages = mean_pool_ages.join(pd.DataFrame({"System Age": system_age_mean}))
mean_ages.to_csv("mean_ages.csv", index=False)

##### age medians #####

# start cumulative mass functions of age
# to have it available allows faster computations,
# since we save tedious numerical integrations
def F0(a):
    res = (
        np.array(LM.a_cum_dist_func(a)).astype(np.float64).reshape((5,)) * start_values
    )
    return res


pool_age_medians = []
system_age_medians = []

pool_age_median = smr.pool_age_distributions_quantiles(
    quantile=0.5, start_age_densities=start_age_densities, F0=F0 #, max_step=0.5
)
system_age_median = smr.system_age_distribution_quantiles(
    quantile=0.5, start_age_densities=start_age_densities #F0=F0
)

median_pool_age = pd.DataFrame(pool_age_median, columns=pool_names)
median_ages = median_pool_age.join(pd.DataFrame({"SystemA Age": system_age_median}))
median_ages.to_csv("median_ages.csv", index=False)


####### forward transit time #####

# years = np.array([1850, 1900, 1950, 2000])
years = np.arange(1851, 2001, 1)
inds = np.array([i for i, t in enumerate(times) if t in years])

ftt_densities = []
ftt_density_func = smr.forward_transit_time_density_func(times=years)
ftt_density = ftt_density_func(ages)

ftt = pd.DataFrame(ftt_density, columns=years)
ftt.to_csv("forward_transit_time.csv", index=False)

# Backward transit times
tol = 1e-1
btt_density = smr.backward_transit_time_density(pool_age_densities)
btt_dens = pd.DataFrame(btt_density, columns=times)
btt_mean = smr.backward_transit_time_moment(1, start_age_moments)
F_btt_sv = smr.cumulative_backward_transit_time_distribution_single_value_func(F0=F0)
btt_median = smr.distribution_quantiles(mr = smr, quantile = 0.5, 
                                            F_sv = F_btt_sv, 
                                            norm_consts = smr.external_output_vector.sum(1), 
                                            start_values = btt_mean, 
                                            method = 'brentq', 
                                            tol = tol)


btt_dens.to_csv("backward_transit_time_densities.csv", index=False)

btt_mean_and_median=pd.DataFrame({'Year':times, 'meanBTT':btt_mean, 'medianBTT':btt_median})
btt_mean_and_median.to_csv("backward_transit_time_mean_and_median.csv", index=False)
