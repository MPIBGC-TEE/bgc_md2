#!/usr/bin/env python
# coding: utf-8

# # Transient age simulations, Emanuel model

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotly.offline import init_notebook_mode, iplot
from sympy import Matrix, symbols, Symbol, Function, latex
from scipy.interpolate import interp1d
from LAPM.linear_autonomous_pool_model import LinearAutonomousPoolModel
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun



########## symbol definitions ##########

# time symbol
time_symbol = symbols('t')

# Terrestrial carbon pools
x_1, x_2, x_3, x_4, x_5 = symbols('x_1 x_2 x_3 x_4 x_5')

# GPP inputs
G = Function('G')(time_symbol)

# Parameters
G_eq = symbols('G_eq')



########## model structure: equilibrium values and fluxes ##########

# equilibrium values
x_1e, x_2e, x_3e, x_4e, x_5e = (37.0, 452.0, 69.0, 81.0, 1121.0)

state_vector = Matrix([x_1, x_2, x_3, x_4, x_5])

B = Matrix([[-(25+21+31)/37,              0,             0,          0, 0],
            [         31/37, -(14+15+2)/452,             0,          0, 0],
            [             0,              0, -(18+12+6)/69,          0, 0],
            [         21/37,         15/452,         12/69, -(45+3)/81, 0],
            [             0,          2/452,          6/69,       3/81, -11/1121]])

u = G * Matrix(5, 1, [77/G_eq, 0, 36/G_eq, 0, 0])


srm=SmoothReservoirModel.from_B_u(state_vector, time_symbol, B, u)

# define the time and age windows of interest
start_year = 1851
end_year   = 2014
max_age    =  250

times = np.arange(start_year, end_year+1, 1)
ages  = np.arange(0, max_age+1, 1)


# Constant GPP
G_emanuel = 113.0

# GPP data
gpp_mean = pd.read_csv('~/cmip6stop/esm-hist/timeSeries/multimodGPPannual.csv')
gpp_anomaly=gpp_mean.GPPmean - gpp_mean.GPPmean[0] # Anomaly with respect to 1850 value
gpp_industrial=G_emanuel + gpp_anomaly

# linear interpolation of the (nonnegative) data points
u_interp = interp1d(gpp_mean.Year, gpp_industrial)

def u_func(t_val):
    # here we could do whatever we want to compute the input function
    return u_interp(t_val)

# define a dictionary to connect the symbols with the according functions
func_set = {G: u_func}


# the system starts in equilibrium
start_values = np.array([x_1e, x_2e, x_3e, x_4e, x_5e])
#start_values = np.zeros(5)

# parameter dictionary
par_dict = {G_eq: G_emanuel}

# create the nonlinear model run and the according cache
smrs = []
smr = SmoothModelRun(srm, par_dict, start_values, times, func_set)
smr.initialize_state_transition_operator_cache(
        lru_maxsize = 5000,
        size        =  100
    )
smrs.append(smr)


soln, _= smr.solve()


##### load linear autonomous pool model in steady state #####

# no fossil fuel inputs
u_eq = Matrix(5, 1, [77, 0, 36, 0, 0])

# force purely numerical treatment of the LAPM
# symbolic treatment would be too slow here
LM = LinearAutonomousPoolModel(u_eq, B, force_numerical=True)

## load equilibrium age densities ##

# the start age densities are given as a function of age that returns
# a vector of mass with that age
def start_age_densities(a):
    # we need to convert from sympy data types to numpy data types
    res =  np.array(LM.a_density(a)).astype(np.float64).reshape((5,)) * start_values
    return res

# get the start mean ages
start_mean_ages = np.array(LM.a_expected_value).astype(np.float64).reshape((5,))
start_age_moments = start_mean_ages.reshape((1,5))



age_densitiess       = []
pool_age_densitiess  = []
system_age_densities = []

p = smr.pool_age_densities_func(start_age_densities)
pool_age_densities = p(ages)
pool_age_densitiess.append(pool_age_densities)
        


system_age_density = smr.system_age_density(pool_age_densities)
system_age_densities.append(system_age_density)

print('Saving age densities', flush = True)
smr.save_pools_and_system_density_csv(
            'Emmanuel_age_dens.csv.gz',
            pool_age_densities,
            system_age_density,
            ages
        )
    
# combine pool and system age densities to one numpy array
age_densities = smr.age_densities(pool_age_densities, system_age_density)
age_densitiess.append(age_densities)
        
print('done', flush = True)


##### mean ages #####

pool_age_means   = []
system_age_means = []

pool_age_mean   = smr.age_moment_vector(1, start_age_moments)
system_age_mean = smr.system_age_moment(1, start_age_moments)
    
pool_age_means.append(pool_age_mean)
system_age_means.append(system_age_mean)


##### age medians #####

# start cumulative mass functions of age
# to have it available allows faster computations,
# since we save tedious numerical integrations
def F0(a):
    res = np.array(LM.a_cum_dist_func(a))        .astype(np.float64).reshape((5,)) * start_values
    return res

pool_age_medians   = []
system_age_medians = []

pool_age_median = smr.pool_age_distributions_quantiles_by_ode(
            0.5,
            start_age_densities,
            F0       = F0,
            max_step = 0.5
        )
system_age_median = smr.system_age_distribution_quantiles_by_ode(
            0.5,
            start_age_densities,
            F0       = F0,
            max_step = 0.5
        )

pool_age_medians.append(pool_age_median)
system_age_medians.append(system_age_median)


####### forward transit time #####

years = np.array([1850, 1900, 1950, 2000])
inds  = np.array([i for i, t in enumerate(times) if t in years])

ftt_densities = []
ftt_density_func = smr.forward_transit_time_density_func(times=years)
ftt_density = ftt_density_func(ages)


