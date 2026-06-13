#!/usr/bin/env python
#pip install -U protobuf


import sys
#import git
sys.path.insert(0, '..')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import json
from scipy.stats import gaussian_kde

from model_specific_helpers import (
    get_example_site_vars,
    make_param_filter_func,
    # make_weighted_cost_func,
    make_param2res,
    # make_param2res_2,
    UnEstimatedParameters,
    EstimatedParameters,
    Observables
)

from general_helpers import (
    make_uniform_proposer,
    make_multivariate_normal_proposer,
    mcmc,
    adaptive_mcmc,
    autostep_mcmc,
    make_feng_cost_func,
    plot_solutions
)

with Path('config.json').open(mode='r') as f:
    conf_dict = json.load(f)

dataPath = Path(conf_dict['dataPath'])
#dataPath = 'E:/BaiduNetdiskDownload/OCN'

npp,  C_leaf, C_wood, C_root, C_litter, C_soil, rh, ra, mrso, tsl = get_example_site_vars(dataPath)


#nyears = 150
nyears = 10  # reduced time span for testing purposes
tot_len = 12 * nyears
obs_tup = Observables(
    C_leaf=C_leaf,
    C_wood=C_wood,
    C_root=C_root,
    C_litter=C_litter,
    C_soil=C_soil,
    rh=rh,
    ra=ra,
)
obs = np.stack(obs_tup, axis=1)[0:tot_len, :]

# save observational data for comparison with model output
pd.DataFrame(obs).to_csv(dataPath.joinpath('obs.csv'), sep=',')

cpa = UnEstimatedParameters(
    C_wood1_0=  0.76,#C_wood1[0],
     C_wood2_0= 0.67,#C_wood2[0],
     C_wood3_0= 0.34,#C_wood3[0],
     C_wood4_0= 0.42,#C_wood4[0],
     C_leaf_0=C_leaf[0],
     C_root_0=C_root[0],
      C_fruit_0= 0.056, #C_fruit[0],
      C_litter1_0= 0.005,#C_litter1[0],
      C_litter2_0= 0.009,#C_litter2[0],
      C_litter3_0= 0.035,#C_litter3[0],
      C_litter4_0= 0.017,#C_litter4[0],
     C_litter5_0= 0.05, #C_litter5[0],
      C_litter6_0= 0.065,#C_litter6[0],
      C_surface_som_0= 0.0044,#C_surface_som[0],
      C_fast_som_0= 0.012,#C_fast_som[0],
      C_slow_som_0= 4,#C_slow_som[0],
     C_pass_som_0= 0.0125,#C_pass_som[0],
       rh_0=rh[0],
       ra_0=ra[0],
      npp=npp,
      number_of_months=tot_len,
     mrso=mrso,
    tsl=tsl,


   # C_wood_0=C_wood[0],
    # C_leaf_0=C_leaf[0],
    # C_root_0=C_root[0],
    # C_litter_0=C_litter[0],
    # C_soil_0=C_soil[0],
    # rh_0=rh[0],
    # ra_0=ra[0],
    # npp=npp,
    #  number_of_months=tot_len,
    #  mrso=mrso,
    #  tsl=tsl,
    beta_fruit = 0.1
)
#print(cpa)
param2res = make_param2res(cpa)

c_min = np.array(
    EstimatedParameters(
        beta_sapwood1=0,
        beta_sapwood2=0,
        beta_leaf=0,
        beta_root=0,
        f_sapwood1_heartwood1=0.001,
        f_sapwood2_heartwood2=0.001,

        f_wood1_liiter1=0.001,  # 2
        f_wood2_liiter2=0.001,  # 2
        f_leaf_liiter3=0.001,  # 2
        f_root_liiter4=0.001,  # 2
        f_fruit_liiter3=0.001,  # 2

        f_litter1_surface_som=0.1,  # 2
        f_litter1_fast_som=0.1,  # 2

        f_litter2_fast_som=0.1,  # 2
        f_litter2_slow_som=0.01,  # 2

        f_litter3_surface_som=0.1,  # 2
        f_litter3_slow_som=0.01,  # 2

        f_litter4_surface_som=0.1,  # 2
        f_litter4_fast_som=0.1,  # 2

        f_litter5_surface_som=0.1,  # 2
        f_litter6_fast_som=0.1,  # 2

        f_surface_som_slow_som=0.001,  # 2

        #  "k_labile",  # 11
        k_leaf=1 / (365 * 2),  # 11
        k_wood1=1 / (365 * 60),  # 12
        k_wood2=1 / (365 * 60),  # 12
        k_wood3=1 / (365 * 60),  # 12
        k_wood4=1 / (365 * 60),  # 12
        k_root=1 / (365 * 30),  # 13
        k_fruit=1 / (365 * 20),
        k_lit1=1 / (365 * 30),  # 14
        k_lit2=1 / (365 * 30),  # 15
        k_lit3=1 / (365 * 30),  # 16        "k_leaf_lit",  # 14
        k_lit4=1 / (365 * 30),  # 15
        k_lit5=1 / (365 * 30),  # 16
        k_lit6=1 / (365 * 30),  # 17
        k_surface_som=1 / (365 * 100),  # 18
        k_fast_som=1 / (365 * 200),  # 18
        k_slow_som=1 / (365 * 500),  # 18
        k_pass_som=1 / (365 * 1000),  # 19
        T_0=-10,
        E=1,
        KM=1
    )
)

c_max = np.array(
    EstimatedParameters(
        beta_sapwood1=0.4,
        beta_sapwood2=0.4,
        beta_leaf=0.4,
        beta_root=0.4,
        f_sapwood1_heartwood1=0.9,
        f_sapwood2_heartwood2=0.9,

        f_wood1_liiter1=0.9,  # 2
        f_wood2_liiter2=0.9,  # 2
        f_leaf_liiter3=0.9,  # 2
        f_root_liiter4=0.9,  # 2
        f_fruit_liiter3=0.9,  # 2

        f_litter1_surface_som=0.9,  # 2
        f_litter1_fast_som=0.9,  # 2

        f_litter2_fast_som=0.9,  # 2
        f_litter2_slow_som=0.9,  # 2

        f_litter3_surface_som=0.9,  # 2
        f_litter3_slow_som=0.9,  # 2

        f_litter4_surface_som=0.9,  # 2
        f_litter4_fast_som=0.9,  # 2

        f_litter5_surface_som=0.9,  # 2
        f_litter6_fast_som=0.9,  # 2

        f_surface_som_slow_som=0.9,  # 2

        #  "k_labile",  # 11
        k_leaf=1 / 30,  # 11
        k_wood1=1 / (365 * 1),  # 12
        k_wood2=1 / (365 * 1),  # 12
        k_wood3=1 / (365 * 1),  # 12
        k_wood4=1 / (365 * 1),  # 12
        k_root=1 / (365 * 0.5),  # 13
        k_fruit=1 / (365 * 2),
        k_lit1=1 / (365 * 0.5),  # 14
        k_lit2=1 / (365 * 0.5),  # 15
        k_lit3=1 / (365 * 0.5),  # 16        "k_leaf_lit",  # 14
        k_lit4=1 / (365 * 0.5),  # 15
        k_lit5=1 / (365 * 0.5),  # 16
        k_lit6=1 / (365 * 0.5),  # 17
        k_surface_som=1 / (365 * 0.5),  # 18
        k_fast_som=1 / (365 * 1),  # 18
        k_slow_som=1 / (365 * 3.5),  # 18
        k_pass_som=1 / (365 * 10),  # 19
        T_0=4,
        E=100,
        KM=100
    )
)
# I commented your original settings and instead used the
# values constructed from the limits (see below)
epa_0 = EstimatedParameters(

    beta_sapwood1=0.15,
    beta_sapwood2=0.15,
    beta_leaf=0.2,
    beta_root=0.2,
    f_sapwood1_heartwood1=0.4,
    f_sapwood2_heartwood2=0.4,

    f_wood1_liiter1=0.5,  # 2
    f_wood2_liiter2=0.5,  # 2
    f_leaf_liiter3=0.4,  # 2
    f_root_liiter4=0.4,  # 2
    f_fruit_liiter3=0.4,  # 2

    f_litter1_surface_som=0.4,  # 2
    f_litter1_fast_som=0.4,  # 2

    f_litter2_fast_som=0.4,  # 2
    f_litter2_slow_som=0.4,  # 2

    f_litter3_surface_som=0.4,  # 2
    f_litter3_slow_som=0.4,  # 2

    f_litter4_surface_som=0.4,  # 2
    f_litter4_fast_som=0.4,  # 2

    f_litter5_surface_som=0.49,  # 2
    f_litter6_fast_som=0.4,  # 2

    f_surface_som_slow_som=0.4,  # 2

    #  "k_labile",  # 11
    k_leaf=1 / (60 * 2),  # 11
    k_wood1=1 / (365 * 30),  # 12
    k_wood2=1 / (365 * 30),  # 12
    k_wood3=1 / (365 * 30),  # 12
    k_wood4=1 / (365 * 30),  # 12
    k_root=1 / (365 * 22),  # 13
    k_fruit=1 / (365 * 4),
    k_lit1=1 / (365 * 3.3),  # 14
    k_lit2=1 / (365 * 3.3),  # 15
    k_lit3=1 / (365 * 3.3),  # 16        "k_leaf_lit",  # 14
    k_lit4=1 / (365 * 11),  # 15
    k_lit5=1 / (365 * 11),  # 16
    k_lit6=1 / (365 * 11),  # 17
    k_surface_som=1 / (365 * 5),  # 18
    k_fast_som=1 / (365 * 18),  # 18
    k_slow_som=1 / (365 * 100),  # 18
    k_pass_som=1 / (365 * 350),  # 19
    T_0=2,
    E=4,
    KM=10

)
# save initial parameters
pd.DataFrame(epa_0).to_csv(dataPath.joinpath('epa_0.csv'), sep=',')

# epa_0 = EstimatedParameters._make(
#         np.concatenate(
#             [
#                 # we don't want to use the average for
#                 # the betas and fs since their sum will immediately
#                 # violate our filter condition 3,4,5,6
#                 np.array((0.25, 0.2, 0.42,0.075, 0.005, 0.35, 0.12,  0.03, 0.37, 0.11, 0.01)),
#                 # but for the rest it se
#                 (c_min+c_max)[11:] / 2.0
#             ]
#         )
# )

isQualified = make_param_filter_func(c_max, c_min)
# check if the current value passes the filter
# to avoid to get stuck in an inescapable series of rejections
print (epa_0)
if not (isQualified(np.array(epa_0))):
    raise ValueError(
        """the current value does not pass filter_func. This is probably due to an initial value chosen outside the permitted range""")

# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa_0,
    filter_func=isQualified,
    param2res=param2res,
    costfunction=make_feng_cost_func(obs),
    nsimu=10000,
    c_max=c_max,
    c_min=c_min,
    acceptance_rate=10,   # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=1   # default value | increase value to reduce initial step size
)
# save the parameters and cost function values for postprocessing
pd.DataFrame(C_autostep).to_csv(dataPath.joinpath('OCN_v2_autostep_da_aa.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(dataPath.joinpath('OCN_v2_autostep_da_j_aa.csv'), sep=',')


# calculate maximum likelihood for each parameter as a peak of posterior distribution
def density_peaks(distr):
    peaks = np.zeros(distr.shape[0])
    for i in range(distr.shape[0]):
        x = distr[i, :]
        pdf = gaussian_kde(x)
        l = np.linspace(np.min(x), np.max(x), len(x))
        density = pdf(l)
        peaks[i] = l[density.argmax()]
    return peaks

max_likelihood_autostep = density_peaks(C_autostep)

# forward run with median, max likelihood, and min J parameter sets
sol_median_autostep = param2res(np.median(C_autostep, axis=1))
sol_likelihood_autostep = param2res(max_likelihood_autostep)
sol_min_J_autostep = param2res(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))])
# export solutions
pd.DataFrame(sol_median_autostep).to_csv(dataPath.joinpath('sol_median_autostep.csv'), sep=',')
pd.DataFrame(sol_likelihood_autostep).to_csv(dataPath.joinpath('sol_likelihood_autostep.csv'), sep=',')
pd.DataFrame(sol_min_J_autostep).to_csv(dataPath.joinpath('sol_min_J_autostep.csv'), sep=',')

# Demo MCMC: with a uniform proposer and fixed step
uniform_prop = make_uniform_proposer(
    c_min,
    c_max,
    D=240,
    filter_func=isQualified
)

C_demo, J_demo = mcmc(
    initial_parameters=epa_0,
    proposer=uniform_prop,
    param2res=param2res,
    # costfunction=make_weighted_cost_func(obs)
    costfunction=make_feng_cost_func(obs),
    nsimu=10000
)
# save the parameters and cost function values for postprocessing
pd.DataFrame(C_demo).to_csv(dataPath.joinpath('OCN_v2_demo_da_aa.csv'), sep=',')
pd.DataFrame(J_demo).to_csv(dataPath.joinpath('OCN_v2_demo_da_j_aa.csv'), sep=',')

max_likelihood_demo = density_peaks(C_demo)

# forward run with median, max likelihood, and min J parameter sets
sol_median_demo = param2res(np.median(C_demo, axis=1))
sol_likelihood_demo = param2res(max_likelihood_demo)
sol_min_J_demo = param2res(C_demo[:, np.where(J_demo[1] == np.min(J_demo[1]))])
# export solutions
pd.DataFrame(sol_median_demo).to_csv(dataPath.joinpath('sol_median_demo.csv'), sep=',')
pd.DataFrame(sol_likelihood_demo).to_csv(dataPath.joinpath('sol_likelihood_demo.csv'), sep=',')
pd.DataFrame(sol_min_J_demo).to_csv(dataPath.joinpath('sol_min_J_demo.csv'), sep=',')

# build a new proposer based on a multivariate_normal distribution using the estimated covariance of the previous run if available
# parameter values of the previous run

# Formal MCMC: with multivariate normal proposer based on adaptive covariance matrix
C_formal, J_formal = adaptive_mcmc(
    initial_parameters=epa_0,
    covv=np.cov(C_demo[:, int(C_demo.shape[1] / 10):]),
    filter_func=isQualified,
    param2res=param2res,
    # costfunction=make_weighted_cost_func(obs)
    costfunction=make_feng_cost_func(obs),
    nsimu=10000,
    sd_controlling_factor=1  # default value | increase value to reduce step size
)
# save the parameters and cost function values for postprocessing
pd.DataFrame(C_formal).to_csv(dataPath.joinpath('OCN_v2_formal_da_aa.csv'), sep=',')
pd.DataFrame(J_formal).to_csv(dataPath.joinpath('OCN_v2_formal_da_j_aa.csv'), sep=',')

max_likelihood_formal = density_peaks(C_formal)
# forward run with median, max likelihood, and min J parameter sets
sol_median_formal = param2res(np.median(C_formal, axis=1))
sol_likelihood_formal = param2res(max_likelihood_formal)
sol_min_J_formal = param2res(C_formal[:, np.where(J_formal[1] == np.min(J_formal[1]))])
# export solutions
pd.DataFrame(sol_median_formal).to_csv(dataPath.joinpath('sol_median_formal.csv'), sep=',')
pd.DataFrame(sol_likelihood_formal).to_csv(dataPath.joinpath('sol_likelihood_formal.csv'), sep=',')
pd.DataFrame(sol_min_J_formal).to_csv(dataPath.joinpath('sol_min_J_formal.csv'), sep=',')

# POSTPROCESSING
#
# The 'solution' of the inverse problem is actually the (joint) posterior
# probability distribution of the parameters, which we approximate by the
# histogram consisting of the mcmc generated samples.
# This joint distribution contains as much information as all its (infinitely
# many) projections to curves through the parameter space combined.
# Unfortunately, for this very reason, a joint distribution of more than two
# parameters is very difficult to visualize in its entirety.
# to do:
#   a) make a movie of color coded samples  of the a priori distribution of the parameters.
#   b) -"-                                  of the a posterior distribution -'-

# Therefore the  following visualizations have to be considered with caution:
# 1.
# The (usual) histograms of the values of a SINGLE parameters can be very
# misleading since e.g. we can not see that certain parameter combination only
# occure together. In fact this decomposition is only appropriate for
# INDEPENDENT distributions of parameters in which case the joint distribution
# would be the product of the distributions of the single parameters.  This is
# however not even to be expected if our prior probability distribution can be
# decomposed in this way. (Due to the fact that the Metropolis Hastings Alg. does not
# produce independent samples )

# df = pd.DataFrame({name: C_formal[:, i] for i, name in enumerate(EstimatedParameters._fields)})
# subplots = df.hist()
# fig = subplots[0, 0].figure
# fig.set_figwidth(15)
# fig.set_figheight(15)
# fig.savefig('histograms.pdf')

# As the next best thing we can create a matrix of plots containing all
# projections to possible  parameter tuples
# (like the pairs plot in the R package FME) but 16x16 plots are too much for one page..
# However the plot shows that we are dealing with a lot of collinearity for this  parameter set

# subplots = pd.plotting.scatter_matrix(df)
# fig = subplots[0, 0].figure
# fig.set_figwidth(15)
# fig.set_figheight(15)
# fig.savefig('scatter_matrix.pdf')

# 2.
# another way to get an idea of the quality of the parameter estimation is
# to plot trajectories.
# A possible aggregation of this histogram to a singe parameter
# vector is the mean which is an estimator of  the expected value of the
# desired distribution.

# fig = plt.figure()
# plot_solutions(
#     fig,
#     times=range(sol_median_formal.shape[0]),
#     var_names=Observables._fields,
#     tup=(sol_median_formal, obs),
#     names=('mean', 'obs')
# )
# fig.savefig('solutions.pdf')
