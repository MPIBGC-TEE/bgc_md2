#!/usr/bin/env python
import sys
sys.path.insert(0,'..')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import json 


from model_specific_helpers import (
    get_example_site_vars,
    make_param_filter_func,
    #make_weighted_cost_func,
    make_param2res,
    #make_param2res_2,
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
    conf_dict=json.load(f) 

dataPath = Path(conf_dict['dataPath'])

npp, C_leaf, C_wood, C_root, C_litter_above, C_litter_below, C_fast_som, C_slow_som, C_pass_som, \
rh, f_veg2litter, f_litter2som, mrso, tsl = get_example_site_vars(dataPath)

nyears=150
#nyears = 10
tot_len = 12*nyears
obs_tup=Observables(
    C_leaf=C_leaf,
    C_wood=C_wood,
    C_root=C_root,
    C_litter_above=C_litter_above,
    C_litter_below=C_litter_below,
    C_fast_som=C_fast_som,
    C_slow_som=C_slow_som,
    C_pass_som=C_pass_som,
    rh=rh,
    f_veg2litter=f_veg2litter,
    f_litter2som=f_litter2som
)
obs = np.stack(obs_tup, axis=1)[0:tot_len,:]

# save observational data for comparison with model output
pd.DataFrame(obs).to_csv(dataPath.joinpath('obs.csv'),sep=',')

cpa = UnEstimatedParameters(
    C_leaf_0=C_leaf[0],
    C_wood_0=C_wood[0],
    C_root_0=C_root[0],
    C_litter_above_0=C_litter_above[0],
    C_litter_below_0=C_litter_below[0],
    C_fast_som_0=C_fast_som[0],
    C_slow_som_0=C_slow_som[0],
    C_pass_som_0=C_pass_som[0],
    rh_0=rh[0],
    f_veg_lit_0=f_veg2litter[0],
    f_lit_soil_0=f_litter2som[0],
    npp=npp,
    number_of_months=tot_len,
    mrso=mrso,
    tsl=tsl
)
param2res = make_param2res(cpa) 

c_min = np.array(
    EstimatedParameters(
        beta_leaf=0,
        beta_wood=0,
        f_leaf_lit2fast_som=0.1,
        f_leaf_lit2slow_som=0.01,
        f_leaf_lit2pass_som=0.001,
        f_wood_lit2fast_som=0.1,
        f_wood_lit2slow_som=0.01,
        f_wood_lit2pass_som=0.001,
        f_root_lit2fast_som=0.1,
        f_root_lit2slow_som=0.01,
        f_root_lit2pass_som=0.001,
        k_leaf=1/(365*2),
        k_wood=1/(365*60),
        k_root=1/(365*30),
        k_leaf_lit=1/(365*60),
        k_wood_lit=1/(365*30),
        k_root_lit=1/(365*30),
        k_fast_som=1/(365*200),
        k_slow_som=1/(365*500),
        k_pass_som=1/(365*1000),
        C_leaf_lit_0=0,
        T_0=-10,
        E=1,
        KM=1
    )
)

c_max = np.array(
    EstimatedParameters(
        beta_leaf=1,
        beta_wood=1,
        f_leaf_lit2fast_som=0.9,
        f_leaf_lit2slow_som=0.9,
        f_leaf_lit2pass_som=0.9,
        f_wood_lit2fast_som=0.9,
        f_wood_lit2slow_som=0.2,
        f_wood_lit2pass_som=0.2,
        f_root_lit2fast_som=0.9,
        f_root_lit2slow_som=0.9,
        f_root_lit2pass_som=0.9,
        k_leaf=1/30,
        k_wood=1/(365*1),
        k_root=1/(365*0.5),
        k_leaf_lit=1/(365*1),
        k_wood_lit=1/(365*0.5),
        k_root_lit=1/(365*0.5),
        k_fast_som=1/(365*1),
        k_slow_som=1/(365*3.5),
        k_pass_som=1/(365*10),
        C_leaf_lit_0=0.4,
        T_0=4,
        E=100,
        KM=100
    )
)
# I commented your original settings and instead used the 
# values constructed from the limits (see below)
epa_0 = EstimatedParameters(
    beta_leaf=0.6,    #  0 (parameters used in original code)
    beta_wood=0.25,    #  1
    f_leaf_lit2fast_som=0.41,  #  2
    f_leaf_lit2slow_som=0.07,#  3
    f_leaf_lit2pass_som=0.02,#  4
    f_wood_lit2fast_som=0.30,  #  5
    f_wood_lit2slow_som=0.12,#  6
    f_wood_lit2pass_som=0.08,#  7
    f_root_lit2fast_som=0.30,  #  8
    f_root_lit2slow_som=0.14,#  9
    f_root_lit2pass_som=0.07,#  10
    k_leaf=1/(60*2),       #  11
    k_wood=1/(365*30),       #  12
    k_root=1/(365*22),       #  13
    k_leaf_lit=1/(365*3.3),	#  14
    k_wood_lit=1/(365*11),	#  15
    k_root_lit=1/(365*11),	#  16
    k_fast_som=1/(365*18),	#  17
    k_slow_som=1/(365*100),	# 18
    k_pass_som=1/(365*350),	# 19
    C_leaf_lit_0=0.3,	# 20
    T_0=2,	# 21
    E=4,	# 22
    KM=10  # 23
)

pd.DataFrame(epa_0).to_csv(dataPath.joinpath('epa_0.csv'),sep=',')

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

isQualified = make_param_filter_func(c_max,c_min)
# check if the current value passes the filter
# to avoid to get stuck in an inescapable series of rejections  
if not(isQualified(np.array(epa_0))):
    raise ValueError("""the current value does not pass filter_func. This is probably due to an initial value chosen outside the permitted range""")

# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
        initial_parameters=epa_0,
        filter_func=isQualified,
        param2res=param2res,
        costfunction=make_feng_cost_func(obs),
        nsimu=10000,
        c_max=c_max,
        c_min=c_min
)
# save the parameters and cost function values for postprocessing
pd.DataFrame(C_autostep).to_csv(dataPath.joinpath('visit_autostep_da_aa.csv'),sep=',')
pd.DataFrame(J_autostep).to_csv(dataPath.joinpath('visit_autostep_da_j_aa.csv'),sep=',')
# forward run with median and min J parameter sets
sol_median_autostep =param2res(np.median(C_autostep,axis=1))
sol_min_J_autostep =param2res(C_autostep[:,np.max(np.where(J_autostep==np.min(J_autostep)))])
# export solutions
pd.DataFrame(sol_median_autostep).to_csv(dataPath.joinpath('sol_median_autostep.csv'),sep=',')
pd.DataFrame(sol_min_J_autostep).to_csv(dataPath.joinpath('sol_min_J_autostep.csv'),sep=',')

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
        #costfunction=make_weighted_cost_func(obs)
        costfunction=make_feng_cost_func(obs),
        nsimu=10000
)
# save the parameters and cost function values for postprocessing
pd.DataFrame(C_demo).to_csv(dataPath.joinpath('visit_demo_da_aa.csv'),sep=',')
pd.DataFrame(J_demo).to_csv(dataPath.joinpath('visit_demo_da_j_aa.csv'),sep=',')
# forward run with median and min J parameter sets
sol_median_demo =param2res(np.median(C_demo,axis=1))
sol_min_J_demo =param2res(C_demo[:,np.max(np.where(J_demo==np.min(J_demo)))])
# export solutions
pd.DataFrame(sol_median_demo).to_csv(dataPath.joinpath('sol_median_demo.csv'),sep=',')
pd.DataFrame(sol_min_J_demo).to_csv(dataPath.joinpath('sol_min_J_demo.csv'),sep=',')

# build a new proposer based on a multivariate_normal distribution using the estimated covariance of the previous run if available
# parameter values of the previous run

# Formal MCMC: with multivariate normal proposer based on adaptive covariance matrix
C_formal, J_formal = adaptive_mcmc(
        initial_parameters=epa_0,
        covv=np.cov(C_demo[:, int(C_demo.shape[1]/10):]),
        filter_func = isQualified, 
        param2res=param2res,
        #costfunction=make_weighted_cost_func(obs)
        costfunction=make_feng_cost_func(obs),
        #nsimu=20000
        nsimu=10000
)
# save the parameters and cost function values for postprocessing
pd.DataFrame(C_formal).to_csv(dataPath.joinpath('visit_formal_da_aa.csv'),sep=',')
pd.DataFrame(J_formal).to_csv(dataPath.joinpath('visit_formal_da_j_aa.csv'),sep=',')
# forward run with median and min J parameter sets
sol_median_formal =param2res(np.median(C_formal,axis=1))
sol_min_J_formal =param2res(C_formal[:,np.max(np.where(J_formal==np.min(J_formal)))])
# export solutions
pd.DataFrame(sol_median_formal).to_csv(dataPath.joinpath('sol_median_formal.csv'),sep=',')
pd.DataFrame(sol_min_J_formal).to_csv(dataPath.joinpath('sol_min_J_formal.csv'),sep=',')

# POSTPROCESSING
#
# The 'solution' of the inverse problem is actually the (joint) posterior
# probability distribution of the parameters, which we approximate by the
# histogram consisting of the mcmc generated samples.  
# This joint distribution contains as much information as all its (infinitly
# many) projections to curves through the parameter space combined.
# Unfortunately, for this very reason, a joint distribution of more than two
# parameters is very difficult to visualize in its entirity. 
# to do: 
#   a) make a movie of color coded samples  of the a priori distribution of the parameters.
#   b) -"-                                  of the a posteriory distribution -'- 

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
df = pd.DataFrame({name :C_formal[:,i] for i,name in enumerate(EstimatedParameters._fields)})
subplots=df.hist()
fig=subplots[0,0].figure
fig.set_figwidth(15)
fig.set_figheight(15)
fig.savefig('histograms.pdf')

# As the next best thing we can create a matrix of plots containing all 
# projections to possible  parameter tuples
# (like the pairs plot in the R package FME) but 16x16 plots are too much for one page..
# However the plot shows that we are dealing with a lot of colinearity for this  parameter set
subplots = pd.plotting.scatter_matrix(df) 
fig=subplots[0,0].figure
fig.set_figwidth(15)
fig.set_figheight(15)
fig.savefig('scatter_matrix.pdf')


# 2.
# another way to get an idea of the quality of the parameter estimation is
# to plot trajectories.
# A possible aggregation of this histogram to a singe parameter
# vector is the mean which is an estimator of  the expected value of the
# desired distribution.

fig = plt.figure()
plot_solutions(
        fig,
        times=range(sol_median_formal.shape[0]),
        var_names=Observables._fields,
        tup=(sol_median_formal, obs),
        names=('mean','obs')
)
fig.savefig('solutions.pdf')

# additional output
#sol_init_par =param2res(epa_0)
#sol_mean =param2res(np.mean(C_formal,axis=1))
#pd.DataFrame(sol_init_par).to_csv(sol_init_par_path,sep=',')
#pd.DataFrame(sol_mean).to_csv(sol_mean_path,sep=',')