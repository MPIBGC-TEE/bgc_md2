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
    make_weighted_cost_func,
    make_param2res,
    make_param2res_2,
    UnEstimatedParameters,
    EstimatedParameters,
    Observables
)

from general_helpers import (
        make_uniform_proposer,
        make_multivariate_normal_proposer,
        mcmc,
        make_feng_cost_func,
        plot_solutions
)

with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

dataPath = Path(conf_dict['dataPath'])

npp, rh, clitter, ccwd, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(dataPath)

#nyears=320
nyears = 10
tot_len = 12*nyears
obs_tup=Observables(
    C_leaf=cleaf,
    C_root=croot,
    C_wood=cwood,
    c_litter=clitter,
    c_soil=csoil,
    respiration=rh
)
obs = np.stack(obs_tup, axis=1)[0:tot_len,:]

cpa = UnEstimatedParameters(
    C_leaf_0=cleaf[0],
    C_root_0=croot[0],
    C_wood_0=cwood[0],
    C_cwd_0=ccwd[0],
    c_litter_0=clitter[0],
    c_soil_0=csoil[0],
    rh_0 = rh[0],
    npp=npp,
    number_of_months=tot_len,
    clay=0.2028,
    silt=0.2808,
    lig_wood=0.4,
    f_wood2CWD=1, 
    f_metlit2mic=0.45
)
param2res = make_param2res(cpa) 


c_min = np.array(
    EstimatedParameters(
        beta_leaf=0.09,
        beta_root=0.09,
        lig_leaf=0.09,
        f_leaf2metlit=0.01,
        f_root2metlit=0.01,
        k_leaf=1/(2*365),
        k_root=1/(365*10),
        k_wood=1/(60*365),
        k_metlit=0.1/(0.1*365),
        k_mic=0.06/(0.137*365),
        k_slowsom=0.06/(5*365),
        k_passsom=0.06/(222.22*365),
        C_metlit_0=clitter[0]/100,
        C_strlit_0=clitter[0]/100,
        C_mic_0=csoil[0]/100,
        C_passom_0=csoil[0]/2
    )
)

c_max = np.array(
    EstimatedParameters(
        beta_leaf=1,               
        beta_root=1,		    
        lig_leaf=0.21,		    
        f_leaf2metlit=1,		    
        f_root2metlit=1,		    
        k_leaf=1/(0.3*365),	    
        k_root=1/(0.8*365),	    
        k_wood=1/365,		    
        k_metlit=1/(365*0.1),	    
        k_mic=0.6/(365*0.137),    
        k_slowsom=0.6/(365*5),	    
        k_passsom=0.6/(222.22*365),   
        C_metlit_0=clitter[0],	    
        C_strlit_0=clitter[0],	    
        C_mic_0=csoil[0]/3,	    
        C_passom_0=csoil[0]/2
    )
)
# I commented your original settings and instead used the 
# values constructed from the limits (see below)
#epa_0 = EstimatedParameters(
#    beta_leaf=0.15,
#    beta_root=0.2,
#    lig_leaf=0.15,
#    f_leaf2metlit=0.28,
#    f_root2metlit=0.6,
#    k_leaf=1/365,
#    k_root=1/(365*5),
#    k_wood=1/(365*40),
#    k_metlit=0.5/(365*0.1),
#    k_mic=0.3/(365*0.137),
#    k_slowsom=0.3/(365*5),
#    k_passsom=0.3/(222.22*365),
#    C_metlit_0=0.05,
#    C_strlit_0=0.1,
#    C_mic_0=1,
#    C_passom_0=10,
#)
epa_0 = EstimatedParameters._make( 
        np.concatenate(
            [   
                # we dont want to use the average for 
                # the betas since their sum will immidiately
                # violate our filter condition 3
                np.array((0.15, 0.2)), 
                # but for the rest it se
                (c_min+c_max)[2:] / 2.0
            ]
        )
)

isQualified = make_param_filter_func(c_max,c_min)
# check if the current value passes the filter
# to avoid to get stuck in an inescapable series of rejections  
if not(isQualified(np.array(epa_0))):
    raise ValueError("""the current value does not pass filter_func. This is probably due to an initial value chosen outside the permitted range""")

uniform_prop = make_uniform_proposer(
    c_min,
    c_max,
    D=30,
    filter_func=isQualified
)

C_demo, J_demo = mcmc(
        initial_parameters=epa_0,
        proposer=uniform_prop,
        param2res=param2res,
        #costfunction=make_weighted_cost_func(obs)
        costfunction=make_feng_cost_func(obs),
        nsimu=100
)
# save the parameters and costfunctionvalues for postprocessing 
pd.DataFrame(C_demo).to_csv(dataPath.joinpath('cable_demo_da_aa.csv'),sep=',')
pd.DataFrame(J_demo).to_csv(dataPath.joinpath('cable_demo_da_j_aa.csv'),sep=',')

# build a new proposer based on a multivariate_normal distribution using the estimated covariance of the previous run if available
# parameter values of the previous run

C_formal, J_formal = adaptive_mcmc(
        initial_parameters=epa_0,
        covv=np.cov(C_demo[:, int(C_demo.shape[9]/10):]),
        filter_func = isQualified, 
        param2res=param2res,
        #costfunction=make_weighted_cost_func(obs)
        costfunction=make_feng_cost_func(obs),
        #nsimu=20000
        nsimu=100
)
formal_aa_path = dataPath.joinpath('cable_formal_da_aa.csv')
formal_aa_j_path = dataPath.joinpath('cable_formal_da_j_aa.csv')
pd.DataFrame(C_formal).to_csv(formal_aa_path,sep=',')
pd.DataFrame(J_formal).to_csv(formal_aa_j_path,sep=',')

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
sol_mean =param2res(np.mean(C_formal,axis=1))

fig = plt.figure()
plot_solutions(
        fig,
        times=range(sol_mean.shape[0]),
        var_names=Observables._fields,
        tup=(sol_mean, obs),
        names=('mean','obs')
)
fig.savefig('solutions.pdf')
