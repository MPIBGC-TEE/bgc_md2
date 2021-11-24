# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # yy_cable

# +
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

import bgc_md2.helper as h
import importlib 
# -

importlib.invalidate_caches()
#mod = importlib.import_module('bgc_md2.models.cable_yuanyuan.source_Luo')
mod = importlib.import_module('bgc_md2.models.cable_yuanyuan.source')
mvs = mod.mvs


h.compartmental_graph(mvs)

# +
from ComputabilityGraphs.TypeSet import TypeSet

mvs.computable_mvar_types()

# -

mvs.get_CompartmentalMatrix()

mvs.get_OutFluxesBySymbol()

# +
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
    make_param2res_sym,
    UnEstimatedParameters,
    EstimatedParameters,
    Observables
)

from general_helpers import (
        make_uniform_proposer,
        make_multivariate_normal_proposer,
        #mcmc,
        make_feng_cost_func,
        mcmc,
        adaptive_mcmc,
        plot_solutions
)

with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

dataPath = Path(conf_dict['dataPath'])

npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(dataPath)

# combine them to a single array which we will later use as input to the costfunction
#nyears=140
nyears = 10
tot_len = 12*nyears
obs = np.stack([cleaf, croot, cwood, clitter, csoil, rh], axis=1)[0:tot_len,:]

c_min=np.array([0.09, 0.09,0.09,0.01,0.01,  1/(2*365), 1/(365*10), 1/(60*365), 0.1/(0.1*365), 0.06/(0.137*365), 0.06/(5*365), 0.06/(222.22*365),clitter[0]/100,clitter[0]/100,csoil[0]/100,csoil[0]/2])
c_max=np.array([1   ,    1,0.21,   1,   1,1/(0.3*365),1/(0.8*365),      1/365,   1/(365*0.1),  0.6/(365*0.137),  0.6/(365*5),  0.6/(222.22*365),    clitter[0],    clitter[0],  csoil[0]/3,csoil[0]])

isQualified = make_param_filter_func(c_max,c_min)
uniform_prop = make_uniform_proposer(
    c_min,
    c_max,
    D=100.0, # this value 
    filter_func=isQualified
)

cpa = UnEstimatedParameters(
    C_leaf_0=cleaf[0],
    C_root_0=croot[0],
    C_wood_0=cwood[0],
    clitter_0=clitter[0],
    csoil_0=csoil[0],
    rh_0 = rh[0],
    npp=npp,
    number_of_months=tot_len,
    clay=0.2028,
    silt=0.2808,
    lig_wood=0.4,
    f_wood2CWD=1, 
    f_metlit2mic=0.45
)
param2res = make_param2res_sym(cpa) 

# +
epa_0 = EstimatedParameters(
    beta_leaf=0.15,
    beta_root=0.2,
    lig_leaf=0.15,
    f_leaf2metlit=0.28,
    f_root2metlit=0.6,
    k_leaf=1/365,
    k_root=1/(365*5),
    k_wood=1/(365*40),
    k_metlit=0.5/(365*0.1),
    k_mic=0.3/(365*0.137),
    k_slowsom=0.3/(365*5),
    k_passsom=0.3/(222.22*365),
    C_metlit_0=0.05,
    C_CWD_0=0.1,
    C_mic_0=1,
    C_passom_0=5,
)
# it is sensible to use the same costfunction for both the demo and
# the formal run so we define it here for both
costfunction=make_feng_cost_func(obs)
#costfunction=make_weighted_cost_fun c(obs)
demo_aa_path = dataPath.joinpath('cable_demo_da_aa.csv')
demo_aa_j_path = dataPath.joinpath('cable_demo_da_j_aa.csv')

C_demo, J_demo = mcmc(
        initial_parameters=epa_0,
        proposer=uniform_prop,
        param2res=param2res,
        costfunction=costfunction,
        nsimu=200
)
#save the parameters and costfunctionvalues for postprocessing 
pd.DataFrame(C_demo).to_csv(demo_aa_path,sep=',')
pd.DataFrame(J_demo).to_csv(demo_aa_j_path,sep=',')

# +
# build a new proposer based on a multivariate_normal distribution using the
# estimated covariance of the previous run if available first we check how many
# accepted parameters we got 
# and then use part of them to compute a covariance matrix for the 
# formal run
covv = np.cov(C_demo[:, int(C_demo.shape[1]/10):]) 


normal_prop = make_multivariate_normal_proposer(
    covv = covv,
    filter_func=isQualified
)
# Look for data from the formal run and use it  for postprocessing 
formal_aa_path = dataPath.joinpath('cable_formal_da_aa.csv')
formal_aa_j_path = dataPath.joinpath('cable_formal_da_j_aa.csv')
C_formal, J_formal = mcmc(
        initial_parameters=epa_0,
        proposer=normal_prop,
        param2res=param2res,
        costfunction=costfunction,
        nsimu=200
)
pd.DataFrame(C_formal).to_csv(formal_aa_path,sep=',')
pd.DataFrame(J_formal).to_csv(formal_aa_j_path,sep=',')

# +
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
# -
from ComputabilityGraphs.helpers import all_mvars
from ComputabilityGraphs.TypeSet import TypeSet
for t in all_mvars(h.bgc_md2_computers()):
    print(t.__name__)

mvs.get_StateVariableTupleTimeDerivative()


