print("zhoujian...")
#!/usr/bin/env python
# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
sys.path.insert(0,'..')

# Throuout the whole script there are '# fixme' comment
# you can search for them to find places where you have to adapt the code
# to your needs. Hopefully I found
# 
#import os
#import glob
#import datetime as dt
#from typing_extensions import Concatenate
import numpy as np
#import math
import pandas as pd
#import calendar
import matplotlib.pyplot as plt
# use platform independent path descriptions so that you can run
# your stuff on windows or linux or mac
from pathlib import Path
import json 

# I use the non mandatory type hints to make the code more readable
# This is escpecially useful for the description of functions that
# are used as arguments for other functions
#from typing import Callable, Tuple 

#from matplotlib import pyplot as plt

# fixme: 
#   The idea is that everybody writes her own version of this 
#   directory 'yy_cable' maybe
#   - 'jon_ybis' or 
#   - 'mm_cable' and so on. 
#   All functions in module model_specific_helpers.py provide model specific results
#   and can not be applied directly but serve as examples.

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
# from general_helpers import mcmc_parallel as mcmc

# fixme: 
#   put the (relative or asolute) location of your data into a small file called 'config.json' and
#   in my case the content looks like this:
#   {"dataPath": "/home/data/yuanyuan"}
#   DO NOT add the file to the repository. It is not only model- but also site specific. 
#   So you are likely to have one for every model on every computer
#   you run this code on.
#   (this example uses an absolute path starting with a '/'
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

dataPath = Path(conf_dict['dataPath'])

# fixme: 
#    Note that the function is imported from 
#    model_specific_helpers which means that you have to provide
#    your version of this fuction which will most likely return different
#    variables 
npp, rh, clitter, csoil, cveg = get_example_site_vars(dataPath)

print('npp.shape=',npp.shape, 'rh.shape=',rh.shape, 'clitter.shape=',clitter.shape, 'csoil.shape=',csoil.shape, 'cveg.shape=',cveg.shape)
# combine them to a single array which we will later use as input to the costfunction
#nyears=140
nyears = 320
tot_len = 12*nyears
rh1 = [0.0 for d in range(nyears)]

# print('len(rh)=',len(rh))
# print('len(rh1)=',len(rh1))

for i in range(len(rh1)):
    #print('i=',i,'(i+1)*12=',(i+1)*12)
    rh1[i] = np.mean(rh[i*12:(i+1)*12])

# Comment by Chenyu Bian due to the cpool and cflux does not have the same shape
obs = np.stack([cveg, clitter, csoil, rh1], axis=1)[0:nyears,:]  # Zhou

# leaf, root , wood, metabolic, structural, CWD, microbial, slow, passive 

# fixme 
# c_min=np.array([0.09, 0.09,0.09,0.01,0.01,  1/(2*365), 1/(365*10), 1/(60*365), 0.1/(0.1*365), 0.06/(0.137*365), 0.06/(5*365), 0.06/(222.22*365),clitter[0]/100,clitter[0]/100,csoil[0]/100,csoil[0]/2])
# c_max=np.array([1   ,    1,0.21,   1,   1,1/(0.3*365),1/(0.8*365),      1/365,   1/(365*0.1),  0.6/(365*0.137),  0.6/(365*5),  0.6/(222.22*365),    clitter[0],    clitter[0],  csoil[0]/3,csoil[0]])

c_min = np.array([0.01, 0.01,
                 0.2, 0.01, 0.01, 0.01, 0.01, 0.01,
                 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
                 # tau_leaf, tau_wood, tau_root, k12, k13, k14, k15
                  0.125, 5.0, 0.3, 1.0 * \
                  10**(-6), 1.0 * 10**(-7), 1.0 * 10**(-8), 1.0 * 10**(-8),
                  0.2, 0.2,  # decompl,decomps,
                  #  1.0 * 10**(-4), 1.0 * 10**(-6), 1.0 * 10**(-6), 0.2, 0.2, #k_leaf, k_wood, k_root, decompl, decompl
                  cveg[0]/10, cveg[0]/10,  # cWood, cRoot
                  # cLMetaLeaf, cLMetaWood, cLStrLeaf, cLStrWood, cLLigLeaf
                  clitter[0]/100, clitter[0]/100, clitter[0]/100,
                   clitter[0]/100, clitter[0]/100,
                  csoil[0]/100, csoil[0]/100, csoil[0]/100, csoil[0]/1000, csoil[0]/50, csoil[0]/50])  # cLMetaRoot, cLStrRoot, cLLigRoot, cSMiroc, cSPort, cSNonprot

c_max = np.array([0.5, 0.99,
                 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1,
                 # tau_leaf, tau_wood, tau_root, k12, k13, k14, k15
                  2.5, 45, 3,  0.2, 0.02, 0.03, 0.8,
                  0.8, 0.9,  # decompl, decomps
                  #   0.2, 0.02, 0.03, 0.8, 0.9,  # k_leaf, k_wood, k_root, decompl, decompl
                  120, 80,  # cWood, cRoot
                  clitter[0]/3,  clitter[0]/2, clitter[0]/3, 
                  clitter[0]/2, clitter[0]/3,
                  csoil[0]/3, csoil[0]/3, csoil[0]/3, csoil[0]/5, csoil[0], csoil[0]])

# fixme
#   this function is model specific: It discards parameter proposals
#   where beta1 and beta2 add up to more than 0.99
# change by zhoujian
cpa = UnEstimatedParameters(
    cveg_0 = cveg[0],
    clitter_0=clitter[0],
    csoil_0=csoil[0],
    rh_0 = rh[0],
    npp=npp,
    lig_frac = 0.5,
    number_of_months=tot_len,
)
isQualified = make_param_filter_func(c_max,c_min,cpa)
uniform_prop = make_uniform_proposer(
    c_min,
    c_max,
    D=100.0, # this value 
    filter_func=isQualified
)

# cpa = UnEstimatedParameters(
#     cveg_0 = cveg[0],
#     clitter_0=clitter[0],
#     csoil_0=csoil[0],
#     rh_0 = rh[0],
#     npp=npp,
#     lig_frac = 0.5,
#     number_of_months=tot_len,
# )
param2res = make_param2res(cpa) 
epa_0 = EstimatedParameters(
    # for ibis
    beta_leaf = 0.2,
    beta_wood = 0.5,
    f_leaf2metLeaflit = 0.4,   #
    f_wood2metWoodlit = 0.3,   #
    f_root2metRootlit = 0.3,   #
    f_leaf2strLeaflit = 0.35,  #
    f_wood2strWoodlit = 0.35,  # 6,  f82
    f_root2strRootlit=0.3,   # 7,  f93
    f_mic2prot=0.5,          # 8,  f14_13
    f_mic2nonprot=0.5,       # 9,  f15_13
    f_prot2mic=0.5,          # 10, f13_14
    f_prot2pass=0.5,         # 11, f16_14
    f_nonprot2mic=0.5,       # 12, f13_15
    f_nonprot2pass=0.5,      # 13, f16_15
    tau_leaf=0.5,            # 14
    tau_wood=15.0,           # 15
    tau_root=0.7,            # 16
    k_mic=0.001,               # 17
    k_protsom=0.0001,          # 18
    k_nonprotsom=0.0001,       # 19
    k_passsom=0.0001,          # 20
    decompl=0.45,              # 21
    decomps=0.5,               # 22
    c_wood0=80,                # 23
    c_root0 = 48,              # 24
    c_LMetaLeaf0=2.3,        # 25
    c_LMetaWood0=2.0,        # 26
    c_LStrLeaf0=1.5,         # 27
    c_LStrWood0=3.0,         # 28
    c_LLigLeaf0=2.0,         # 29
    c_LMetaRoot0=2.0,        # 30
    c_LStrRoot0=4.0,         # 31
    c_LLigRoot0=8.0,         # 32
    c_SMiroc0=2.0,           # 33
    c_SProt0=20,             # 34
    c_Snonprot0=100         # 35
)
# it is sensible to use the same costfunction for both the demo and
# the formal run so we define it here for both
# costfunction=make_feng_cost_func(obs)

# comment this line by bian
# costfunction=make_weighted_cost_func(obs)
costfunction = make_weighted_cost_func(cveg, clitter, csoil, rh)


# Look for data from the demo run and use it to compute the covariance matrix if necessarry
demo_aa_path = dataPath.joinpath('ibis_demo_da_aa.csv')
demo_aa_j_path = dataPath.joinpath('ibis_demo_da_j_aa.csv')
if not demo_aa_path.exists():

    print("Did not find demo run results. Will perform  demo run")
    C_demo, J_demo = mcmc(
            initial_parameters=c_min+0.01*(c_max-c_min),#epa_0,
            proposer=uniform_prop,
            param2res=param2res,
            costfunction=costfunction,
            nsimu=1000
    )
    # print('C_demo=',C_demo[0:5])
    # print('J_demo=',J_demo[0:5])

    # save the parameters and costfunctionvalues for postprocessing 
    pd.DataFrame(C_demo).to_csv(demo_aa_path,sep=',')
    pd.DataFrame(J_demo).to_csv(demo_aa_j_path,sep=',')
    print("output demo file ...")
else:
    print("""Found {p} from a previous demo run. 
    If you also want to recreate the demo output then move the file!
    """.format(p = demo_aa_path)) 
    C_demo = pd.read_csv(demo_aa_path).to_numpy()
    J_demo = pd.read_csv(demo_aa_j_path).to_numpy() 

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
formal_aa_path = dataPath.joinpath('ibis_formal_da_aa.csv')
formal_aa_j_path = dataPath.joinpath('ibis_formal_da_j_aa.csv')
if not formal_aa_path.exists():
    print("Did not find results. Will perform formal run")
    C_formal, J_formal = mcmc(
            initial_parameters=c_min+0.01*(c_max-c_min),  # epa_0,
            proposer=uniform_prop,  # normal_prop,
            param2res=param2res,
            costfunction=costfunction,
            nsimu=1000
    )
    pd.DataFrame(C_formal).to_csv(formal_aa_path,sep=',')
    pd.DataFrame(J_formal).to_csv(formal_aa_j_path,sep=',')
    print("output formal file")

else:
    print("""Found {p} from a previous demo run. 
If you also want recreate the output then move the file!
""".format(p = formal_aa_path)) 
    C_formal = pd.read_csv(formal_aa_path).to_numpy()
    J_formal = pd.read_csv(formal_aa_j_path).to_numpy() 

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
# print('type(sol_mean)=', type(sol_mean))
# print('sol_mean',sol_mean)

sol_mean_array = np.asarray(sol_mean)
# print('type(sol_mean_array)=', type(sol_mean_array))
# print('sol_mean_array.shape=', sol_mean_array.shape)
# print('sol_mean_array.size=',sol_mean_array.size)
# print('line 315')
fig = plt.figure()

newSimuRh = np.zeros((nyears,1))
for i in range(nyears):
    newSimuRh[i] = np.mean(sol_mean[3][i*12:(i+1)*12])
sol_mean_new = np.squeeze(np.stack([sol_mean[0], sol_mean[1], sol_mean[2], newSimuRh], axis=1))

# print("zhoujian: shapebet", sol_mean_new.shape, obs.shape)
plot_solutions(
        fig,
        times=range(sol_mean_new.shape[0]),
        var_names=Observables._fields,
        tup=(sol_mean_new, obs),
        names=('mean','obs')
)
fig.savefig('solutions.pdf')


