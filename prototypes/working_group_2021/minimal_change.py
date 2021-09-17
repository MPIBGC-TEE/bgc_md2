#!/usr/bin/env python

# Throuout the whole script there are '# fixme' comment
# you can search for them to find places where you have to adapt the code
# to your needs. Hopefully I found
# 
import os
import sys
import glob
import datetime as dt
from typing_extensions import Concatenate
import numpy as np
import math
import pandas as pd
import calendar
import matplotlib 
# use platform independent path descriptions so that you can run
# your stuff on windows or linux or mac
import pathlib

# I use the non mandatory type hints to make the code more readable
# This is escpecially useful for the description of functions that
# are used as arguments for other functions
from typing import Callable, Tuple 

from matplotlib import pyplot as plt

# fixme: 
#   The idea is that everybody writes her own version of this 
#   module 'yy_cable_specific_helpers.py' maybe
#   - 'jon_ybis_specific_helpers.py' or 
#   - 'mm_cable_specific helpers' and so on. 
#   All functions in module yy_cable_specific_helpers provide model specific results
#   and can not be applied directly but serve as examples.
#   If you use the same naming convention inside the module you only have to
#   change the name of the module in the next line

from yy_cable_specific_helpers import (
    get_example_site_vars,
    make_param_filter_func,
    make_weighted_cost_func,
    make_param2res,
    make_param2res_2,
    UnEstimatedParameters,
    EstimatedParameters
)

from general_helpers import (
        make_uniform_proposer,
        make_multivariate_normal_proposer,
        mcmc,
        make_feng_cost_func
)

# fixme: 
#   put the (relative or asolute) location of your data here
#   (this example uses an absolute path starting with a '/'
dataPath = pathlib.Path('/home/data/yuanyuan') 

# fixme: 
#    Note that the function is imported from 
#    yy_cable_specific_helpers which means that you have to provide
#    your version of this fuction which will most likely return different
#    variables 
npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(dataPath)

# combine them to a single array which we will later use as input to the costfunction
#nyears=140
nyears = 10
tot_len = 12*nyears
obs = np.stack([cleaf, croot, cwood, clitter, csoil, rh], axis=1)[0:tot_len,:]


# leaf, root , wood, metabolic, structural, CWD, microbial, slow, passive 



# fixme 
c_min=np.array([0.09, 0.09,0.09,0.01,0.01,  1/(2*365), 1/(365*10), 1/(60*365), 0.1/(0.1*365), 0.06/(0.137*365), 0.06/(5*365), 0.06/(222.22*365),clitter[0]/100,clitter[0]/100,csoil[0]/100,csoil[0]/2])
c_max=np.array([1   ,    1,0.21,   1,   1,1/(0.3*365),1/(0.8*365),      1/365,   1/(365*0.1),  0.6/(365*0.137),  0.6/(365*5),  0.6/(222.22*365),    clitter[0],    clitter[0],  csoil[0]/3,csoil[0]])

# fixme
#   this function is model specific: It discards parameter proposals
#   where beta1 and beta2 add up to more than 0.99
isQualified = make_param_filter_func(c_max,c_min)
uniform_prop = make_uniform_proposer(
    c_min,
    c_max,
    D=10.0,
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
param2res = make_param2res(cpa) #pa=[beta1,beta2, lig_leaf, f41,f42, kleaf,kroot,kwood,kmet,kmic, kslow,kpass, cmet_init, cstr_init, cmic_init, cpassive_init ]
pa=            [0.15,  0.2,0.15,0.28, 0.6,      1/365,  1/(365*5), 1/(365*40), 0.5/(365*0.1),  0.3/(365*0.137),  0.3/(365*5),  0.3/(222.22*365),          0.05,           0.1,           1,         5]
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
    CWD_0=0.1,
    C_mic_0=1,
    C_passom_0=5,
)
nsimu_demo = 2000    
C_demo, J_demo = mcmc(
        initial_parameters=epa_0,
        proposer=uniform_prop,
        param2res=param2res,
        #costfunction=make_weighted_cost_func(obs)
        costfunction=make_feng_cost_func(obs),
        nsimu=nsimu_demo
)
# save the parameters and costfunctionvalues for postprocessing 
pd.DataFrame(C_demo).to_csv(dataPath.joinpath('cable_demo_da_aa.csv'),sep=',')
pd.DataFrame(J_demo).to_csv(dataPath.joinpath('cable_demo_da_j_aa.csv'),sep=',')

# build a new proposer based on a multivariate_normal distribution using the estimated covariance of the previous run if available
# parameter values of the previous run

#from IPython import embed; embed()
normal_prop = make_multivariate_normal_proposer(
    covv = np.cov(C_demo[:, int(nsimu_demo/10):]), # the part of the demo run samples to use (here the last 90%) 
    filter_func=isQualified
)
C_formal, J_formal = mcmc(
        initial_parameters=epa_0,
        proposer=normal_prop,
        param2res=param2res,
        #costfunction=make_weighted_cost_func(obs)
        costfunction=make_feng_cost_func(obs),
        #nsimu=20000
        nsimu=2000
)

