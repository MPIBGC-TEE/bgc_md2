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
#   and can not be applied generally 
from yy_cable_specific_helpers import (
    get_example_site_vars,
    make_param_filter_func,
    make_weighted_cost_func,
    make_matrix_simu
)

from general_helpers import (
        make_uniform_proposer,
        mcmc
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
# from IPython import embed; embed()
obs = np.stack([cleaf, croot, cwood, clitter, csoil, rh], axis=1)



# leaf, root , wood, metabolic, structural, CWD, microbial, slow, passive 




c_min=np.array([0.09, 0.09,0.09,0.01,0.01,  1/(2*365), 1/(365*10), 1/(60*365), 0.1/(0.1*365), 0.06/(0.137*365), 0.06/(5*365), 0.06/(222.22*365),clitter[0]/100,clitter[0]/100,csoil[0]/100,csoil[0]/2])
c_max=np.array([1   ,    1,0.21,   1,   1,1/(0.3*365),1/(0.8*365),      1/365,   1/(365*0.1),  0.6/(365*0.137),  0.6/(365*5),  0.6/(222.22*365),    clitter[0],    clitter[0],  csoil[0]/3,csoil[0]])

# fixme
#   this function is model specific: It discards parameter proposals
#   where beta1 and beta2 
isQualified = make_param_filter_func(c_max,c_min)
GenerateParamValues = make_uniform_proposer(c_min, c_max,filter_func=isQualified)

matrix_simu = make_matrix_simu(
    cleaf_0=cleaf[0],
    croot_0=croot[0],
    cwood_0=cwood[0],
    clitter_0=clitter[0],
    csoil_0=csoil[0],
    npp=npp,
    nyears=140,
    clay=0.2028,
    silt=0.2808,
    lig_wood=0.4
) 
#pa=[beta1,beta2, lig_leaf, f41,f42, kleaf,kroot,kwood,kmet,kmic, kslow,kpass, cmet_init, cstr_init, cmic_init, cpassive_init ]
pa=            [0.15,  0.2,0.15,0.28, 0.6,      1/365,  1/(365*5), 1/(365*40), 0.5/(365*0.1),  0.3/(365*0.137),  0.3/(365*5),  0.3/(222.22*365),          0.05,           0.1,           1,         5]
    
C_upgraded, J_upgraded = mcmc(
        initial_parameters=pa,
        proposer=GenerateParamValues,
        forward_simulation=matrix_simu,
        costfunction=make_weighted_cost_func(obs)
)


df=pd.DataFrame(C_upgraded)
df_j=pd.DataFrame(J_upgraded)
#df.columns = ["std_tomcat","r_tomcat","std_lmdz","r_lmdz","std_jamstec","r_jamstec"]
#df.index = ['all', 'ha', 'ar', 'sa','ds','hu'] 
 
df.to_csv('/home/599/yh4968/cable_demo_da_aa.csv',sep=',')
df_j.to_csv('/home/599/yh4968/cable_demo_da_j_aa.csv',sep=',')
