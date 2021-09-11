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
    make_weighted_cost_func
)

from general_helpers import make_uniform_proposer

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

days=[31,28,31,30,31,30,31,31,30,31,30,31]
nyears =140 
tot_len = nyears*12
clay = 0.2028
silt= 0.2808
lig_wood=0.4

def matrix_simu(pa):
    # Construct B matrix 
    beta1=pa[0]; beta2=pa[1]; beta3= 1- beta1- beta2
    B = np.array([beta1, beta2, beta3, 0, 0, 0, 0,0,0]).reshape([9,1])   # allocation
    # Now construct A matrix
    lig_leaf = pa[2]

    f41 = pa[3]; f42 = pa[4]; f51 = 1-f41; f52 = 1-f42; f63 = 1;
    f74 = 0.45; f75 = 0.45*(1-lig_leaf); 
    f85 = 0.7*lig_leaf; f86 = 0.4*(1-lig_wood);
    f96 = 0.7*lig_wood;  
    f87=(0.85 - 0.68 * (clay+silt))* (0.997 - 0.032*clay)
    f97=(0.85 - 0.68 * (clay+silt))* (0.003 + 0.032*clay)
    f98=0.45 * (0.003 + 0.009 *clay)

    A = np.array([  -1,   0,   0,   0,   0,   0,   0,   0,   0,
                     0,  -1,   0,   0,   0,   0,   0,   0,   0,
                     0,   0,  -1,   0,   0,   0,   0,   0,   0,
                   f41, f42,   0,  -1,   0,   0,   0,   0,   0,
                   f51, f52,   0,   0,  -1,   0,   0,   0,   0,
                     0,   0, f63,   0,   0,  -1,   0,   0,   0,
                     0,   0,   0, f74, f75,   0,  -1,   0,   0,
                     0,   0,   0,   0, f85, f86, f87,  -1,   0,
                     0,   0,   0,   0,   0, f96, f97, f98,  -1 ]).reshape([9,9])   # tranfer

    #turnover rate per day of pools: 
    temp = [pa[5],pa[6],pa[7], pa[8],pa[8]/(5.75*np.exp(-3*pa[2])), pa[8]/20.6, pa[9],pa[10], pa[11]]
    K = np.zeros(81).reshape([9, 9])
    for i in range(0, 9):
        K[i][i] = temp[i]
      
    x_fin=np.zeros((tot_len,9))
    rh_fin=np.zeros((tot_len,1))
    # leaf, root , wood, metabolic, structural, CWD, microbial, slow, passive 
    x_init = np.array([cleaf[0],croot[0],cwood[0],pa[12],pa[13],clitter[0]-pa[12]-pa[13],pa[14],csoil[0]- pa[14] - pa[15], pa[15]]).reshape([9,1])   # Initial carbon pool size
    X=x_init   # initialize carbon pools 
    jj=0
    for y in np.arange(0,nyears):
       for m in np.arange(0,12):
         npp_in = npp[jj] 
         co2_rh = 0  
         for d in np.arange(0,days[m]):
             X=X + B*npp_in + np.array(A@K@X).reshape([9,1])
             co2_rate = [0,0,0, (1-f74)*K[3,3],(1-f75-f85)*K[4,4],(1-f86-f96)*K[5,5], (1- f87 - f97)*K[6,6], (1-f98)*K[7,7], K[8,8] ]
             co2=np.sum(co2_rate*X.reshape(1,9))
             co2_rh = co2_rh + co2/days[m]   # monthly average rh 
         x_fin[jj,:]=X.reshape(1,9)
         rh_fin[jj,0]=co2_rh
         jj= jj+1
         
    return x_fin, rh_fin     


c_min=np.array([0.09, 0.09,0.09,0.01,0.01,  1/(2*365), 1/(365*10), 1/(60*365), 0.1/(0.1*365), 0.06/(0.137*365), 0.06/(5*365), 0.06/(222.22*365),clitter[0]/100,clitter[0]/100,csoil[0]/100,csoil[0]/2])
c_max=np.array([1   ,    1,0.21,   1,   1,1/(0.3*365),1/(0.8*365),      1/365,   1/(365*0.1),  0.6/(365*0.137),  0.6/(365*5),  0.6/(222.22*365),    clitter[0],    clitter[0],  csoil[0]/3,csoil[0]])

# fixme
#   this function is model specific: It discards parameter proposals
#   where beta1 and beta2 
isQualified = make_param_filter_func(c_max,c_min)
GenerateParamValues = make_uniform_proposer(c_min, c_max,filter_func=isQualified)
def mcmc(
        forward_simulation: Callable[[Tuple], Tuple[np.ndarray, np.ndarray]],
        costfunction: Callable[[np.ndarray],np.float64]
    ) -> Tuple:
    ##==== MCMC =======================================================================================================
    #pa=[beta1,beta2, lig_leaf, f41,f42, kleaf,kroot,kwood,kmet,kmic, kslow,kpass, cmet_init, cstr_init, cmic_init, cpassive_init ]
    #===
    
    pa=            [0.15,  0.2,0.15,0.28, 0.6,      1/365,  1/(365*5), 1/(365*40), 0.5/(365*0.1),  0.3/(365*0.137),  0.3/(365*5),  0.3/(222.22*365),          0.05,           0.1,           1,         5]
    
    np.random.seed(seed=10)
    
    paramNum=len(pa)
    nsimu    = 20000
    
    upgraded=0;
    J_last = 400
    C_op = pa
    
    C_upgraded = np.zeros((paramNum, nsimu))
    J_upgraded = np.zeros((1, nsimu))
    
    for simu in range(nsimu):
        c_new = GenerateParamValues(C_op)
    
        x_simu,rh_simu = forward_simulation(c_new)
        out_simu = np.concatenate([x_simu,rh_simu],axis=1)
        J_new = costfunction(out_simu)
    
        delta_J =  J_last - J_new;
        
        randNum = np.random.uniform(0, 1)
        if (min(1.0, np.exp(delta_J)) > randNum):
                C_op=c_new;
                J_last=J_new;
                C_upgraded[:,upgraded]=C_op;
                J_upgraded[:,upgraded]=J_last; 
                upgraded=upgraded+1;
    
    return C_upgraded, J_upgraded

# main
C_upgraded, J_upgraded = mcmc(
        forward_simulation=matrix_simu,
        costfunction=make_weighted_cost_func(obs)
)


df=pd.DataFrame(C_upgraded)
df_j=pd.DataFrame(J_upgraded)
#df.columns = ["std_tomcat","r_tomcat","std_lmdz","r_lmdz","std_jamstec","r_jamstec"]
#df.index = ['all', 'ha', 'ar', 'sa','ds','hu'] 
 
df.to_csv('/home/599/yh4968/cable_demo_da_aa.csv',sep=',')
df_j.to_csv('/home/599/yh4968/cable_demo_da_j_aa.csv',sep=',')
