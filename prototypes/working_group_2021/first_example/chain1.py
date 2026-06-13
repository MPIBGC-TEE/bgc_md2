#!/usr/bin/env python
import os
import sys
import glob
import netCDF4 as nc
import datetime as dt
import numpy as np
import math
import pandas as pd
import calendar
import matplotlib 

from matplotlib import pyplot as plt


# Read NetCDF data  ******************************************************************************************************************************
file= "/g/data/p66/yh4968/cmip6_1pct/access/npp_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
ds = nc.Dataset(file)
var_npp = ds.variables['npp'][:,:,:]
ds.close()

file= "/g/data/p66/yh4968/cmip6_1pct/access/rh_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
ds = nc.Dataset(file)
var_rh = ds.variables['rh'][:,:,:]
ds.close()

file= "/g/data/p66/yh4968/cmip6_1pct/access/cLeaf_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
ds = nc.Dataset(file)
var_cleaf = ds.variables['cLeaf'][:,:,:]
ds.close()


file= "/g/data/p66/yh4968/cmip6_1pct/access/cRoot_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
ds = nc.Dataset(file)
var_croot = ds.variables['cRoot'][:,:,:]
ds.close()


file= "/g/data/p66/yh4968/cmip6_1pct/access/cVeg_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
ds = nc.Dataset(file)
var_cveg = ds.variables['cVeg'][:,:,:]
ds.close()


file= "/g/data/p66/yh4968/cmip6_1pct/access/cSoil_Emon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
ds = nc.Dataset(file)
var_csoil = ds.variables['cSoil'][:,:,:]
lat_data=ds.variables['lat'][:].data
ds.close()

file= "/g/data/p66/yh4968/cmip6_1pct/access/cLitter_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
ds = nc.Dataset(file)
var_clitter = ds.variables['cLitter'][:,:,:]
ds.close()

# pick up 1 site 
npp= var_npp[:,58,159]* 86400   #   kg/m2/s kg/m2/day; 
rh= var_rh[:,58,159]*86400;   # per s to per day  
clitter= var_clitter[:,58,159]; csoil= var_csoil[:,58,159]; 
cveg= var_cveg[:,58,159]; cleaf= var_cleaf[:,58,159]; croot= var_croot[:,58,159]; cwood = cveg - cleaf - croot; 

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
    x_init = np.array([cleaf[0],croot[0],cwood[0],pa[12],pa[13],clitter[0]-pa[12]-pa[13],pa[14], csoil[0]- pa[14] - pa[15],pa[15]]).reshape([9,1])   # Initial carbon pool size
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




def GenerateParamValues_init(c_op):
   flag = True
   while (flag):
      c_new = c_op + (np.random.random((paramNum)) - 0.5)*(c_max - c_min)/10.0
      #c_new = c_op + (np.random.normal(0, 1, paramNum))*(c_max - c_min)/15.0
      if (isQualified(c_new)):
         flag = False
   return c_new


def GenerateParamValues(c_op,covv):
   flag = True
   while (flag):
      c_new = c_op + np.random.multivariate_normal(np.zeros(len(c_op)), covv)
      #c_new = c_op + (np.random.normal(0, 1, paramNum))*(c_max - c_min)/15.0
      if (isQualified(c_new)):
         flag = False
   return c_new





def isQualified(c):
   flag = True
   for i in range(paramNum):
      if(c[i] > c_max[i] or c[i] < c_min[i]):
         flag = False
         break
      if(c[0] + c[1] > 0.99):
         flag = False
         break
   return flag


##==== MCMC =======================================================================================================
#pa=[beta1,beta2, lig_leaf, f41,f42, kleaf,kroot,kwood,kmet,kmic, kslow,kpass, cmet_init, cstr_init, cmic_init, cslow_init ]
#===
c_min=np.array([0.09, 0.09,0.09,0.01,0.01,  1/(2*365), 1/(365*10), 1/(60*365), 0.1/(0.1*365), 0.06/(0.137*365), 0.06/(5*365), 0.06/(222.22*365),clitter[0]/100,clitter[0]/100,csoil[0]/100,csoil[0]/2])
pa=            [0.15,  0.2,0.15,0.28, 0.6,      1/365,  1/(365*5), 1/(365*40), 0.5/(365*0.1),  0.3/(365*0.137),  0.3/(365*5),  0.3/(222.22*365),          0.05,           0.1,           1,         5]
c_max=np.array([1   ,    1,0.21,   1,   1,1/(0.3*365),1/(0.8*365),      1/365,   1/(365*0.1),  0.6/(365*0.137),  0.6/(365*5),  0.6/(222.22*365),    clitter[0],    clitter[0],  csoil[0]/3,csoil[0]])

np.random.seed(seed=10)

paramNum=len(pa)
nsimu    = 20000

upgraded=0;
J_last = 400
C_op = pa

C_upgraded = np.zeros((paramNum, nsimu))
J_upgraded = np.zeros((1, nsimu))


coeff=pd.read_csv("/home/599/yh4968/cable_demo_da_aa.csv")
aa=np.array(coeff)
bb=aa[0:19,1000:9300]
covv=np.cov(bb)




for simu in range(nsimu):
    c_new = GenerateParamValues(C_op,covv)

    x_simu,rh_simu = matrix_simu(c_new)

    J_obj1 = np.mean (( x_simu[:,0] - cleaf[0:tot_len] )**2)/(2*np.var(cleaf[0:tot_len]))
    J_obj2 = np.mean (( x_simu[:,1] - croot[0:tot_len] )**2)/(2*np.var(croot[0:tot_len]))
    J_obj3 = np.mean (( x_simu[:,2] - cwood[0:tot_len] )**2)/(2*np.var(cwood[0:tot_len]))
    J_obj4 = np.mean (( np.sum(x_simu[:,3:6],axis=1)- clitter[0:tot_len] )**2)/(2*np.var(clitter[0:tot_len]))
    J_obj5 = np.mean (( np.sum(x_simu[:,6:9],axis=1)- csoil[0:tot_len] )**2)/(2*np.var(csoil[0:tot_len]))
    
    J_obj6 = np.mean (( rh_simu[:,0] - rh[0:tot_len] )**2)/(2*np.var(rh[0:tot_len]))
    
    J_new= (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 )/200+ J_obj6/4

    delta_J =  J_last - J_new;
    
    randNum = np.random.uniform(0, 1)
    if (min(1.0, np.exp(delta_J)) > randNum):
            C_op=c_new;
            J_last=J_new;
            C_upgraded[:,upgraded]=C_op;
            J_upgraded[:,upgraded]=J_last; 
            upgraded=upgraded+1;

df=pd.DataFrame(C_upgraded)
df_j=pd.DataFrame(J_upgraded)
#df.columns = ["std_tomcat","r_tomcat","std_lmdz","r_lmdz","std_jamstec","r_jamstec"]
#df.index = ['all', 'ha', 'ar', 'sa','ds','hu'] 
 
df.to_csv('/home/599/yh4968/cable_demo_da_chain1.csv',sep=',')
df_j.to_csv('/home/599/yh4968/cable_demo_da_j_chain1.csv',sep=',')


