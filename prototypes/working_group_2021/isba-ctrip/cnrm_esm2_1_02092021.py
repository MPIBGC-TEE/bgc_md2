#coding for CNRM_ESM2_1 CIMP6 output with Matrix Equation
"""
Created on Thu Sep  2 20:36:42 2021

@author: Qing-Fang Bi
"""

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
file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/npp_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_npp = ds.variables['npp'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/rh_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_rh = ds.variables['rh'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/fLitterSoil_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_flitterSoil = ds.variables['fLitterSoil'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/fVegLitter_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_fveglitter = ds.variables['fVegLitter'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cLeaf_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_cleaf = ds.variables['cLeaf'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cStem_Emon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_cstem = ds.variables['cStem'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cWood_Emon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_cwood = ds.variables['cWood'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cLitterSurf_Emon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_clitterabove = ds.variables['cLitterSurf'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cLitterBelow_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_clitterbelow = ds.variables['cLitterBelow'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cRoot_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_croot = ds.variables['cRoot'][:,:,:]
ds.close()


file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cVeg_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_cveg = ds.variables['cVeg'][:,:,:]
ds.close()


file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cSoil_Emon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_csoil = ds.variables['cSoil'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cSoilFast_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_csoilfast = ds.variables['cSoilFast'][:,:,:]
ds.close()

file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cSoilMedium_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_csoilmedium = ds.variables['cSoilMedium'][:,:,:]
ds.close()


file= "/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/cSoilSlow_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"
ds = nc.Dataset(file)
var_csoilslow = ds.variables['cSoilSlow'][:,:,:]
lat_data=ds.variables['lat'][:].data
lon_data=ds.variables['lon'][:].data
ds.close()


# pick up 1 site   62.8125 W, 17.5S
npp= var_npp[:,99,24]* 86400   #   kg/m2/s to kg/m2/day; 
rh= var_rh[:,99,24]*86400;   # per s to per day  
cstem= var_cstem[:,99,24]; csoil= var_csoil[:,99,24]; 
cveg= var_cveg[:,99,24]; cleaf= var_cleaf[:,99,24]; croot= var_croot[:,99,24]; cwood= var_cwood[:,99,24];
clitterabove= var_clitterabove[:,99,24];clitterbelow= var_clitterbelow[:,99,24];
csoilfast= var_csoilfast[:,99,24]; csoilmedium= var_csoilmedium[:,99,24]; csoilslow= var_csoilslow[:,99,24]; 

clay=0.4
silt=0.2 #??

days=[31,28,31,30,31,30,31,31,30,31,30,31]
nyears =140 
tot_len = nyears*12


def matrix_simu(pa):
    # Construct B matrix 
    beta1=pa[0]
    B = np.array([beta1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).reshape([13,1])   # allocation
    # Now construct A matrix
    f21 = pa[1]; f32 = pa[2];f41 = pa[3];f42 = pa[4];f43 = pa[5];f53 = pa[6];f64 = pa[7];
    lig_aboveleaf= pa[8]; #lignin fraction of above-ground structural litter
    lig_belowleaf = pa[9];#lignin fraction of below-ground structural litter
    f71 = 0.3084; f72 = 0.402; f73 = 0.402; f75 = 0.402; f81 = 0.6916; f82 = 0.598; f83 = 0.598; f85 = 0.598; 
    f94 = 0.402;  f96 = 0.402; f104 = 0.598; f106 = 0.598; f117 = 0.55*(1-lig_aboveleaf); f118 = 0.45;
    f119 = 0.45*(1-lig_belowleaf); f1110 = 0.45; f1112 = 0.42; f1113 = 0.45;  f127 = 0.7 * lig_aboveleaf; f129 = 0.7 * lig_belowleaf; f1211 = 1 - 0.004 - (0.85- 0.68 * (clay+silt));
    f1311 = 0.004; f1312 = 0.03
    
    A = np.array([ -1,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                   f21,-1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                   0,  f32, -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                   f41,f42, f43, -1,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
                   0,   0,  f53,  0,  -1,   0,   0,   0,   0,   0,   0,   0,   0, 
                   0,   0,   0,  f64,  0,  -1,   0,   0,   0,   0,   0,   0,   0, 
                   f71, f72, f73, 0,  f75,  0,  -1,   0,   0,   0,   0,   0,   0,
                   f81, f82, f83, 0,  f85,  0,   0,  -1,   0,   0,   0,   0,   0,
                   0,   0,   0,  f94,  0, f96,   0,   0,  -1,   0,   0,   0,   0,
                   0,   0,   0,  f104, 0,f106,   0,   0,   0,   -1,  0,   0,   0,
                   0,   0,   0,   0,   0,   0, f117, f118,f119,f1110,-1,f1112,f1113,
                   0,   0,   0,   0,   0,   0, f127,  0,  f129, 0, f1211, -1,  0,
                   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, f1311,f1312,-1]).reshape([13,13])   # tranfer
    
    #turnover rate per day of pools: 
    temp = [pa[10],pa[11],pa[12], pa[13],pa[14], pa[15], 1/(0.245*np.exp(-3*pa[8])),1/0.066, 1/(0.245*np.exp(-3*pa[9])), 1/0.066, 1/0.149 * (1 - 0.75 * (clay+silt)), 1/5.37, 1/241]
    K = np.zeros(169).reshape([13, 13]) 
    for i in range(0, 13):
        K[i][i] = temp[i]
      
    x_fin=np.zeros((tot_len,13))
    rh_fin=np.zeros((tot_len,1))
    # leaf, active structural biomass,passive structural biomass, below ground structural biomass, above ground woody biomass, below ground woody biomass,Aboveground structure litter,Aboveground metabolic litter,Belowground structure litte,Belowground metabolic litter, active soil, slow soil, passive soil
    x_init = np.array([cleaf[0],pa[16],pa[17],pa[18],cwood[0]+cstem[0],croot[0]-pa[18],pa[19],clitterabove[0]-pa[19],pa[20],clitterbelow[0]-pa[20],csoilfast[0],csoilmedium[0],csoilslow[0]]).reshape([13,1])   # Initial carbon pool size
    X=x_init   # initialize carbon pools 
    jj=0
    for y in np.arange(0,nyears):
       for m in np.arange(0,12):
         npp_in = npp[jj] 
         co2_rh = 0  
         for d in np.arange(0,days[m]):
             X=X + B*npp_in + np.array(A@K@X).reshape([13,1])
             co2_rate = [0,0,0,0,0,0,(1-f117-f127)*K[6,6],(1-f118)*K[7,7],(1-f119-f129)*K[8,8],(1-f1110)*K[9,9],(1-f1211-f1311)*K[10,10],(1-f1112-f1312)*K[11,11],(1-f1113)*K[12,12]]
             co2=np.sum(co2_rate*X.reshape(1,13))
             co2_rh = co2_rh + co2/days[m]   # monthly average rh 
         x_fin[jj,:]=X.reshape(1,13)
         rh_fin[jj,0]=co2_rh
         jj= jj+1
         
    return x_fin, rh_fin     

 
 
def GenerateParamValues(c_op):
   flag = True
   while (flag):
      c_new = c_op + (np.random.random((paramNum)) - 0.5)*(c_max - c_min)/10.0
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
#pa=[beta1,beta2, lig_leaf, f41,f42, kleaf,kroot,kwood,kmet,kmic, kslow,kpass, cmet_init, cstr_init, cmic_init, cpassive_init ]
#===
    
c_min=np.array([0.09,  0.09, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2,   0.01, 0.01,    1/(2*365), 1/(365*100), 1/(365*100), 0.1/(365*0.5),  0.03/(365*3),  0.03/(365*2),  cleaf[0]/100, 0.0001, croot[0]/100,    clitterabove[0]/100,   clitterbelow[0]/100])

pa=            [0.4,  0.2, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2,   0.6916, 0.3084,    1/365, 1/(365*10), 1/(365*5), 0.5/(365*0.5),  0.3/(365*3),  0.3/(365*2),  0.003, 0.001, 0.005,    0.01,   0.01]

c_max=np.array([1,   1, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2,    0.9, 0.6,     1/(0.3*365), 1/(365*0.8), 1/(365*0.8), 1/(365*0.5),  1/(365*3),  0.6/(365*2),  cleaf[0], 0.1, croot[0],    clitterabove[0],   clitterbelow[0]])

np.random.seed(seed=10)

paramNum=len(pa)
nsimu    = 1000

upgraded=0;
J_last = 400
C_op = pa

C_upgraded = np.zeros((paramNum, nsimu))
J_upgraded = np.zeros((1, nsimu))

for simu in range(nsimu):
    c_new = GenerateParamValues(C_op)
    x_simu, rh_simu = matrix_simu(c_new)
    J_obj1 = np.mean (( x_simu[:,0] - cleaf[0:tot_len] )**2)/(2*np.var(cleaf[0:tot_len]))
    J_obj2 = np.mean (( np.sum(x_simu[:,6:7],axis=1)- clitterabove[0:tot_len] )**2)/(2*np.var(clitterabove[0:tot_len]))
    J_obj3 = np.mean (( np.sum(x_simu[:,8:9],axis=1)- clitterbelow[0:tot_len] )**2)/(2*np.var(clitterbelow[0:tot_len]))
    J_obj4 = np.mean (( x_simu[:,10] - csoilfast[0:tot_len] )**2)/(2*np.var(csoilfast[0:tot_len]))
    J_obj5 = np.mean (( x_simu[:,11] - csoilmedium[0:tot_len] )**2)/(2*np.var(csoilmedium[0:tot_len]))
    J_obj6 = np.mean (( x_simu[:,12] - csoilslow[0:tot_len] )**2)/(2*np.var(csoilslow[0:tot_len])) 
    
    J_obj7 = np.mean (( rh_simu[:,0] - rh[0:tot_len] )**2)/(2*np.var(rh[0:tot_len]))
    
    J_new= (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 + J_obj6 )/200+ J_obj7/4
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

"/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc/npp_Lmon_CNRM-ESM2-1_1pctCO2-bgc_r1i1p1f2_gr_185001-198912.nc"

df.to_csv('/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/CNRM_ESM2_1_da_aa.csv',sep=',')
df_j.to_csv('/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/CNRM_ESM2_1_da_j_aa.csv',sep=',')
df.to_csv("/C/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/CNRM_ESM2_1_da_aa.csv',sep=','")




