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
from pathlib import Path
from matplotlib import pyplot as plt

# Read NetCDF data  ******************************************************************************************************************************
#pull in amazon cell; TRENDY is -90:90 N and -180:180 E
#Not sure exactly but ~ -62.5 E and -17 N is lat cell 74, lon cell 118
ex_lat = 74
ex_lon = 118

dataPath = Path('/home/data/jon')
file= dataPath.joinpath("YIBs_S0_Monthly_npp.nc")
ds = nc.Dataset(file)
var_npp = ds.variables['npp'][:,ex_lat,ex_lon] 
ds.close()
#convert kg/m2/s to kg/m2/day
npp = var_npp*86400 


file= dataPath.joinpath("YIBs_S0_Monthly_rh.nc")
ds = nc.Dataset(file)
var_rh = ds.variables['rh'][:,ex_lat,ex_lon]
ds.close()
#convert per s to per day
rh = var_rh*86400 

file= dataPath.joinpath("YIBs_S0_Monthly_ra.nc")
ds = nc.Dataset(file)
var_ra = ds.variables['ra'][:,ex_lat,ex_lon]
ds.close()
#convert per s to per day
ra = var_ra*86400 

file= dataPath.joinpath("YIBs_S0_Annual_cSoil.nc")
ds = nc.Dataset(file)
var_cSoil = ds.variables['cSoil'][:,ex_lat,ex_lon]
ds.close()
csoil = var_cSoil

file= dataPath.joinpath("YIBs_S0_Annual_cVeg.nc")
ds = nc.Dataset(file)
var_cVeg = ds.variables['cVeg'][:,ex_lat,ex_lon]
ds.close()
cveg = var_cVeg

#Still have these variables from TRENDY compared to Yuanyuan CMIP6 example:
#npp, rh, csoil, cveg

#I additionally have the following from TRENDY:
#lai, gpp, ra

#I no longer have the following data streams from TRENDY:
#clitter, cleaf, croot, cwood, cwd

#assign time values
days=[31,28,31,30,31,30,31,31,30,31,30,31]
nyears =320 
tot_len = nyears*12
#assign constants
clay = 0.2028
silt= 0.2808

#define matrix simulation - CASA example
def matrix_simu(pa):
    #define matrix dimensions
    abh = 12  #height of A and B matrix
    aw = 12   #width of A matrix
    bw = 1    #width of B matrix
    
    #assign betas (allocation) from input parameters
    beta_leaf=pa[0]
    beta_root=pa[1]
    beta_wood=1-beta_leaf-beta_root
    
    #create B matrix
    B = np.zeros(abh).reshape([abh,bw])  #empty matrix of zeros
    B[0] = beta_leaf   #fill matrix with beta values
    B[1] = beta_root   #remember python n starts at zero
    B[2] = beta_wood
    
    #define f values
    #leaf(p1),root(p2),wood(p3),cwd(p4),samet(p5),sastr(p6),samic(p7),
    #slmet(p8),slstr(p9),slmic(p10),slow(p11),arm(p12)
    fsamet_leaf = pa[14]             	#optimized leaf -> surface metabolic
    fsastr_leaf = (1-fsamet_leaf-0.2) #remainder to surface str minus 20% loss
    fslmet_root = pa[15]              #optimized root -> soil metabolic
    fslstr_root = (1-fslmet_root-0.2) #remainder to soil str minus 20% loss
    fcwd_wood = 0.4                   #fixed - made up values so far
    fsamic_cwd = pa[16]               #optimized cwd ->  surface microbial
    fslow_cwd = (1-fsamic_cwd-0.2)    #remainder to slow soil minus 20% loss
    fsamic_samet = 0.3                #fixed - made up values so far
    fsamic_sastr = 0.1                #fixed - made up values so far
    fslow_sastr = 0.1                 #fixed - made up values so far
    fslow_samic = 0.1                 #fixed - made up values so far
    fslmic_slmet = 0.4                #fixed - made up values so far
    fslmic_slstr = 0.3                #fixed - made up values so far
    fslow_slstr = 0.2                 #fixed - made up values so far
    fslow_slmic = 0.2                 #fixed - made up values so far
    farm_slmic = 0.1                  #fixed - made up values so far
    fslmic_slow = 0.15                #fixed - made up values so far
    farm_slow = 0.45*(0.003+0.009*clay)	#copied CABLES slow to passive
    fslmic_arm = 0.15                 #fixed - made up values so far
    
    #Create A matrix
    A = np.zeros(abh*aw).reshape([abh,aw])   #create empty A matrix
    np.fill_diagonal(A,-1)  	#replace diagonal with -1s
    #leaf(p1),root(p2),wood(p3),cwd(p4),samet(p5),sastr(p6),samic(p7),
    #slmet(p8),slstr(p9),slmic(p10),slow(p11),arm(p12)
    A[4,0] = fsamet_leaf   	#f5_1 #add in transfers - A is a zero ordered matrix
    A[5,0] = fsastr_leaf   	#f6_1 #minus one from matrix position for [row,col]
    A[7,1] = fslmet_root   	#f8_2 #f8_2, to pool8 from pool2 matrix pos [7,1]
    A[8,1] = fslstr_root   	#f9_2
    A[3,2] = fcwd_wood      #f4_3
    A[6,3] = fsamic_cwd     #f7_4
    A[10,3] = fslow_cwd    	#f11_4
    A[6,4] = fsamic_samet   #f7_5
    A[6,5] = fsamic_sastr   #f7_6
    A[10,5] = fslow_sastr   #f11_6
    A[10,6] = fslow_samic   #f11_7
    A[9,7] = fslmic_slmet   #f10_8
    A[9,8] = fslmic_slstr   #f10_9
    A[10,8] = fslow_slstr   #f11_9
    A[10,9] = fslow_slmic   #f11_10
    A[11,9] = farm_slmic   	#f12_10
    A[9,10] = fslmic_slow   #f10_11
    A[11,10] = farm_slow  	#f12_11
    A[9,11] =fslmic_arm    	#f10_12
    
    #turnover rate per day of pools:
    #leaf(p1),root(p2),wood(p3),cwd(p4),samet(p5),sastr(p6),samic(p7),
    #slmet(p8),slstr(p9),slmic(p10),slow(p11),arm(p12) 
    k_val = np.array(
        [
            pa[2],
            pa[3],
            pa[4], 
            pa[5],
            pa[6], 
            pa[7], 
            pa[8],
            pa[9], 
            pa[10], 
            pa[11], 
            pa[12], 
            pa[13]
        ]
    )
    K = np.zeros(abh*aw).reshape([12,12])
    np.fill_diagonal(K, k_val)
      
    x_fin=np.zeros((tot_len,12))
    rh_fin=np.zeros((tot_len,1))
    ra_fin=np.zeros((tot_len,1))
    cveg_m = np.zeros((tot_len,1))
    csoil_m = np.zeros((tot_len,1))
    #leaf(p1),root(p2),wood(p3),cwd(p4),samet(p5),sastr(p6),samic(p7),
    #slmet(p8),slstr(p9),slmic(p10),slow(p11),arm(p12)
    #x_init = np.array([cleaf[0],croot[0],cwood[0],cwd[0],pa[17],pa[18],
    #clitter[0]-pa[17]-pa[18],pa[19],pa[20],pa[21],pa[22],
    #csoil[0]-pa[19]-pa[20]-pa[21]-pa[22]]).reshape([12,1])
    x_init = np.array(
        [
            pa[17],
            pa[18],
            cveg[1]-pa[17]-pa[18],
            pa[19],
            pa[20],
            pa[21],
            pa[22],
            pa[23],
            pa[24],
            pa[25],
            pa[26],
            csoil[1]-pa[23]-pa[24]-pa[25]-pa[26]
        ]
    ).reshape([12,1])
    X=x_init   # initialize carbon pools 
    jj=0
    for y in np.arange(0,nyears):
        for m in np.arange(0,12):
            npp_in = npp[jj]
            co2_rh = 0
            co2_ra = 0
            for d in np.arange(0,days[m]):
                X=X + B*npp_in + np.array(A@K@X).reshape([12,1])
                #leaf(p1),root(p2),wood(p3),cwd(p4),samet(p5),sastr(p6),samic(p7),slmet(p8),slstr(p9),slmic(p10),slow(p11),arm(p12)
                co2_hrate = [
                    0,
                    0,
                    0, 
                    (1-fsamic_cwd-fslow_cwd)*K[3,3], 
                    (1-fsamic_samet)*K[4,4], 
                    (1-fsamic_sastr-fslow_sastr)*K[5,5], 
                    (1-fslow_samic)*K[6,6], 
                    (1-fslmic_slmet)*K[7,7], 
                    (1-fslmic_slstr-fslow_slstr)*K[8,8], 
                    (1-fslow_slmic-farm_slmic)*K[9,9], 
                    (1-fslmic_slow-farm_slow)*K[10,10], 
                    (1-fslmic_arm)*K[11,11]
                ]
                co2h=np.sum(co2_hrate*X.reshape(1,12))
                co2_rh = co2_rh + co2h/days[m]   # monthly average rh 
                #similarly calculate autotrophic respiration ra
                co2_arate = [
                    (1-fsamet_leaf-fsastr_leaf)*K[0,0], 
                    (1-fslmet_root-fslstr_root)*K[1,1],
                    (1-fcwd_wood)*K[2,2],
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0
                ]
                co2a=np.sum(co2_arate*X.reshape(1,12))
                co2_ra = co2_ra + co2a/days[m]
            x_fin[jj,:]=X.reshape(1,12)
            rh_fin[jj,0]=co2_rh
            ra_fin[jj,0]=co2_ra
            #one option make the annual cveg/csoil into repeated monthly values
            cveg_m[jj,0] =cveg[y]
            csoil_m[jj,0]=csoil[y]
            jj= jj+1
        
    return x_fin,rh_fin,ra_fin,cveg_m,csoil_m     

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

c_min=np.array(
    [
        0.09,               # 1  - beta1
        0.09,               # 2  - beta2
        1/(2*365),          # 3  - kleaf
        1/(365*10),         # 4  - kroot
        1/(60*365),         # 5  - kwood
        1/(365*10),         # 6  - kcwd
        0.1/(0.1*365),      # 7  - ksamet 
        0.1/(0.1*365),      # 8  - ksastr 
        0.06/(0.137*365),   # 9  - ksamic
        0.06/(5*365),       # 10 - kslmet
        0.06/(5*365),       # 11 - kslstr
        0.06/(5*365),       # 12 - kslmic
        0.06/(5*365),       # 13 - kslow
        0.06/(222.22*365),  # 14 - karm
        0,                  # 15 - f5_1
        0,                  # 16 - f9_2 
        0,                  # 17 - f11_4
        cveg[0]/100,        # 18 - leaf_i
        cveg[0]/100,        # 19 - root_i
        cveg[0]/100,        # 20 - cwd_i
        cveg[0]/10000,      # 21 - samet_i
        cveg[0]/10000,      # 22 - sastr_i
        cveg[0]/10000,      # 23 - samic_i
        csoil[0]/100,       # 24 - slmet_i
        csoil[0]/100,       # 25 - slstr_i
        csoil[0]/100,       # 26 - slmic_i
        csoil[0]/100        # 27 - slow_i
    ]
)

pa=[
    0.15,                   # 1  - beta1
    0.2,                    # 2  - beta2
    1/365,                  # 3  - kleaf        
    1/(365*5),              # 4  - kroot
    1/(365*40),             # 5  - kwood
    1/(365*5),              # 6  - kcwd 
    0.5/(365*0.1),          # 7  - ksamet
    0.5/(365*0.1),          # 8  - ksastr
    0.3/(365*0.137),        # 9  - ksamic 
    0.3/(365),              # 10 - kslmet
    0.3/(365),              # 11 - kslstr 
    0.3/(365),              # 12 - kslmic
    0.3/(365*5),            # 13 - kslow
    0.3/(222.22*365),       # 14 - karm
    0.3,                    # 15 - f5_1
    0.3,                    # 16 - f9_2
    0.3,                    # 17 - f11_4
    cveg[0]/3,              # 18 - leaf_i
    cveg[0]/3,              # 19 - root_i
    cveg[0]/50,             # 20 - cwd_i
    cveg[0]/300,            # 21 - samet_i
    cveg[0]/300,            # 22 - sastr_i
    cveg[0]/500,            # 23 - samic_i
    csoil[0]/10,            # 24 - slmet_i
    csoil[0]/10,            # 25 - slstr_i
    csoil[0]/10,            # 26 - slmic_i
    csoil[0]/10             # 27 - slow_i 
]

c_max=np.array(
    [
        1,                  # 1  - beta1
        1,                  # 2  - beta2
        1/(0.3*365),        # 3  - kleaf
        1/(0.8*365),        # 4  - kroot
        1/365,              # 5  - kwood
        1/(0.8*365),        # 6  - kcwd 
        1/(365*0.1),        # 7  - ksamet  
        1/(365*0.1),        # 8  - ksastr
        0.6/(365*0.137),    # 9  - ksamic
        0.6/(365),          # 10 - kslmet
        0.6/(365),          # 11 - kslstr
        0.6/(365),          # 12 - kslmic
        0.6/(365*5),        # 13 - kslow
        0.6/(222.22*365),   # 14 - karm
        0.8,                # 15 - f5_1
        0.8,                # 16 - f9_2
        0.8,                # 17 - f11_4 
        cveg[0]/2,          # 18 - leaf_i
        cveg[0]/2,          # 19 - root_i
        cveg[0]/10,         # 20 - cwd_i
        cveg[0]/100,        # 21 - samet_i
        cveg[0]/100,        # 22 - sastr_i
        cveg[0]/100,        # 23 - samic_i
        csoil[0]/2,         # 24 - slmet_i
        csoil[0]/2,         # 25 - slstr_i
        csoil[0]/2,         # 26 - slmic_i
        csoil[0]/2          # 27 - slow_i
    ]
)

np.random.seed(seed=10)

paramNum=len(pa)
nsimu   =10

###############################
#initial covariance run
###############################

upgraded=0;
J_last = 4000000
C_op = pa

C_upgraded = np.zeros((paramNum, nsimu))
J_upgraded = np.zeros((1, nsimu))

for simu in range(nsimu):
    c_new = GenerateParamValues_init(C_op)
    
    x_simu,rh_simu,ra_simu,cveg_m,csoil_m = matrix_simu(c_new)

    J_obj4 = np.mean((np.sum(x_simu[:,0:2],axis=1)-cveg_m[0:tot_len])**2)/(2*np.var(cveg_m[0:tot_len]))
    J_obj5 = np.mean((np.sum(x_simu[:,7:11],axis=1)-csoil_m[0:tot_len])**2)/(2*np.var(csoil_m[0:tot_len]))
    J_obj6 = np.mean((rh_simu[:,0]-rh[0:tot_len])**2)/(2*np.var(rh[0:tot_len]))
    J_obj7 = np.mean((ra_simu[:,0]-ra[0:tot_len])**2)/(2*np.var(ra[0:tot_len]))
    #print("j values: ")
    #print(J_obj4)
    #print(J_obj5)
    #print(J_obj6)
    #print(J_obj7)
    J_new= J_obj4 + J_obj5 + J_obj6 + J_obj7

    #print("J_new: " + str(J_new))
    delta_J =  J_last - J_new
    #print("delta_J: " + str(delta_J))
    
    randNum = np.random.uniform(0, 1)
    if (min(1.0, np.exp(delta_J)) > randNum):
        C_op=c_new
        J_last=J_new
        C_upgraded[:,upgraded]=C_op
        J_upgraded[:,upgraded]=J_last 
        upgraded=upgraded+1

#output parameter and J values
df=pd.DataFrame(C_upgraded)
df_j=pd.DataFrame(J_upgraded)
df.to_csv('./YIBs_covar.csv',sep=',')
df_j.to_csv('./YIBs_j_covar.csv',sep=',')

###########################################
#final da run with covariance matrix
###########################################

#how many initial MCMC iterations do we toss from covar run (first 10%)
burn_len = nsimu*0.1

upgraded=0
J_last = 400
C_op = pa

coeff=df
aa=np.array(coeff)
bb=aa[0:paramNum,int(burn_len):(np.count_nonzero(aa[0]))]
covv=np.cov(bb)

C_upgraded = np.zeros((paramNum, nsimu))
J_upgraded = np.zeros((1, nsimu))

for simu in range(nsimu):
    c_new = GenerateParamValues(C_op,covv)
    
    x_simu,rh_simu,ra_simu,cveg_m,csoil_m = matrix_simu(c_new)

    J_obj4 = np.mean (( np.sum(x_simu[:,0:2],axis=1)- cveg_m[0:tot_len] )**2)/(2*np.var(cveg_m[0:tot_len]))
    J_obj5 = np.mean (( np.sum(x_simu[:,7:11],axis=1)- csoil_m[0:tot_len] )**2)/(2*np.var(csoil_m[0:tot_len]))
    J_obj6 = np.mean (( rh_simu[:,0] - rh[0:tot_len] )**2)/(2*np.var(rh[0:tot_len]))
    J_obj7 = np.mean (( ra_simu[:,0] - ra[0:tot_len] )**2)/(2*np.var(ra[0:tot_len]))
    
    J_new= (J_obj4 + J_obj5)/1000 + J_obj6 + J_obj7
    
    delta_J =  J_last - J_new
    
    randNum = np.random.uniform(0, 1)
    if (min(1.0, np.exp(delta_J)) > randNum):
        C_op=c_new
        J_last=J_new
        C_upgraded[:,upgraded]=C_op
        J_upgraded[:,upgraded]=J_last 
        upgraded=upgraded+1

#output parameter and J values
df=pd.DataFrame(C_upgraded)
df_j=pd.DataFrame(J_upgraded)
df.to_csv('./YIBs_chain1.csv',sep=',')
df_j.to_csv('./YIBs_j_chain1.csv',sep=',')

#remove zero values from C_ and J_upgraded
df = df.loc[:, (df != 0).any(axis=0)]
df_j = df_j.loc[:, (df_j != 0).any(axis=0)]

#select best parameter values
bestpar_index = int(df_j.idxmin(axis=1)[0])
print(bestpar_index)
bestpar = df.iloc[:,bestpar_index]
print(bestpar)

#output observation data to graph
x_simu,rh_simu,ra_simu,cveg_m,csoil_m = matrix_simu(bestpar)

#output simulation data to graph
df=pd.DataFrame(x_simu)
df.to_csv('./simu_best_x.csv', sep=',')
df=pd.DataFrame(rh_simu)
df.to_csv('./simu_best_rh.csv', sep=',')
df=pd.DataFrame(ra_simu)
df.to_csv('./simu_best_ra.csv', sep=',')

#numpy plotting example
xs = np.linspace(0, 2 * np.pi, 100)
ys1 = np.sin(xs)
ys2 = np.sin(xs + np.pi / 16)

fig, ax = plt.subplots(figsize=(9.2, 5))
ax.set_title("blue chasing red")
ax.plot(xs,ys1,color='red', label='sin(x)')
ax.plot(xs,ys2,color='blue', label='sin(x+pi/16)')
ax.legend()
fig.savefig('myfigure.pdf')
