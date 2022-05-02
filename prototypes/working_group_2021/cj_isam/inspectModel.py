# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
# load HTML to adjust jupyter settings
from IPython.display import HTML
# #%load_ext autoreload
# #%autoreload 2
# adjust jupyter display to full screen width
display(HTML("<style>.container { width:100% !important; }</style>"))


from source import mvs

# -

mvs.get_CompartmentalMatrix()

# we can also print the whole mass balance equation
import bgc_md2.display_helpers as dh
dh.mass_balance_equation(mvs)

# we can also plot a picture
import bgc_md2.helper as h
h.compartmental_graph(mvs)

# ### Connect to the form ispecific to the tracability analysis.
#
# We are aiming for a decomposition of the compartmental matrix $B$ into three factors 
# $B(t)= \xi(t)  A K $ 
# with $ \xi$ and $K$ diagonal. 
# `bgc_md2` can automatically compute a decomposition into $ B=A N(t)$ where $N(t)=\xi(t)K$ but
# which part of $N(t)$ should go into $\xi(t)$ is easier to specify manually. 
#
# We will first compute the $B=A N$ decomposition and then specify $\xi$.
#
#

# Copy the content of the above cell into a file `model_specific_helpers_2.py` 
# and then import it and call the function to check that it works.

from pathlib import Path
import json
import model_specific_helpers_2 as msh
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 
msh.download_my_TRENDY_output(conf_dict)

import model_specific_helpers_2 as msh
#svs,dvs=msh.get_example_site_vars(dataPath=Path(conf_dict["dataPath"])) #also works
svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
svs_0 = msh.Observables(*map(lambda v: v[0],svs))
dvs_0 = msh.Drivers(*map(lambda v: v[0],dvs))


dvs.npp.shape

from sympy import Symbol
par_dict={
    Symbol(k):v for k,v in 
    {
        "beta_NWT": 0.3,
        "beta_AGWT": 0.3,
        "beta_TR": 0.10000000000000003,
        "beta_GVF": 0.15,
        "beta_GVR": 0.15,
        "r_C_NWT_rh": 0,
        "r_C_AGWT_rh": 0,
        "r_C_TR_rh": 0,
        "r_C_GVF_rh": 0,
        "r_C_GVR_rh": 0,
        "r_C_AGML_rh": 0.00678082191780822,
        "r_C_AGSL_rh": 0.0354794520547945,
        "r_C_AGMS_rh": 0.00800000000000000,
        "r_C_YHMS_rh": 0.00246575342465753,
        "r_C_BGDL_rh": 0.0200000000000000,
        "r_C_BGRL_rh": 0.000600000000000000,
        "r_C_BGMS_rh": 0.00132000000000000,
        "r_C_SHMS_rh": 4.00000000000000e-5,
        "r_C_NWT_2_C_AGML": 0.00116438356164384,
        "r_C_NWT_2_C_AGSL": 0.000205479452054795,
        "r_C_AGWT_2_C_AGSL": 9.13242009132420e-5,
        "r_C_TR_2_C_BGDL": 9.21544209215442e-5,
        "r_C_TR_2_C_BGRL": 3.23785803237858e-5,
        "r_C_GVF_2_C_AGML": 7.76255707762557e-5,
        "r_C_GVF_2_C_AGSL": 1.36986301369863e-5,
        "r_C_GVR_2_C_BGDL": 9.21544209215442e-5,
        "r_C_GVR_2_C_BGRL": 3.23785803237858e-5,
        "r_C_AGML_2_C_AGMS": 0.00554794520547945,
        "r_C_AGSL_2_C_AGMS": 0.00760273972602740,
        "r_C_AGSL_2_C_YHMS": 0.00760273972602740,
        "r_C_AGMS_2_C_YHMS": 0.0120000000000000,
        "r_C_YHMS_2_C_AGMS": 0.00271232876712329,
        "r_C_YHMS_2_C_SHMS": 0.000301369863013699,
        "r_C_BGDL_2_C_SHMS": 0.00739726027397260,
        "r_C_BGRL_2_C_BGMS": 0.000110958904109589,
        "r_C_BGRL_2_C_SHMS": 0.000110958904109589,
        "r_C_BGMS_2_C_SHMS": 0.000488219178082192,
        "r_C_SHMS_2_C_BGMS": 1.47945205479452e-5 
    }.items()    
}

# To be able to run the model forward we not only have to replace parameter symbols by values but symbolic functions by normal python functions.
# In our case the functions for $NPP$ and $\xi$ have to be provided. NPP_fun will interpolate the NPP for the day in question from the data. Which we have to load. 
# We will later store these functions in  `model_specific_helpers.py` which resides in the same folder as this notebook. You will have to adapt them to your data set. 

# +
import numpy as np 
#svs_0=msh.Observables(*map(lambda v: v[0],svs))

X_0= np.array((
    svs_0.cVeg/5,
    svs_0.cVeg/5,
    svs_0.cVeg/5,
    svs_0.cVeg/5,
    svs_0.cVeg/5,
    svs_0.cLitter/4,
    svs_0.cLitter/4,
    svs_0.cLitter/4,
    svs_0.cLitter/4,
    svs_0.cSoil/4,
    svs_0.cSoil/4,
    svs_0.cSoil/4,
    svs_0.cSoil/4,
))#.reshape(9,)
# -

# ## Data assimilation
# Until now we have used only the initial values of the observations. 
# The next step is to decide which parameters we want to consider fixed and which to be estimated.
# This distinction helps, to keep the to create a function which only takes the estimated parameters and thus can be used by a generalized mcmc as will become clear.
#
# We can change which parameters we fix and which we estimate later or can have several approaches for the same symbolic model.
# The distinction is not model inherent but just a reflection of our choice for data assimilation.
# The more parameter values we can find out from the literature the fewer values we have to estimate.  

cpa=msh.Constants(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 npp_0=dvs.npp[0] * 86400,   # kg/m2/s kg/m2/day
 rh_0=svs_0.rh * 86400,   # kg/m2/s kg/m2/day
 ra_0=svs_0.ra * 86400,   # kg/m2/s kg/m2/day
 r_C_NWT_rh=0,
 r_C_AGWT_rh=0,
 r_C_TR_rh=0,
 r_C_GVF_rh=0,
 r_C_GVR_rh=0,
 r_C_AGML_rh=0.55*4.5/365,
 r_C_AGSL_rh=0.7*18.5/365,
 r_C_AGMS_rh=0.4*7.3/365,
 r_C_YHMS_rh=0.45*2.0/365,
 k_C_BGDL=10/365,
 k_C_BGRL=0.3/365,
 k_C_BGMS=0.66/365,
 k_C_SHMS=0.02/365,
 r_C_AGML_2_C_AGMS=0.45*4.5/365,
 r_C_AGMS_2_C_YHMS=0.6*7.3/365,
 r_C_YHMS_2_C_AGMS=0.9*0.55*2.0/365,
 r_C_YHMS_2_C_SHMS=0.1*0.55*2.0/365,
 #number_of_months=len(svs.rh)
 number_of_months=3840 # for testing and tuning mcmc
)

# ### Finding better start values for the data assimilation
# You don't have to do this. It's a heuristic approach to find a better starting position.

# +
import sys
sys.path.insert(0,'..') # necessary to import general_helpers
import general_helpers as gh
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2

def make_steady_state_iterator_sym(
        mvs,
        V_init,
        par_dict,
        func_dict
    ):
    B_func, u_func = gh.make_B_u_funcs_2(mvs,par_dict,func_dict)  
    def f(it,tup):
        X,_,_=tup
        b = u_func(it,X)
        B = B_func(it,X)
        return (X,b,B)
  
    return TimeStepIterator2(
        initial_values=V_init,
        f=f)
# calculate steady state
func_dict=msh.make_func_dict(svs,dvs)
B_func, u_func = gh.make_B_u_funcs_2(mvs,par_dict,func_dict)  


# +
import math
it_sym = make_steady_state_iterator_sym(
    mvs,
    V_init=(X_0,u_func(0,X_0),B_func(0,X_0)),
    par_dict=par_dict,
    func_dict=func_dict
)

Bs=[]
bs=[]
max_day=(cpa.number_of_months-1)*30
print(max_day,gh.day_2_month_index(max_day))

for i in range(max_day):
    res=it_sym.__next__()
    bs.append(res[1])
    Bs.append(res[2])
B_mean=np.stack(Bs).mean(axis=0)
b_mean=np.stack(bs).mean(axis=0)
B_mean,b_mean
np.linalg.inv(Bs[0])


# +
# calculate pseudo steady state
X_ss = np.linalg.solve(B_mean, (-b_mean))

steady_state_dict={str(name): X_ss[i,0] for i,name in enumerate(mvs.get_StateVariableTuple())}
# -

# create a start parameter tuple for the mcmc. The order has to be the same as when you created the namedtupl3 
# If you don't you get a "TypeError". 
epa_0=msh.EstimatedParameters(
    fwt=0.5907770914828289, 
    fgv=0.10708374044873868, 
    fco=0.9502719613629499, 
    fml=0.6985590765466911, 
    fd=0.8108017779961694, 
    k_C_NWT=0.0018600810916478165, 
    k_C_AGWT=0.00017354142452106252, 
    k_C_TR=0.00016065843641210772, 
    k_C_GVF=0.00022102017216433633, 
    k_C_GVR=0.00017926856125131916, 
    f_C_AGSL_2_C_AGMS=0.20853426509202325, 
    f_C_BGRL_2_C_SHMS=0.24638112975102788, 
    C_NWT_0=0.39641121927763323, 
    C_AGWT_0=1.0098899271611432, 
    C_GVF_0=0.1784893310039542, 
    C_GVR_0=2.1680315400436174, 
    C_AGML_0=0.1251689278629053, 
    C_AGSL_0=0.005800531050824444, 
    C_BGDL_0=0.0484130929152639, 
    C_AGMS_0=0.10074331291791151, 
    C_YHMS_0=0.5036084965444287, 
    C_SHMS_0=8.080067914918454,
    
    #fwt=0.22720510630531535, 
    #fgv=0.13339749844173146, 
    #fco=0.48468873886205804, 
    #fml=0.6625400558788637, 
    #fd=0.609731866092099, 
    #k_C_NWT=0.0015651547605874845, 
    #k_C_AGWT=0.00020579549570896914, 
    #k_C_TR=0.00021655197553090258,
    #k_C_GVF=0.00028880916566528537, 
    #k_C_GVR=0.0001574617151611989, 
    #f_C_AGSL_2_C_AGMS=0.16522686655577523, 
    #f_C_BGRL_2_C_SHMS=0.5461047001295201,
    #C_NWT_0= svs_0.cVeg * 0.06,#0.1674883218647344,
    #C_AGWT_0= svs_0.cVeg * 0.095,#0.29255052631195533,
    #C_GVF_0= svs_0.cVeg * 0.1,#0.3359968317564281,
    #C_GVR_0= svs_0.cVeg * 0.65,#2.3328171507922555,
    #C_AGML_0= svs_0.cLitter * 0.6,#0.02963751068264173,
    #C_AGSL_0= svs_0.cLitter * 0.1,#0.0004373128929123478,
    #C_BGDL_0= svs_0.cLitter * 0.1,#0.0005250285545731843,
    #C_AGMS_0= svs_0.cSoil * 0.05,#0.07115408727464849,
    #C_YHMS_0= svs_0.cSoil * 0.1,#0.17191831703645521,
    #C_SHMS_0= svs_0.cSoil * 0.7,#8.037110582162018
)




# +
# now test it 
import matplotlib.pyplot as plt
from general_helpers import plot_solutions

param2res_sym = msh.make_param2res_sym(mvs,cpa,dvs)
xs= param2res_sym(epa_0)
#obs=np.column_stack([ np.array(v) for v in svs])
obs=np.column_stack((np.repeat(svs.cVeg, 12),np.repeat(svs.cLitter, 12),np.repeat(svs.cSoil, 12),svs.rh,svs.ra))
#xs.shape

# +
obs=obs[0:cpa.number_of_months,:] #cut 
obs[:,3:4]=obs[:,3:4]
n=cpa.number_of_months

# convert to yearly output if necessary
#obs_yr=np.zeros(int(cpa.number_of_months/12)*obs.shape[1]).refwt=0.538141061754704, fgv=0.09965430660231399, fco=0.899657000739055, fml=0.7055535672710194, fd=0.8300632677297293, k_C_NWT=0.0017297080940248545, k_C_AGWT=0.00017468511195159253, k_C_TR=0.00016223935511977146, k_C_GVF=0.0002115629437369954, k_C_GVR=0.00016807958377661582, f_C_AGSL_2_C_AGMS=0.23986217759314643, f_C_BGRL_2_C_SHMS=0.17612859148192067, C_NWT_0=0.24871325996896232, C_AGWT_0=1.7557640443395466, C_GVF_0=0.20250040970682928, C_GVR_0=2.4078877379683274, C_AGML_0=0.13706684206636224, C_AGSL_0=0.03301021821996376, C_BGDL_0=0.07767639566417399, C_AGMS_0=0.5340042592206394, C_YHMS_0=0.1376951143246505, C_SHMS_0=7.329872222951247fwt=0.538141061754704, fgv=0.09965430660231399, fco=0.899657000739055, fml=0.7055535672710194, fd=0.8300632677297293, k_C_NWT=0.0017297080940248545, k_C_AGWT=0.00017468511195159253, k_C_TR=0.00016223935511977146, k_C_GVF=0.0002115629437369954, k_C_GVR=0.00016807958377661582, f_C_AGSL_2_C_AGMS=0.23986217759314643, f_C_BGRL_2_C_SHMS=0.17612859148192067, C_NWT_0=0.24871325996896232, C_AGWT_0=1.7557640443395466, C_GVF_0=0.20250040970682928, C_GVR_0=2.4078877379683274, C_AGML_0=0.13706684206636224, C_AGSL_0=0.03301021821996376, C_BGDL_0=0.07767639566417399, C_AGMS_0=0.5340042592206394, C_YHMS_0=0.1376951143246505, C_SHMS_0=7.329872222951247shape([int(cpa.number_of_months/12),obs.shape[1]])  
#for i in range(obs.shape[1]):
#    obs_yr[:,i]=monthly_to_yearly(obs[:,i])
#obs=obs_yr
#n=int(cpa.number_of_months/12)

print("Forward run with initial parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(12, 4), dpi=80)
plot_solutions(
        fig,
        times=range(n),
        var_names=msh.Observables._fields,
        tup=(xs,obs)
        #tup=(obs,)
)
fig.savefig('solutions.pdf')

# +
import matplotlib.pyplot as plt

print("Forward run with initial parameters")
plt.figure(figsize=(12,10), dpi=80)
plt.figure(1)

ax0 = plt.subplot(221)
ax0.plot(xs[:,0],label='Simulation',color="blue")
ax0.plot(obs[:,0],label='TRENDY',color="red")
plt.xlabel("Months since 1700",size=13)
plt.ylabel("Vegetation C (kg m-2)",size=13)
#ax0.legend(loc='best')
ax0.tick_params(labelsize=12)

ax1 = plt.subplot(222)
ax1.plot(xs[:,1],label='Simulation',color="blue")
ax1.plot(obs[:,1],label='TRENDY',color="red")
plt.xlabel("Months since 1700",size=13)
plt.ylabel("Litter C (kg m-2)",size=13)
ax1.legend(loc='best')
ax1.tick_params(labelsize=12)

ax2 = plt.subplot(223)
ax2.plot(xs[:,2],label='Simulation',color="blue")
ax2.plot(obs[:,2],label='TRENDY',color="red")
plt.xlabel("Months since 1700",size=13)
plt.ylabel("Soil C (kg m-2)",size=13)
ax2.legend(loc='best')
ax2.tick_params(labelsize=12)

ax3 = plt.subplot(224)
ax3.plot(xs[:,3],label='Simulation',color="blue")
ax3.plot(obs[:,3],label='TRENDY',color="red")
plt.xlabel("Months since 1700",size=13)
plt.ylabel("Heterotrophic Respiration (kg m-2 s-1)",size=13)
ax3.tick_params(labelsize=12)

plt.savefig('ISAM_test.pdf')



# +
epa_min=msh.EstimatedParameters(
    fwt=0.3,
    fgv=0.05,
    fco=0.7,
    fml=0.5,
    fd=0.6,
    k_C_NWT=1/(365*10),
    k_C_AGWT=1/(365*30),
    k_C_TR=1/(365*30),
    k_C_GVF=1/(365*30),
    k_C_GVR=1/(365*40),
    f_C_AGSL_2_C_AGMS=0.2*0.3,
    f_C_BGRL_2_C_SHMS=0.1,
    C_NWT_0=0,
    C_AGWT_0=0,
    C_GVF_0=0,
    C_GVR_0=0,
    C_AGML_0=0,
    C_AGSL_0=0,
    C_BGDL_0=0,
    C_AGMS_0=0,
    C_YHMS_0=0,
    C_SHMS_0=0,
)
   
epa_max=msh.EstimatedParameters(
    fwt=0.8,
    fgv=0.5,
    fco=0.98,
    fml=0.95,
    fd=0.95,
    k_C_NWT=1/(365*0.8),
    k_C_AGWT=1/(365*8),
    k_C_TR=1/(365*8),
    k_C_GVF=1/(365*8),
    k_C_GVR=1/(365*8),
    f_C_AGSL_2_C_AGMS=0.9*0.3,
    f_C_BGRL_2_C_SHMS=0.7,
    C_NWT_0=svs_0.cVeg,
    C_AGWT_0=svs_0.cVeg,
    C_GVF_0=svs_0.cVeg,
    C_GVR_0=svs_0.cVeg,
    C_AGML_0=svs_0.cLitter,
    C_AGSL_0=svs_0.cLitter,
    C_BGDL_0=svs_0.cLitter,
    C_AGMS_0=svs_0.cSoil,
    C_YHMS_0=svs_0.cSoil,
    C_SHMS_0=svs_0.cSoil,
)

# +
#import test_helpers as th
#ta=th.make_test_args(conf_dict,msh,mvs)
# -

# ### mcmc to optimize parameters 
#

# +
from general_helpers import autostep_mcmc, make_param_filter_func, make_feng_cost_func

isQualified = make_param_filter_func(epa_max, epa_min)
param2res = msh.make_param2res_sym(mvs,cpa,dvs)

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa_0,
    filter_func=isQualified,
    param2res=param2res,
    costfunction=make_feng_cost_func(obs),
    nsimu=500, # for testing and tuning mcmc
    #nsimu=20000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=15,   # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=1,   # default value | increase value to reduce initial step size
    K=2 # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

# +
# optimized parameter set (lowest cost function)
par_opt=np.min(C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields),1),axis=1)
epa_opt=msh.EstimatedParameters(*par_opt)
mod_opt = param2res(epa_opt)  

print("Forward run with optimized parameters (blue) vs TRENDY output (orange)")
fig = plt.figure(figsize=(12, 4), dpi=80)
plot_solutions(
        fig,
        #times=range(cpa.number_of_months),
        times=range(int(cpa.number_of_months)), # for yearly output
        var_names=msh.Observables._fields,
        tup=(mod_opt,obs)
)

fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath=Path(conf_dict["dataPath"]) # save output to data directory (or change it)

import pandas as pd
pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('ISAM_da_aa.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('ISAM_da_j_aa.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('ISAM_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('ISAM_optimized_solutions.csv'), sep=',')

# +
import matplotlib.pyplot as plt

print("Forward run with initial parameters")
plt.figure(figsize=(12,10), dpi=80)
plt.figure(1)

ax0 = plt.subplot(221)
ax0.plot(obs[:,0],label='TRENDY',color="red")
ax0.plot(mod_opt[:,0],label='Simulation',color="blue")
plt.xlabel("Months since 1700",size=13)
plt.ylabel("Vegetation C (kg m-2)",size=13)
#ax0.legend(loc='best')
ax0.tick_params(labelsize=12)

ax1 = plt.subplot(222)
ax1.plot(mod_opt[:,1],label='Simulation',color="blue")
ax1.plot(obs[:,1],label='TRENDY',color="red")
plt.xlabel("Months since 1700",size=13)
plt.ylabel("Litter C (kg m-2)",size=13)
ax1.legend(loc='best')
ax1.tick_params(labelsize=12)

ax2 = plt.subplot(223)
ax2.plot(obs[:,2],label='TRENDY',color="red")
ax2.plot(mod_opt[:,2],label='Simulation',color="blue")
plt.xlabel("Months since 1700",size=13)
plt.ylabel("Soil C (kg m-2)",size=13)
ax2.legend(loc='best')
ax2.tick_params(labelsize=12)

ax3 = plt.subplot(224)
ax3.plot(mod_opt[100:300,3],label='Simulation',color="blue")
ax3.plot(obs[100:300,3],label='TRENDY',color="red")
plt.xlabel("Months since 1700",size=13)
plt.ylabel("Heterotrophic Respiration (kg m-2 s-1)",size=13)
ax3.tick_params(labelsize=12)

plt.savefig('ISAM_opt.pdf')


# -



print("Optimized parameters: ", epa_opt)
#par_dict_opt={
#    beta_leaf: epa_opt.beta_leaf,
#    beta_wood: epa_opt.beta_wood,
#    T_0: epa_opt.T_0,
#    E: epa_opt.E,
#    KM: epa_opt.KM,
    #r_C_leaf_rh: 0,
    #r_C_wood_rh: 0,
    #r_C_root_rh: 0,
#    r_C_leaf_litter_rh: epa_opt.r_C_leaf_litter_rh,
#    r_C_wood_litter_rh: epa_opt.r_C_wood_litter_rh,
#    r_C_root_litter_rh: epa_opt.r_C_root_litter_rh,
#    r_C_soil_fast_rh: epa_opt.r_C_soil_fast_rh,
#    r_C_soil_slow_rh: epa_opt.r_C_soil_slow_rh,
#    r_C_soil_passive_rh: epa_opt.r_C_soil_passive_rh,
#    r_C_leaf_2_C_leaf_litter: epa_opt.r_C_leaf_2_C_leaf_litter,
#    r_C_wood_2_C_wood_litter: epa_opt.r_C_wood_2_C_wood_litter,
#    r_C_root_2_C_root_litter: epa_opt.r_C_root_2_C_root_litter,
#    r_C_leaf_litter_2_C_soil_fast: epa_opt.r_C_leaf_litter_2_C_soil_fast,
#    r_C_leaf_litter_2_C_soil_slow: epa_opt.r_C_leaf_litter_2_C_soil_slow,
#    r_C_leaf_litter_2_C_soil_passive: epa_opt.r_C_leaf_litter_2_C_soil_passive,
#    r_C_wood_litter_2_C_soil_fast: epa_opt.r_C_wood_litter_2_C_soil_fast,
#    r_C_wood_litter_2_C_soil_slow: epa_opt.r_C_wood_litter_2_C_soil_slow,
#    r_C_wood_litter_2_C_soil_passive: epa_opt.r_C_wood_litter_2_C_soil_passive,
#    r_C_root_litter_2_C_soil_fast: epa_opt.r_C_root_litter_2_C_soil_fast,
#    r_C_root_litter_2_C_soil_slow: epa_opt.r_C_root_litter_2_C_soil_slow,
#    r_C_root_litter_2_C_soil_passive: epa_opt.r_C_root_litter_2_C_soil_passive 
#}
#print("Optimized parameters dictionary: ", par_dict_opt)

# ### Traceability analysis  
#
# #### Outline
# The traceability analysis defines several diagnostic variables using as much algebraic structure of the mass balance equation as is available.
# Not all diagnostic variables are possible for all compartmental models. 
#
# We chose here to introduce the diagnostic variables not all at once but rather in the order of decreasing generality.
#
# The first diagnostic variables are available for all compartmental models and need no additional assumptions. 
# In the later parts of this section we then assume to be able to identify more and more specific terms in the mass balance equation and use those to derive and trace ever more specific diagnostics.
# Thus the very first part is valid for all models but how many of the later parts are applicable to a specific model  depends on how much we know about it.  
#
#
# #### Derivation of the matrix decomposition 
# Compartmental models (well mixed mass balanced) can be written in as an ordinary differential equation in matrix form that relates the momentary value of the (time) derivative $\frac{d X}{d t}$ of an yet unknown function $X$ to the momentary value of $X$ itself.   
# $$
# \frac{d X}{d t}= M(X,t) X + I(X,t) \quad (1)   
# $$ 
# where $X$ is the statevector representing the pool contents, $M$ the "Compartmental matrix" and $I$ the input vector.
# Together with a startvalue $X_0$ it constitutes an "initial value problem" (ivp) which can be solved numerically by moving step by step forward in time.
#
# Note: 
#
# It is mathematical standard notation to use $X$ in the *formulation* of the ivp (representing the momentary value) althoug *after we have solved it* the solution is expressed as function of time $X(t)$. This avoids confusion since everything appering with arguments is recognizable as explicitly calculable *before* we have solved the ivp.
#
# The system "nonautonomous" (if they depend on time $t$) and "nonlinear" if the dependent on $X$.
# It is always possible to factorize $M(X,t)$ into a product $M=A(X,t) K(X,t)$ where $K$ is a  diagonal matrix.
# and $I=B(t)*u(t)$ where $u$ is a scalar.
# Using these we arrive at 
# $$
# \frac{d X}{d t}= A(X,t) K(X,t) X + B(X,t) u(X,t)  
# $$
# ##### Linearity assumption
# If we assume the model to be linear and nonautonomous the dependency on $X$ vanishes and we have
#
# $$
# \frac{d X}{d t}= A(t) K(t) X + B(t) u(t) . 
# $$
#
# ##### Factorizability  assumption
# Although this is not possible in general in many published models the nonautonous part  can be further localized into a diagonal matrix $\xi(t)$ so that we can achieve constant $A$ and $K$ which allows more specific interpretation.
#
# $$
# \frac{d X}{d t}= A \xi(t) K X + B(t)u(t)
# $$
#
# ##### Factorizability of $\xi$ assumption 
# In some cases we can resolve $\xi$ further.
# $$
# \frac{d X}{d t}= A \xi_temp(t) \xi_mois(t) K X + B(t)u(t)
# $$
#
# #### Definition of diagnostic variables
#
# ##### Storage capacity $X_c$ and storage potential $X_p$
# These variables can be defined for any compartmental system and do not require either linearity nor factorizability. 
# We can rearrange eq. $(1)$ and give names to the two summands. 
# $$
# X = M^{-1}(X,t) \left( \frac{d X}{d t}-I(X,t) \right) \\ 
#   = \underbrace{M^{-1}(X,t) \frac{d X}{d t}}_{X_c} - \underbrace{M^{-1}(X,t)I(X,t)}_{X_p} \\
#   = X_c - X_p
# $$
# Note:
# This is not to be read as a recipe to compute $X$.
# The equation becomes a bit clearer if we adapt the nomenclature to express that we *have solved the ivp* and know its solution $X(t)$  
# <!-- and therefore also  the derivative $\frac{d X}{d t}=M(X(t),t) X(t) + I(X(t),t)=\prime{X}(t)$ -->
# By substituting the solution $X(t)$ we get the recipes to compute:
# $$
# X_p(t) = M^{-1}(X(t),t)I(X(t),t)  \\ 
# X_c(t) = X(t)-X_p(t) \\ 
# $$
# we see that all the ingredients become explicit functions of time.   
# Since all values are taken at the same time $t$ we can drop the time dependence
# in the notation and write an equation we can use in the iterator.
# $$
# X_p = M^{-1}I(X,t)  \\ 
# X_c = X + X_p \\ 
# $$
#
# ##### Residence time
# The influx $I$ can always be written as $I=b u$ where the scalar $u=\sum_{k=1\dots n} I_k$  and the dimensionless vector $b=I/u$ where $\sum_{k=1\dots n} b_k =1$.
# Assumimg that the pool contents (the components of $X$)  have dimension $mass$ we can infer from eq. (1) that $M$ has dimension $\frac{1}{time}$.
# The components of the (inverse) matrix $M^{-1}$ have therefore dimension $time$. Accordingly the product $RT= M^{-1} b$ is a vector of the same shape as $X$  whose components have dimesion $time$.
# In the context of the Traceability Framework $RT$ is therefore called *residence time*.
#
# Notes on nomenclature: 
# 1. The term *residence time* is not universally used with the same connotation outside the context of the *Traceability Analysis*.
#
# 1. It is not *the time of residence* of the particles in the system for the following reasons:
#     1. In well mixed systems particles can reside in a pool for different times from zero to infinity.
#     1. You could compute the mean of these times over all particles exiting a pool, but even then the result is in general not equal to the above mentioned $rt$.
#     1. The mean residence time would only coincide with the definition above if the system was in equilibrium (which it clearly is not as e.g $NPP(t)$ shows.)
#     1. The origin of the term is probably most easily understood as the generalization of a one dimensional rate equation $\frac{d}{dt} x = m x + u$ 
#        If $r$ and $u$ are constant then the mean residence time is $rt= m^{-1}$. If we start with the rate as property of the model the *residence time* 
#        can be defined as the inverse of this rate. The above definition is the generalization of this simple relationship to matrices and vectors.
#        The matrix $M^{-1}$ takes the role of the number $\frac{1}{m}$ . In the context of the *Traceability Analysis* $M^{-1}$ is called *Chasing Time*. 
#

# +
it_sym_trace = msh.make_traceability_iterator(mvs,dvs,cpa,epa_opt)
ns=10*360 #1500
StartVectorTrace=gh.make_StartVectorTrace(mvs)
nv=len(StartVectorTrace._fields)
res_trace= np.zeros((ns,nv))
for i in range(ns):
    res_trace[i,:]=it_sym_trace.__next__().reshape(nv)
#res_trace

import matplotlib.pyplot as plt
n=len(mvs.get_StateVariableTuple())
fig=plt.figure(figsize=(20,(n+1)*10), dpi=80)
axs=fig.subplots(n+1,2)
days=list(range(ns))


for i in range(n):
    
    ax = axs[i,0]
    #  the solution
    pos=i
    ax.plot(
        days,
        res_trace[:,i],
        label=StartVectorTrace._fields[pos],
        color='blue'
    )
    # X_p
    pos=i+n
    ax.plot(
        days,
        res_trace[:,pos],
        label=StartVectorTrace._fields[pos],
        color='red'
    )
    # X_c
    pos=i+2*n
    ax.plot(
        days,
        res_trace[:,pos],
        label=StartVectorTrace._fields[pos],
        color='yellow'
    )
    ax.legend()
    
    ax = axs[i,1]
    # RT
    pos=i+3*n
    ax.plot(
        days,
        res_trace[:,pos],
        label=StartVectorTrace._fields[pos],
        color='black'
    )
    ax.legend()
    
axs[n,0].plot(
    days,
    [msh.make_npp_func(dvs)(d) for d in days],
    label='NPP',
    color='green'
)
axs[n,0].legend()
# -

