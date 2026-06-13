#!/usr/bin/env python
import sys
sys.path.insert(0,'..')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import json 

#parallel chains -  JW
from os import environ
from dask_mpi import initialize
initialize()
from distributed import Client
client=Client()

#import model specific functions
from model_specific_helpers import (
    EstimatedParameters, 
    make_param_filter_func,
    UnEstimatedParameters, 
    Parameters, 
    StateVariables,
    ModelParameters,
    Observables,
    monthly_to_yearly,
    pseudo_daily_to_yearly,
    get_example_site_vars,
    get_variables_from_files,
    make_param2res,
)

#import general functions
from general_helpers import (
        make_uniform_proposer,
        make_multivariate_normal_proposer,
        mcmc,
        autostep_mcmc,
        make_feng_cost_func,
        make_jon_cost_func, 
        plot_solutions,
	autostep_mcmc
)

#load file path from json file
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 
dataPath = Path(conf_dict['dataPath'])

#get data streams
npp, rh, ra, cveg, csoil = get_example_site_vars(Path(conf_dict['dataPath']))
nyears = 320
obs_tup=Observables(
    h_respiration=monthly_to_yearly(rh),
    #a_respiration=monthly_to_yearly(ra),
    c_veg=cveg,
    c_soil=csoil
)
obs = np.stack(obs_tup, axis=1)[0:nyears,:]

#set estimated parameters
epa0 = EstimatedParameters(
    beta_leaf=0.21,
    beta_root=0.27,
    k_leaf=0.014,
    k_root=0.022,
    k_wood=0.003,
    k_cwd=0.005,
    k_samet=0.01,
    k_sastr=0.001,
    k_samic=0.05,
    k_slmet=0.040,
    k_slstr=0.0039,
    k_slmic=0.005,
    k_slow=0.00001,
    k_arm=3.27E-06,
    f_samet_leaf=0.28,
    f_slmet_root=0.34,
    f_samic_cwd=0.29,
    C_leaf_0=0.06,
    C_root_0=0.08,
    C_cwd_0=0.5,
    C_samet_0=0.4,
    C_sastr_0=1.07,
    C_samic_0=0.49,
    C_slmet_0=0.83,
    C_slstr_0=2.07,
    C_slmic_0=1.04,
    C_slow_0=5.66
)

#set fixed parameters
cpa = UnEstimatedParameters(
    npp=npp,
    rh_0 = monthly_to_yearly(rh)[0],
    ra_0 = monthly_to_yearly(ra)[0],    
    C_veg_0=cveg[0],
    C_soil_0=csoil[0],
    clay=0.2028,
    silt=0.2808,
    nyears=320
)
#make model function
param2res = make_param2res(cpa)

# set max/min parameters limits 
c_min=np.array(epa0)*0.001
c_max=np.array(epa0)*1000
c_max[0]=1
c_max[1]=1
c_max[14]=1
c_max[15]=1
c_max[16]=1
c_max[17]=cveg[0]
c_max[18]=cveg[0]
c_max[19]=csoil[0]
c_max[20]=csoil[0]
c_max[21]=csoil[0]
c_max[22]=csoil[0]
c_max[23]=csoil[0]
c_max[24]=csoil[0]
c_max[25]=csoil[0]
c_max[26]=csoil[0]

#   this function is model specific: It discards parameter proposals
#   where beta1 and beta2 add up to more than 0.99
isQualified = make_param_filter_func(c_max,c_min,cveg[0],csoil[0])
uniform_prop = make_uniform_proposer(
    c_min,
    c_max,
    D=15, # this value 
    filter_func=isQualified
)

#set cost function 
costfunction=make_jon_cost_func(obs)

#define uniform parallel mcmc wrapper function
def uniform_parallel_mcmc(_):
    #calculate length of parameters
    par_len = len(epa0)
    #randomly perturb parameters for each chain by up to +- 50%
    flag=True
    while(flag):
        pertb = epa0*np.random.uniform(low=0.85, high=1.15, size=(par_len,)).astype(float)
        if(isQualified(pertb)):
            flag=False
    return(
        autostep_mcmc(
            initial_parameters=pertb,
            filter_func=isQualified,
            param2res=param2res,
            costfunction=costfunction,
            nsimu=50000,
            c_max=c_max,
            c_min=c_min
	    )
	)
        #mcmc(
        #    initial_parameters=pertb,
        #    proposer=uniform_prop,
        #    param2res=param2res,
        #    costfunction=costfunction,
        #    nsimu=20000
        #)

# set file names for saving data
uni_c_path = dataPath.joinpath('yibs_pmcmc_uniform_c.csv')
uni_j_path = dataPath.joinpath('yibs_pmcmc_uniform_j.csv')

# Parallel uniform distribution run 
print("starting parallel run")
[
    [c_uni1,j_uni1],
    [c_uni2,j_uni2],
    [c_uni3,j_uni3],
    [c_uni4,j_uni4],
    [c_uni5,j_uni5],
    [c_uni6,j_uni6],
    [c_uni7,j_uni7],
    [c_uni8,j_uni8],
    [c_uni9,j_uni9],
    [c_uni10,j_uni10]
] = client.gather(
        client.map(
            uniform_parallel_mcmc, 
            range(0,10)
        )
    )

#concatenate chains
C_cat = np.concatenate(
    (
        c_uni1,
        c_uni2,
        c_uni3,
        c_uni4,
        c_uni5,
        c_uni6,
        c_uni7,
        c_uni8,
        c_uni9,
        c_uni10 
    ), axis=1
)

#concatenate cost function
J_cat = np.concatenate(
    (
        j_uni1,
        j_uni2,
        j_uni3,
        j_uni4,
        j_uni5,
        j_uni6,
        j_uni7,
        j_uni8,
        j_uni9,
        j_uni10 
    ), axis=1
)

# +
# save the parameters and costfunctionvalues for postprocessing 
pd.DataFrame(C_cat).to_csv(uni_c_path,sep=',')
pd.DataFrame(J_cat).to_csv(uni_j_path,sep=',')

#subset to best pars
best_pars = C_cat[:,J_cat[1,:].argmin()]
# -

# #sort lowest to highest
# #indx = np.argsort(J_cat) 
# #C_demo = C_cat[np.arange(C_cat.shape[0])[:,None], indx]
# #J_demo = J_cat[np.arange(J_cat.shape[0])[:,None], indx]
#
# # formal run using normal distribution and cov matrix from uniform run
# #covv = np.cov(C_demo[:, 0:int(C_demo.shape[1]*0.4)]) #lowest 10% by cost 
# #normal_prop = make_multivariate_normal_proposer(
# #   covv = covv,
# #  filter_func=isQualified
# #)
#
# #define normal parallel mcmc wrapper
# #def normal_parallel_mcmc(_):
# #   return(
# #      adaptive_mcmc(
#  #         initial_parameters=C.demo[:,0],
#    #        covv=covv,
#     #       filter_func=isQualified,
#     #       param2res=param2res,
#     #       costfunction=costfunction,
#     #       nsimu=2000
#     #   )
#    #)
#
# # formal run 
# #[
# #   [c_form1,j_form1],
#  #  [c_form2,j_form2],
# #   [c_form3,j_form3],
# #   [c_form4,j_form4],
# #   [c_form5,j_form5],
# #   [c_form6,j_form6],
# #   [c_form7,j_form7],
# #   [c_form8,j_form8],
# #   [c_form9,j_form9],
# #   [c_form10,j_form10]
# #] = client.gather(
# #       client.map(
# #           normal_parallel_mcmc, 
# #           range(0,10)
#  #      )
#  #  )
#
# #concatenate chains
# #C_cat = np.concatenate(
# #   (
# #       c_form1,
#  #      c_form2,
#   #     c_form3,
#   #     c_form4,
#   #     c_form5,
#   #     c_form6,
#   #     c_form7,
#   #     c_form8,
#   #     c_form9,
#   #     c_form10 
#   # ), axis=1
# #)
#
# #concatenate cost function
# #J_cat = np.concatenate(
# #   (
# #       j_form1,
# #       j_form2,
# #       j_form3,
# #       j_form4,
# #       j_form5,
# #       j_form6,
# #       j_form7,
# #       j_form8,
# #       j_form9,
# #       j_form10 
# #   ), axis=1
# #)
#
# #sort lowest to highest
# #indx = np.argsort(J_cat) 
# #C_demo = C_cat[np.arange(C_cat.shape[0])[:,None], indx]
# #J_demo = J_cat[np.arange(J_cat.shape[0])[:,None], indx]
#
# #print chain5 output as test
# #formal_c_path = dataPath.joinpath('yibs_pmcmc_normal_c.csv')
# #formal_j_path = dataPath.joinpath('yibs_pmcmc_normal_j.csv')
# #pd.DataFrame(C_demo).to_csv(formal_c_path,sep=',')
# #pd.DataFrame(J_demo).to_csv(formal_j_path,sep=',')
#    
# #use output csv file for post processing
# #C_formal = pd.read_csv(formal_c_path).to_numpy()
# #J_formal = pd.read_csv(formal_j_path).to_numpy()
#
# #subset to lowest cost subset of mulitple chains (lowest 10%)
# #C_formal = C_formal[:, :int(C_formal.shape[1]*0.1)]
#

# # POSTPROCESSING 
# #
# # The 'solution' of the inverse problem is actually the (joint) posterior
# # probability distribution of the parameters, which we approximate by the
# # histogram consisting of the mcmc generated samples.  
# # This joint distribution contains as much information as all its (infinitly
# # many) projections to curves through the parameter space combined.
# # Unfortunately, for this very reason, a joint distribution of more than two
# # parameters is very difficult to visualize in its entirity. 
# # to do: 
# #   a) make a movie of color coded samples  of the a priori distribution of the parameters.
# #   b) -"-                                  of the a posteriory distribution -'- 
#
# # Therefore the  following visualizations have to be considered with caution:
# # 1.
# # The (usual) histograms of the values of a SINGLE parameters can be very
# # misleading since e.g. we can not see that certain parameter combination only
# # occure together. In fact this decomposition is only appropriate for
# # INDEPENDENT distributions of parameters in which case the joint distribution
# # would be the product of the distributions of the single parameters.  This is
# # however not even to be expected if our prior probability distribution can be
# # decomposed in this way. (Due to the fact that the Metropolis Hastings Alg. does not
# # produce independent samples ) 
# #df = pd.DataFrame({name :C_formal[:,i] for i,name in enumerate(EstimatedParameters._fields)})
# #subplots=df.hist()
# #fig=subplots[0,0].figure
# #fig.set_figwidth(15)
# #fig.set_figheight(15)
# #fig.savefig('histograms.pdf')
#
# # As the next best thing we can create a matrix of plots containing all 
# # projections to possible  parameter tuples
# # (like the pairs plot in the R package FME) but 16x16 plots are too much for one page..
# # However the plot shows that we are dealing with a lot of colinearity for this  parameter set
# #subplots = pd.plotting.scatter_matrix(df) 
# #fig=subplots[0,0].figure
# #fig.set_figwidth(15)
# #fig.set_figheight(15)
# #fig.savefig('scatter_matrix.pdf')
#
#
# # 2.
# # another way to get an idea of the quality of the parameter estimation is
# # to plot trajectories.
# # A possible aggregation of this histogram to a singe parameter
# # vector is the mean which is an estimator of  the expected value of the
# # desired distribution.
#
sol_mean =param2res(best_pars)

fig = plt.figure()
plot_solutions(
        fig,
        times=np.array(range(nyears)),
        var_names=Observables._fields,
        tup=(sol_mean, obs),
        names=('best','obs')
)
fig.savefig('solutions.pdf')
