

# identify script location (if using R-studio). If not working: input path manually into setwd()
library(Metrics)
library(rstudioapi)  
script_location<-normalizePath(rstudioapi::getActiveDocumentContext()$path)
script_location<-substr(script_location, 1, regexpr("CCNRM_ESM2_1_main.R",script_location)-1)
setwd("C:/Users/qbi/Desktop/Modeltraining/CMIP6/R/CCNRM_ESM2_1_pct1CO2bgc/try")


# Throughout the whole script there are '# fixme' comment
# you can search for them to find places where you have to adapt the code
# to your needs.

# use platform independent path descriptions so that you can run
# your stuff on windows or linux or mac

# I use the non mandatory type hints to make the code more readable
# This is especially useful for the description of functions that
# are used as arguments for other functions

# fixme: 
#   The idea is that everybody writes her own version of this 
#   directory 'yy_cable' maybe
#   - 'jon_ybis' or 
#   - 'mm_cable' and so on. 
#   All functions in module model_specifiC_helpers.py provide model specific results
#   and can not be applied directly but serve as examples.

source ("CCNRM_ESM2_1_model_specific_helpers.R")
source ("general_helpers.R") 

# fixme: 
#   put the (relative or asolute) location of your data into a small file called 'config.json' and
#   in my case the content looks like this:
#   {"dataPath": "/home/data/yuanyuan"}
#   DO NOT add the file to the repository. It is not only model- but also site specific. 
#   So you are likely to have one for every model on every computer
#   you run this code on.
#   (this example uses an absolute path starting with a '/'
library(rjson)
dataPath <- fromJSON(file = "config.json")[[1]]

# fixme: 
#    Note that the function is imported from 
#    model_specifiC_helpers which means that you have to provide
#    your version of this fuction which will most likely return different
#    variables 
# get data for a selected site: data path, longitude, latitude
lon=33.3
lat=50.0
f_clay_silt=0.4
# read data from NetCDF if using the script for he 1st time
# dat<-get_example_site_vars(dataPath, lon, lat)

# read data from a csv file if previously saved
dat<-get_data_from_file(dataPath)

# combine them to a single array which we will later use as input to the costfunction
nyears=130
nyears =10
tot_len = 12*nyears
library(dplyr)
obs = dplyr::select(dat, 
     C_leaf,
     C_stem,
     C_wood,
     C_root,
     C_litter,
     C_litter_above,
     C_litter_below,    
     C_fastsom,
     C_mediumsom,
     C_slowsom,
     C_soil,
     rh,
     f_veg2litter,
     f_litter2som,
)
obs<-obs[1:tot_len,]

# fixme 

C_min=c(0,0, 0.1,   0.001, 0.01,0.1, 0.01, 0.01,0.01, 0.1, 0.001, 0.01,0.01,  1/(365*2),1/(365*60),1/(365*60), 1/(365*30),       1/(365*60),1/(365*10),   1/(365*30),  1/(365*30),     1/(365*20), 1/(365*50), 1/(365*500), 0.05, 0.1, 0.01,0.01,  0.1,1,1)
C_max=c(0.9,0.9,0.9,  0.2, 0.8, 0.9,  0.8,0.9,0.8,  0.9,  0.9, 0.9, 0.6,   1/30, 1/365, 1/365,1/(365*0.5),    1/365,  1/(365*0.5),    1/(365*0.5),  1/(365*0.5),    1/(365*1),1/(365*3.5),1/(365*20),  4, 4, 1.9,1.9,  4, 100, 100)

C_min=c(0,0, 0, 0,0,0,0, 0, 0,0,0,0, 0, 0,0,0,0, 0, 0,0,0,0, 0, 0,0,0,0, 0, 0,0,0)
C_max=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

C_min=c(
  0.09,
  0.09,
  0.09,
  0.01,  #  4
  0.01,#  5
  0.01,#  6
  0.01,  # 7
  0.01,#  8
  0.01,#  9
  0.01,  #  10
  0.01,#  11
  0.09,#  12
  0.09,     
  1/(2*365),       #  14
  1/(60*365),       #  15
  1/(60*365),       #  16
  1/(365*10),       #  17
  1/(365*60),	#  18
  0.1/(0.1*365),	#  19
  1/(365*200),	#  20
  1/(365*60), #21
  1/(365*20),	#  22
  0.06/(5*365),	# 23
  0.06/(222.22*365),	# 24
  cpa$C_litter_above_0/100, # 25
  cpa$C_litter_above_0/100, #- 26
  cpa$C_litter_below_0/100, #- 27
  cpa$C_litter_below_0/100, #- 28
  0.1,
  1,
  1
)

C_max=c(
  1,
  1,       
  1,
  0.5,  #  4
  0.8,#  5
  0.95,#  6
  0.95,  # 7
  0.9,#  8
  0.5,#  9
  0.9,  #  10
  0.9,#  11
  0.3,#  12
  0.3,
  1/(0.3*365),       #  14
  1/365,       #  15
  1/365,       #  16
  1/(0.8*365),       #  17
  1/365,	#  18
  1/(365*0.1),	#  19
  1/(365*0.5),	#  20
  1/(365*0.1), #21
  0.6/(365*0.137),	#  22
  0.6/(365*5),	# 23
  0.6/(222.22*365),	# 24
  cpa$C_litter_above_0, # 25
  cpa$C_litter_above_0, #- 26
  cpa$C_litter_below_0, #- 27
  cpa$C_litter_below_0, #- 28
  4, 
  100, 
  100
)





# fixme
#   this function is model specific: It discards parameter proposals
#   where beta1 and beta2 add up to more than 0.99
isQualified = make_param_filter_func(C_max,C_min)
uniform_prop = make_uniform_proposer(
  C_min,
  C_max,
  D=10.0,
  filter_func=isQualified
)

cpa = list(
  C_leaf_0=dat$C_leaf[1],
  C_stem_0=dat$C_stem[1],
  C_wood_0=dat$C_wood[1],
  C_root_0=dat$C_root[1],
  C_litter_above_0=dat$C_litter_above[1],
  C_litter_below_0=dat$C_litter_below[1],
  C_fastsom_0=dat$C_fastsom[1],
  C_mediumsom_0=dat$C_mediumsom[1],
  C_slowsom_0=dat$C_slowsom[1],
  npp=dat$npp,
  number_of_months=tot_len,
  tsl=dat$tsl, # soil temperature - for dynamic turnover rates
  mrso=dat$mrso, # soil moisture - for dynamic turnover rates
  ts=dat$ts # air temperature - for dynamic turnover rates
)

param2res = make_param2res(cpa) #pa=[beta1,beta2, lig_leaf, f41,f42, kleaf,kroot,kwood,kmet,kmic, kslow,kpass, cmet_init, cstr_init, cmiC_init, cpassive_init ]

epa_0 = list(
  beta_leaf=0.33,    #  1 beta 1 - allocation to leaf 
  beta_stem=0.25,    #  2 beta 2 - allocation to stem
  beta_wood=0.35,    #  3 beta 3 - allocation to wood
  f_leaf2abovestruclitter=0.15,  #  4
  f_stem2abovestruclitter=0.4,#  5
  f_wood2abovestruclitter=0.2,#  6
  f_leaf2abovemetablitter=0.8,  # 7
  f_stem2abovemetablitter=0.5,#  8
  f_wood2abovemetablitter=0.05,#  9
  f_root2belowstruclitter=0.2,  #  10
  f_root2belowmetablitter=0.7,#  11
  f_aboveligninfrac=0.8,#  12
  f_belowligninfrac=0.8,       #  13
  k_leaf=7/365,       #  14
  k_stem=1/365,       #  15
  k_wood=1/(365*40),       #  16
  k_root=1/(365*28),       #  17
  k_above_structurelit=1/(365*4),	#  18
  k_above_metaboliclit=1/(365*2),	#  19
  k_below_structurelit=1/(365*200),	#  20
  k_below_metaboliclit=1/(365*50), #21
  k_fastsom=1/(365*1.5),	#  22
  k_mediumsom=0.6/(365*22),	# 23
  k_slowsom=0.3/(365*450),	# 24
  C_above_structurelit_0=0.8, # 25
  C_above_metaboliclit_0=0.75, #- 26
  C_below_structurelit_0=0.4, #- 27
  C_below_metaboliclit_0=0.4, #- 28
  T0=2,	# 22
  E=4,	# 23
  KM=10  # 24
)

########################## this is test of forward run and visualisation of initial fit ################################3
test = param2res(epa_0)
summary(as.data.frame(test))
{
  par(mfrow=c(4, 4)) # make 3x4 plots in 1 window
  
  for (i in 1:length(names(test))) {
    plot(test[[i]], type="l", col="red", xlab="month",
         ylim=c(min(min(test[[i]]),min(obs[[i]])),max(max(test[[i]]),max(obs[[i]]))), 
         ylab=names(test)[i], main=names(test)[i])
    lines(obs[[i]], col="blue")
  }
  plot(2, xlim=c(0,1), ylim=c(0,0.9), axes = F, main="legend", ylab="")
  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Modelled"),
         col=c("blue", "red"), lty=1, cex=1)
  
  par(mfrow=c(1, 1)) # return to single plot mode
}
#################################################################################################
# MCMC demo run

nsimu_demo = 2000  

mcmc_demo = mcmc(
  initial_parameters=epa_0,
  proposer=uniform_prop,
  param2res=param2res,
  costfunction=make_weighted_cost_func(obs),
  #costfunction=make_feng_cost_func(obs),
  nsimu=nsimu_demo,
  K_accept=0.5 # modifier to reduce acceptance rate
)




initial_parameters=epa_0
proposer=uniform_prop
param2res=param2res
#costfunction=make_weighted_cost_func(obs),
costfunction=make_feng_cost_func(obs)
nsimu=nsimu_demo
K_accept=0.5 # modifier to reduce acceptance rate









  # """
  # perform the Markov chain Monte Carlo simulation an returns a tuple of the array of sampled parameter(tuples) with shape (len(initial_parameters),nsimu) and the array of costfunction values with shape (q,nsimu)
  # 
  # :param initial_parameters: The initial guess for the parameter (tuple) to be estimated
  # :param proposer: A function that proposes a new parameter(tuple) from a given parameter (tuple).
  # :param param2res: A function that given a parameter(tuple) returns
  # the model output, which has to be an array of the same shape as the observations used to
  # build the costfunction.
  # :param costfunction: A function that given a model output returns a real number. It is assumed to be created for a specific set of observations, which is why they do not appear as an argument.
  # :param nsimu: The length of the chain
  # """
  set.seed(10)
  
  paramNum=length(initial_parameters)
  upgraded=0
  
  C_op = initial_parameters
  first_out = param2res(C_op)
  J_last = costfunction(first_out)
  #J_last = 400 # original code
  
  C_upgraded = rep(0,paramNum*nsimu)
  C_upgraded = matrix(C_upgraded, nrow = nsimu, byrow = TRUE)
  J_upgraded = rep(0, nsimu)
  
  simu = 1
  for (simu in 1:nsimu) {
    if (simu%%100==0) {print (paste0("simulation ",simu, " out of ", nsimu))} 
    
    
    
    paramNum = length(C_op)
    flag = T
    while (flag) {
      C_new = as.data.frame(C_op) + (runif((paramNum)) - 0.5)*(C_max - C_min)/10.0
      #C_new = as.data.frame(C_op) + (runif((paramNum)) - 0.5)*(C_max - C_min)/15.0
      if (filter_func(C_new)){flag = F}
    }
    names(C_new)=names(C_op)
    C_new=as.list(C_new)
    
    
    
    
    C_new = proposer(C_op)
    
    out_simu = param2res(C_new)
    J_new = costfunction(out_simu)
    
    delta_J =  J_last - J_new;
    
    randNum = runif(1)
    
    if (min(1.0, exp(delta_J)) > randNum/K_accept) {
      C_op=C_new;
      J_last=J_new;
      C_upgraded[upgraded,]=unlist(C_op, use.names=FALSE);
      J_upgraded[upgraded]=J_last; 
      upgraded=upgraded+1 
    }
  }
  # select only rows with upgraded parameters, discard 0s at the end
  C_upgraded<-C_upgraded[1:upgraded,]
  J_upgraded<-J_upgraded[1:upgraded] 
  acceptance_rate<-upgraded/nsimu
  


















# save demo parameters and costfunction values for postprocessing 

df=data.frame(mcmc_demo[[1]])
df_j=data.frame(mcmc_demo[[2]])

write.csv(df,paste0(dataPath,'/visit_demo_da_aa.csv'))
write.csv(df_j,paste0(dataPath,'/visit_demo_da_j_aa.csv'))

names(df)<-names(epa_0)

# visualize parameter distribution
{
  par(mfrow=c(4, 6)) # make 4x6 plots in 1 window - good for 24 parameters
  for (i in 1:length(df)) {hist(df[[i]], breaks=20, main=names(df)[i])}
  par(mfrow=c(1, 1)) # return to single plot mode
}
# build a new proposer based on a multivariate_normal distribution using the estimated covariance of the previous run if available
# parameter values of the previous run

normal_prop = make_multivariate_normal_proposer(
  
  covv = cov(df[as.integer(length(df)*0.9):length(df),]),  # the part of the demo run samples to use (here the last 90%)
  filter_func=isQualified
)
##############  MCMC formal run ###############
nsimu_formal = 5000
mcmc_formal = mcmc(
  initial_parameters=epa_0,
  proposer=normal_prop,
  param2res=param2res,
  #costfunction=make_weighted_cost_func(obs),
  costfunction=make_feng_cost_func(obs),
  nsimu=nsimu_formal,
  K_accept=0.5 # modifier to reduce acceptance rate
)

# save the parameters and costfunction values for postprocessing 

df=data.frame(mcmc_formal[[1]])
df_j=data.frame(mcmc_formal[[2]])

write.csv(df,paste0(dataPath,'/visit_formal_da_aa.csv'))
write.csv(df_j,paste0(dataPath,'/visit_formal_da_j_aa.csv'))

######################################## explore MCMC results ###################################################

# calculate median of each parameter distribution
epa_median=rep(0,length(epa_0))
for (i in 1:length(epa_0)) {epa_median[i]=median(df[[i]])}
names(epa_median)<-names(epa_0)
epa_median<-as.list(epa_median)
# calculate mode of each parameter distribution
epa_mode=rep(0,length(epa_0))
for (i in 1:length(epa_0)) {
  uniqv <- unique(df[[i]])
  epa_mode[i] <- uniqv[which.max(tabulate(match(df[[i]], uniqv)))]}
names(epa_mode)<-names(epa_0)
epa_mode<-as.list(epa_mode)
# select parameters from a set with minimal cost function distribution
epa_min_J<-as.list(df[df_j==min(df_j[df_j!=0])])
names(epa_min_J)<-names(epa_0)

# visualize parameter distributions with median, mode and min cost function
{
  par(mfrow=c(4, 6)) # make 4x6 plots in 1 window
  for (i in 1:length(df)) {
    hist(df[[i]], breaks=100, col="lightgray", border="lightgray", probability=T, main=names(epa_0)[i])
    lines(stats::density(df[[i]]), col="black", lwd=2)
    abline(v=epa_median[[i]],col="red",lwd=2)
    abline(v=epa_mode[[i]],col="orange",lwd=2)
    abline(v=epa_min_J[[i]],col="green",lwd=2)
  }
  par(mfrow=c(1, 1)) # return to single plot mode
}

# run the model with optimized parameters
optimized_median = param2res(epa_median)
optimized_mode = param2res(epa_mode)
optimized_min_J = param2res(epa_min_J)
{ # plot model output with optimized parameters
  par(mfrow=c(3, 4)) # make 3x4 plots in 1 window
  N=tot_len # number of years to plot
  
  # change number of years to plot for very dynamic pools and fluxes (optional)
  for (i in 1:length(names(optimized_median))) {
    if (names(optimized_median)[i]=="C_leaf" || names(optimized_median)[i]=="rh" || names(optimized_median)[i]=="f_veg2litter"
        || names(optimized_median)[i]=="f_litter2som") {N=tot_len%/%10} else (N=tot_len)
    
    plot(optimized_median[[i]][1:N], type="l", col="red", xlab="month",
         ylim=c(min(min(optimized_median[[i]]),min(obs[[i]])),max(max(optimized_median[[i]]),max(obs[[i]]))), 
         ylab=names(optimized_median)[i], main=names(optimized_median)[i])
    lines(obs[[i]], col="blue")
    lines(optimized_mode[[i]], col="orange")
    lines(optimized_min_J[[i]], col="green")
  }
  plot(2, xlim=c(0,1), ylim=c(0,0.9), axes = F, main="legend", ylab="")
  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Median", "Mode", "Min Cost Function"),
         col=c("blue", "red", "orange", "green"), lty=1, cex=1)
  
  par(mfrow=c(1, 1)) # return to single plot mode
}
print(as.data.frame(epa_0))
print(as.data.frame(epa_median))
print(as.data.frame(epa_mode))
print(as.data.frame(epa_min_J))
print(paste0("Acceptance rate (demo): ",mcmc_demo[[3]]))
print(paste0("Acceptance rate (formal): ",mcmc_formal[[3]]))
# calculate RMSE for each observed/predicted variable
RMSE_median<-data.frame (param=names(obs), RMSE=0)
RMSE_mode<-data.frame (param=names(obs), RMSE=0)
RMSE_min_J<-data.frame (param=names(obs), RMSE=0)
for (i in 1:length(names(obs))){
  RMSE_median[i,2]<-rmse(obs[[i]], optimized_median[[i]])
  RMSE_mode[i,2]<-rmse(obs[[i]], optimized_mode[[i]])
  RMSE_min_J[i,2]<-rmse(obs[[i]], optimized_min_J[[i]])
}
# best parameter values based on minimal RMSE
for (i in 1:length(names(obs))){
  print(names(obs)[i])
  if (RMSE_median[i,2] < RMSE_mode[i,2] && RMSE_median[i,2] < RMSE_min_J[i,2])  {print("best value - median")}
  if (RMSE_mode[i,2] < RMSE_median[i,2] && RMSE_mode[i,2] < RMSE_min_J[i,2])   {print("best value - mode")}
  if (RMSE_min_J[i,2] < RMSE_mode[i,2] && RMSE_min_J[i,2] < RMSE_median[i,2])   {print("best value - min J")}
}
