library(Metrics)
# identify script location (if using R-studio). If not working: input path manually into setwd()
library(rstudioapi)  
script_location<-normalizePath(rstudioapi::getActiveDocumentContext()$path)
script_location<-substr(script_location, 1, regexpr("main.R",script_location)-1)
setwd(script_location)

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

source ("model_specifiC_helpers.R")
source ("../general_helpers.R") 

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

# read data from NetCDF if using the script for he 1st time
#dat<-get_example_site_vars(dataPath, lon, lat)

# read data from a csv file if previously saved
dat<-get_data_from_file(dataPath)
    
# combine them to a single array which we will later use as input to the costfunction
nyears=150
#nyears = 20
tot_len = 12*nyears
library(dplyr)
obs = dplyr::select(dat, 
    C_leaf,
    C_wood,
    C_root,
    C_litter_above,
    C_litter_below,    
    C_fastsom,
    C_slowsom,
    C_passsom,
        rh,
    f_veg2litter,
    f_litter2som,
)
obs<-obs[1:tot_len,]

# fixme 

 
C_min=c(0,  0,   0.1,  0.01,  0.001,  0.1,  0.01,  0.01, 0.1,  0.1,   0.001,  1/(365*2),  1/(365*60), 1/(365*30), 1/(365*60),   1/(365*10), 1/(365*30),    1/(365*20), 1/(365*50), 1/(365*500),    0, 0.1,    1,     1)
#pa=c(0.25, 0.2, 0.42,0.075, 0.005, 0.35, 0.12,  0.03, 0.37, 0.11, 0.01, 1/60,  1/(365*12),   1/(365*5),  1/(365*2),    1/(365*6),    1/(365*2.8),  1/(365*4.5), 1/(365*25),    1/(365*325),  0.3, 2,      4,     10)
C_max=c(1,  1,   0.9,   0.9,  0.9,  0.9,  0.2,  0.2,  0.9,  0.9,  0.9,   1/30,    1/365,   1/(365*0.5),   1/365,     1/(365*0.5),  1/(365*0.5),    1/(365*1),  1/(365*3.5),    1/(365*20),   0.4,  4,  100,   100)

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
    C_wood_0=dat$C_wood[1],
    C_root_0=dat$C_root[1],
    C_litter_above_0=dat$C_litter_above[1],
    C_litter_below_0=dat$C_litter_below[1],
    C_fastsom_0=dat$C_fastsom[1],
    C_slowsom_0=dat$C_slowsom[1],
    C_passsom_0=dat$C_passsom[1],
    npp=dat$npp,
    number_of_months=tot_len,
    tsl=dat$tsl, # soil temperature - for dynamic turnover rates
    mrso=dat$mrso, # soil moisture - for dynamic turnover rates
    ts=dat$ts # air temperature - for dynamic turnover rates
)

param2res = make_param2res(cpa) #pa=[beta1,beta2, lig_leaf, f41,f42, kleaf,kroot,kwood,kmet,kmic, kslow,kpass, cmet_init, cstr_init, cmiC_init, cpassive_init ]

epa_0 = list(
    beta_leaf=0.25,    #  1 (parameters used in original code) 
    beta_wood=0.2,    #  2
    f_leaflit2fastsom=0.41,  #  3
    f_leaflit2slowsom=0.07,#  4
    f_leaflit2passsom=0.02,#  5
    f_woodlit2fastsom=0.30,  #  6
    f_woodlit2slowsom=0.12,#  7
    f_woodlit2passsom=0.08,#  8
    f_rootlit2fastsom=0.30,  #  9
    f_rootlit2slowsom=0.14,#  10
    f_rootlit2passsom=0.07,#  11
    k_leaf=1/60,       #  12
    k_wood=1/(365*12),       #  13
    k_root=1/(365*5),       #  14
    k_leaflit=1/(365*2),	#  15
    k_woodlit=1/(365*6.7),	#  16
    k_rootlit=1/(365*3.5),	#  17
    k_fastsom=1/(365*6.7),	#  18
    k_slowsom=1/(365*28),	# 19
    k_passsom=1/(365*87),	# 20
    C_leaflit_0=0.3,	# 21
    T0=2,	# 22
    E=4,	# 23
    KM=10  # 24
)

########################## this is test of forward run and visualization of initial fit ################################3
test = param2res(epa_0)
summary(as.data.frame(test))
{
par(mfrow=c(3, 4)) # make 3x4 plots in 1 window

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
#################################    MCMC   ##############################################################

# MCMC demo run

nsimu_demo = 2000  
mcmc_demo = mcmc(
        initial_parameters=epa_0,
        proposer=uniform_prop,
        param2res=param2res,
        #costfunction=make_weighted_cost_func(obs),
        costfunction=make_feng_cost_func(obs),
        nsimu=nsimu_demo,
        K_accept=0.5 # modifier to reduce acceptance rate
)
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

####################################################################################################