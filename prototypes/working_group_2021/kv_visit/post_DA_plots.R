# identify script location (if using R-studio). If not working: input path manually into setwd()
library(rstudioapi)
script_location<-normalizePath(rstudioapi::getActiveDocumentContext()$path)
script_location<-substr(script_location, 1, regexpr("post_DA_plots.R",script_location)-1)
setwd(script_location)
model_name="VISIT"
# read data path from config.json
library(rjson)
dataPath <- fromJSON(file = "config.json")[[1]]
time_stamp<-substr(Sys.time(),1,regexpr(":",Sys.time())-1)

############################################ Plot posterior (formal) distributions of parameters ##########################33

# load posterior distributions after the formal run
df<-read.csv(paste0(dataPath,"/visit_formal_da_aa.csv"))
df<-df[,2:dim(df)[2]]
df<-t(df)

df<-as.data.frame(df)

df_j<-read.csv(paste0(dataPath,"/visit_formal_da_j_aa.csv"))
df_j<-df_j[-c(1)]
df_j<-t(df_j)
df_j<-as.data.frame(df_j)

par_names<-c("beta_leaf","beta_wood","f_leaflit2fastsom", "f_leaflit2slowsom", "f_leaflit2passsom", "f_woodlit2fastsom", "f_woodlit2slowsom", "f_woodlit2passsom",
             "f_rootlit2fastsom", "f_rootlit2slowsom", "f_rootlit2passsom", "k_leaf",   "k_wood",   "k_root",   "k_leaflit","k_woodlit","k_rootlit","k_fastsom",
             "k_slowsom","k_passsom","C_leaflit_0","T_0","E","KM" )

# load initial parameters
epa_0<-read.csv(paste0(dataPath,"/epa_0.csv"))
epa_0<-epa_0[,2:dim(epa_0)[2]]
epa_0<-t(epa_0)

epa_0<-as.data.frame(epa_0)

# calculate median of each parameter distribution
epa_median=rep(0,length(df))
for (i in 1:length(df)) {epa_median[i]=median(df[[i]])}
epa_median<-as.list(epa_median)
# calculate mean of each parameter distribution
epa_mean=rep(0,length(df))
for (i in 1:length(df)) {epa_mean[i]=mean(df[[i]])}
epa_mean<-as.list(epa_mean)
# select parameters from a set with minimal cost function distribution
epa_min_J<-as.list(df[df_j[2]==min(df_j[2]),])
# calculate max likelihood of each parameter distribution
epa_likelihood=rep(0,length(df))
for (i in 1:length(df)) {epa_likelihood[i]<-density(df[[i]])$x[which.max(density(df[[i]])$y)]}

# plot as a .jpeg file
jpeg(paste0(dataPath,'/posterior_distr_formal_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
  par(mfrow=c(4, 6)) # make 4x6 plots in 1 window
  for (i in 1:length(df)) {
    hist(df[[i]], breaks=100, col="lightgray", border="lightgray", probability=T, main=par_names[i])
    lines(stats::density(df[[i]]), col="black", lwd=2)
    abline(v=epa_min_J[[i]],col="green",lwd=3)
    abline(v=epa_likelihood[[i]],col="orange",lwd=2)
  }

  par(mfrow=c(1, 1))
}
dev.off()

############################################ Plot posterior (autostep) distributions of parameters ##########################33

# load posterior distributions after the demo run

df<-read.csv(paste0(dataPath,"/visit_autostep_da_aa.csv"))
df<-df[,2:dim(df)[2]]
df<-t(df)

df<-as.data.frame(df)

df_j<-read.csv(paste0(dataPath,"/visit_autostep_da_j_aa.csv"))
df_j<-df_j[-c(1)]
df_j<-t(df_j)
df_j<-as.data.frame(df_j)

par_names<-c("beta_leaf","beta_wood","f_leaflit2fastsom", "f_leaflit2slowsom", "f_leaflit2passsom", "f_woodlit2fastsom", "f_woodlit2slowsom", "f_woodlit2passsom",
             "f_rootlit2fastsom", "f_rootlit2slowsom", "f_rootlit2passsom", "k_leaf",   "k_wood",   "k_root",   "k_leaflit","k_woodlit","k_rootlit","k_fastsom",
             "k_slowsom","k_passsom","C_leaflit_0","T_0","E","KM" )

epa_0<-read.csv(paste0(dataPath,"/epa_0.csv"))
epa_0<-epa_0[,2:dim(epa_0)[2]]
epa_0<-t(epa_0)

epa_0<-as.data.frame(epa_0)

# calculate median of each parameter distribution
epa_median=rep(0,length(df))
for (i in 1:length(df)) {epa_median[i]=median(df[[i]])}
epa_median<-as.list(epa_median)
# calculate mean of each parameter distribution
epa_mean=rep(0,length(df))
for (i in 1:length(df)) {epa_mean[i]=mean(df[[i]])}
epa_mean<-as.list(epa_mean)
# select parameters from a set with minimal cost function distribution
epa_min_J<-as.list(df[df_j[2]==min(df_j[2]),])
# calculate max likelihood of each parameter distribution
epa_likelihood=rep(0,length(df))
for (i in 1:length(df)) {epa_likelihood[i]<-density(df[[i]])$x[which.max(density(df[[i]])$y)]}

# plot as a .jpeg file
jpeg(paste0(dataPath,'/posterior_distr_autostep_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
  par(mfrow=c(4, 6)) # make 4x6 plots in 1 window
  for (i in 1:length(df)) {
    hist(df[[i]], breaks=100, col="lightgray", border="lightgray", probability=T, main=par_names[i])
    lines(stats::density(df[[i]]), col="black", lwd=2)
    abline(v=epa_min_J[[i]],col="green",lwd=3)
    abline(v=epa_likelihood[[i]],col="orange",lwd=2)
  }

  par(mfrow=c(1, 1))
}
dev.off()

########################## Plot cost function evolution ############################################################

# load cost function data from demo and formal run
df_j<-read.csv(paste0(dataPath,"/visit_formal_da_j_aa.csv"))
df_j_demo<-read.csv(paste0(dataPath,"/visit_demo_da_j_aa.csv"))
df_j_autostep<-read.csv(paste0(dataPath,"/visit_autostep_da_j_aa.csv"))

# transform the data for plotting
df_j<-df_j[-c(1)]
df_j_demo<-df_j_demo[-c(1)]
df_j_autostep<-df_j_autostep[-c(1)]

df_j<-t(df_j)
df_j<-as.data.frame(df_j)
df_j_demo<-t(df_j_demo)
df_j_demo<-as.data.frame(df_j_demo)
df_j_autostep<-t(df_j_autostep)
df_j_autostep<-as.data.frame(df_j_autostep)

# plot on a log scale
jpeg(paste0(dataPath,'/cost_functions_log_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
  plot(df_j, type="l", col="green", xlab="iterations", ylab="log cost value", main="Cost functions over time (log scale)",
       xlim=c(0,max(max(df_j[1]),max(df_j_demo[1]), max(df_j_autostep[1]))), 
       ylim=c(min(min(df_j[2]),min(df_j_demo[2]),min(df_j_autostep[2])),max(max(df_j[2]),max(df_j_demo[2],max(df_j_autostep[2])))),
       log="y") 
  lines(df_j_demo, col="red")
  lines(df_j_autostep, col="orange")
  legend(x="topright", legend=c("demo run", "formal run", "demo auto-step"),
         col=c("red", "green", "orange"), lty=1, cex=1)
}
dev.off()

# Same plot with only last 50% iterations autostep and formal runs

df_j1<-df_j[df_j[1]>max(df_j[1]*0.5),]
df_j_autostep1<-df_j_autostep[df_j_autostep[1]>max(df_j[1]*0.5),]

jpeg(paste0(dataPath,'/cost_functions_last_50_percent_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
  plot(df_j1, type="l", col="green", xlab="iterations", ylab="log cost value", main="Last 50% of accepted cost functions",
       xlim=c(min(min(df_j1[1]),min(df_j_autostep1[1])),max(max(df_j1[1]),max(df_j_autostep1[1]))), 
       ylim=c(min(min(df_j1[,2]),min(df_j_autostep1[2])),max(max(df_j1[2]),max(df_j_autostep1[2]))),
       log="y") 
  lines(df_j_autostep1, col="orange")
  legend(x="topright", legend=c("formal run", "demo auto-step"),
         col=c("green", "orange"), lty=1, cex=1)
}
dev.off()
############################# Plot forward run results ########################################################

# load forward simulation results
sol_median_demo<-read.csv(paste0(dataPath,"/sol_median_demo.csv")) # using medians of posterior distributions
sol_likelihood_demo<-read.csv(paste0(dataPath,"/sol_likelihood_demo.csv")) # using max likelihood parameter set
sol_min_J_demo<-read.csv(paste0(dataPath,"/sol_min_J_demo.csv")) # using parameter set with the lowest cost function 
sol_median_formal<-read.csv(paste0(dataPath,"/sol_median_formal.csv")) # using medians of posterior distributions
sol_likelihood_formal<-read.csv(paste0(dataPath,"/sol_likelihood_formal.csv")) # using max likelihood parameter set
sol_min_J_formal<-read.csv(paste0(dataPath,"/sol_min_J_formal.csv")) # using parameter set with the lowest cost function 
sol_median_autostep<-read.csv(paste0(dataPath,"/sol_median_autostep.csv")) # using medians of posterior distributions
sol_likelihood_autostep<-read.csv(paste0(dataPath,"/sol_likelihood_autostep.csv")) # using max likelihood parameter set
sol_min_J_autostep<-read.csv(paste0(dataPath,"/sol_min_J_autostep.csv")) # using parameter set with the lowest cost function 

# transform data
sol_median_demo<-sol_median_demo[,2:length(sol_median_demo)]
sol_likelihood_demo<-sol_likelihood_demo[,2:length(sol_likelihood_demo)]
sol_min_J_demo<-sol_min_J_demo[,2:length(sol_min_J_demo)]
sol_median_formal<-sol_median_formal[,2:length(sol_median_formal)]
sol_likelihood_formal<-sol_likelihood_formal[,2:length(sol_likelihood_formal)]
sol_min_J_formal<-sol_min_J_formal[,2:length(sol_min_J_formal)]
sol_median_autostep<-sol_median_autostep[,2:length(sol_median_autostep)]
sol_likelihood_autostep<-sol_likelihood_autostep[,2:length(sol_likelihood_autostep)]
sol_min_J_autostep<-sol_min_J_autostep[,2:length(sol_min_J_autostep)]


tot_len<-dim(sol_median_formal)[1]

# load observations
obs<-read.csv(paste0(dataPath,"/obs.csv")) 
obs<-obs[2:length(obs)]
names(obs)<-c(    "C_leaf",
                  "C_wood",
                  "C_root",
                  "C_litter_above",
                  "C_litter_below",    
                  "C_fastsom",
                  "C_slowsom",
                  "C_passsom",
                  "rh",
                  "f_veg2litter",
                  'f_litter2som')

names(sol_median_formal)<-names(obs)

# plot as a .jpeg file - MAXIMUM LIKELIHOOD (distribution peak)
jpeg(paste0(dataPath,'/post_da_simu_likelihood', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{ # plot model output with optimized parameters
  par(mfrow=c(3, 4)) # make 3x4 plots in 1 window
  N=tot_len # number of years to plot
  
  for (i in 1:length(names(sol_median_formal))) {
    # # change number of years to plot for very dynamic pools and fluxes (optional)
    #if (names(sol_min_J)[i]=="C_leaf" || names(sol_min_J)[i]=="rh" || names(sol_min_J)[i]=="f_veg2litter"
    #    || names(sol_min_J)[i]=="f_litter2som") {N=tot_len%/%10} else (N=tot_len)
    
    plot(sol_likelihood_formal[[i]][1:N], type="l", col="green", xlab="month",
         ylim=c(min(min(sol_likelihood_formal[[i]]),min(obs[[i]])),max(max(sol_likelihood_formal[[i]]),max(obs[[i]]))), 
         ylab=names(sol_median_formal)[i], main=names(sol_median_formal)[i])
    lines(obs[[i]], col="blue")
    lines(sol_likelihood_autostep[[i]], col="orange")
  }
  plot(0, type = 'n', axes = FALSE, ann = FALSE, main="legend")
  legend(x = "topleft", legend=c("CMIP-6 Output", "Formal run", "Demo-autostep"),
         col=c("blue", "green", "orange"), lty=1, cex=1)
  
  par(mfrow=c(1, 1)) # return to single plot mode
}
dev.off()

# plot as a .jpeg file - MIN_J
jpeg(paste0(dataPath,'/post_da_simu_min_J_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{ # plot model output with optimized parameters
  par(mfrow=c(3, 4)) # make 3x4 plots in 1 window
  N=tot_len # number of years to plot
  
  for (i in 1:length(names(sol_median_formal))) {
    # # change number of years to plot for very dynamic pools and fluxes (optional)
    #if (names(sol_min_J)[i]=="C_leaf" || names(sol_min_J)[i]=="rh" || names(sol_min_J)[i]=="f_veg2litter"
    #    || names(sol_min_J)[i]=="f_litter2som") {N=tot_len%/%10} else (N=tot_len)
    
    plot(sol_min_J_formal[[i]][1:N], type="l", col="green", xlab="month",
         ylim=c(min(min(sol_min_J_formal[[i]]),min(obs[[i]])),max(max(sol_min_J_formal[[i]]),max(obs[[i]]))), 
         ylab=names(sol_median_formal)[i], main=names(sol_median_formal)[i])
    lines(obs[[i]], col="blue")
    lines(sol_min_J_autostep[[i]], col="orange")
  }
  plot(0, type = 'n', axes = FALSE, ann = FALSE, main="legend")
  legend(x = "topleft", legend=c("CMIP-6 Output", "Formal run", "Demo-autostep"),
         col=c("blue", "green", "orange"), lty=1, cex=1)
  
  par(mfrow=c(1, 1)) # return to single plot mode
}
dev.off()

# plot as a .jpeg file - Different solutions of formal mcmc
jpeg(paste0(dataPath,'/post_da_simu_formal_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{ # plot model output with optimized parameters
  par(mfrow=c(3, 4)) # make 3x4 plots in 1 window
  N=tot_len # number of years to plot
  
  for (i in 1:length(names(sol_median_formal))) {
    # # change number of years to plot for very dynamic pools and fluxes (optional)
    #if (names(sol_min_J)[i]=="C_leaf" || names(sol_min_J)[i]=="rh" || names(sol_min_J)[i]=="f_veg2litter"
    #    || names(sol_min_J)[i]=="f_litter2som") {N=tot_len%/%10} else (N=tot_len)
    
    plot(sol_min_J_formal[[i]][1:N], type="l", col="green", xlab="month",
         ylim=c(min(min(sol_min_J_formal[[i]]),min(obs[[i]])),max(max(sol_min_J_formal[[i]]),max(obs[[i]]))), 
         ylab=names(sol_median_formal)[i], main=names(sol_median_formal)[i])
    lines(obs[[i]], col="blue")
    lines(sol_min_J_formal[[i]], col="green")
    lines(sol_likelihood_formal[[i]], col="orange")
  }
  plot(0, type = 'n', axes = FALSE, ann = FALSE, main="legend")
  legend(x = "topleft", legend=c("CMIP-6 Output", "Min Cost", "Max Likelihood"),
         col=c("blue", "green", "orange"), lty=1, cex=1)
  
  par(mfrow=c(1, 1)) # return to single plot mode
}
dev.off()

# plot as a .jpeg file - Different solutions of autostep mcmc
jpeg(paste0(dataPath,'/post_da_simu_autostep_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{ # plot model output with optimized parameters
  par(mfrow=c(3, 4)) # make 3x4 plots in 1 window
  N=tot_len # number of years to plot
  
  for (i in 1:length(names(sol_median_formal))) {
    # # change number of years to plot for very dynamic pools and fluxes (optional)
    #if (names(sol_min_J)[i]=="C_leaf" || names(sol_min_J)[i]=="rh" || names(sol_min_J)[i]=="f_veg2litter"
    #    || names(sol_min_J)[i]=="f_litter2som") {N=tot_len%/%10} else (N=tot_len)
    
    plot(sol_min_J_autostep[[i]][1:N], type="l", col="green", xlab="month",
         ylim=c(min(min(sol_min_J_autostep[[i]]),min(obs[[i]])),max(max(sol_min_J_autostep[[i]]),max(obs[[i]]))), 
         ylab=names(sol_median_formal)[i], main=names(sol_median_formal)[i])
    lines(obs[[i]], col="blue")
    lines(sol_min_J_autostep[[i]], col="green")
    lines(sol_likelihood_autostep[[i]], col="orange")
  }
  plot(0, type = 'n', axes = FALSE, ann = FALSE, main="legend")
  legend(x = "topleft", legend=c("CMIP-6 Output", "Min Cost", "Max Likelihood"),
         col=c("blue", "green", "orange"), lty=1, cex=1)
  
  par(mfrow=c(1, 1)) # return to single plot mode
}
dev.off()

{
print("Minimized Cost Functions:")
print(paste0("Demo: ", round(min(df_j_demo[2]),2)," at iteration ", df_j_demo$V1[which(df_j_demo[2]==min(df_j_demo[2]))] ))
print(paste0("Formal: ", round(min(df_j[2]),2)," at iteration ", df_j$V1[which(df_j[2]==min(df_j[2]))] ))
print(paste0("Autostep: ", round(min(df_j_autostep[2]),2)," at iteration ", df_j_autostep$V1[which(df_j_autostep[2]==min(df_j_autostep[2]))] ))
}
