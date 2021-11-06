# # identify script location (if using R-studio). If not working: input path manually into setwd()
# library(rstudioapi)  
# script_location<-normalizePath(rstudioapi::getActiveDocumentContext()$path)
# script_location<-substr(script_location, 1, regexpr("post_DA_plots.R",script_location)-1)
# setwd(script_location)
model_name="VISIT"
# read data path from config.json
library(rjson)
dataPath <- fromJSON(file = "config.json")[[1]]
time_stamp<-substr(Sys.time(),1,regexpr(":",Sys.time())-1)

############################################ Plot posterior (demo) distributions of parameters ##########################33

# load posterior distributions after the demo run

df<-read.csv(paste0(dataPath,"/visit_demo_da_aa.csv"))
df<-df[,2:dim(df)[2]]
df<-t(df)

df<-as.data.frame(df)

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
epa_min_J<-as.list(df[df_j==min(df_j),])

# plot as a .jpeg file
jpeg(paste0(dataPath,'/posterior_distr_demo_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
  par(mfrow=c(4, 6)) # make 4x6 plots in 1 window
  for (i in 1:length(df)) {
    hist(df[[i]], breaks=100, col="lightgray", border="lightgray", probability=T, main=par_names[i])
    lines(stats::density(df[[i]]), col="black", lwd=2)
    abline(v=epa_min_J[[i]],col="green",lwd=3)
    abline(v=epa_median[[i]],col="red",lwd=2)
    #abline(v=epa_mean[[i]],col="orange",lwd=2)

    #abline(v=epa_0[[i]],col="purple",lwd=2)
  }
  #  plot(2, xlim=c(0,1), ylim=c(0,0.9), axes = F, main="legend", ylab="")
  #  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Median", "Min Cost Function", "Initial parameters"),
  #         col=c("blue", "red", "green", "purple"), lty=1, bty="n")
  
  par(mfrow=c(1, 1))
}
dev.off()

############################################ Plot posterior (formal) distributions of parameters ##########################33

# load posterior distributions after the formal run
df<-read.csv(paste0(dataPath,"/visit_formal_da_aa.csv"))
df<-df[,2:dim(df)[2]]
df<-t(df)

df<-as.data.frame(df)

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
epa_min_J<-as.list(df[df_j==min(df_j),])

# plot as a .jpeg file
jpeg(paste0(dataPath,'/posterior_distr_formal_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
  par(mfrow=c(4, 6)) # make 4x6 plots in 1 window
  for (i in 1:length(df)) {
    hist(df[[i]], breaks=100, col="lightgray", border="lightgray", probability=T, main=par_names[i])
    lines(stats::density(df[[i]]), col="black", lwd=2)
    abline(v=epa_median[[i]],col="red",lwd=2)
    #abline(v=epa_mean[[i]],col="orange",lwd=2)
    abline(v=epa_min_J[[i]],col="green",lwd=3)
    #abline(v=epa_0[[i]],col="purple",lwd=2)
  }
#  plot(2, xlim=c(0,1), ylim=c(0,0.9), axes = F, main="legend", ylab="")
#  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Median", "Min Cost Function", "Initial parameters"),
#         col=c("blue", "red", "green", "purple"), lty=1, bty="n")
  
  par(mfrow=c(1, 1))
}
dev.off()

############################################ Plot posterior (autostep) distributions of parameters ##########################33

# load posterior distributions after the demo run

df<-read.csv(paste0(dataPath,"/visit_autostep_da_aa.csv"))
df<-df[,2:dim(df)[2]]
df<-t(df)

df<-as.data.frame(df)

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
epa_min_J<-as.list(df[df_j==min(df_j),])

# plot as a .jpeg file
jpeg(paste0(dataPath,'/posterior_distr_autostep_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
  par(mfrow=c(4, 6)) # make 4x6 plots in 1 window
  for (i in 1:length(df)) {
    hist(df[[i]], breaks=100, col="lightgray", border="lightgray", probability=T, main=par_names[i])
    lines(stats::density(df[[i]]), col="black", lwd=2)
    abline(v=epa_min_J[[i]],col="green",lwd=3)
    abline(v=epa_median[[i]],col="red",lwd=2)
    #abline(v=epa_mean[[i]],col="orange",lwd=2)
    
    #abline(v=epa_0[[i]],col="purple",lwd=2)
  }
  #  plot(2, xlim=c(0,1), ylim=c(0,0.9), axes = F, main="legend", ylab="")
  #  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Median", "Min Cost Function", "Initial parameters"),
  #         col=c("blue", "red", "green", "purple"), lty=1, bty="n")
  
  par(mfrow=c(1, 1))
}
dev.off()

########################## Plot cost function evolution ############################################################

# load cost function data from demo and formal run
df_j<-read.csv(paste0(dataPath,"/visit_formal_da_j_aa.csv"))
df_j_demo<-read.csv(paste0(dataPath,"/visit_demo_da_j_aa.csv"))
df_j_autostep<-read.csv(paste0(dataPath,"/visit_autostep_da_j_aa.csv"))

# transform the data for plotting
df_j<-df_j[df_j!=0]
df_j_demo<-df_j_demo[df_j_demo!=0]
df_j_autostep<-df_j_autostep[df_j_autostep!=0]
df_j<-rbind(c(1:length(df_j)),df_j)
df_j_demo<-rbind(c(1:length(df_j_demo)),df_j_demo)
df_j_autostep<-rbind(c(1:length(df_j_autostep)),df_j_autostep)
df_j<-t(df_j)
df_j<-as.data.frame(df_j)
df_j_demo<-t(df_j_demo)
df_j_demo<-as.data.frame(df_j_demo)
df_j_autostep<-t(df_j_autostep)
df_j_autostep<-as.data.frame(df_j_autostep)


# plot as a .jpeg file
jpeg(paste0(dataPath,'/cost_functions_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
plot(df_j, type="l", col="green", xlab="accepted cost functions", ylab="cost value", main="Cost functions over time",
     xlim=c(0,max(dim(df_j)[1],dim(df_j_demo)[1])), 
     ylim=c(min(min(df_j),min(df_j_demo)),max(max(df_j[,2]),max(df_j_demo[,2])))) 
lines(df_j_demo, col="red")
lines(df_j_autostep, col="orange")

legend(max(dim(df_j)[1],dim(df_j_demo)[1])*0.7, max(max(df_j),max(df_j_demo))*0.9, legend=c("demo run", "formal run", "demo auto-step"),
       col=c("red", "green", "orange"), lty=1, cex=1)
}
dev.off()

# Same plot on a log scale
jpeg(paste0(dataPath,'/cost_functions_log_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
  plot(df_j, type="l", col="green", xlab="accepted cost functions", ylab="log cost value", main="Cost functions over time (log scale)",
       xlim=c(0,max(dim(df_j)[1],dim(df_j_demo)[1])), 
       ylim=c(min(min(df_j),min(df_j_demo)),max(max(df_j[,2]),max(df_j_demo[,2]))),
       log="y") 
  lines(df_j_demo, col="red")
  lines(df_j_autostep, col="orange")
  legend(max(dim(df_j)[1],dim(df_j_demo)[1])*0.7, max(max(df_j),max(df_j_demo))*0.9, legend=c("demo run", "formal run", "demo auto-step"),
         col=c("red", "green", "orange"), lty=1, cex=1)
}
dev.off()

# Same plot with only last 10% accepted cost functions

df_j1<-df_j[df_j$V1>length(df_j[[1]])*0.1,]
df_j1$V1<-1:length(df_j1$V1)
df_j_demo1<-df_j_demo[df_j_demo$V1>length(df_j_demo[[1]])*0.1,]
df_j_demo1$V1<-1:length(df_j_demo1$V1)
df_j_autostep1<-df_j_autostep[df_j_autostep$V1>length(df_j_autostep[[1]])*0.1,]
df_j_autostep1$V1<-1:length(df_j_autostep1$V1)

jpeg(paste0(dataPath,'/cost_functions_last_10_percent_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{
  plot(df_j1, type="l", col="green", xlab="accepted cost functions", ylab="cost value", main="Last 90% of accepted cost functions",
       xlim=c(0,max(dim(df_j1)[1],dim(df_j_demo1)[1])), 
       ylim=c(min(min(df_j1),min(df_j_demo1)),max(max(df_j1[,2]),max(df_j_demo1[,2])))) 
  lines(df_j_demo1, col="red")
  lines(df_j_autostep1, col="orange")
  legend(max(dim(df_j1)[1],dim(df_j_demo1)[1])*0.7, max(max(df_j1),max(df_j_demo1))*0.9, legend=c("demo run", "formal run", "demo auto-step"),
         col=c("red", "green", "orange"), lty=1, cex=1)
}
dev.off()
############################# Plot forward run results ########################################################

# load forward simulation results
#sol_init_par<-read.csv(paste0(dataPath,"/sol_init_par.csv")) # using initial parameters
#sol_mean<-read.csv(paste0(dataPath,"/sol_mean.csv")) # using means of posterior distributions
sol_median_demo<-read.csv(paste0(dataPath,"/sol_median_demo.csv")) # using medians of posterior distributions
sol_min_J_demo<-read.csv(paste0(dataPath,"/sol_min_J_demo.csv")) # using parameter set with the lowest cost function 
sol_median_formal<-read.csv(paste0(dataPath,"/sol_median_formal.csv")) # using medians of posterior distributions
sol_min_J_formal<-read.csv(paste0(dataPath,"/sol_min_J_formal.csv")) # using parameter set with the lowest cost function 
sol_median_autostep<-read.csv(paste0(dataPath,"/sol_median_autostep.csv")) # using medians of posterior distributions
sol_min_J_autostep<-read.csv(paste0(dataPath,"/sol_min_J_autostep.csv")) # using parameter set with the lowest cost function 

# transform data
#sol_init_par<-sol_init_par[,2:length(sol_init_par)]
#sol_mean<-sol_mean[,2:length(sol_mean)]
sol_median_demo<-sol_median_demo[,2:length(sol_median_demo)]
sol_min_J_demo<-sol_min_J_demo[,2:length(sol_min_J_demo)]
sol_median_formal<-sol_median_formal[,2:length(sol_median_formal)]
sol_min_J_formal<-sol_min_J_formal[,2:length(sol_min_J_formal)]
sol_median_autostep<-sol_median_autostep[,2:length(sol_median_autostep)]
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

names(sol_median)<-names(obs)

# plot as a .jpeg file - MEDIANS
jpeg(paste0(dataPath,'/post_da_simu_median_', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{ # plot model output with optimized parameters
  par(mfrow=c(3, 4)) # make 3x4 plots in 1 window
  N=tot_len # number of years to plot
  
  for (i in 1:length(names(sol_median))) {
    # # change number of years to plot for very dynamic pools and fluxes (optional)
    #if (names(sol_median)[i]=="C_leaf" || names(sol_median)[i]=="rh" || names(sol_median)[i]=="f_veg2litter"
    #    || names(sol_median)[i]=="f_litter2som") {N=tot_len%/%10} else (N=tot_len)
    
    plot(sol_median_demo[[i]][1:N], type="l", col="red", xlab="month",
         ylim=c(min(min(sol_median_formal[[i]]),min(obs[[i]])),max(max(sol_median_formal[[i]]),max(obs[[i]]))), 
         ylab=names(sol_median)[i], main=names(sol_median)[i])
    lines(obs[[i]], col="blue")
    #lines(sol_mean[[i]], col="orange")
    lines(sol_median_formal[[i]], col="green")
    lines(sol_median_autostep[[i]], col="orange")
  }
  plot(2, xlim=c(0,1), ylim=c(0,0.9), axes = F, main="legend", ylab="")
  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Demo run", "Formal run", "Demo-autostep"),
         col=c("blue", "red", "green", "orange"), lty=1, cex=1)
  
  par(mfrow=c(1, 1)) # return to single plot mode
}
dev.off()

# plot as a .jpeg file - MIN_J
jpeg(paste0(dataPath,'/post_da_simu_min_J', model_name,'_',time_stamp,'.jpeg'), width = 15, height = 8, units = 'in', res = 300)
{ # plot model output with optimized parameters
  par(mfrow=c(3, 4)) # make 3x4 plots in 1 window
  N=tot_len # number of years to plot
  
  for (i in 1:length(names(sol_min_J))) {
    # # change number of years to plot for very dynamic pools and fluxes (optional)
    #if (names(sol_min_J)[i]=="C_leaf" || names(sol_min_J)[i]=="rh" || names(sol_min_J)[i]=="f_veg2litter"
    #    || names(sol_min_J)[i]=="f_litter2som") {N=tot_len%/%10} else (N=tot_len)
    
    plot(sol_min_J_demo[[i]][1:N], type="l", col="red", xlab="month",
         ylim=c(min(min(sol_min_J_formal[[i]]),min(obs[[i]])),max(max(sol_min_J_formal[[i]]),max(obs[[i]]))), 
         ylab=names(sol_min_J)[i], main=names(sol_min_J)[i])
    lines(obs[[i]], col="blue")
    #lines(sol_mean[[i]], col="orange")
    lines(sol_min_J_formal[[i]], col="green")
    lines(sol_min_J_autostep[[i]], col="orange")
  }
  plot(2, xlim=c(0,1), ylim=c(0,0.9), axes = F, main="legend", ylab="")
  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Demo run", "Formal run", "Demo-autostep"),
         col=c("blue", "red", "green", "orange"), lty=1, cex=1)
  
  par(mfrow=c(1, 1)) # return to single plot mode
}
dev.off()
