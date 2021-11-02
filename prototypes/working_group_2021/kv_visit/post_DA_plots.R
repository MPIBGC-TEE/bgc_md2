# identify script location (if using R-studio). If not working: input path manually into setwd()
library(rstudioapi)  
script_location<-normalizePath(rstudioapi::getActiveDocumentContext()$path)
script_location<-substr(script_location, 1, regexpr("post_DA_plots.R",script_location)-1)
setwd(script_location)

# read data path from config.json
library(rjson)
dataPath <- fromJSON(file = "config.json")[[1]]

############################################ Plot posterior distributions of parameters ##########################33

df<-read.csv(paste0(dataPath,"/visit_formal_da_aa.csv"))
df<-t(df)
df<-as.data.frame(df)

par_names<-c("beta_leaf","beta_wood","f_leaflit2fastsom", "f_leaflit2slowsom", "f_leaflit2passsom", "f_woodlit2fastsom", "f_woodlit2slowsom", "f_woodlit2passsom",
             "f_rootlit2fastsom", "f_rootlit2slowsom", "f_rootlit2passsom", "k_leaf",   "k_wood",   "k_root",   "k_leaflit","k_woodlit","k_rootlit","k_fastsom",
             "k_slowsom","k_passsom","C_leaflit_0","T_0","E","KM" )

# calculate median of each parameter distribution
epa_median=rep(0,length(df))
for (i in 1:length(df)) {epa_median[i]=median(df[[i]])}
epa_median<-as.list(epa_median)
# calculate mean of each parameter distribution
epa_mean=rep(0,length(df))
for (i in 1:length(df)) {epa_mean[i]=mean(df[[i]])}
epa_mean<-as.list(epa_mean)
# select parameters from a set with minimal cost function distribution
epa_min_J<-as.list(df[df_j==min(df_j[df_j!=0])])

# visualize parameter distributions with median, mode and min cost function
jpeg('posterior_distr_visit.jpeg', width = 20, height = 10, units = 'in', res = 300)
{
  par(mfrow=c(4, 6)) # make 4x6 plots in 1 window
  for (i in 1:length(df)) {
    hist(df[[i]], breaks=100, col="lightgray", border="lightgray", probability=T, main=par_names[i])
    lines(stats::density(df[[i]]), col="black", lwd=2)
    abline(v=epa_median[[i]],col="red",lwd=2)
    abline(v=epa_mean[[i]],col="orange",lwd=2)
    #abline(v=epa_min_J[[i]],col="green",lwd=2)
  }
}
dev.off()
########################## Plot cost function evolution #########################################################3333

df_j<-read.csv(paste0(dataPath,"/visit_formal_da_j_aa.csv"))
df_j_demo<-read.csv(paste0(dataPath,"/visit_demo_da_j_aa.csv"))
df_j<-df_j[df_j!=0]
df_j_demo<-df_j_demo[df_j_demo!=0]
df_j<-rbind(c(1:length(df_j)),df_j)
df_j_demo<-rbind(c(1:length(df_j_demo)),df_j_demo)

df_j<-t(df_j)
df_j<-as.data.frame(df_j)
df_j_demo<-t(df_j_demo)
df_j_demo<-as.data.frame(df_j_demo)

jpeg('cost_functions_visit.jpeg', width = 20, height = 10, units = 'in', res = 300)
{
plot(df_j, type="l", col="blue", xlab="accepted cost functions", ylab="cost value", main="Cost functions over time",
     xlim=c(0,max(dim(df_j)[1],dim(df_j_demo)[1])), 
     ylim=c(min(min(df_j),min(df_j_demo)),max(max(df_j),max(df_j_demo)))) 
lines(df_j_demo, col="red")
legend(max(dim(df_j)[1],dim(df_j_demo)[1])*0.7, max(max(df_j),max(df_j_demo))*0.9, legend=c("demo run", "formal run"),
       col=c("red", "blue"), lty=1, cex=1)
}
dev.off()
############################# Plot forward run results ########################################################
sol_mean<-read.csv(paste0(dataPath,"/sol_mean.csv"))
sol_median<-read.csv(paste0(dataPath,"/sol_median.csv"))
sol_min_J<-read.csv(paste0(dataPath,"/sol_min_J.csv"))

sol_mean<-sol_mean[,2:length(sol_mean)]
sol_median<-sol_median[,2:length(sol_median)]
sol_min_J<-sol_min_J[,2:length(sol_min_J)]

tot_len<-dim(sol_mean)[1]

dat<-read.csv(paste0(dataPath,"/dat.csv")) 
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
names(sol_median)<-c(    "C_leaf",
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
jpeg('post_da_simu_visit.jpeg', width = 20, height = 10, units = 'in', res = 300)
{ # plot model output with optimized parameters
  par(mfrow=c(3, 4)) # make 3x4 plots in 1 window
  N=tot_len # number of years to plot
  
  # change number of years to plot for very dynamic pools and fluxes (optional)

  for (i in 1:length(names(sol_median))) {
    #if (names(sol_median)[i]=="C_leaf" || names(sol_median)[i]=="rh" || names(sol_median)[i]=="f_veg2litter"
    #    || names(sol_median)[i]=="f_litter2som") {N=tot_len%/%10} else (N=tot_len)
    
    plot(sol_median[[i]][1:N], type="l", col="red", xlab="month",
         ylim=c(min(min(sol_median[[i]]),min(obs[[i]])),max(max(sol_median[[i]]),max(obs[[i]]))), 
         ylab=names(sol_median)[i], main=names(sol_median)[i])
    lines(obs[[i]], col="blue")
    lines(sol_mean[[i]], col="orange")
    lines(sol_min_J[[i]], col="green")
  }
  plot(2, xlim=c(0,1), ylim=c(0,0.9), axes = F, main="legend", ylab="")
  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Median", "Mean", "Min Cost Function"),
         col=c("blue", "red", "orange", "green"), lty=1, cex=1)
  
  par(mfrow=c(1, 1)) # return to single plot mode
}
dev.off()
