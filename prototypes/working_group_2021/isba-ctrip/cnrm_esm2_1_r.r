################################# Script for running draft CNRM_ESM2_1 matrix function #################################


library (raster)
library(ncdf4)
library(sp)
library(gridExtra)
library(ggplot2)
setwd("C:/Users/qbi/Desktop/Modeltraining/CMIP6/ISBA_CTRIP/output/pct1CO2bgc")

######################################## Import CMIP6 outputs ######################################################

# pick up 1 site: 33.3 E, 50.0 N
point<-SpatialPoints(cbind(33.3,50.0), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))

# load all netCDF files in the folder
files<-list.files(pattern="..nc") # list all netCDF files in the folder
(var_names<-substr(files, 0 ,regexpr("_",files)-1)) # derive variable names from file names
dat<-data.frame(point)
for (i in 1:length(var_names)){
  r<-stack(files[i])
  r_dat<-extract(r,point)
  dat<-cbind(dat,r_dat[1,])
}
(names(dat)<-c("lon","lat", var_names)) # assign variable names


# correct fluxes from per second to per day
dat$gpp<-dat$gpp*86400
dat$npp<-dat$npp*86400
dat$fVegLitter<-dat$fVegLitter*86400
dat$fLitterSoil<-dat$fLitterSoil*86400
dat$rh<-dat$rh*86400

# explore the data
summary(dat)

#save and load the data
write.csv(dat,"dat.csv")
dat<-read.csv("dat.csv")

########################################## Define matrix simulation function ########################################

# define number of days, years and months
days = c(31,28,31,30,31,30,31,31,30,31,30,31)
nyears = 150 
tot_len = nyears*12

########################################## Define matrix simulation function ########################################

# define number of days, years and months
days = c(31,28,31,30,31,30,31,31,30,31,30,31)
nyears = 150 
tot_len = nyears*12
clay=0.4  
silt=0.2 #??

## Define matrix simulation function
matrix_simu<-function(pa) {
  # B vector
  beta1=pa[1]
  B = c(beta1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)   # allocation
  # A matrix
  f21 = pa[2]; f32 = pa[3];f41 = pa[4];f42 = pa[5];f43 = pa[6];f53 = pa[7];f64 = pa[8];
  lig_aboveleaf= pa[9]; #lignin fraction of above-ground structural litter
  lig_belowleaf = pa[10];#lignin fraction of below-ground structural litter
  f71 = 0.3084; f72 = 0.402; f73 = 0.402; f75 = 0.402; 
  f81 = 0.6916; f82 = 0.598; f83 = 0.598; f85 = 0.598;f94 = 0.402;  
  f96 = 0.402; f104 = 0.598; f106 = 0.598; f117 = 0.55*(1-lig_aboveleaf); f118 = 0.45;
  f119 = 0.45*(1-lig_belowleaf); f1110 = 0.45; f1112 = 0.42; f1113 = 0.45;  
  f127 = 0.7*lig_aboveleaf; f129 = 0.7*lig_belowleaf; f1211 = 1 - 0.004 - (0.85- 0.68*(clay+silt));
  f1311 = 0.004; f1312 = 0.03
  A = c(-1,  0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        f21, -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,  f32,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
        f41,f42, f43,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
        0,   0,  f53,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0, 
        0,   0,   0,  f64,   0,  -1,   0,   0,   0,   0,   0,   0,   0, 
        f71, f72, f73,  0, f75,   0,  -1,   0,   0,   0,   0,   0,   0,
        f81, f82, f83,  0, f85,   0,   0,  -1,   0,   0,   0,   0,   0,
        0,   0,   0,  f94,   0,  f96,  0,   0,  -1,   0,   0,   0,   0,
        0,   0,   0,  f104,  0,  f106, 0,   0,   0,   -1,  0,   0,   0,
        0,   0,   0,   0,   0,   0, f117, f118,f119,f1110,-1,f1112,f1113,
        0,   0,   0,   0,   0,   0, f127,  0,  f129, 0, f1211, -1,  0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0, f1311,f1312,-1)
  A = matrix(A, nrow = 13, byrow = TRUE)
  
  # K matrix 
  temp = c(pa[11],pa[12],pa[13], pa[14],pa[15], pa[16], 1/0.245*exp(-3*pa[9]), 1/0.066, 1/0.245*exp(-3*pa[10]), 1/0.066, 1/0.149*(1 - 0.75*(clay+silt)), 1/5.37, 1/241) 
  K = rep(0,169)
  K = matrix(K, nrow = 13, byrow = TRUE) 
  for (i in 1:13) { K[i,i] = temp[i] }
  # X vector
  x=rep(0,tot_len)
  x_fin=data.frame(x,x,x,x,x,x,x,x,x,x,x,x,x); names(x_fin)=c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13")
  rh_fin=rep(0,tot_len)
  # Initial carbon pool size derived from 1st year outputs where possible
  x_init = c(dat$cLeaf[1],  #leaf biomass - unknown
             pa[17], #active structural biomass, linked to leaf biomass through nitrogen dilution  - unknown
             pa[18], #passive structural biomass - unknown
             pa[19], #below ground structural biomass  - unknown
             dat$cWood[1] + dat$cStem[1], #above ground woody biomass (for trees) - unknown
             dat$cRoot[1]-pa[19], #below ground woody biomass (for trees)  - unknown
             pa[20], #Aboveground structure litter - unknown
             pa[21], #Aboveground metabolic litter - unknown
             pa[22], #Belowground structure litter - unknown
             dat$cLitterBelow[1] - pa[22], #Belowground metabolic litter  - unknown
             dat$cSoilFast[1], # soil fast
             dat$cSoilMedium[1], # soil medium
             dat$cSoilSlow[1])  # soil slow  
  X=x_init   # initialize carbon pools
  
  jj=1
  for (y in 1:nyears){
    for (m in 1:12){
      npp_in = dat$npp[jj] 
      co2_rh = 0; f_veg_lit = 0; f_lit_soil = 0  
      # soil moisture factor ksi: monthly soil moisture data "mrso" multiplied by sensitivity parameters for different pools
      ksi=c(1,1,1,1,1,1,pa[23]*dat$mrso[jj],pa[24]*dat$mrso[jj],pa[25]*dat$mrso[jj],pa[26]*dat$mrso[jj],pa[27]*dat$mrso[jj],pa[28]*dat$mrso[jj],pa[29]*dat$mrso[jj])
      
      for (d in 1:days[m]) {
        # matrix equation with ksi
        X = X + B * npp_in + A %*% K %*% X * ksi
        # deriving rh from each litter and soil pool as 1 - sum of transfer coefficients 
        co2_rate = c(0,0,0,0,0,0,(1-f117-f127)*K[7,7],(1-f118)*K[8,8],(1-f119-f129)*K[9,9],(1-f1110)*K[10,10],(1-f1211-f1311)*K[11,11],(1-f1112-f1312)*K[12,12],(1-f1113)*K[13,13])
        co2=sum(co2_rate*X)
        co2_rh = co2_rh + co2/days[m]   # monthly average rh
      }
      x_fin[jj,]=X
      rh_fin[jj]=co2_rh
      jj= jj+1
    }
  } 
  # outputs: C pools, heterotrophic respiration, litterfall, humus formation
  result<-list(x_fin, rh_fin)
  return(result)
}

################################################### Run matrix model ################################################
paramNames=c("beta1","f21","f32","f41","f42","f43","f53","f64",
             "lig_aboveleaf","lig_belowleaf", 
             "k1", "k2", "k3", "k4", "k5", "k6", 
             "x2_init", "x3_init","x4_init","x7_init","x8_init","x9_init",
             "sens_moist_x7","sens_moist_x8","sens_moist_x9","sens_moist_x10","sens_moist_x11","sens_moist_x12","sens_moist_x13")
# example parameters for running the model
pa<-c(
  0.4, # beta 1 - allocation to leaf 
  0.2, #f21
  0.05, # f32
  0.2, # f41 
  0.2, #f
  0.2, # f
  0.2, #f
  0.2,  # f
  0.6916, # LX1 - lignin fraction of X1
  0.3084, # LX3 - lignin fraction of X3
  1/365, # k1 - turnover rate of leaf biomass 
  1/(365*10), # k2 - turnover rate of active structural biomass
  1/(365*5), # k3 - turnover rate of passive structural biomass
  1/(365*0.5), # k4 - turnover rate of below ground structural biomass 
  1/(365*3), # k5 - turnover rate of above ground woody biomass 
  1/(365*2), # k6 - turnover rate of below ground woody biomass 
  0.003, # X2_init - initial size of active structural biomass
  0.001, # X3_init - initial size of passive structural biomass
  0.005, # X4_init - initial size of below ground structural biomass
  0.01, # X7_init - initial size of Aboveground structure litter 
  0,01,
  0.01,# X9_init - initial size of Belowground structure litter
  0.0003, # sensitivity to moisture of aboveground structural litter
  0.0003, # sensitivity to moisture of aboveground metabolic litter
  0.0003, # sensitivity to moisture of belowground structural litter
  0.0003, # sensitivity to moisture of belowground metabolic litter
  0.003, # sensitivity to moisture of fast soil
  0.003, # sensitivity to moisture of medium soil
  0.003  # sensitivity to moisture of slow soil
)

# test-run the model
test<-matrix_simu(pa)

# view results
View(test[[1]]) # monthly pool sizes
View(test[[2]]) # monthly rh
#View(test[[3]]) # monthly litterfall
#View(test[[4]]) # monthly humus formation

# Compare modeled results with CMIP6 output for 1st 10 years (blue = modeled, red = CMIP6 output)
par(mfrow=c(3, 4)) 
{
  plot(test[[1]]$x1[1:240], type="l", col="red",  xlab="month", 
       ylim=c(min(min(test[[1]]$x1[1:240]),min(dat$cLeaf[1:240])),max(max(test[[1]]$x1[1:240]),max(dat$cLeaf[1:240]))), 
       xaxp  = c(0, 120, 10), ylab="cLeaf", main="Leaf Pool modelled vs CMIP6")
  lines(dat$cLeaf[1:240], col="blue")
  
  plot(test[[1]]$x4[1:240]+test[[1]]$x6[1:240], type="l", col="red", xlab="month", 
       ylim=c(min(min(test[[1]]$x4[1:240]+test[[1]]$x6[1:240]),min(dat$cRoot[1:240])),max(max(test[[1]]$x4[1:240]+test[[1]]$x6[1:240]),max(dat$cRoot[1:240]))), 
       xaxp  = c(0, 120, 10),ylab="cRoot", main="Root Pool modelled vs CMIP6")
  lines(dat$cRoot[1:240], col="blue")
  
  plot(test[[1]]$x7[1:240]+test[[1]]$x8[1:240], type="l", col="red", 
       ylim=c(min(min(test[[1]]$x7[1:240]+test[[1]]$x8[1:240]),min(dat$cLitterSurf[1:240])),max(max(test[[1]]$x4[1:240]+test[[1]]$x5[1:240]),max(dat$cLitterSurf[1:240]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cLitterAbove", main="Above-ground Litter (leaf+wood) modelled vs CMIP6")
  lines(dat$cLitterSurf[1:240], col="blue")
  
  plot(test[[1]]$x9[1:240]+test[[1]]$x10[1:240], type="l", col="red", 
       ylim=c(min(min(test[[1]]$x9[1:240]+test[[1]]$x10[1:240]),min(dat$cLitterBelow[1:240])),max(max(test[[1]]$x9[1:240]+test[[1]]$x10[1:240]),max(dat$cLitterBelow[1:240]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cLitterBelow", main="Below-ground Litter (root) modelled vs CMIP6")
  lines(dat$cLitterBelow[1:240], col="blue")
  
  plot(test[[1]]$x11[1:240], type="l", col="red", 
       ylim=c(min(min(test[[1]]$x11[1:240]),min(dat$cSoilFast[1:240])),max(max(test[[1]]$x11[1:240]),max(dat$cSoilFast[1:240]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilFast", main="Fast SOM pool modelled vs CMIP6")
  lines(dat$cSoilFast[1:240], col="blue")
  
  plot(test[[1]]$x12[1:240], type="l", col="red",  
       ylim=c(min(min(test[[1]]$x12[1:240]),min(dat$cSoilMedium[1:240])),max(max(test[[1]]$x12[1:240]),max(dat$cSoilMedium[1:240]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilMEdium", main="Medium, SOM pool modelled vs CMIP6")
  lines(dat$cSoilMedium[1:240], col="blue")
  
  plot(test[[1]]$x13[1:240], type="l", col="red",  
       ylim=c(min(min(test[[1]]$x13[1:240]),min(dat$cSoilSlow[1:240])),max(max(test[[1]]$x13[1:240]),max(dat$cSoilSlow[1:240]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilSlow", main="Slow SOM pool modelled vs CMIP6")
  lines(dat$cSoilSlow[1:240], col="blue")
  
  plot(test[[2]][1:240], type="l", col="red",  
       ylim=c(min(min(test[[2]][1:240]),min(dat$rh[1:240])),max(max(test[[2]][1:240]),max(dat$rh[1:240]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="Rh", main="Heterotrophic respiration modelled vs CMIP6")
  lines(dat$rh[1:240], col="blue")

  plot(2, axes = F,xlim=c(0,1), ylim=c(0,1),  main="legend", ylab="")
  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Modelled"),
         col=c("blue", "red"), lty=1, cex=1)
}

#error happens about ylim in the above code, so I used the previous one

{plot(test[[1]]$x1[1:120], type="l", col="blue", ylim=c(0,0.2), xlab="month", ylab="cLeaf", main="Leaf Pool modelled vs CMIP6")
  lines(dat$cLeaf[1:120], col="red")}

{plot(test[[1]]$x4[1:240]+test[[1]]$x6[1:240], type="l", col="blue", ylim=c(0,1), xlab="month", ylab="cRoot", main="Root Pool modelled vs CMIP6")
  lines(dat$cRoot[1:120], col="red")}

{plot(test[[1]]$x7[1:240]+test[[1]]$x8[1:240], type="l", col="blue", ylim=c(0,2), xlab="month", ylab="cLitterAbove", main="Above-ground Litter modelled vs CMIP6")
  lines(dat$cLitterSurf[1:120], col="red")}

{plot(test[[1]]$x9[1:240]+test[[1]]$x10[1:240], type="l", col="blue", ylim=c(0,1), xlab="month", ylab="cLitterBelow", main="Below-ground Litter (root) modelled vs CMIP6")
  lines(dat$cLitterBelow[1:120], col="red")}

{plot(test[[1]]$x11[1:120], type="l", col="blue", ylim=c(0,1), xlab="month", ylab="cSoilFast", main="Fast SOM pool modelled vs CMIP6")
  lines(dat$cSoilFast[1:120], col="red")}

{plot(test[[1]]$x12[1:120], type="l", col="blue", ylim=c(0,10), xlab="month", ylab="cSoilMEdium", main="Medium, SOM pool modelled vs CMIP6")
  lines(dat$cSoilMedium[1:120], col="red")}

{plot(test[[1]]$x13[1:120], type="l", col="blue", ylim=c(0,10), xlab="month", ylab="cSoilSlow", main="Slow SOM pool modelled vs CMIP6")
  lines(dat$cSoilSlow[1:120], col="red")}

{plot(test[[2]][1:120], type="l", col="blue", ylim=c(0,0.002), xlab="month", ylab="Rh", main="Heterotrophic respiration modelled vs CMIP6")
  lines(dat$rh[1:120], col="red")}

####################################################### MCMC##############################################################
# assign min-max values to all parameters

c_min=c(0.09, 0.09, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2,   0.01, 0.01,    1/(2*365), 1/(365*100), 1/(365*100), 0.1/(365*0.5),  0.03/(365*3),  0.03/(365*2),  dat$cLeaf[1]/100, 0.0001, dat$croot[1]/100,dat$cLitterAbove[1]/100, dat$cLitterAbove[1]/100, dat$cLitterBelow[1]/100,   0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,0.0001)

pa=c(0.4,  0.2, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2,   0.6916, 0.3084,    1/365, 1/(365*10), 1/(365*5), 0.5/(365*0.5),  0.3/(365*3),  0.3/(365*2),  0.003, 0.001, 0.005,    0.01,  0,01, 0.01, 0.0003, 0.0003, 0.0003, 0.0003, 0.003, 0.003,0,003 )

c_max=c(1,  1, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2,    0.9, 0.6,     1/(0.3*365), 1/(365*0.8), 1/(365*0.8), 1/(365*0.5),  1/(365*3),  0.6/(365*2),  dat$cLeaf[0], 0.1, dat$cRoot[0],    dat$cLitterAbove[0],    dat$cLitterAbove[0],  dat$cLitterBelow[0], 0.003,  0.003,  0.03,  0.03,  0.03, 0.03,0,03)

# set random numbers table to get reproducible results
set.seed(8192021)

paramNum=length(pa)
# number of simulations
nsimu    = 2000

upgraded=0
J_last = 400
C_op = pa

C_upgraded = rep(0,paramNum*nsimu)
C_upgraded = matrix(C_upgraded, nrow = paramNum, byrow = TRUE)
J_upgraded = rep(0, nsimu)

isQualified=function(c) {
  flag = T
  for (i in 1:paramNum){
    # check for min-max values
    if (c[i] > c_max[i] || c[i] < c_min[i]){
      flag = F
      break}
    # check for total outgoing flux fraction for each pool not exceeding 1
    if(c[1] + c[2] > 1){
      flag = F
      break}
    if(c[3] + c[8] > 1){
      flag = F
      break}
    if(c[4] + c[9] > 1){
      flag = F
      break}
    if(c[5] + c[10] > 1){
      flag = F
      break}
    if(c[6] + c[13] > 1){
      flag = F
      break}
    if(c[7] + c[11] > 1){
      flag = F
      break}
  }
  return (flag)
}

GenerateParamValues=function(c_op){
  flag = T
  while (flag){
    c_new = c_op + (runif((paramNum)) - 0.5)*(c_max - c_min)/15.0
    if (isQualified(c_new)){  flag = F }
  }
  return (c_new)
}

for (simu in 1:nsimu){
  print (paste0("simulation ",simu))
  c_new = GenerateParamValues(C_op)
  x_simu = matrix_simu(c_new)[[1]]
  rh_simu = matrix_simu(c_new)[[2]]
  litterfall_simu = matrix_simu(c_new)[[3]]
  humification_simu = matrix_simu(c_new)[[4]]
  # cost functions
  J_obj1 = mean (( x_simu$x1 - dat$cLeaf )**2)/(2*var(dat$cLeaf))
  J_obj2 = mean (( x_simu$x7+x_simu$x8 - dat$cLitterSurf)**2)/(2*var(dat$cLitterSurf))
  J_obj3 = mean (( x_simu$x9+x_simu$x10 - dat$cLitterBelow)**2)/(2*var(dat$cLitterBelow))
  J_obj4 = mean (( x_simu$x11 - dat$cSoilFast )**2)/(2*var(dat$cSoilFast))
  J_obj5 = mean (( x_simu$x12 - dat$cSoilMedium )**2)/(2*var(dat$cSoilMedium))
  J_obj6 = mean (( x_simu$x13 - dat$cSoilSlow )**2)/(2*var(dat$cSoilSlow))
  J_obj7 = mean (( rh_simu - dat$rh )**2)/(2*var(dat$rh))
  # total cost function
  J_new= (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 + J_obj6 + J_obj7)/200
  
  delta_J =  J_last - J_new;
  
  randNum = runif(1)
  if (min(1.0, exp(delta_J)) > randNum) {
    C_op=c_new;
    J_last=J_new;
    C_upgraded[,upgraded]=C_op;
    J_upgraded[upgraded]=J_last; 
    upgraded=upgraded+1 }
}

# results: distributions of parameters and cost functions
df=data.frame(C_upgraded)
df_j=data.frame(J_upgraded)

write.csv(df,'visit_demo_da_aa.csv')
write.csv(df_j,'visit_demo_da_j_aa.csv')
################################################### Run matrix model ################################################
paramNames=c("beta1","f21","f32","f41","f42","f43","f53","f64",
             "lig_aboveleaf","lig_belowleaf", 
             "k1", "k2", "k3", "k4", "k5", "k6", 
             "x2_init", "x3_init","x4_init","x7_init","x8_init","x9_init",
             "sens_moist_x7","sens_moist_x8","sens_moist_x9","sens_moist_x10","sens_moist_x11","sens_moist_x12","sens_moist_x13")

# plot distribution of optimized parameters
distr<-t(df) #transpose optimized parameter distributions
distr<-as.data.frame(distr)
names(distr)=paramNames
{ # create density plots of parameter distributions using ggplot2
  p1=ggplot(distr, aes(x=beta1)) + geom_density() + geom_vline(aes(xintercept=median(beta1)),color="blue", linetype="dashed", size=1)
  p2=ggplot(distr, aes(x=f21)) + geom_density() + geom_vline(aes(xintercept=median(f21)),color="blue", linetype="dashed", size=1)
  p3=ggplot(distr, aes(x=f32)) + geom_density() + geom_vline(aes(xintercept=median(f32)),color="blue", linetype="dashed", size=1)
  p4=ggplot(distr, aes(x=f41)) + geom_density() + geom_vline(aes(xintercept=median(f41)),color="blue", linetype="dashed", size=1)
  p5=ggplot(distr, aes(x=f42)) + geom_density() + geom_vline(aes(xintercept=median(f42)),color="blue", linetype="dashed", size=1)
  p6=ggplot(distr, aes(x=f43)) + geom_density() + geom_vline(aes(xintercept=median(f43)),color="blue", linetype="dashed", size=1)
  p7=ggplot(distr, aes(x=f53)) + geom_density() + geom_vline(aes(xintercept=median(f53)),color="blue", linetype="dashed", size=1)
  p8=ggplot(distr, aes(x=f64)) + geom_density() + geom_vline(aes(xintercept=median(f64)),color="blue", linetype="dashed", size=1)
  p9=ggplot(distr, aes(x=lig_aboveleaf)) + geom_density() + geom_vline(aes(xintercept=median(lig_aboveleaf)),color="blue", linetype="dashed", size=1)
  p10=ggplot(distr, aes(x=lig_belowleaf)) + geom_density() + geom_vline(aes(xintercept=median(flig_belowleaf)),color="blue", linetype="dashed", size=1)
  p11=ggplot(distr, aes(x=k1)) + geom_density() + geom_vline(aes(xintercept=median(k1)),color="blue", linetype="dashed", size=1)
  p12=ggplot(distr, aes(x=k2)) + geom_density() + geom_vline(aes(xintercept=median(k2)),color="blue", linetype="dashed", size=1)
  p13=ggplot(distr, aes(x=k3)) + geom_density() + geom_vline(aes(xintercept=median(k3)),color="blue", linetype="dashed", size=1)
  p14=ggplot(distr, aes(x=k4)) + geom_density() + geom_vline(aes(xintercept=median(k4)),color="blue", linetype="dashed", size=1)
  p15=ggplot(distr, aes(x=k5)) + geom_density() + geom_vline(aes(xintercept=median(k5)),color="blue", linetype="dashed", size=1)
  p16=ggplot(distr, aes(x=k6)) + geom_density() + geom_vline(aes(xintercept=median(k6)),color="blue", linetype="dashed", size=1)
  p17=ggplot(distr, aes(x=k4)) + geom_density() + geom_vline(aes(xintercept=median(k4)),color="blue", linetype="dashed", size=1)
  p18=ggplot(distr, aes(x=k5)) + geom_density() + geom_vline(aes(xintercept=median(k5)),color="blue", linetype="dashed", size=1)
  p19=ggplot(distr, aes(x=k6)) + geom_density() + geom_vline(aes(xintercept=median(k6)),color="blue", linetype="dashed", size=1)
  p20=ggplot(distr, aes(x=x2_init)) + geom_density() + geom_vline(aes(xintercept=median(x2_init)),color="blue", linetype="dashed", size=1)
  p21=ggplot(distr, aes(x=x3_init)) + geom_density() + geom_vline(aes(xintercept=median(x3_init)),color="blue", linetype="dashed", size=1)
  p22=ggplot(distr, aes(x=x4_init)) + geom_density() + geom_vline(aes(xintercept=median(x4_init)),color="blue", linetype="dashed", size=1)
  p23=ggplot(distr, aes(x=x7_init)) + geom_density() + geom_vline(aes(xintercept=median(x7_init)),color="blue", linetype="dashed", size=1)
  p24=ggplot(distr, aes(x=x8_init)) + geom_density() + geom_vline(aes(xintercept=median(x8_init)),color="blue", linetype="dashed", size=1)
  p25=ggplot(distr, aes(x=x9_init)) + geom_density() + geom_vline(aes(xintercept=median(x9_init)),color="blue", linetype="dashed", size=1)
  p26=ggplot(distr, aes(x=sens_moist_x7)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x7)),color="blue", linetype="dashed", size=1)
  p27=ggplot(distr, aes(x=sens_moist_x8)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x8)),color="blue", linetype="dashed", size=1)
  p28=ggplot(distr, aes(x=sens_moist_x9)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x9)),color="blue", linetype="dashed", size=1)
  p29=ggplot(distr, aes(sens_moist_x10)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x10)),color="blue", linetype="dashed", size=1)
  p30=ggplot(distr, aes(sens_moist_x11)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x11)),color="blue", linetype="dashed", size=1)
  p31=ggplot(distr, aes(sens_moist_x12)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x12)),color="blue", linetype="dashed", size=1)
  p32=ggplot(distr, aes(sens_moist_x13)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x13)),color="blue", linetype="dashed", size=1)
  # arrange density plots in a 5x4 grid
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32, nrow = 6)
}
# set optimized parameters as median values of the distribution
pa_final=rep(0,paramNum)
for (i in 1:paramNum) {pa_final[i]=median(distr[[i]])}
# compare original and optimized parameters
pa
pa_final
# run the simulation with optimized parameters
optimized = matrix_simu(pa_final)

#compare results
par(mfrow=c(3, 4)) #set plot window to 3x4 grid
{
  plot(test[[1]]$x1[1:120], type="l", col="red",  xlab="month", 
       ylim=c(min(min(test[[1]]$x1[1:120]),min(dat$cLeaf[1:120])),max(max(test[[1]]$x1[1:120]),max(dat$cLeaf[1:120]))), 
       xaxp  = c(0, 120, 10), ylab="cLeaf", main="Leaf Pool modelled vs CMIP6")
  lines(dat$cLeaf[1:120], col="blue")
  lines(optimized[[1]]$x1[1:120], col="green")
  
  plot(test[[1]]$x2[1:120], type="l", col="red", xlab="month", 
       ylim=c(min(min(test[[1]]$x2[1:120]),min(dat$cWood[1:120])),max(max(test[[1]]$x2[1:120]),max(dat$cWood[1:120]))), 
       xaxp  = c(0, 120, 10),ylab="cWood", main="Wood Pool modelled vs CMIP6")
  lines(dat$cWood[1:120], col="blue")
  lines(optimized[[1]]$x2[1:120], col="green")
  
  plot(test[[1]]$x3[1:120], type="l", col="red", xlab="month", 
       ylim=c(min(min(test[[1]]$x3[1:120]),min(dat$cRoot[1:120])),max(max(test[[1]]$x3[1:120]),max(dat$cRoot[1:120]))), 
       xaxp  = c(0, 120, 10),ylab="cRoot", main="Root Pool modelled vs CMIP6")
  lines(dat$cRoot[1:120], col="blue")
  lines(optimized[[1]]$x3[1:120], col="green")
  
  plot(test[[1]]$x4[1:120]+test[[1]]$x5[1:120], type="l", col="red", 
       ylim=c(min(min(test[[1]]$x4[1:120]+test[[1]]$x5[1:120]),min(dat$cLitterAbove[1:120])),max(max(test[[1]]$x4[1:120]+test[[1]]$x5[1:120]),max(dat$cLitterAbove[1:120]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cLitterAbove", main="Above-ground Litter (leaf+wood) modelled vs CMIP6")
  lines(dat$cLitterAbove[1:120], col="blue")
  lines(optimized[[1]]$x4[1:120]+optimized[[1]]$x5[1:120], col="green")
  
  plot(test[[1]]$x6[1:120], type="l", col="red", 
       ylim=c(min(min(test[[1]]$x6[1:120]),min(dat$cLitterBelow[1:120])),max(max(test[[1]]$x6[1:120]),max(dat$cLitterBelow[1:120]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cLitterBelow", main="Below-ground Litter (root) modelled vs CMIP6")
  lines(dat$cLitterBelow[1:120], col="blue")
  lines(optimized[[1]]$x6[1:120], col="green")
  
  plot(test[[1]]$x7[1:120], type="l", col="red", 
       ylim=c(min(min(test[[1]]$x7[1:120]),min(dat$cSoilFast[1:120])),max(max(test[[1]]$x7[1:120]),max(dat$cSoilFast[1:120]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilFast", main="Fast SOM pool modelled vs CMIP6")
  lines(dat$cSoilFast[1:120], col="blue")
  lines(optimized[[1]]$x7[1:120], col="green")
  
  plot(test[[1]]$x8[1:120], type="l", col="red",  
       ylim=c(min(min(test[[1]]$x8[1:120]),min(dat$cSoilMedium[1:120])),max(max(test[[1]]$x8[1:120]),max(dat$cSoilMedium[1:120]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilMEdium", main="Medium, SOM pool modelled vs CMIP6")
  lines(dat$cSoilMedium[1:120], col="blue")
  lines(optimized[[1]]$x8[1:120], col="green")
  
  plot(test[[1]]$x9[1:120], type="l", col="red",  
       ylim=c(min(min(test[[1]]$x9[1:120]),min(dat$cSoilMedium[1:120])),max(max(test[[1]]$x9[1:120]),max(dat$cSoilMedium[1:120]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilSlow", main="Slow SOM pool modelled vs CMIP6")
  lines(dat$cSoilSlow[1:120], col="blue")
  lines(optimized[[1]]$x9[1:120], col="green")
  
  plot(test[[2]][1:120], type="l", col="red",  
       ylim=c(min(min(test[[2]][1:120]),min(dat$rh[1:120])),max(max(test[[2]][1:120]),max(dat$rh[1:120]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="Rh", main="Heterotrophic respiration modelled vs CMIP6")
  lines(dat$rh[1:120], col="blue")
  lines(optimized[[2]][1:120], col="green")
  
  plot(test[[3]][1:120], type="l", col="red", 
       ylim=c(min(min(test[[3]][1:120]),min(dat$fVegLitter[1:120])),max(max(test[[3]][1:120]),max(dat$fVegLitter[1:120]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="fVegLitter", main="Litterfall modelled vs CMIP6")
  lines(dat$fVegLitter[1:120], col="blue")
  lines(optimized[[3]][1:120], col="green")
  
  plot(test[[4]][1:120], type="l", col="red", 
       ylim=c(min(min(test[[4]][1:120]),min(dat$fLitterSoil[1:120])),max(max(test[[4]][1:120]),max(dat$fLitterSoil[1:120]))), 
       xlab="month", xaxp  = c(0, 120, 10),ylab="fLitterSoil", main="Humus formation modelled vs CMIP6")
  lines(dat$fLitterSoil[1:120], col="blue")
  lines(optimized[[4]][1:120], col="green")
  
  plot(2, xlim=c(0,1), ylim=c(0,1), axes = F, main="legend", ylab="")
  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Modelled Non-optimized", "Modelled Optimized"),
         col=c("blue", "red", "green"), lty=1, cex=1)
  
}
#set plot window back to 1x1
par(mfrow=c(1, 1))