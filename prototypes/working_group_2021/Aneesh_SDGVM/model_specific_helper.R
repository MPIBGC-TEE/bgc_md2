## Aneesh SDGVM model specific helper file

rm(list=ls())
library(raster)
library(sp)
library(ncdf4)
library(dplyr)
library(mvtnorm)
source ("./general_helpers.R") 


# fixme:
# Your parameters will most likely differ but you can still use the
# distinctions between different sets of parameters. The aim is to make
# the different tasks in the code more obvious. In principal you could
# have a lot of overlapping sets and just have to keep them consistent. 
# 'namedtuples' can be used just as normal tuples by functions
# that are not aware of the names. They can still use the positions like 
# in the original code

# @Kostia and the 'R'tists: 
# It is not necessary to replicate the complete functionality of the #
# namedtuple classes. A simple 'R' approximation is a list with named entries 
# pa=list(C_leaf_0=2,...)

# This set is used by the functions that produce the 
# specific ingredients (functions) that will be run by
# mcmc alg.

UnEstimatedParameters = list(
  C_leaf_wood_0=0,
  C_root_0=0,
  C_soil_passive_0=0,
  npp=0,
  number_of_months=0
)

# This set is used (as argument) by the functions that are called
# inside the mcmc
EstimatedParameters = list(
  beta1=0,
  f31=0,
  f52=0,
  f73=0,
  f74=0,
  f85=0,
  f95=0,
  f86=0,
  f97=0,
  f98=0,
  f108=0,
  f89=0,
  f109=0,
  f810=0,
  k1=0,
  k2=0,
  k3=0,
  k4=0,
  k5=0,
  k6=0,
  k7=0,
  k8=0,
  k9=0,
  k10=0,
  C_agb_stru_litter_0=0,
  C_agb_met_litter_0=0,
  C_blw_stru_litter_0=0,
  C_blw_met_litter_0=0,
  #C_surface_microbe_0=0,
  C_soil_microbe_0=0,
  C_slow_soil_0=0
  #C_soil_passive_0=0
)

# This is the set off all 
Parameters=c(EstimatedParameters,UnEstimatedParameters)

# This set defines the order of the c pools
# The order is crucial for the compatibility
# with the matrices (B and b) If you change it
# the matrix changes

#Carbon pools in the model available from trendy output
Statevariables = list(
  C_leaf_wood=0,
  C_root=0,
  C_litter=0,
  C_soil=0
)

#carbon pools in your model
Observables = list(
  C_leaf_wood=0,
  C_root=0,
  C_agb_stru_litter=0,
  C_agb_met_litter=0,
  C_blw_stru_litter=0,
  C_blw_met_litter=0,
  C_soil_microbe=0,
  C_slow_soil=0,
  rh=0
)

# We define another set of parameters which describes
# the parameters of the matrices A,K and the vector b
# and drivers like npp (in form of arrays)
# but does not include start values and hyperparameters like the 'number_of_months'
# This distinction is helpful for the forward simulation where the
# distinction between estimated and constant is irrelevant.

ModelParameters = Parameters[!( names(Parameters) %in% list(
  'C_leaf_wood_0',
  'C_root_0',
  'C_agb_stru_litter_0',
  'C_agb_met_litter_0',
  'C_blw_stru_litter_0',
  'C_blw_met_litter_0',
  'C_soil_microbe_0',
  'C_slow_soil_0',
  'rh_0',
  'number_of_months')
)]
# to do
# 1.) make a namedtuple for the yycable data or use xarray to create a multifile dataset

#    # Read NetCDF data  ******************************************************************************************************************************

# get_example_site_vars<-function(dataPath, lon, lat){
#   
#   # pick up 1 site
#   point<-SpatialPoints(cbind(lon,lat), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
#   
#   # load all netCDF files in the folder
#   files<-list.files(dataPath, pattern="..nc") # list all netCDF files in the folder
#   #(var_names<-substr(files, 0 ,regexpr("_",files)-1)) # derive variable names from file names
#   (var_names=files)
#   dat<-data.frame(point)
#   
#   for (i in 1:length(var_names)){
#     r<-stack(paste0(dataPath,"/",files[i]))
#     r_dat<-extract(r,point)
#     dat<-cbind(dat,r_dat[1,])
#   }
#   names(dat)<-c("lon","lat", var_names) # assign variable names

# setwd("/Users/aneeshchandel/Desktop/Working_group/S1")  
npp<-stack("npp.nc")
clitter <- stack("cLitter.nc")
cRoot<-stack("cRoot.nc")
cSoil <- stack("cSoil.nc")
cVeg<-stack("cVeg.nc")
rh<-stack("rh.nc")

point<-SpatialPoints(cbind(33.3,50.0), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))

# extract data for each variable
npp_data<-extract(npp, point)
clitter_data <- extract(clitter, point)
cRoot_data<-extract(cRoot, point)
cSoil_data <- extract(cSoil, point)
cVeg_data<-extract(cVeg, point)
rh_data<-extract(rh, point)


dat <- data.frame(
  npp=npp_data[1,c(1:1400)],    #for reduced datapoints to make it uniform with other Trendy outputs from 1900 to 2014
  rh = rh_data[1,c(1:1400)]     #npp and rh are monthly
)

dat2 <- data.frame(
  litter  =clitter_data[1,c(204:320)],   #from 1900 to 2014
  root    =cRoot_data[1,c(204:320)],     #yearly
  soil    =cSoil_data[1,c(204:320)],
  Veg     =cVeg_data[1,c(204:320)]
)

leaf_wood = dat2$Veg - dat2$root
dat2$leaf_wood <- leaf_wood

  # correct fluxes from per second to per day
dat$npp<-dat$npp*86400
dat$rh<-dat$rh*86400
  

# function to check if mcmc-generated parameters can be accepted
make_param_filter_func<-function(C_max, C_min) {
  
  isQualified<-function(c){
    # fixme
    #   this function is model specific: It discards parameter proposals
    #   where beta1 + beta3 are > 1 and 
    #   where sum of transfer coefficients from the same pool > 1
    paramNum = length(c)
    flag = T
    for (i in 1:paramNum){
      if(c[i] > C_max[i] || c[i] < C_min[i]) {flag = F; break}
      if(c[6] + c[7] > 1){flag = F; break}
      if(c[10] + c[11] > 1){flag = F; break}
      if(c[12] + c[13] > 1){flag = F; break}
      #if(c[9] + c[10] + c[11] > 1){flag = F; break}
    }
    return (flag)
  }
  return (isQualified)
}

make_weighted_cost_func<-function(obs) {
  
  costfunction<-function(out_simu){
    # fixme
    #   as indicated by the fact that the function lives in this
    #   model-specific module it is not apropriate for (all) other models.
    #   There are model specific properties:
    #   1.) The weight for the different observation streams
    #
    
    # we assume the model output to be in the same shape and order
    # as the observation
    # this convention has to be honored by the forward_simulation as well
    # which in this instance already compresses the 3 different litter pools
    # to C_litter and the 3 different soil pools to one
    
    #this J_obj should represent the exact output pools from the Trendy output
    J_obj1 = mean (( out_simu$pools$C_leaf_wood - obs$C_leaf_wood)**2)/(2*var(obs$C_leaf_wood))
    J_obj3 = mean (( out_simu$pools$C_root - obs$C_root )**2)/(2*var(obs$C_root))
    J_obj4 = mean (( out_simu$pools$C_litter - obs$C_litter )**2)/(2*var(obs$C_litter))
    J_obj5 = mean (( out_simu$rh - obs$rh )**2)/(2*var(obs$rh))
     
    J_new= (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5)/200
    
    return (J_new)
  }
  return (costfunction)
}

make_param2res<-function(cpa){
  # """The returned function 'param2res' is used by the metropolis hastings algorithm
  #   to produce a result from a proposed parameter(tuple).
  #   'param2res' will be called only with the parameters to be estimated 'pe'.
  #   The constant parameters are the parameters of the present function 'pc' and are automatically 
  #   available in the returned 'param2res' (closure).
  # 
  #   In the case of this model 'param2res' has to perform the following
  #   tasks.
  #   -   use both sets of parameters 'pc' and 'pe' to build a forward model
  #   -   run the forward model to produce output for the times
  #       when observations are available.(Here monthly)
  #   -   project the output of the forward model to the observed variables.
  #       (here summation of all different soil-pools to compare to the single
  #       observed soil-pool and the same for the litter pools)
  # 
  #   In this version of the function all tasks are performed at once as in
  #   the original script.
  #   See the alternative implementation to see how all the tasks can be
  #   delegated to smaller functions.
  #   """
  # fixme
  # This function is model-specific in several ways:
  # 0. Which are the fixed and which are the estimated variables
  # 1. The matrices for the forward simulation
  # 2. the driver (here monthly npp)
  # 3. the projection of the simulated pool values to observable variables
  #    (here summation of
  
  param2res<-function(epa){
    # pa is a numpy array when pa comes from the predictor
    # so we transform it to be able to use names instead of positions
    #epa=EstimatedParameters
    
    #SDGVM records observation on ervery 30 days so all months are considered to be 30 days
    days = c(30,30,30,30,30,30,30,30,30,30,30,30)# Construct b vector
    
    # leaf, root, wood
    beta1=epa$beta1; beta2=1 -beta1 
    
    B = c(beta1, beta2, 0, 0, 0, 0, 0, 0, 0, 0)   # allocation
    # transfer coefficients
    # f41=1; f52=1; f63=1
    # f74 = epa$f_leaflit2fastsom; f84=epa$f_leaflit2slowsom;  f94=epa$f_leaflit2passsom;
    # f75=epa$f_woodlit2fastsom; f85=epa$f_woodlit2slowsom;  f95=epa$f_woodlit2passsom;
    # f76=epa$f_rootlit2fastsom; f86=epa$f_rootlit2slowsom; f96=epa$f_rootlit2passsom
    
    f31=epa$f31;f41=1-epa$f31; f52=epa$f52; f62=1 - epa$f52; f73 = epa$f73;
    f74=epa$f74; f85=epa$f85; f95=epa$f95; f86=epa$f86; f97=epa$f97; f98=epa$f98;
    f108=epa$f108; f89=epa$f89; f109=epa$f109; f810=epa$f810
    # A matrix
    A = c(-1,   0,   0,   0,   0,   0,     0,   0,     0,    0,
           0,  -1,   0,   0,   0,   0,     0,   0,     0,    0,
           f31, 0,  -1,   0,   0,   0,     0,   0,     0,    0,
           f41, 0,   0,  -1,   0,   0,     0,   0,     0,    0,
           0,  f52,  0,   0,  -1,   0,     0,   0,     0,    0,
           0,  f62,  0,   0,   0,  -1,     0,   0,     0,    0,
           0,   0,  f73, f74,  0,   0,    -1,   0,     0,    0,
           0,   0,   0,   0,  f85,  f86,   0,  -1,     f89,  f810, 
           0,   0,   0,   0,  f95,   0,   f97,  f98,  -1,     0,
           0,   0,   0,   0,   0,   0,     0,   f108,  f109,  -1)
    A = matrix(A, nrow = 10, byrow = TRUE)
    #turnover rate per day of pools:
    temp = c(epa$k1,epa$k2,epa$k3,epa$k4,epa$k5,epa$k6,epa$k7,epa$k8,epa$k9,epa$k10 ) 
    K = rep(0,100)
    #temp = c(epa$k_leaf,epa$k_wood,epa$k_root, epa$k_leaflit, epa$k_woodlit, epa$k_rootlit, epa$k_fastsom, epa$k_slowsom, epa$k_passsom)
    # K matrix
    #K = rep(0,81)
    K = matrix(K, nrow = 10, byrow = TRUE)
    for (i in 1:10) { K[i,i] = temp[i] }
    
    # X vector
    x=rep(0,cpa$number_of_months)
    x_fin=data.frame(x,x,x,x,x,x,x,x,x,x); 
    names(x_fin)=c("leaf_wood","root","leaf_wood_str_lit","leaf_wood_met_lit","root_str_lit","root_met_lit","surface_microbe", "fastsom","slowsom","passsom")
    rh_fin=rep(0,cpa$number_of_months);   f_veg_lit_fin=rep(0,cpa$number_of_months);   f_lit_soil_fin=rep(0,cpa$number_of_months)
    
    # "leaf_wood","root","leaf_wood_str_lit","leaf_wood_met_lit","root_str_lit","root_met_lit","surface_microbe", "fastsom","slowsom","passsom"
    #this represent the number of pools in your model i.e SDGVM in my case
    x_init = c(cpa$C_leaf_wood_0, # leaf & wood
               cpa$C_root_0, # root
               epa$C_agb_stru_litter_0, # leaf and wood structural litter
               epa$C_agb_met_litter_0, # leaf and wood metabolic litter
               epa$C_blw_stru_litter_0, #root structural litter
               epa$C_blw_met_litter_0, #root metabolic litter
               epa$C_surface_microbe_0,
               epa$C_soil_microbe_0, # soil microbe
               epa$C_slow_soil_0, # slow soil pool
               cpa$C_soil_passive_0)  # soil passive  
    X=x_init   # initialize carbon pools 
    
    # initialize first respiration value
    #co2_rh=cpa.rh_0
    # fixme:
    # slight change to the original
    # I would like to start the solution with the initial values
    # m=0 means after 0 moths = in the initial step
    #B=A@K
    #pa=Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)
    #B=make_compartmental_matrix_func(pa)(0,X)
    
    ############################ Additional functions specific for VISITe ##############
    
    # modifier for leaf pool turnover rate for deciduous forest depending on air temp
    # leaf_fall<-function (TA)
    # {
    #   lf=0
    #   if (TA<281) {lf<-1}
    #   return (lf)
    # }
    # # modifier for soil pools turnover rate to implement time-dependent soil respiration
    # Rh_calc<-function (TS, M, T0, E, KM)
    # {
    #   TS<-TS-273.15
    #   if (TS>T0) {rh_out=exp(E*(1/(10-T0)-1/(TS-T0))) *M/(KM+M)}  else {rh_out=0}
    #   return(rh_out)
    # }
    ######################################################################################
    
    for (m in 1:cpa$number_of_months){
      npp_in = cpa$npp[m] 
      co2_rh = 0; #f_veg_lit = 0; f_lit_soil = 0  
      
      # environmental factor ksi: modifies leaf fall and turnover rates / respiration 
      # depending on monthly surface t "ts", soil t "tsl" and soil moisture data "mrso"
      # rh_modifier - VISIT-specific - remove if not relevant for your model
      #rh_modifier=Rh_calc(cpa$tsl[m], cpa$mrso[m], epa$T0, epa$E, epa$KM)
      # ksi vector (remove if no environmental modifiers)
      #ksi=c(leaf_fall(cpa$ts[m]), # leaf pool
      #      1, # wood pool
      #      1, # root pool
      #      rh_modifier, # leaf liter
      #      rh_modifier, # wood litter
      #      rh_modifier, # root litter
      #      rh_modifier, # fast soil
      #      rh_modifier, # slow soil
      #      rh_modifier) # passive soil
      
      for (d in 1:30) {
        
        # matrix equation with ksi (remove ksi if no environmental modifiers)
        X = X + B * npp_in + A %*% K %*% X #* ksi
        # deriving rh from each litter and soil pool as 1 - sum of transfer coefficients 
        co2_rate = c(0,
                     0, 
                     (1-f73)*K[3,3],#*ksi[4],
                     (1-f74)*K[4,4],#*ksi[5],
                     (1-f85-f95)*K[5,5],#*ksi[6], 
                     (1-f86)*K[6,6],#*ksi[7], 
                     (1-f97)*K[7,7],#*ksi[8], 
                     (1-f98-f108)*K[8,8],
                     (1-f89-f109)*K[9,9],
                     (1-f810)*K[10,10])#*ksi[9])
        co2=sum(co2_rate*X)
        #co2_rh = co2_rh + co2/days[m%%12+1]   # monthly average rh
        co2_rh = co2_rh + co2/30   # monthly average rh
        # deriving litterfall
        
        litterfall_rate = c((f31 + f41)*K[1,1],(f52+f62)*K[2,2],0,0,0,0,0,0,0,0)
        litterfall=sum(litterfall_rate*X)
        #f_veg_lit=f_veg_lit+litterfall/days[m%%12]
        #f_veg_lit=f_veg_lit+litterfall/30
        # deriving humus formation
        litter_to_soil_rate = c(0,0,f73*K[3,3],f74*K[4,4],(f85+f95)*K[5,5],f86*K[6,6],0,0,0,0)
        litter_to_soil=sum(litter_to_soil_rate*X)
        f_lit_soil=f_lit_soil+litter_to_soil/30
      }
      x_fin[m,]=X
      rh_fin[m]=co2_rh
      #f_veg_lit_fin[m]=f_veg_lit
      f_lit_soil_fin[m]=f_lit_soil
    }
    
    # We create an output that has the same shape
    # as the observations to make the costfunctions
    # easier.
    # To this end we project our 10 output variables of the matrix simulation
    # onto the 6 data streams by summing up the 3 litter pools into one
    # and also the 3 soil pools into one
    # out_simu = list(
    #   C_leaf=x_fin$leaf,
    #   C_wood=x_fin$wood,
    #   C_root=x_fin$root,
    #   C_litter_above=x_fin$leaflit+x_fin$woodlit,
    #   C_litter_below=x_fin$rootlit,
    #   C_fastsom=x_fin$fastsom,
    #   C_slowsom=x_fin$slowsom,
    #   C_passsom=x_fin$passsom,
    #   rh=rh_fin,
    #   f_veg2litter=f_veg_lit_fin,
    #   f_litter2som=f_lit_soil_fin
    # )
    names(x_fin)=c("leaf_wood","root","leaf_wood_str_lit","leaf_wood_met_lit","root_str_lit","root_met_lit","surface_microbe", "fastsom","slowsom","passsom")
    pools = list(
      C_leaf_wood=x_fin$leaf_wood,
      C_root=x_fin$root,
      C_litter = x_fin$leaf_wood_str_lit + x_fin$leaf_wood_met_lit + x_fin$root_str_lit +
        x_fin$root_met_lit + x_fin$surface_microbe,
      C_soil = x_fin$fastsom + x_fin$slowsom + x_fin$passsom,
      f_veg2litter=f_veg_lit_fin,
      f_litter2som=f_lit_soil_fin
    )
    rh = rh_fin
    out_simu = list (pools, rh)
    
    return (out_simu)
  }
  return (param2res)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# 
# # alternative implementation not yet translated to R