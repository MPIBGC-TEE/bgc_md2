library(raster)
library(sp)
library(ncdf4)
library(dplyr)
source ("general_helpers.R") 

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
  C_leaf_0=0,
  C_stem_0=0,
  C_wood_0=0,
  C_root_0=0,
  C_litter_0=0,
  C_litter_above_0=0,
  C_litter_below_0=0,
  C_fastsom_0=0,
  C_mediumsom_0=0,
  C_slowsom_0=0,
  npp=0,
  number_of_months=0,
  tsl=0, # soil temperature - for dynamic turnover rates
  mrso=0, # soil moisture - for dynamic turnover rates
  ts=0 # air temperature - for dynamic turnover rates
)

# This set is used (as argument) by the functions that are called
# inside the mcmc
EstimatedParameters = list(
  beta_leaf=0,    #  1 (indices uses in original code) 
  beta_stem=0,    #  2
  beta_wood=0,    #  3
  f_leaf2abovestruclitter=0,  #  4
  f_stem2abovestruclitter=0,#  5
  f_wood2abovestruclitter=0,#  6
  f_leaf2abovemetablitter=0,  # 7
  f_stem2abovemetablitter=0,#  8
  f_wood2abovemetablitter=0,#  9
  f_root2belowstruclitter=0,  #  10
  f_root2belowmetablitter=0,#  11
  f_aboveligninfrac=0,#  12
  f_belowligninfrac=0,       #  13
  k_leaf=0,       #  14
  k_stem=0,       #  15
  k_wood=0,       #  16
  k_root=0,       #  17
  k_above_structurelit=0,	#  18
  k_above_metaboliclit=0,	#  19
  k_below_structurelit=0,	#  20
  k_below_metaboliclit=0, #21
  k_fastsom=0,	#  22
  k_mediumsom=0,	# 23
  k_slowsom=0,	# 24
  C_above_structurelit_0=0, # 25
  C_above_metaboliclit_0=0, #- 26
  C_below_structurelit_0=0, #- 27
  C_below_metaboliclit_0=0, #- 28
  T0=0,	# 22
  E=0,	# 23
  KM=0  # 24
)


# This is the set off all 
Parameters=c(EstimatedParameters,UnEstimatedParameters)

# This set defines the order of the c pools
# The order is crucial for the compatibility
# with the matrices (B and b) If you change ist
# the matrix changes
StateVariables = list(
  C_leaf=0,
  C_stem=0,
  C_wood=0,
  C_root=0,
  C_above_structurelit=0,
  C_above_metaboliclit=0, 
  C_below_structurelit=0, 
  C_below_metaboliclit=0,
  C_fastsom=0,
  C_mediumsom=0,
  C_slowsom=0
)
Observables = list(
  C_leaf=0,
  C_stem=0,
  C_wood=0,
  C_root=0,
  C_litter=0,
  C_litter_above=0,
  C_litter_below=0,    
  C_fastsom=0,
  C_mediumsom=0,
  C_slowsom=0,
  rh=0,
  f_veg2litter=0,
  f_litter2som=0
)

# We define another set of parameters which describes
# the parameters of the matrices A,K and the vector b
# and drivers like npp (in form of arrays)
# but does not include start values and hyperparameters like the 'number_of_months'
# This distinction is helpful for the forward simulation where the
# distinction between estimated and constant is irrelevant.
ModelParameters = Parameters[!( names(Parameters) %in% list(
  'C_leaf_0',
  'C_stem_0',
  'C_root_0',
  'C_wood_0',
  'C_litter_above_0',
  'C_litter_below_0',
  'C_fastsom_0',
  'C_mediumsom_0',
  'C_slowsom_0',
  'rh_0',
  'C_above_structurelit_0',
  'C_above_metaboliclit_0',
  'C_below_structurelit_0', 
  'C_below_metaboliclit_0',
  'number_of_months')                    
)] 

# to do
# 1.) make a namedtuple for the yycable data or use xarray to create a multifile dataset

#    # Read NetCDF data  ******************************************************************************************************************************
get_example_site_vars<-function(dataPath, lon, lat){
  
  # pick up 1 site
  point<-SpatialPoints(cbind(lon,lat), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  
  # load all netCDF files in the folder
  files<-list.files(dataPath, pattern="..nc") # list all netCDF files in the folder
  (var_names<-substr(files, 0 ,regexpr("_",files)-1)) # derive variable names from file names
  dat<-data.frame(point)
  
  for (i in 1:length(var_names)){
    r<-stack(paste0(dataPath,"/",files[i]))
    r_dat<-extract(r,point)
    dat<-cbind(dat,r_dat[1,])
  }
  names(dat)<-c("lon","lat", var_names) # assign variable names
  
  # correct fluxes from per second to per day
  dat$gpp<-dat$gpp*86400
  dat$npp<-dat$npp*86400
  dat$fVegLitter<-dat$fVegLitter*86400
  dat$fLitterSoil<-dat$fLitterSoil*86400
  dat$rh<-dat$rh*86400
  
  # explore the data
  print("Importing data:")
  names(dat)
  #save the data to a file
  
  output=data.frame(
    C_leaf=dat$cLeaf,
    C_stem=dat$cStem,
    C_wood=dat$cWood,
    C_root=dat$cRoot,
    c_litter=dat$cLitter,
    C_litter_above=dat$cLitterSurf,
    C_litter_below=dat$cLitterBelow,
    C_fastsom=dat$cSoilFast,
    C_mediumsom=dat$cSoilMedium,
    C_slowsom=dat$cSoilSlow,
    C_soil=dat$cSoil,
    npp=dat$npp,
    gpp=dat$gpp,
    rh=dat$rh,
    f_veg2litter=dat$fVegLitter,
    f_litter2som=dat$fLitterSoil,
    tsl=dat$tsl, # soil temperature - for dynamic turnover rates
    mrso=dat$mrsos, # soil moisture - for dynamic turnover rates
    ts=dat$ts # air temperature - for dynamic turnover rates
  )
  write.csv(output,paste0(dataPath,"/dat.csv"))
  return(output)
  
}  

# read data from a csv file if previously saved  
get_data_from_file<-function(dataPath){
  dat<-read.csv(paste0(dataPath,"/dat.csv")) 
  #output=list(dat$npp,dat$nppLeaf,dat$nppWood,dat$nppRoot,dat$fVegLitter,dat$fLitterSoil,dat$fLitterSoil, dat$rh,
  #            dat$cLeaf, dat$cVeg, dat$cRoot, dat$cLitterAbove, dat$cLitterBelow, dat$cSoilFast, dat$cSoilMedium, dat$cSoilSlow, 
  #            dat$tsl, dat$mrso, dat$ts)
  return(dat)
}

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
      if(c[1] + c[2] > 1){flag = F; break}
      if(c[3] + c[4] + c[5] > 1){flag = F; break}
      if(c[6] + c[7] + c[8] > 1){flag = F; break}
      if(c[9] + c[10] + c[11] > 1){flag = F; break}
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
    
    J_obj1 = mean (( out_simu$C_leaf - obs$C_leaf)**2)/(2*var(obs$C_leaf))
    J_obj2 = mean (( out_simu$C_stem - obs$C_stem )**2)/(2*var(obs$C_stem))
    J_obj3 = mean (( out_simu$C_wood - obs$C_wood )**2)/(2*var(obs$C_wood))
    J_obj4 = mean (( out_simu$C_root - obs$C_root )**2)/(2*var(obs$C_root))
    #J_obj5 = mean (( out_simu$C_litter - obs$C_litter )**2)/(2*var(obs$C_litter))
    #J_obj6 = mean (( out_simu$C_litter_below - obs$C_litter_below )**2)/(2*var(obs$C_litter_below))
    J_obj5 = mean (( out_simu$C_soil - obs$C_soil )**2)/(2*var(obs$C_soil))
    #J_obj7 = mean (( out_simu$C_fastsom - obs$C_fastsom )**2)/(2*var(obs$C_fastsom))
    # J_obj8 = mean (( out_simu$C_mediumsom - obs$C_mediumsom )**2)/(2*var(obs$C_mediumsom))
    # J_obj9 = mean (( out_simu$C_slowsom - obs$C_slowsom )**2)/(2*var(obs$C_slowsom))
    J_obj6 = mean (( out_simu$rh - obs$rh )**2)/(2*var(obs$rh))

    
    J_new= (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 + J_obj6)/200
    
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
    days = c(31,28,31,30,31,30,31,31,30,31,30,31)# Construct b vector
    
    # leaf, root, wood
    beta1=epa$beta_leaf; beta2=epa$beta_stem; beta3=epa$beta_wood 
    beta4 = 1- beta1- beta2- beta3
    B = c(beta1, beta2, beta3, beta4, 0, 0, 0,0,0,0,0)   # allocation
    # transfer coefficients
    f51 =epa$f_leaf2abovestruclitter; f52 =epa$f_stem2abovestruclitter;f53 =epa$f_wood2abovestruclitter;  #pa[4]+pa[7]<=1;pa[5]+pa[8]<=1;pa[6]+pa[9]<=1;
    f61 =epa$f_leaf2abovemetablitter; f62 =epa$f_stem2abovemetablitter;f63 =epa$f_wood2abovemetablitter; 
    f74 =epa$f_root2belowstruclitter; 
    f84 =epa$f_root2belowmetablitter; #pa[10]+pa[11]<=1;
    lig_aboveleaf=epa$f_aboveligninfrac; #lignin fraction of above-ground structural litter  pa[12]<=3
    lig_belowleaf =epa$f_belowligninfrac;#lignin fraction of below-ground structural litte   pa[13]<=2.2
    f95 = 0.55*(1-lig_aboveleaf); f96 = 0.45; f97 = 0.45*(1-lig_belowleaf); f98 = 0.45; f910 = 0.42; f911 = 0.55; 
    f105 = 0.7*lig_aboveleaf; f107 =0.7*lig_belowleaf;f109 = 0.146+0.68*f_clay_silt;  #clay+silt<=1.25
    f119 = 0.004; f1110 = 0.03
    # A matrix
    A = c(-1,  0,    0,   0,   0,   0,   0,   0,   0,   0,   0,
          0,  -1,    0,   0,   0,   0,   0,   0,   0,   0,   0,
          0,   0,   -1,   0,   0,   0,   0,   0,   0,   0,   0,
          0,   0,    0,  -1,   0,   0,   0,   0,   0,   0,   0,
          f51, f52, f53,  0,  -1,   0,   0,   0,   0,   0,   0,
          f61, f62, f63,  0,   0,  -1,   0,   0,   0,   0,   0,
          0,   0,    0, f74,   0,   0,  -1,   0,   0,   0,   0,
          0,   0,    0, f84,   0,   0,   0,  -1,   0,   0,   0,
          0,   0,   0,    0,  f95, f96, f97,f98,  -1,f910, f911,
          0,   0,   0,    0,  f105, 0,f107,   0, f109,  -1,  0,
          0,   0,   0,   0,    0,   0,   0,   0,f119,f1110, -1)
    A = matrix(A, nrow = 11, byrow = TRUE)
    #turnover rate per day of pools:
    temp = c(epa$k_leaf,epa$k_stem,epa$k_wood,epa$k_root, epa$k_above_structurelit, epa$k_above_metaboliclit,epa$k_below_structurelit, epa$k_below_metaboliclit,epa$k_fastsom, epa$k_mediumsom, epa$k_slowsom)
    # K matrix
    K = rep(0,121)
    K = matrix(K, nrow = 11, byrow = TRUE)
    for (i in 1:11) { K[i,i] = temp[i] }
    
    # X vector
    x=rep(0,cpa$number_of_months)
    x_fin=data.frame(x,x,x,x,x,x,x,x,x,x,x); 
    names(x_fin)=c("leaf","stem","wood","root","above_structurelit","above_metaboliclit","below_structurelit","below_metaboliclit","fastsom","mediumsom","slowsom")
    rh_fin=rep(0,cpa$number_of_months);   f_veg_lit_fin=rep(0,cpa$number_of_months);   f_lit_soil_fin=rep(0,cpa$number_of_months)
    
    # "leaf","stem","wood","root","abovstruclit","abovmetalit","beltruclit","belmetalit","soil_fast","soil_medium","soil_slow"
    x_init = c(cpa$C_leaf_0, # leaf
               cpa$C_stem_0, # stem
               cpa$C_wood_0, # stem
               cpa$C_root_0, # root
               epa$C_above_structurelit_0, #- unknown
               epa$C_above_metaboliclit_0, #- unknown
               epa$C_below_structurelit_0, #- unknown
               epa$C_below_metaboliclit_0, #- unknown
               cpa$C_fastsom_0, # soil active
               cpa$C_mediumsom_0, # soil intermediate
               cpa$C_slowsom_0)  # soil passive  
    X=x_init   # initialize carbon pools 
  
    # initialize first respiration value
    #co2_rh=cpa$rh_0
    # fixme:
    # slight change to the original
    # I would like to start the solution with the initial values
    # m=0 means after 0 moths = in the initial step
    #B=A@K
    #pa=Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)
    #B=make_compartmental_matrix_func(pa)(0,X)
    
    ############################ Additional functions specific for VISITe ##############
    
    # modifier for leaf pool turnover rate for deciduous forest depending on air temp
    leaf_fall<-function (TA)
       {
         lf=0
         if (TA<281) {lf<-1}
         return (lf)
        }
    # modifier for soil pools turnover rate to implement time-dependent soil respiration
    Rh_calc<-function (TS, M, T0, E, KM)
        {
          TS<-TS-273.15
         if (TS>T0) {rh_out=exp(E*(1/(10-T0)-1/(TS-T0))) *M/(KM+M)}  else {rh_out=0}
         return(rh_out)
        }
    ######################################################################################
    
    for (m in 1:cpa$number_of_months){
      npp_in = cpa$npp[m] 
      co2_rh = 0; f_veg_lit = 0; f_lit_soil = 0  
      
      # environmental factor ksi: modifies leaf fall and turnover rates / respiration 
      # depending on monthly surface t "ts", soil t "tsl" and soil moisture data "mrso"
      # rh_modifier - VISIT-specific - remove if not relevant for your model
      rh_modifier=Rh_calc(cpa$tsl[m], cpa$mrso[m], epa$T0, epa$E, epa$KM)
      #ksi vector (remove if no environmental modifiers)
      ksi=c(leaf_fall(cpa$ts[m]), # leaf pool
             1, # stem pool
             1, # wood pool
             1, # root pool
             rh_modifier, # above structure litter
             rh_modifier, # above metabolic litter
             rh_modifier, # below structure litter
             rh_modifier, # below metabolic litter
             rh_modifier, # fast soil
             rh_modifier, # medium soil
             rh_modifier) # slow soil
      
      for (d in 1:days[m%%12+1]) {
        
        # matrix equation with ksi (remove ksi if no environmental modifiers)
        X = X + B * npp_in +  A %*% K %*% X * ksi
        # deriving rh from each litter and soil pool as 1 - sum of transfer coefficients 
        co2_rate = c(0,0,0,0, 
                     (1-f95-f105)*K[5,5]*ksi[5],
                     (1-f96)*K[6,6]*ksi[6],
                     (1-f97-f107)*K[7,7]*ksi[7], 
                     (1-f98)*K[8,8]*ksi[8], 
                     (1-f109-f119)*K[9,9]*ksi[9], 
                     (1-f910-f1110)*K[10,10]*ksi[10], 
                     (1-f911)*K[11,11]*ksi[11])
        co2=sum(co2_rate*X)
        co2_rh = co2_rh + co2/days[m%%12+1]   # monthly average rh
        # deriving litterfall
        litterfall_rate = c((f51+f61)*K[1,1]*ksi[1],(f52+f62)*K[2,2]*ksi[2],(f53+f63)*K[3,3]*ksi[3],(f74+f84)*K[4,4]*ksi[4],0,0,0,0,0,0,0)
        litterfall=sum(litterfall_rate*X)
        f_veg_lit=f_veg_lit+litterfall/days[m%%12+1]
        # deriving humus formation
        litter_to_soil_rate = c(0,0,0,0,f95*K[5,5]*ksi[5]+f105*K[5,5]*ksi[5],
                                f96*K[6,6]*ksi[6],
                                f97*K[7,7]*ksi[7]+f107*K[7,7]*ksi[7],
                                f98*K[8,8]*ksi[8],
                                0,0,0)
        litter_to_soil=sum(litter_to_soil_rate*X)
        f_lit_soil=f_lit_soil+litter_to_soil/days[m%%12+1]
      }
      x_fin[m,]=X
      rh_fin[m]=co2_rh
      f_veg_lit_fin[m]=f_veg_lit
      f_lit_soil_fin[m]=f_lit_soil
    }
    
    # We create an output that has the same shape
    # as the observations to make the costfunctions
    # easier.
    # To this end we project our 10 output variables of the matrix simulation
    # onto the 6 data streams by summing up the 3 litter pools into one
    # and also the 3 soil pools into one
    out_simu = list(
      C_leaf=x_fin$leaf,
      C_stem=x_fin$stem,
      C_wood=x_fin$wood,
      C_root=x_fin$root,
      C_litter=x_fin$above_structurelit+x_fin$above_metaboliclit+x_fin$below_structurelit+x_fin$below_metaboliclit,
      # C_litter_above=x_fin$above_structurelit+x_fin$above_metaboliclit,
      # C_litter_below=x_fin$below_structurelit+x_fin$below_metaboliclit,
      C_soil=x_fin$fastsom+x_fin$mediumsom+x_fin$slowsom,
      # C_fastsom=x_fin$fastsom,
      # C_mediumsom=x_fin$mediumsom,
      # C_slowsom=x_fin$slowsom,
      rh=rh_fin 
      # f_veg2litter=f_veg_lit_fin,
      # f_litter2som=f_lit_soil_fin
    )
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