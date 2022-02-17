library (mvtnorm)
make_uniform_proposer<-function(C_max, C_min, D, filter_func) {
  # '''Returns a function that will be used by the mcmc algorithm to propose
  # a new parameter value tuple based on a given one. 
  # The two arrays C_max and C_min define the boundaries
  # of the n-dimensional rectangular domain for the parameters and must be of
  # the same shape.  After a possible parameter value has been sampled the
  # filter_func will be applied to it to either accept or discard it.  So
  # filter func must accept parameter array and return either True or False'''
  
  GenerateParamValues<-function(C_op){
    paramNum = length(C_op)
    flag = T
    while (flag) {
      C_new = as.data.frame(C_op) + (runif((paramNum)) - 0.5)*(C_max - C_min)/D
      #C_new = as.data.frame(C_op) + (runif((paramNum)) - 0.5)*(C_max - C_min)/15.0
      if (filter_func(C_new)){flag = F}
    }
    names(C_new)=names(C_op)
    C_new=as.list(C_new)
    return (C_new)
  }
  return (GenerateParamValues)
}
make_multivariate_normal_proposer<-function(covv, filter_func){
  # """Returns a function that will be used by mcmc algorithm to propose
  # a new parameter(tuple) based on a given one.
  # param: covv: The covariance matrix (usually estimated from a previously run chain)
  # """
  GenerateParamValues<-function(C_op){
    flag = T
    #C_op2<-unlist(C_op, use.names=FALSE)
    while (flag){
      C_new = as.data.frame(C_op) + rmvnorm(1,mean = rep(0, length(C_op)),sigma=covv)
      if (filter_func(C_new)) {flag = F}
    }
    names(C_new)=names(C_op)
    C_new=as.list(C_new)
    return (C_new)
  }
  return (GenerateParamValues)
}

mcmc<-function(initial_parameters, proposer, param2res, costfunction, nsimu, K_accept=1) { # K_accept - coefficient to modify acceptance rate
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
  
  for (simu in 1:nsimu) {
    if (simu%%100==0) {print (paste0("simulation ",simu, " out of ", nsimu))} 
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
  
  output<-list(C_upgraded, J_upgraded, acceptance_rate)
  return (output)
}
make_feng_cost_func<-function(obs){
  # first unpack the observation array into its parts and estimate the mean
  # and standard deviation of each observable which has to be done only once
  # for a set of observables, hence we do it outside the actual costfunction
  # which will be used
  time_dim_ind = 0
  means = colMeans(obs)
  sigmas = sqrt(base::colSums((obs-means)^2))
  
  costfunction<-function (out_simu){
    out_simu<-as.data.frame(out_simu)
    output=sum( ((out_simu - means)/sigmas - (obs-means)/sigmas)^2) 
    return (output)
  }    
  return (costfunction)     
}