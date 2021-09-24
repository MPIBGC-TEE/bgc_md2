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
        C_op2<-unlist(C_op, use.names=FALSE)
        while (flag) {
            C_new = C_op2 + (runif((paramNum)) - 0.5)*(C_max - C_min)/10.0
            #C_new = C_op + (runif((paramNum)) - 0.5)*(C_max - C_min)/15.0
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
        C_op2<-unlist(C_op, use.names=FALSE)
        while (flag){
            C_new = C_op2 + rmvnorm(1,mean = rep(0, length(C_op2)),sigma=covv)
            if (filter_func(C_new)) {flag = F}
        }
        names(C_new)=names(C_op)
        C_new=as.list(C_new)
        return (C_new)
    }
return (GenerateParamValues)
}

mcmc<-function(initial_parameters, proposer, param2res, costfunction, nsimu) {
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

    upgraded=0;
    C_op = initial_parameters
    first_out = param2res(C_op)
    J_last = costfunction(first_out)
    #J_last = 400 # original code

    C_upgraded = rep(0,paramNum*nsimu)
    C_upgraded = matrix(C_upgraded, nrow = nsimu, byrow = TRUE)
    J_upgraded = rep(0, nsimu)

    for (simu in 1:nsimu) {
        print (paste0("simulation ",simu, " out of ", nsimu)) 
        C_new = proposer(C_op)

        out_simu = param2res(C_new)
        J_new = costfunction(out_simu)

        delta_J =  J_last - J_new;

        randNum = runif(1)
        if (min(1.0, exp(delta_J)) > randNum) {
            C_op=C_new;
            J_last=J_new;
            C_upgraded[upgraded,]=unlist(C_op, use.names=FALSE);
            J_upgraded[upgraded]=J_last; 
            upgraded=upgraded+1 }
        }
        output<-list(C_upgraded, J_upgraded)
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

#' def day_2_month_index(d):
#'     return months_by_day_arr()[(d%365)]
#' 
#' @lru_cache
#' def months_by_day_arr():
#'     days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]
#' return  np.concatenate(
#'     tuple(
#'         map(
#'             lambda m: m*np.ones(
#'                 days_per_month[m],
#'                 dtype=np.int64
#'             ),
#'             range(12)
#'         )
#'     )
#' )
#' }
#' 
#' def month_2_day_index(ns):
#'     """ computes the index of the day at the end of the month n in ns
#'     this works on vectors and is faster than a recursive version working
#'     on a single index (since the smaller indices are handled anyway)
#'     """
#' 
#' # We first compute the sequence of day indices up to the highest month in ns
#' # and then select from this sequence the day indices for the months in ns
#' days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]
#' dpm =  (days_per_month[i%len(days_per_month)] for i in range(max(ns)))
#' # compute indices for which we want to store the results which is the
#' # list of partial sums of the above list  (repeated)
#' 
#' def f(acc,el):
#'     if len(acc) <1:
#'     res = (el,)
#' else:
#'     last = acc[-1]
#' res = acc + (el+last,)
#' return res
#' day_indices_for_continuous_moths = reduce(
#'     f,
#'     dpm,
#'     (0,) 
#' )
#' day_indices = reduce( 
#'     lambda acc,n: acc + [day_indices_for_continuous_moths[n]], #for n=0 we want 0
#'     ns,
#'     []
#' )
#' return day_indices
#' 
#' class TimeStepIterator2():
#'     """iterator for looping forward over the results of a difference equation
#'     X_{i+1}=f(X_{i},i)"""
#' 
#' def __init__(
#'     self,
#'     initial_values, # a tupel of values that will be
#'     f, # the function to compute the next ts
#'     max_it = False
#' ):
#'     self.initial_values = initial_values
#' self.f= f
#' self.reset()
#' self.max_it = max_it
#' 
#' def reset(self):
#'     self.i = 0
#' self.ts = self.initial_values
#' 
#' def __iter__(self):
#'     self.reset()
#' return(self)
#' 
#' def __next__(self):
#'     if self.max_it:
#'     if self.i == self.max_it:
#'     raise StopIteration
#' 
#' ts = copy(self.ts)
#' ts_new = self.f(self.i, ts)
#' self.ts = ts_new
#' self.i += 1
#' return ts
#' 
#' def values(self,day_indices):
#'     # we traverse the iterator to the highest index and
#'     # collect the results we want to keep in a list (acc)
#'     tsi=copy(self)
#' tsi.reset()
#' def g(acc, i):
#'     v = tsi.__next__()
#' if i in day_indices:
#'     acc += [v] 
#' return acc
#' xs =  reduce(g,range(max(day_indices)+1),([]))
#' return xs
#' 
#' 
#' def respiration_from_compartmental_matrix(B,X):
#'     """This function computes the combined respiration from all pools"""
#' return -np.sum(B@X) 
#' 
#' 
#' def plot_solutions(
#'     fig,
#'     times, 
#'     var_names,
#'     tup,
#'     names=None
#' ):
#'     if names is None:
#'     names = tuple(str(i) for i in range(len(tup)))
#' 
#' assert(all([tup[0].shape == el.shape for el in tup]))
#' 
#' if tup[0].ndim == 1:
#'     n_times = tup[0].shape[0]
#' ax = fig.subplots(1,1)
#' for i,sol in enumerate(tup):
#'     ax.plot(
#'         np.array(times).reshape(n_times,), 
#'         sol,
#'         marker="o",
#'         label=names[i]
#'     )
#' ax.set_title(var_names[0])
#' ax.legend()
#' else:
#'     n_times, n_vars = tup[0].shape
#' 
#' fig.set_figheight(n_vars*fig.get_figwidth())
#' axs = fig.subplots(n_vars,1)
#' colors =('red','blue','green','orange')
#' for j in range(n_vars):
#'     for i,sol in enumerate(tup):
#'     axs[j].plot(
#'         np.array(times).reshape(n_times,), 
#'         sol[:, j],
#'         marker="o",
#'         label=names[i]
#'     )
#' axs[j].set_title(var_names[j])
#' axs[j].legend()
