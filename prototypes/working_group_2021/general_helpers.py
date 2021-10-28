import numpy as np
from tqdm import tqdm
from typing import Callable, Tuple, Iterable
from functools import reduce, lru_cache
from copy import copy
from time import time

def make_uniform_proposer(
        c_max: Iterable,
        c_min: Iterable,
        D: float, 
        filter_func: Callable[[np.ndarray], bool],
    ) -> Callable[[Iterable], Iterable]:
    '''Returns a function that will be used by the mcmc algorithm to propose
    a new parameter value tuple based on a given one. 
    The two arrays c_max and c_min define the boundaries
    of the n-dimensional rectengular domain for the parameters and must be of
    the same shape.  After a possible parameter value has been sampled the
    filter_func will be applied to it to either accept or discard it.  So
    filter func must accept parameter array and return either True or False'''
    
    g = np.random.default_rng()
    def GenerateParamValues(c_op):
        paramNum = len(c_op)
        keep_searching = True
        while (keep_searching):
            c_new = c_op + (g.random((paramNum)) - 0.5)*(c_max - c_min)/D
            if (filter_func(c_new)):
                keep_searching = False
        return c_new

    return GenerateParamValues

def make_multivariate_normal_proposer(
        covv: np.ndarray,
        filter_func: Callable[[Iterable], bool],
    ) -> Callable[[Iterable], Iterable]:
    """Returns a function that will be used by mcmc algorithm to propose
    a new parameter(tuple) based on a given one.
    param: covv: The covariance matrix (usually estimated from a previously run chain)
    """

    def GenerateParamValues(c_op):
        flag = True
        while (flag):
           c_new = c_op + np.random.multivariate_normal(np.zeros(len(c_op)), covv)
           if (filter_func(c_new)):
              flag = False
        return c_new

    return GenerateParamValues


def mcmc(
        initial_parameters: Iterable,
        proposer: Callable[[Iterable],Iterable], 
        param2res: Callable[[np.ndarray], np.ndarray],
        costfunction: Callable[[np.ndarray],np.float64],
        nsimu: int 
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    performs the Markov chain Monte Carlo simulation an returns a tuple of the array of sampled parameter(tuples) with shape (len(initial_parameters),nsimu) and the array of costfunction values with shape (q,nsimu)

    :param initial_parameters: The initial guess for the parameter (tuple) to be estimated
    :param proposer: A function that proposes a new parameter(tuple) from a given parameter (tuple).
    :param param2res: A function that given a parameter(tuple) returns
    the model output, which has to be an array of the same shape as the observations used to
    build the costfunction.
    :param costfunction: A function that given a model output returns a real number. It is assumed to be created for a specific set of observations, which is why they do not appear as an argument.
    :param nsimu: The length of the chain
    """
    np.random.seed(seed=10)
    
    paramNum=len(initial_parameters)
    
    upgraded=0;
    C_op = initial_parameters
    tb=time()
    first_out = param2res(C_op)
    # 
    #print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' upgrade2 iworker ' + str(iworker) +\
    #					': ' + str(upgrade2) + ' out of ' + str(simu) + ' cost: ' + str(np.round(cost_new, 2)) +\
    #					' r2: soc ' + str(np.round(r2_soc, 2)) + ' litter ' + str(np.round(r2_litter, 2)) +\
    #					' cwd ' + str(np.round(r2_cwd, 2)) + ' hr ' + str(np.round(r2_hr, 2)) + \
    #					' bgoc ' + str(np.round(r2_bgoc, 2)))
    #			
    J_last = costfunction(first_out)
    print('first_iteration done after' + str(time()-tb))
    #J_last = 400 # original code
    
    # intialize the result arrays to the maximum length
    # Depending on many of the parameters will be accepted only 
    # a part of them will be filled with real values
    C_upgraded = np.zeros((paramNum, nsimu))
    J_upgraded = np.zeros((1, nsimu))
    
    #for simu in tqdm(range(nsimu)):
    st =time() 
    for simu in range(nsimu):
        c_new = proposer(C_op)
        out_simu = param2res(c_new)
        J_new = costfunction(out_simu)
    
        delta_J =  J_last - J_new;
        
        randNum = np.random.uniform(0, 1)
        if (min(1.0, np.exp(delta_J)) > randNum):
                C_op=c_new;
                J_last=J_new;
                C_upgraded[:,upgraded]=C_op;
                J_upgraded[:,upgraded]=J_last; 
                upgraded=upgraded+1;
        # print some metadata 
        # (This could be added to the outputfile later)
        
        if simu%10==0 or simu == (nsimu-1):
            print(
 """ 
#(upgraded): {n}
over all acceptance ratio till now: {r}% 
progress: {simu:05d}/{nsimu:05d} {pbs} {p:02d}%
time elapsed: {minutes:02d}:{sec:02d}.
""".format(
                n=upgraded,
                r=int(upgraded/(simu+1)*100),
                simu=simu,
                nsimu=nsimu,
                pbs='|'+int(50*simu/(nsimu-1))*'#'+int((1-simu/(nsimu-1))*50)*' '+'|',
                p=int(simu/(nsimu-1)*100),
                minutes=int((time()-st)/60),
                sec=int((time()-st)%60)
            ),
            end='\033[5A' # print alway on the same spot of the screen...
            )

    #remove the part of the arryas that is still filled with zeros
    useful_slice = slice(0,upgraded)
    return C_upgraded[:,useful_slice], J_upgraded[:,useful_slice]

def make_feng_cost_func(
        obs: np.ndarray,
    ) -> Callable[[np.ndarray],np.float64]:
    # Note:
    # in our code the dimension 0 is the time
    # and dimension 1 the pool index
    n = obs.shape[1]
    means = obs.mean(axis=0)
    mean_centered_obs = obs-means
    # now we compute a scaling factor per observable stream
    # fixme mm 10-28-2021
    #   The dominators in this case are actually the TEMPORAL variances of the data streams
    #   which i find weierd
    denominators=np.sum(mean_centered_obs**2,axis=0) 
    #   The desired effect of automatically adjusting weight could be achieved
    #   by the mean itself.
    # dominators = means
    def costfunction(mod: np.ndarray) ->np.float64:
        cost= np.mean(
                np.sum((obs - mod)**2,axis=0)/denominators*100
        ) 
        # this is equivalent to the following slightly less compact but slightly more readable code
        #cost = np.mean(
        #        [
        #            np.sum((obs[:,i] - mod[:,i])**2)/ denominators[i]*100 
        #            for i in range(n)
        #        ]
        #)
        return cost
    
    return costfunction     

def make_jon_cost_func(
        obs: np.ndarray,
    ) -> Callable[[np.ndarray],np.float64]:
    # Note:
    # in our code the dimension 0 is the time
    # and dimension 1 the pool index
    n = obs.shape[1]
    means = obs.mean(axis=0)
    denominators = means**2
    mean_centered_obs = obs-means

    def costfunction(mod: np.ndarray) ->np.float64:
        cost= np.mean(
                np.sum((obs - mod)**2,axis=0)/denominators 
        ) 
        return cost
    
    return costfunction     


def day_2_month_index(d):
    return months_by_day_arr()[(d%365)]

@lru_cache
def months_by_day_arr():
    days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    return  np.concatenate(
        tuple(
            map(
                lambda m: m*np.ones(
                    days_per_month[m],
                    dtype=np.int64
                ),
                range(12)
            )
        )
    )


def month_2_day_index(ns):
    """ computes the index of the day at the end of the month n in ns
    this works on vectors and is faster than a recursive version working
    on a single index (since the smaller indices are handled anyway)
    """

    # We first compute the sequence of day indices up to the highest month in ns
    # and then select from this sequence the day indices for the months in ns
    days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    dpm =  (days_per_month[i%len(days_per_month)] for i in range(max(ns)))
    # compute indices for which we want to store the results which is the
    # list of partial sums of the above list  (repeated)

    def f(acc,el):
        if len(acc) <1:
            res = (el,)
        else:
            last = acc[-1]
            res = acc + (el+last,)
        return res
    day_indices_for_continuous_moths = reduce(
        f,
        dpm,
        (0,) 
    )
    day_indices = reduce( 
        lambda acc,n: acc + [day_indices_for_continuous_moths[n]], #for n=0 we want 0
        ns,
        []
    )
    return day_indices

class TimeStepIterator2():
    """iterator for looping forward over the results of a difference equation
    X_{i+1}=f(X_{i},i)"""

    def __init__(
        self,
        initial_values, # a tupel of values that will be
        f, # the function to compute the next ts
        max_it = False
    ):
        self.initial_values = initial_values
        self.f= f
        self.reset()
        self.max_it = max_it

    def reset(self):
        self.i = 0
        self.ts = self.initial_values

    def __iter__(self):
        self.reset()
        return(self)

    def __next__(self):
        if self.max_it:
            if self.i == self.max_it:
                raise StopIteration

        ts = copy(self.ts)
        ts_new = self.f(self.i, ts)
        self.ts = ts_new
        self.i += 1
        return ts

    def values(self,day_indices):
        # we traverse the iterator to the highest index and
        # collect the results we want to keep in a list (acc)
        tsi=copy(self)
        tsi.reset()
        def g(acc, i):
            v = tsi.__next__()
            if i in day_indices:
                acc += [v] 
            return acc
        xs =  reduce(g,range(max(day_indices)+1),([]))
        return xs

    
def respiration_from_compartmental_matrix(B,X):
    """This function computes the combined respiration from all pools"""
    return -np.sum(B@X) 


def plot_solutions(
        fig,
        times, 
        var_names,
        tup,
        names=None
    ):
    if names is None:
        names = tuple(str(i) for i in range(len(tup)))

    assert(all([tup[0].shape == el.shape for el in tup]))

    if tup[0].ndim == 1:
        n_times = tup[0].shape[0]
        ax = fig.subplots(1,1)
        for i,sol in enumerate(tup):
            ax.plot(
                np.array(times).reshape(n_times,), 
                sol,
                marker="o",
                label=names[i]
            )
            ax.set_title(var_names[0])
            ax.legend()
    else:
        n_times, n_vars = tup[0].shape

        fig.set_figheight(n_vars*fig.get_figwidth())
        axs = fig.subplots(n_vars,1)
        colors =('red','blue','green','organge')
        for j in range(n_vars):
            for i,sol in enumerate(tup):
                axs[j].plot(
                    np.array(times).reshape(n_times,), 
                    sol[:, j],
                    marker="+",
                    label=names[i]
                )
                axs[j].set_title(var_names[j])
                axs[j].legend()
