import numpy as np
from typing import Callable, Tuple

def make_uniform_proposer(
        c_max: np.ndarray,
        c_min: np.ndarray,
        filter_func: Callable[[np.ndarray], bool],
    ) -> Callable[[np.ndarray],np.ndarray]:
    '''Returns a function that will be used by the mcmc algorithm to propose
    new parameter values The two arrays c_max and c_min define the boundaries
    of the n-dimensional rectengular domain for the parameters and must be of
    the same shape.  After a possible parameter value has been sampled the
    filter_func will be applied to it to either accept or discard it.  So
    filter func must accept parameter array and return either True or False'''
    

    def GenerateParamValues(c_op):
        paramNum = len(c_op)
        flag = True
        while (flag):
           c_new = c_op + (np.random.random((paramNum)) - 0.5)*(c_max - c_min)/10.0
           #c_new = c_op + (np.random.normal(0, 1, paramNum))*(c_max - c_min)/15.0
           if (filter_func(c_new)):
              flag = False
        return c_new

    return GenerateParamValues

def mcmc(
        initial_parameters,
        proposer,
        forward_simulation: Callable[[Tuple], Tuple[np.ndarray, np.ndarray]],
        costfunction: Callable[[np.ndarray],np.float64]
    ) -> Tuple:
    
    np.random.seed(seed=10)
    
    paramNum=len(initial_parameters)
    nsimu    = 20000
    
    upgraded=0;
    J_last = 400
    C_op = initial_parameters
    
    C_upgraded = np.zeros((paramNum, nsimu))
    J_upgraded = np.zeros((1, nsimu))
    
    for simu in range(nsimu):
        c_new = proposer(C_op)
    
        out_simu = forward_simulation(c_new)
        J_new = costfunction(out_simu)
    
        delta_J =  J_last - J_new;
        
        randNum = np.random.uniform(0, 1)
        if (min(1.0, np.exp(delta_J)) > randNum):
                C_op=c_new;
                J_last=J_new;
                C_upgraded[:,upgraded]=C_op;
                J_upgraded[:,upgraded]=J_last; 
                upgraded=upgraded+1;
    
    return C_upgraded, J_upgraded

