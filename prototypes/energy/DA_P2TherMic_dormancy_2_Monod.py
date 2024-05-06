from sympy import Symbol, symbols, Function, Tuple, lambdify, exp, Abs, Matrix, diff
#from P2TherMic_dormancy_2_abstract import state_var_tuple
from collections import namedtuple
import numpy as np
import torch
from sbi import utils as utils
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)
from scipy.integrate import solve_ivp
import P2TherMic_dormancy_2_Monod as mod 

estimated_parameters=namedtuple(
    "estimated_parameters",
    [
        "S",
        "X",
        "mu_max",
        "m_max",
       # "K_S",
       # "V_h",
       # "Y_X_S",
       # "Y_X_N",
       # "rho_B",
       # "Theta",
       # "alpha",
       # "S_T", 
       # "Delta_C", 
       # "H_X", 
       # "H_S", 
       # "H_NH3", 
    ]
)
consts=namedtuple(
    "consts",
    [ 
        "K_S",
        "R",
		"XN",
		"CO2",
		"q",
        "Y_X_S",
        "Y_X_N",
        "rho_B",
        "Theta",
        "alpha",
        "S_T", 
        "Delta_C", 
        "H_X", 
        "H_S", 
        "H_NH3", 
    ]
)
#all_parameter={**estimated_parameters._
## create pardict
#par_dict={
#            Symbol(k)
def simulation_wrapper(params: np.array):
    """
    Returns summary statistics from conductance values in `params`.
    
    Summarizes the output of the HH simulator and converts it to `torch.Tensor`.
    """
    cp=consts(
        **{
            "K_S": 8,
            "R":0,
		    "XN":0,
		    "CO2":0,
		    "q":0,
            "Y_X_S": 0.59,
            "Y_X_N": 0.2,
            "rho_B": 1.2,
            "Theta": 0.3,
            "alpha": 0.1,
            "S_T": 1, 
            "Delta_C": 1, # fixme mm": unknown
            "H_X": 1, # fixme mm": unknown
            "H_S": 1, # fixme mm": unknown
            "H_NH3": 1, # fixme mm": unknown
        }
    
        )
    #params=torch.reshape(
    #    params,
    #    (len(estimated_parameters._fields),)
    #)
    #from IPython import embed;embed()
    ep=estimated_parameters(
        *params
    )
    all_params={
        **cp._asdict(),
        **ep._asdict()
    }
    par_dict={
        Symbol(k): v 
        for k,v in all_params.items()
        if Symbol(k) in mod.par_dict.keys()
    }

    y0=mod.state_var_tuple.subs(all_params)
    sol=solve_ivp(
        fun=lambdify(
            (
                mod.t,
                mod.state_var_tuple
            ),
            mod.rhs_sym.subs(par_dict),"numpy"),
        t_span=(0,500),
        y0=y0,
        dense_output=True
    ).sol
    
    plot_times=np.linspace(sol.t_min,sol.t_max,100)
    values=sol(plot_times)

    #obs = run_HH_model(params)
    #summstats = torch.as_tensor(calculate_summary_statistics(obs))
    #return summstats
    return torch.as_tensor(values.reshape(-1))


prior_ranges= {
    "S": (0.5, 100),
    "X": (0.1, 1),
    "mu_max": (0.05, 0.5),
    "m_max": (0.005, 0.05),
    #"K_S": (1, 100),
    #"V_h": (0.1, 100),
    #"Y_X_S": (0.2, 0.8),
    #"Y_X_N": (0.05, 0.8),
    #"rho_B": (1.1, 1.6),
    #"Theta": (0.05, 0.5),
    #"S_T": (0.01, 100), 
    #"Delta_C" : (0,2) , # fixme mm: unknown 
    #"H_X": (0,2) , # fixme mm: unknown  
    #"H_S": (0,2) , # fixme mm: unknown  
    #"H_NH3": (0,2) , # fixme mm: unknown  
}
prior_min = np.array(
    estimated_parameters(
        **{ 
            k: min
            for k,( min, max) in prior_ranges.items()
        }
    )
)
prior_max = np.array(
    estimated_parameters(
        **{
            k: max
            for k, (min, max) in prior_ranges.items()
        }
    )
)
prior=utils.torchutils.BoxUniform(
        low=torch.as_tensor(prior_min),
        high=torch.as_tensor(prior_max)
)
