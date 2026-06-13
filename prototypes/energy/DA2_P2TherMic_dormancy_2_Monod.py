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
#from P2TherMic_dormancy_2_Monod import state_var_tuple

estimated_parameters=namedtuple(
    "estimated_parameters",
    [
        "S",
        "X",
        "mu_max",
        "m_max",
        "K_S",
        "V_h",
        "Y_X_S",
        "Y_X_N",
        "rho_B",
        "Theta",
        "S_T", 
        # "Delta_C", 
        # "H_X", 
        # "H_S", 
        # "H_NH3", 
    ]
)
consts=namedtuple(
    "consts",
    [ 
        #"S", 
        #"X",
        #"mu_max",
        #"m_max",
        #"K_S",
        "R",
		"XN",
		"CO2",
		"q",
        #"Y_X_S",
        #"Y_X_N",
        #"rho_B",
        #"Theta",
        #"S_T", 
        "alpha",
        "Delta_C", 
        "H_X", 
        "H_S", 
        "H_NH3", 
    ]
)


def simulation_times(obs):
    # at the moment fake 
    t_min=0
    t_max=100
    n=50
    return np.linspace(t_min,t_max,n)

def simulation_wrapper(params: np.array):
    """
    Returns summary statistics from conductance values in `params`.
    
    Summarizes the output of the HH simulator and converts it to `torch.Tensor`.
    """
    cp=consts(
        **{
            #"K_S": 8,
            "R":0,
		    "XN":0,
		    "CO2":0,
		    "q":0,
            #"Y_X_S": 0.59,
            #"Y_X_N": 0.2,
            #"rho_B": 1.2,
            #"Theta": 0.3,
            #"S_T": 1, 
            "alpha": 0.1,
            "Delta_C": 1, # fixme mm": unknown
            "H_X": 1, # fixme mm": unknown
            "H_S": 1, # fixme mm": unknown
            "H_NH3": 1, # fixme mm": unknown
        }
    )
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
        t_span=(0,250),
        y0=y0,
        dense_output=True
    ).sol
    
    plot_times=simulation_times("fake")
    values=sol(plot_times)
    # sbi can only handle one dimensional results 
    # so we have to flatten the result.
    flat_values = values.reshape(-1)
    return torch.as_tensor(flat_values)


prior_ranges= {
    "S": (0.5, 100),
    "X": (0.1, 1),
    "mu_max": (0.05, 0.5),
    "m_max": (0.005, 0.05),
    "K_S": (1, 100),
    "V_h": (0.1, 100),
    "Y_X_S": (0.2, 0.8),
    "Y_X_N": (0.05, 0.8),
    "rho_B": (1.1, 1.6),
    "Theta": (0.05, 0.5),
    "S_T": (0.01, 100), 
    #"Delta_C" : (0,2) , # fixme mm: unknown 
    #"H_X": (0,2) , # fixme mm: unknown  
    #"H_S": (0,2) , # fixme mm: unknown  
    #"H_NH3": (0,2) , # fixme mm: unknown  
} 
ep0 = np.array( 
    estimated_parameters(
        **{
            k: {
                    **mod.start_value_dict,
                    **mod.par_dict
            }[Symbol(k)] 
            for k in estimated_parameters._fields
        }
    )
)
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
