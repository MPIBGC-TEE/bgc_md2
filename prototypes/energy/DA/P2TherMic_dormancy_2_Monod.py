from sympy import Symbol, symbols, Function, Tuple, lambdify, exp, Abs, Matrix, diff
import P2TherMic_dormancy_2_abstract as mod 
#from P2TherMic_dormancy_2_abstract import state_var_tuple
from collections import namedtuple
import numpy as np

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
        "alpha",
        "S_T", 
        "Delta_C", 
        "H_X", 
        "H_S", 
        "H_NH3", 
    ]
)
consts=namedtuple(
    "consts",
    [ "R","XN","CO2","q" ]
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
    cp=const(**{"R":0,"XN":0,"CO2":0,"q":0})
    ep=estimated_parameters(*params)
    all_params={
        **cp._asdict(),
        **estimated_parameters._asdict()
    }
    y0=mod.state_var_tuple.subs(all_params)
    sol=solve_ivp(
        fun=lambdify(
            (
                mod.t,
                mod.state_var_tuple
            ),
            mod.rhs_sym.subs(mod.par_dict),"numpy"),
        t_span=(0,500),
        y0=mod.state_var_tuple.subs(mod.start_value_dict),
        dense_output=True
    ).sol
    
    plot_times=np.linspace(sol.t_min,sol.t_max,100)
    values=sol(plot_times)

    #obs = run_HH_model(params)
    #summstats = torch.as_tensor(calculate_summary_statistics(obs))
    #return summstats
    return 5

#start_value_dict={
#    R: 0,
#    S: 2.0,
#    X: 0.5,
#    XN: 0,
#    CO2: 0,
#    q: 0,
#}
#par_dict = {
#        mu_max: 0.2,
#        m_max: 0.02,
#        K_S: 8,
#        V_h: 8,
#        Y_X_S: 0.59,
#        Y_X_N: 0.2,
#        rho_B: 1.2,
#        Theta: 0.3,
#        alpha: 0.1,
#        S_T: 1, 
#        Delta_C: 1, # fixme mm: unknown
#        H_X: 1, # fixme mm: unknown
#        H_S: 1, # fixme mm: unknown
#        H_NH3: 1, # fixme mm: unknown
#}        
#rhs_sym=mod.rhs_sym.subs(sd_Monod)
