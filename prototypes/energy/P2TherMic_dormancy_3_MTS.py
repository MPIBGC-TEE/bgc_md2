from sympy import Symbol, symbols, Function, Tuple, lambdify, exp, Abs, Matrix, diff
import P2TherMic_dormancy_3_abstract as mod 
from P2TherMic_dormancy_3_abstract import state_var_tuple


sym_dict={ 
    **mod.sym_dict,
    "mu_max": "",
    "rho_B" : "",
    "Theta": "",
    "V_h": "",
    "m_max": "",
    "S_T": "",
    "alpha":"",

}
# Make symbols from  the strings that we can later use in expressions  
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

sd_MTS ={
    r_X: X*mu_max*exp(-Theta/(S*V_h*Y_X_S*rho_B)),
    r_M: X*m_max*exp(-Theta/(S*V_h*Y_X_S*rho_B)),
    r_N: X*m_max*exp(-Theta/(S*V_h*Y_X_S*rho_B)),
    r_M_max: X*m_max,
    Phi: 1/(1 + exp(-S*rho_B/Theta + S_T)/(S_T*alpha))
}
start_value_dict={
    R: 0,
    S: 2.0,
    S_S: 2.0, #fixme mm ??
    X: 0.5,
    XN: 0,
    CO2: 0,
    q: 0,
}
par_dict = {
        mu_max: 0.2,
        m_max: 0.02,
        V_h: 8,
        Y_X_S: 0.59,
        Y_X_N: 0.2,
        rho_B: 1.2,
        Theta: 0.3,
        alpha: 0.1,
        S_T: 1, 
        Delta_C: 1, # fixme mm: unknown
        H_X: 1, # fixme mm: unknown
        H_S: 1, # fixme mm: unknown
        H_NH3: 1, # fixme mm: unknown
        f_d: 2,
}        
rhs_sym=mod.rhs_sym.subs(sd_MTS)
