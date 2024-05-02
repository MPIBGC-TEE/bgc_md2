from sympy import Symbol, symbols, Function, Tuple, lambdify, exp, Abs, Matrix, diff
import P2TherMic_dormancy_2_abstract as mod 
from P2TherMic_dormancy_2_abstract import state_var_tuple


sym_dict={ 
    **mod.sym_dict,
    "mu_max": "",
    "rho_B" : "",
    "K_S": "",
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

sd_Monod={
        r_X: (mu_max * S * rho_B / Theta) / (K_S + S * rho_B / Theta) * X,
        r_M: (m_max * S * rho_B / Theta) / (K_S + S * rho_B / Theta) * X,
        r_N: (m_max * S * rho_B / Theta) / (K_S + S * rho_B / Theta) * X,  # fixme mm: r_N was not defined in the document just used r_M
        r_M_max: m_max * X,
        Phi: 1 / (
            (exp(S_T - S * rho_B/Theta)/(alpha * S_T)) + 1
        )
}
start_value_dict={
    R: 0,
    S: 2.0,
    X: 0.5,
    XN: 0,
    CO2: 0,
    q: 0,
}
par_dict = {
        mu_max: 0.2,
        m_max: 0.02,
        K_S: 8,
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
}        
rhs_sym=mod.rhs_sym.subs(sd_Monod)
