from pathlib import Path
from sympy import Symbol, symbols, Function, Tuple, lambdify, exp, Abs, Matrix, diff
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

sym_dict={
    "t": "time",    
    "S": "concentration of substrate [mol g^{-1}]",
    "X": "concentration of Biomass [mol g^{-1}]",
    "R": "Physiological state [1]",
    "XN": "concentration of Microbial necromass [mol g^{-1}]",
    "CO2": "concentration of CO2 [mol g^{-1}]",
    "q": "Heat {kJ  g^{-1}",
    "r_X": "microbial growth |Monot [mol g^{-1} s^{-1}]",
    "mu_max": "",
    "rho_B" : "",
    "K_S": "",
    "Theta": "",
    "Phi": "",
    "Y_X_S": "",
    "Y_X_N": "", 
    "V_h": "",
    "r_M": "maintenance flux [mol g^{-1} s^{-1}]",
    "r_M_max": "maximum maintenance flux [mol g^{-1} s^{-1}]",
    "r_N": "", # fixme mm: r_N was not defined in the document just used r_M
    "m_max": "",
    "rho_s": "",
    "S_T": "",
    "alpha":"",
    "Delta_C": "",
    "H_S":  "",
    "H_X":  "",
    "H_NH3":  "",

}
# Make symbols from  the strings that we can later use in expressions  
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

#sym_func_dict={
#}
# Make symbols from  the strings that we can later use in expressions  
#for k in sym_func_dict.keys():
#    code=k+" = Function('{0}')".format(k)
#    exec(code)

state_var_tuple=Tuple(
        S,
        X,
        R,
        XN,
        CO2,
        q,
)    
rhs_sym=Tuple(
        - R*r_X - r_M,
        (Y_X_S*r_X + Y_X_N*r_N)*R-(r_M_max-r_M),
        Y_X_S*r_X*(Phi-R),
        -r_N*R+(r_M_max - r_M),
        ((1-Y_X_S)*r_X + (1-Y_X_N)*r_N) * R + r_M,
        (r_X*(Y_X_S*Delta_C*H_X - 0.143* Delta_C*H_NH3) - r_N * (1-Y_X_N)*Delta_C*H_X) * R -r_M*Delta_C * H_S
)       
sd_Monod={
        r_X: (mu_max * S * rho_B / Theta) / (K_S + S * rho_B / Theta) * X,
        r_M: (m_max * S * rho_B / Theta) / (K_S + S * rho_B / Theta) * X,
        r_N: (m_max * S * rho_B / Theta) / (K_S + S * rho_B / Theta) * X,  # fixme mm: r_N was not defined in the document just used r_M
        r_M_max: m_max * X,
        Phi: 1 / (
            (exp(S_T - S * rho_B/Theta)/(alpha * S_T)) + 1
        )
}
rhs_sym_Monod = rhs_sym.subs(sd_Monod)

start_value_dict={
    R: 0,
    S: 2.0,
    X: 0.5,
    XN: 0,
    CO2: 0,
    q: 0,
}
y0=state_var_tuple.subs(start_value_dict)
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

sol_Monod=solve_ivp(
    fun=lambdify((t,state_var_tuple),rhs_sym_Monod.subs(par_dict),"numpy"),
    t_span=(0,500),
    y0=y0,
    dense_output=True
).sol

plot_times=np.linspace(sol_Monod.t_min,sol_Monod.t_max,100)
fig=plt.figure()
ax=fig.add_subplot()
values=sol_Monod(plot_times)
for i in range(values.shape[0]):
    ax.plot(
        plot_times,
        values[i,:],
        label=str(state_var_tuple[i])
    ) 
ax.legend()

fig.savefig(f"{Path(__file__).stem}.pdf")
