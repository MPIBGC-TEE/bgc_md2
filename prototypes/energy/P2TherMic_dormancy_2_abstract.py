from sympy import Symbol, symbols, Function, Tuple, lambdify, exp, Abs, Matrix, diff

sym_dict={
    "t": "time",    
    "S": "concentration of substrate [mol g^{-1}]",
    "X": "concentration of Biomass [mol g^{-1}]",
    "R": "Physiological state [1]",
    "XN": "concentration of Microbial necromass [mol g^{-1}]",
    "CO2": "concentration of CO2 [mol g^{-1}]",
    "q": "Heat {kJ  g^{-1}",
    "r_X": "microbial growth [mol g^{-1} s^{-1}]",
    "Phi": "",
    "Y_X_S": "",
    "Y_X_N": "", 
    "r_M": "maintenance flux [mol g^{-1} s^{-1}]",
    "r_M_max": "maximum maintenance flux [mol g^{-1} s^{-1}]",
    "r_N": "", # fixme mm: r_N was not defined in the document just used r_M
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
