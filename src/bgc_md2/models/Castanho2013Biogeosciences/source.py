from sympy import var, symbols, Symbol, ImmutableMatrix, diag, exp, Rational
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
)
from bgc_md2.resolve.MVarSet import MVarSet

sym_dict={
        'C_f':      'Carbon in foliage',
        'C_r':      'Carbon in roots',
        'C_w':      'Carbon in woody tissue',
        'NPP_i': 'Net Primary Production for PFT$_i$',
        'a_il':  'Fraction of annual NPP allocated to leaves for PFT$_i$',
        'a_is':  'Fraction of annual NPP allocated to stem for PFT$_i$',
        'a_ir':  'Fraction of annual NPP allocated to roots for PFT$_i$',
        'tau_il': 'Residence time of carbon in leaves for PFT$_i$  ',
        'tau_is': 'Residence time of carbon in stem for PFT$_i$  ',
        'tau_ir': 'Residence time of carbon in roots for PFT$_i$ ',
}
for name in sym_dict.keys():
    var(name)

x = StateVariableTuple((C_f, C_r, C_w ))
u = NPP_i
b = (a_il, a_is, a_ir)
Input = InputTuple(u*ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-1/tau_il, 1/-tau_is, 1/-tau_ir)
)
t = TimeSymbol("t")

mvs = MVarSet({
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the the complete system
})
