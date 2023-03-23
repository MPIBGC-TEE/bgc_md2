from sympy import Symbol, Function, exp, Piecewise, diff
from . import source_1 as s1
sym_dict={
    'T_0': 'critical temperature',
    'E': 'activation energy',
    'KM': 'critical moisture',
}

for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

func_dict={
    'TAS': 'air temperature',
    'mrso': 'soil moisture',
}

for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=s1.mvs.get_TimeSymbol()
#TAS_C=TAS(t)-273.15
TAS_C=s1.TAS-273.15

# approximate soil T at 20cm from air T (from https://doi.org/10.1155/2019/6927045)
TS = 0.748*TAS_C + 6.181 

xi = exp(E*(1/(10-T_0)-1/(TS-T_0))) * mrso(t)/(KM+mrso(t))
d_xi=diff(xi,s1.TAS)
##xi=Piecewise(
##    (exp(E*(1/(10-T_0)-1/(TS-T_0))) * mrso_f(t)/(KM+mrso_f(t)),TS > T_0),
##    (0,True)
##)
#
sym2func={s1.TAS: TAS(t)}

# the main part of this module is this little dictionary
# that makes it possible to exchange the symbolic xi with 
# the actual one.
subs_dict={
    s1.xi: xi.subs(sym2func),
    diff(s1.xi,s1.TAS): d_xi.subs(sym2func)
}
