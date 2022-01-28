from sympy import symbols, Function, exp, var, Piecewise
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers
from .source_minus_1 import mvs as mvs_base 

Bib = mvs_base.get_BibInfo()
for name in Bib.sym_dict.keys():
    var(name)
mrso=Function('mrso')
tsl =Function('tsl')
NPP = Function('NPP')
subs_dict ={xi: Piecewise(
    (exp(E * (1 / (10 - T0) - 1 / (tsl(t) - T0))) * mrso(t) / (KM + mrso(t)),tsl(t)>T0),
    #(0,tsl(t)<T0)
    (0,True) #catch all
    )
}

s = {
        mvs_base.get_InternalFluxesBySymbol().subs(subs_dict),
        mvs_base.get_OutFluxesBySymbol().subs(subs_dict)
}
mvs=mvs_base.update(s)
