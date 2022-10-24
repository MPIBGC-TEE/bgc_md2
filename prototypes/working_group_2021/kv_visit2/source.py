from sympy import Symbol, Function, exp, Piecewise
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.models.BibInfo import BibInfo
import source_1 as s1
from subs_1 import subs_dict

def subs_xi(var):
    return var.subs(subs_dict)

mvs=CMTVS(
    {
        s1.mvs.get_TimeSymbol(),
        s1.mvs.get_Temperature(),
        s1.mvs.get_StateVariableTuple(),
        s1.mvs.get_CarbonStateVariableTuple(),
        s1.mvs.get_VegetationCarbonStateVariableTuple(),
        s1.mvs.get_SoilCarbonStateVariableTuple(),
        subs_xi(s1.mvs.get_InFluxesBySymbol()),
        subs_xi(s1.mvs.get_OutFluxesBySymbol()),
        subs_xi(s1.mvs.get_InternalFluxesBySymbol()),
    },
    computers=s1.mvs.computers
)    
