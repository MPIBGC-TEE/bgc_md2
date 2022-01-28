from sympy import var, Symbol
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers    
from .source_Luo import mvs


from bgc_md2.resolve.mvars import (
        OutFluxesBySymbol,
        InternalFluxesBySymbol
)
mvs = CMTVS(
    {
        mvs.get_InFluxesBySymbol(),
        mvs.get_TimeSymbol(),
        mvs.get_StateVariableTuple(),
        OutFluxesBySymbol(
            {
                dp :
                Symbol("k_"+ str(dp)+ "_to_out")*dp
                for dp in mvs.get_OutFluxesBySymbol().keys()
            }
        ),        
        InternalFluxesBySymbol(
            {
                (dp,rp) : Symbol("k_"+ str(dp)+ "_to_" + str(rp))*dp
                for dp,rp in mvs.get_InternalFluxesBySymbol().keys()
            }
        )
    },
    bgc_md2_computers()
)
