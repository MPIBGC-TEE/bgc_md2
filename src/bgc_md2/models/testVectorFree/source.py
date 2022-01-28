from sympy import symbols
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)

I_vl, I_vw = symbols("I_vl I_vw")
vl, vw = symbols("vl vw")
k_vl, k_vw = symbols("k_vl k_vw")


mvs = CMTVS(
    {
        StateVariableTuple((vl, vw)),
        TimeSymbol("t"),
        InFluxesBySymbol({vl: I_vl, vw: I_vw}),
        OutFluxesBySymbol({vl: k_vl * vl, vw: k_vw * vw}),
        InternalFluxesBySymbol({(vl, vw): k_vl * vl, (vw, vl): k_vw * vw}),
    },
    bgc_md2_computers()

)
