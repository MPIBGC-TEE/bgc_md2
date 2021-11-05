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

# from sympy.vector import CoordSysND,express
# fixme mm:
# add this boilerplatecode automatically
# from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

# C=get_CooordSystem()
I_vl, I_vw = symbols("I_vl I_vw")
vl, vw = symbols("vl vw")
k_vl, k_vw = symbols("k_vl k_vw")


# the keys of the internal flux dictionary are tuples (source_pool,target_pool)

# srm:SmoothReservoirModel
# srm=SmoothReservoirModel.from_state_variable_indexed_fluxes(
#     in_fluxes
#    ,out_fluxes
#    ,internal_fluxes
# )


# specialVars = {
mvs = CMTVS(
    {
        InFluxesBySymbol({vl: I_vl, vw: I_vw}),
        OutFluxesBySymbol({vl: k_vl * vl, vw: k_vw * vw}),
        InternalFluxesBySymbol({(vl, vw): k_vl * vl, (vw, vl): k_vw * vw}),
        TimeSymbol("t"),
        StateVariableTuple((vl, vw))
        #srm
    },
    bgc_md2_computers()

)
