from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (InFluxesBySymbol, InternalFluxesBySymbol,
                                   OutFluxesBySymbol, StateVariableTuple,
                                   TimeSymbol)
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from ComputabilityGraphs.CMTVS import CMTVS
from sympy import Matrix

from . import MR_C  # watch out! maintenance resp correction term
from . import (B_L, B_OH, B_OS, B_R, B_TH, B_TS, C_L, C_R, C_S, G_TS, GL, GPP,
               GR, ML, MR, MS, E, G_OS_from_CS, G_OS_from_E, t)
from .leaves.source import leaves
from .other.source import other
from .roots.source import roots
from .trunk.source import trunk

# system leaves + roots
intersect = ({E: GPP}, {E: ML + GL + MR + GR})
LR = hr.combine(leaves, roots, {}, {}, intersect)

# other + trunk
intersect = ({E: GPP, C_S: 0}, {E: MS + G_OS_from_E + G_TS, C_S: G_OS_from_CS})

OT = hr.combine(other, trunk, {}, {}, intersect)

# complete model
intersect = (
    {E: GPP, C_S: 0},
    {
        E: ML + GL + MR + GR + MS + G_OS_from_E + G_TS + MR_C,
        C_S: G_OS_from_CS,
    },  # MR_C = maintenance correction term
)

LROT = hr.combine(LR, OT, {}, {}, intersect)

srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    Matrix([E, B_L, C_L, B_OS, B_OH, C_S, B_TH, B_TS, C_R, B_R]),
    t,
    LROT[1],
    LROT[2],
    LROT[3],
)

mvs = CMTVS(
    {
        InFluxesBySymbol(LROT[1]),
        OutFluxesBySymbol(LROT[2]),
        InternalFluxesBySymbol(LROT[3]),
        TimeSymbol("t"),
        StateVariableTuple((E, B_L, C_L, B_OS, B_OH, C_S, B_TH, B_TS, C_R, B_R)),
    },
    bgc_md2_computers(),
)
