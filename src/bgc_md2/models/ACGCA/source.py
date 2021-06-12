from sympy import Matrix

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems import helpers_reservoir as hr

from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)

from bgc_md2.models.ACGCA.leaves.source import leaves
from bgc_md2.models.ACGCA.roots.source import roots
from bgc_md2.models.ACGCA.other.source import other
from bgc_md2.models.ACGCA.trunk.source import trunk

from bgc_md2.models.ACGCA.__init__ import (
    t,
    GPP,
    E, B_L, C_L, B_OS, B_OH, C_S, B_TH, B_TS, C_R, B_R,
    ML, GL, MR, GR,
    MS, G_OS, G_TS
)

# system leaves + roots
intersect = (
    {E: GPP},
    {E: ML + GL + MR + GR}
)
LR = hr.combine(leaves, roots, {}, {}, intersect)

# other + trunk
intersect = (
    {E: GPP, C_S: 0},
    {E: MS + G_OS + G_TS, C_S: 0}
)

OT = hr.combine(other, trunk, {}, {}, intersect)

# complete model
intersect = (
    {E: GPP, C_S: 0},
    {E: ML + GL + MR + GR + MS + G_OS + G_TS, C_S: 0}
)

LROT = hr.combine(LR, OT, {}, {}, intersect)

srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    Matrix([E, B_L, C_L, B_OS, B_OH, C_S, B_TH, B_TS, C_R, B_R]),
    t,
    LROT[1],
    LROT[2],
    LROT[3]
)

mvs = MVarSet({
    InFluxesBySymbol(LROT[1]),
    OutFluxesBySymbol(LROT[2]),
    InternalFluxesBySymbol(LROT[3]),
    TimeSymbol("t"),
    StateVariableTuple((E, B_L, C_L, B_OS, B_OH, C_S, B_TH, B_TS, C_R, B_R))
})

