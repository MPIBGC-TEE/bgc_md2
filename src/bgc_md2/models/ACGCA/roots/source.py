from sympy import Matrix

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)

from bgc_md2.models.ACGCA.__init__ import (
    t,
    GPP,
    E, B_R, C_R,
    MR, GR,
    S_R, f_R, C_gR, delta_R, eta_R
)


roots_sv_set = set([
    E,
    B_R, C_R,
])

roots_in_fluxes = {E: GPP}

roots_out_fluxes = {
    E: MR + GR,
    B_R: S_R * B_R,
}

roots_internal_fluxes = {
    (E, B_R): f_R * C_gR/(C_gR+delta_R) * eta_R * E,
    (E, C_R): f_R * delta_R/(C_gR+delta_R) * E,
    (C_R, E): S_R * C_R,
}

roots = (
    roots_sv_set,
    roots_in_fluxes,
    roots_out_fluxes,
    roots_internal_fluxes
)

srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    Matrix([E, B_R, C_R]),
    t,
    roots_in_fluxes,
    roots_out_fluxes,
    roots_internal_fluxes
)

mvs = CMTVS(
    {
        InFluxesBySymbol(roots[1]),
        OutFluxesBySymbol(roots[2]),
        InternalFluxesBySymbol(roots[3]),
        TimeSymbol("t"),
        StateVariableTuple((E, B_R, C_R))
    },
    bgc_md2_computers()
)


