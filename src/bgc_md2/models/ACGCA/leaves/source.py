from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)

from bgc_md2.models.ACGCA.__init__ import (
    GPP,
    E, B_L, C_L,
    ML, GL,
    S_L, f_L, C_gL, delta_L, eta_L
)


leaves_sv_set = set([
    E,
    B_L, C_L,
])

leaves_in_fluxes = {E: GPP}

leaves_out_fluxes = {
    E: ML + GL,
    B_L: S_L * B_L,
}

leaves_internal_fluxes = {
    (E, B_L): f_L * C_gL/(C_gL+delta_L) * eta_L * E,
    (E, C_L): f_L * delta_L/(C_gL+delta_L) * E,
    (C_L, E): S_L * C_L,
}

leaves = (
    leaves_sv_set,
    leaves_in_fluxes,
    leaves_out_fluxes,
    leaves_internal_fluxes
)

mvs = MVarSet({
    InFluxesBySymbol(leaves[1]),
    OutFluxesBySymbol(leaves[2]),
    InternalFluxesBySymbol(leaves[3]),
    TimeSymbol("t"),
    StateVariableTuple((E, B_L, C_L))
})

