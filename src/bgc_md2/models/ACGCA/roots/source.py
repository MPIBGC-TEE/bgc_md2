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

mvs = MVarSet({
    InFluxesBySymbol(roots[1]),
    OutFluxesBySymbol(roots[2]),
    InternalFluxesBySymbol(roots[3]),
    TimeSymbol("t"),
    StateVariableTuple((E, B_R, C_R))
})


