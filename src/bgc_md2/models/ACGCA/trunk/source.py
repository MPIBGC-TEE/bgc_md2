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
    E, C_S, B_TH, B_TS, B_OS,
    MS, GS_T,
    f_T, C_gW, delta_W, eta_W, v_T,
    C_gHW,
    zeta_dw, zeta_gluc
)

# trunk
trunk_sv_set = set([
    E,
    C_S,
    B_TS, B_TH
])

trunk_in_fluxes = {E: GPP}
trunk_out_fluxes = {E: MS + GS_T}

trunk_internal_fluxes = {
    # sapwood production
    (E, C_S): f_T * delta_W/(C_gW+delta_W) * E,
    (E, B_TS): f_T * C_gW/(C_gW+delta_W) * eta_W * E,
    
    # trunk heartwood production
    (B_TS, B_TH): v_T * B_TS,
    (C_S, B_TH): v_T * 1/C_gHW * B_TS/(B_OS+B_TS) * zeta_dw/zeta_gluc * C_S
}

trunk = (
    trunk_sv_set,
    trunk_in_fluxes,
    trunk_out_fluxes,
    trunk_internal_fluxes
)

mvs = MVarSet({
    InFluxesBySymbol(trunk[1]),
    OutFluxesBySymbol(trunk[2]),
    InternalFluxesBySymbol(trunk[3]),
    TimeSymbol("t"),
    StateVariableTuple((E, C_S, B_TH, B_TS))
})



