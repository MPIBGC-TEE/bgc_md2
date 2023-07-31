from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (InFluxesBySymbol, InternalFluxesBySymbol,
                                   OutFluxesBySymbol, StateVariableTuple,
                                   TimeSymbol)
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from ComputabilityGraphs.CMTVS import CMTVS
from sympy import Matrix

from .. import (B_OS, B_TH, B_TS, C_S, G_TS, GPP, MS, C_gHW, C_gW, E, delta_W,
                eta_W, f_T, t, v_T, zeta_dw, zeta_gluc)

# trunk
trunk_sv_set = set([E, C_S, B_TS, B_TH])

trunk_in_fluxes = {E: GPP}
trunk_out_fluxes = {E: MS + G_TS}

trunk_internal_fluxes = {
    # sapwood production
    (E, C_S): f_T * delta_W / (C_gW + delta_W) * E,
    (E, B_TS): f_T * C_gW / (C_gW + delta_W) * eta_W * E,
    # trunk heartwood production
    (B_TS, B_TH): v_T * B_TS,
    (C_S, B_TH): v_T * 1 / C_gHW * B_TS / (B_OS + B_TS) * zeta_dw / zeta_gluc * C_S,
}

trunk = (trunk_sv_set, trunk_in_fluxes, trunk_out_fluxes, trunk_internal_fluxes)

srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    Matrix([E, C_S, B_TH, B_TS]),
    t,
    trunk_in_fluxes,
    trunk_out_fluxes,
    trunk_internal_fluxes,
)

mvs = CMTVS(
    {
        InFluxesBySymbol(trunk[1]),
        OutFluxesBySymbol(trunk[2]),
        InternalFluxesBySymbol(trunk[3]),
        TimeSymbol("t"),
        StateVariableTuple((E, C_S, B_TH, B_TS)),
    },
    bgc_md2_computers(),
)
