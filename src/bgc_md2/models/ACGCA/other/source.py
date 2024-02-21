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
    E, C_S, B_OH, B_OS, B_TS,
    MS,
    G_OS_from_E, G_OS_from_CS,
    S_O, f_O, C_gW, delta_W, eta_W, v_O,
    C_gHW,
    zeta_dw, zeta_gluc,
    f_CS
)

# "other" means branches + coarse roots
other_sv_set = set([
    E,
    C_S,
    B_OS, B_OH,
    G_OS_from_CS
])

other_in_fluxes = {E: GPP}

other_out_fluxes = {
    E: MS + G_OS_from_E,
    B_OS: S_O * B_OS,
    B_OH: S_O * B_OH,
    C_S: G_OS_from_CS
}

other_internal_fluxes = {
    # sapwood production and senesence
    (E, C_S): f_O * delta_W/(C_gW+delta_W) * E,
    (C_S, E): S_O * B_OS/(B_OS+B_TS) * C_S,
    (E, B_OS): f_O * C_gW/(C_gW+delta_W) * eta_W * E,

    # use of reserves if necessary
    (C_S, B_OS): f_CS * eta_W * C_S,

    # coarse roots and branches heartwood production
    (B_OS, B_OH): v_O * B_OS,
    (C_S, B_OH): v_O * 1/C_gHW * B_OS/(B_OS+B_TS) * zeta_dw/zeta_gluc * C_S,
}

other = (
    other_sv_set,
    other_in_fluxes,
    other_out_fluxes,
    other_internal_fluxes
)

srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    Matrix([E, C_S, B_OH, B_OS]),
    t,
    other_in_fluxes,
    other_out_fluxes,
    other_internal_fluxes
)

mvs = CMTVS(
    {
        InFluxesBySymbol(other[1]),
        OutFluxesBySymbol(other[2]),
        InternalFluxesBySymbol(other[3]),
        TimeSymbol("t"),
        StateVariableTuple((E, C_S, B_OH, B_OS))
    },
    bgc_md2_computers()
)


