from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (InFluxesBySymbol, InternalFluxesBySymbol,
                                   OutFluxesBySymbol, StateVariableTuple,
                                   TimeSymbol)
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from ComputabilityGraphs.CMTVS import CMTVS
from sympy import Matrix

from .. import B_L, C_L, GL, GPP, ML, S_L, C_gL, E, delta_L, eta_L, f_L, t

leaves_sv_set = set(
    [
        E,
        B_L,
        C_L,
    ]
)

leaves_in_fluxes = {E: GPP}

leaves_out_fluxes = {
    E: ML + GL,
    B_L: S_L * B_L,
}

leaves_internal_fluxes = {
    (E, B_L): f_L * C_gL / (C_gL + delta_L) * eta_L * E,
    (E, C_L): f_L * delta_L / (C_gL + delta_L) * E,
    (C_L, E): S_L * C_L,
}

leaves = (leaves_sv_set, leaves_in_fluxes, leaves_out_fluxes, leaves_internal_fluxes)

srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    Matrix([E, B_L, C_L]),
    t,
    leaves_in_fluxes,
    leaves_out_fluxes,
    leaves_internal_fluxes,
)

mvs = CMTVS(
    {
        InFluxesBySymbol(leaves[1]),
        OutFluxesBySymbol(leaves[2]),
        InternalFluxesBySymbol(leaves[3]),
        TimeSymbol("t"),
        StateVariableTuple((E, B_L, C_L)),
    },
    bgc_md2_computers(),
)
