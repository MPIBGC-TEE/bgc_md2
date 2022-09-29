from sympy import Matrix, symbols

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

# one wood product pool
WP_S, WP_L = symbols("WP_S WP_L")

state_vector = Matrix([WP_S, WP_L])

WP_S_input, WP_L_input = symbols("WP_S_input WP_L_input")
input_fluxes = {
    WP_S: WP_S_input,
    WP_L: WP_L_input
}

internal_fluxes = dict()

k_S, k_L = symbols("k_S k_L")
output_fluxes = {
    WP_S: k_S * WP_S,
    WP_L: k_L * WP_L
}

t = symbols("t")
srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    state_vector,
    t,
    input_fluxes,
    output_fluxes,
    internal_fluxes
)

mvs = CMTVS(
    {
        InFluxesBySymbol(input_fluxes),
        OutFluxesBySymbol(output_fluxes),
        InternalFluxesBySymbol(internal_fluxes),
        TimeSymbol(t.name),
        StateVariableTuple(state_vector)
    },
    bgc_md2_computers()
)

