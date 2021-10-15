from sympy import Matrix, symbols

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

from bgc_md2.resolve.MVarSet import MVarSet
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

r_S, r_L = symbols("r_S r_L")
output_fluxes = {
    WP_S: r_S * WP_S,
    WP_L: r_L * WP_L
}

t = symbols("t")
srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    state_vector,
    t,
    input_fluxes,
    output_fluxes,
    internal_fluxes
)

mvs = MVarSet({
    InFluxesBySymbol(input_fluxes),
    OutFluxesBySymbol(output_fluxes),
    InternalFluxesBySymbol(internal_fluxes),
    TimeSymbol(t.name),
    StateVariableTuple(state_vector)
})

