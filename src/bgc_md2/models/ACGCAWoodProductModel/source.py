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
Wood_products = symbols("Wood-products")

state_vector = Matrix([Wood_products])

WP_input = symbols("WP_input")
input_fluxes = {
    Wood_products: WP_input
}

internal_fluxes = dict()

alpha_WP = symbols("alpha_WP")
output_fluxes = {
    Wood_products: alpha_WP * Wood_products
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

