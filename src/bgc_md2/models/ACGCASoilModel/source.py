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

# two litter pools (maybe fast and slow)
Litter_1, Litter_2 = symbols("Litter_1 Litter_2")

# two SOM pools (maybe fast and slow)
Soil_1, Soil_2 = symbols("Soil_1 Soil_2")

state_vector = Matrix([Litter_1, Litter_2, Soil_1, Soil_2])

L1, L2 = symbols("L1 L2")
input_fluxes = {
    Litter_1: L1,
    Litter_2: L2
}

alpha_L1S1, alpha_L1S2 = symbols("alpha_L1S1, alpha_L1S2")
alpha_L2S1, alpha_L2S2 = symbols("alpha_L2S1, alpha_L2S2")
alpha_S1S2, alpha_S2S1 = symbols("alpha_S1S2, alpha_S2S1")
internal_fluxes = {
    (Litter_1, Soil_1): alpha_L1S1 * Litter_1,
#    (Litter_1, Soil_2): alpha_L1S2 * Litter_2,
    (Litter_2, Soil_1): alpha_L2S1 * Litter_1,
#    (Litter_2, Soil_2): alpha_L2S2 * Litter_2,
    (Soil_1, Soil_2): alpha_S1S2 * Soil_1,
#    (Soil_2, Soil_1): alpha_S2S1 * Soil_2
}

alpha_L1, alpha_L2 = symbols("alpha_L1, alpha_L2")
alpha_S1, alpha_S2 = symbols("alpha_S1, alpha_S2")
output_fluxes = {
    Litter_1: alpha_L1,
    Litter_2: alpha_L2,
    Soil_1: alpha_S1,
    Soil_2: alpha_S2
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

