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

# two litter pools (maybe fast and slow)
Litter, CWD = symbols("Litter CWD")

# one SOM pool
SOM = symbols("SOM")

state_vector = Matrix([Litter, CWD, SOM])

L1, L2 = symbols("L1 L2")
input_fluxes = {
    Litter: L1,
    CWD: L2
}

alpha_1, alpha_2 = symbols("alpha_1 alpha_2")
internal_fluxes = {
    (Litter, SOM): alpha_1 * Litter,
    (CWD, SOM): alpha_2 * CWD
}

r_Litter, r_CWD, r_SOM = symbols("r_Litter r_CWD r_SOM ")
output_fluxes = {
    Litter: r_Litter * Litter,
    CWD: r_CWD * CWD,
    SOM: r_SOM * SOM
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

