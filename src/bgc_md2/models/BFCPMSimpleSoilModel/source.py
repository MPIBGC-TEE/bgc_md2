from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (InFluxesBySymbol, InternalFluxesBySymbol,
                                   OutFluxesBySymbol, StateVariableTuple,
                                   TimeSymbol)
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from ComputabilityGraphs.CMTVS import CMTVS
from sympy import Matrix, symbols

# two litter pools (maybe fast and slow)
Litter, CWD = symbols("Litter CWD")

# one SOC pool
SOC = symbols("SOC")

state_vector = Matrix([Litter, CWD, SOC])

u_L, u_CWD = symbols("u_L u_CWD")
input_fluxes = {Litter: u_L, CWD: u_CWD}

# total turnover rate
k_Litter, k_CWD, k_SOC = symbols("k_Litter k_CWD k_SOC")

# respiration fractions
f_Litter, f_CWD = symbols("f_Litter f_CWD")

# transition rates
b_Litter = (1 - f_Litter) * k_Litter
b_CWD = (1 - f_CWD) * k_CWD

internal_fluxes = {(Litter, SOC): b_Litter * Litter, (CWD, SOC): b_CWD * CWD}

r_Litter = f_Litter * k_Litter
r_CWD = f_CWD * k_CWD
r_SOC = 1.0 * k_SOC

output_fluxes = {Litter: r_Litter * Litter, CWD: r_CWD * CWD, SOC: r_SOC * SOC}

t = symbols("t")
srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    state_vector, t, input_fluxes, output_fluxes, internal_fluxes
)

mvs = CMTVS(
    {
        InFluxesBySymbol(input_fluxes),
        OutFluxesBySymbol(output_fluxes),
        InternalFluxesBySymbol(internal_fluxes),
        TimeSymbol(t.name),
        StateVariableTuple(state_vector),
    },
    bgc_md2_computers(),
)
