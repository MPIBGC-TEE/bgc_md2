from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (InFluxesBySymbol, InternalFluxesBySymbol,
                                   OutFluxesBySymbol, StateVariableTuple,
                                   TimeSymbol)
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from ComputabilityGraphs.CMTVS import CMTVS
from sympy import Matrix, symbols

# litter and coarse woody debris pools
Litter = symbols("Litter")

# dead-wood classes
DWC_1, DWC_2, DWC_3, DWC_4, DWC_5, DWC_6 = symbols(
    "DWC_1 DWC_2 DWC_3 DWC_4 DWC_5 DWC_6"
)

# one SOC pool
SOC = symbols("SOC")

state_vector = Matrix([Litter, DWC_1, DWC_2, DWC_3, DWC_4, DWC_5, DWC_6, SOC])

u_Litter, u_1, u_2, u_3, u_4, u_5, u_6 = symbols("u_Litter u_1 u_2 u_3 u_4 u_5 u_6")
input_fluxes = {
    Litter: u_Litter,
    DWC_1: u_1,
    DWC_2: u_2,
    DWC_3: u_3,
    DWC_4: u_4,
    DWC_5: u_5,
    DWC_6: u_6,
}

# total turnover rate
k_Litter, k_SOC = symbols("k_Litter k_SOC")
k_1, k_2, k_3, k_4, k_5, k_6 = symbols("k_1 k_2 k_3 k_4 k_5 k_6")

# respiration fractions
f_Litter = symbols("f_Litter")
f_1, f_2, f_3, f_4, f_5, f_6 = symbols("f_1 f_2 f_3 f_4 f_5 f_6")

# transition rates
b_Litter = (1 - f_Litter) * k_Litter
b_1 = (1 - f_1) * k_1
b_2 = (1 - f_2) * k_2
b_3 = (1 - f_3) * k_3
b_4 = (1 - f_4) * k_4
b_5 = (1 - f_4) * k_5
b_6 = (1 - f_4) * k_6

internal_fluxes = {
    (Litter, SOC): b_Litter * Litter,
    (DWC_1, SOC): b_1 * DWC_1,
    (DWC_2, SOC): b_2 * DWC_2,
    (DWC_3, SOC): b_3 * DWC_3,
    (DWC_4, SOC): b_4 * DWC_4,
    (DWC_5, SOC): b_5 * DWC_5,
    (DWC_6, SOC): b_6 * DWC_6,
}

# respiration rates
r_Litter = f_Litter * k_Litter
r_SOC = 1.0 * k_SOC
r_1 = f_1 * k_1
r_2 = f_2 * k_2
r_3 = f_3 * k_3
r_4 = f_4 * k_4
r_5 = f_5 * k_5
r_6 = f_6 * k_6

output_fluxes = {
    Litter: r_Litter * Litter,
    DWC_1: r_1 * DWC_1,
    DWC_2: r_2 * DWC_2,
    DWC_3: r_3 * DWC_3,
    DWC_4: r_4 * DWC_4,
    DWC_5: r_5 * DWC_5,
    DWC_6: r_6 * DWC_6,
    SOC: r_SOC * SOC,
}

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
