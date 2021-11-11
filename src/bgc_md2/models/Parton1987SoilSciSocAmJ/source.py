import numpy as np
from sympy import var, ImmutableMatrix, exp
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    NumericParameterization,
    NumericStartValueDict,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict = {
        'C_1': 'structural soil surface litter pool'
        ,'C_2': 'metabolic soil surface litter pool'
        ,'C_3': 'structural soil litter pool'
        ,'C_4': 'metabolic soil litter pool'
        ,'C_5': 'active soil organic matter pool'
        ,'C_6': 'slow soil organic matter pool'
        ,'C_7': 'passive soil organic matter pool'
        ,'K_1': 'maximum decomposition rate of structural soil surface litter'
        ,'K_2': 'maximum decomposition rate of metabolic soil surface litter'
        ,'K_3': 'maximum decomposition rate of structural soil litter'
        ,'K_4': 'maximum decomposition rate of metabolic soil litter'
        ,'K_5': 'maximum decomposition rate of active soil organic matter'
        ,'K_6': 'maximum decomposition rate of slow soil organic matter'
        ,'K_7': 'maximum decomposition rate of passive organic matter'
        ,'k_1': 'decomposition rate of structural soil surface litter'
        ,'k_3': 'decomposition rate of structural soil litter'
        ,'k_5': 'decomposition rate of active soil organic matter'
        ,'LN': 'lignin-to-nitrogen ratio'
        ,'Ls': 'fraction of structural material that is lignin'
        ,'Tx': 'silt and clay fraction of the soil'
        ,'A_l': 'lignin fraction that is composed in structural pools (equals Ls?)'
        ,'E_s': '"fraction of carbon lost as CO$_2$ when active soil organic matter is decomposed and stabilized into slow organic matter"'
        ,'F_m': 'fraction of incoming metabolic litter'
        ,'F_s': 'fraction of incoming structural litter'
        ,'J_1': 'organic matter input to surface'
        ,'J_2': 'organic matter input to soil'
        ,'alpha_51': 'flux coefficient from strucutral soil surface litter pool to active soil organic matter pool'
        ,'alpha_53': 'flux coefficient from strucutral soil litter pool to active soil organic matter pool'
        ,'alpha_61': 'flux coefficient from strucutral soil surface litter pool to slow soil organic matter pool'
        ,'alpha_63': 'flux coefficient from strucutral soil litter pool to slow soil organic matter pool'
        ,'alpha_65': 'flux coefficient from strucutral soil surface litter pool to slow soil organic matter pool'
        ,'f_T': 'function of temperature'
        ,'f_W': 'function of soil moisture'
}

for name in sym_dict.keys():
    var(name)
k_1 = K_1 * exp(-3 * Ls)
k_3 = K_3 * exp(-3 * Ls)
k_5 = K_5 * (1 - 0.75 * Tx)
E_s = 0.85 - 0.68 * Tx
F_m = 0.85 - 0.018 * LN
F_s = 1 - F_m
alpha_51 = 0.55 * (1 - A_l)
alpha_53 = 0.45 * (1 - A_l)
alpha_61 = 0.7 * A_l
alpha_63 = alpha_61
alpha_65 = 1 - E_s - 0.004
t = TimeSymbol("t") # unit: "year" #??? monthly
x = StateVariableTuple((C_1, C_2, C_3, C_4, C_5, C_6, C_7))
u = InputTuple((J_1*F_s, J_1*F_m, J_2*F_s, J_2*F_m, 0, 0, 0))
xi = f_T * f_W #environmental effects multiplier (DEFAG)
A = ([        -k_1,        0,            0,        0,            0,        0,        0],
     [           0,     -K_2,            0,        0,            0,        0,        0],
     [           0,        0,         -k_3,        0,            0,        0,        0],
     [           0,        0,            0,     -K_4,            0,        0,        0],
     [alpha_51*k_1, 0.45*K_2, alpha_53*k_3, 0.45*K_4,         -k_5, 0.42*K_6, 0.45*K_7],
     [alpha_61*k_1,        0, alpha_63*k_3,        0, alpha_65*k_5,     -K_6,        0],
     [           0,        0,            0,        0,    0.004*k_5, 0.03*K_6,     -K_7])
B = CompartmentalMatrix(xi * ImmutableMatrix(A))
#Original values without effects of temperature and soil moisture:
np1 = NumericParameterization(
    par_dict={
K_1: 0.076, K_2: 0.28, K_3: 0.094, K_4: 0.35, K_5: 0.14, K_6: 0.0038, K_7: 0.00013, xi: 1},
    func_dict=frozendict({})
)
nsv1 = NumericStartValueDict({
C_1: 100, C_2: 200, C_3: 00, C_4: 0, C_5: 0, C_6: 0, C_7: 0}) #faked values by Markus to make it run
#ntimes = NumericSimulationTimes(np.arange())

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="Century",
            longName="", 
            version="1",
            entryAuthor="Holger Metzler",
            entryAuthorOrcid="0000-0002-8239-1601",
            entryCreationDate="10/03/2016",
            doi="10.2136/sssaj1987.03615995005100050015x",
            sym_dict=sym_dict
        ),
        B,  # the overall compartmental matrix
        u,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
        np1,
        nsv1,
    #    ntimes
    },
    bgc_md2_computers()
)
