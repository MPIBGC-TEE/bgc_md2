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
        'C_ssl': 'structural soil surface litter pool'
        ,'C_msl': 'metabolic soil surface litter pool'
        ,'C_sl': 'structural soil litter pool'
        ,'C_ml': 'metabolic soil litter pool'
        ,'C_asom': 'active soil organic matter pool'
        ,'C_ssom': 'slow soil organic matter pool'
        ,'C_7': 'passive soil organic matter pool'
        ,'K_ssl': 'maximum decomposition rate of structural soil surface litter'
        ,'K_msl': 'maximum decomposition rate of metabolic soil surface litter'
        ,'K_sl': 'maximum decomposition rate of structural soil litter'
        ,'K_ml': 'maximum decomposition rate of metabolic soil litter'
        ,'K_asom': 'maximum decomposition rate of active soil organic matter'
        ,'K_ssom': 'maximum decomposition rate of slow soil organic matter'
        ,'K_7': 'maximum decomposition rate of passive organic matter'
        ,'k_ssl': 'decomposition rate of structural soil surface litter'
        ,'k_sl': 'decomposition rate of structural soil litter'
        ,'k_asom': 'decomposition rate of active soil organic matter'
        ,'LN': 'lignin-to-nitrogen ratio'
        ,'Ls': 'fraction of structural material that is lignin'
        ,'Tx': 'silt and clay fraction of the soil'
        ,'A_l': 'lignin fraction that is composed in structural pools (equals Ls?)'
        ,'E_s': '"fraction of carbon lost as CO$_2$ when active soil organic matter is decomposed and stabilized into slow organic matter"'
        ,'F_m': 'fraction of incoming metabolic litter'
        ,'F_s': 'fraction of incoming structural litter'
        ,'J_1': 'organic matter input to surface'
        ,'J_2': 'organic matter input to soil'
        ,'alpha_asom_from_ssl': 'flux coefficient from strucutral soil surface litter pool to active soil organic matter pool'
        ,'alpha_asom_from_sl': 'flux coefficient from strucutral soil litter pool to active soil organic matter pool'
        ,'alpha_ssom1': 'flux coefficient from strucutral soil surface litter pool to slow soil organic matter pool'
        ,'alpha_ssom3': 'flux coefficient from strucutral soil litter pool to slow soil organic matter pool'
        ,'alpha_ssom5': 'flux coefficient from strucutral soil surface litter pool to slow soil organic matter pool'
        ,'f_T': 'function of temperature'
        ,'f_W': 'function of soil moisture'
}

for name in sym_dict.keys():
    var(name)
k_ssl = K_ssl * exp(-3 * Ls)
k_sl = K_sl * exp(-3 * Ls)
k_asom = K_asom * (1 - 0.75 * Tx)
E_s = 0.85 - 0.68 * Tx
F_m = 0.85 - 0.018 * LN
F_s = 1 - F_m
alpha_asom_from_ssl = 0.55 * (1 - A_l)
alpha_asom_from_sl = 0.45 * (1 - A_l)
alpha_ssom1 = 0.7 * A_l
alpha_ssom3 = alpha_ssom1
alpha_ssom5 = 1 - E_s - 0.004
t = TimeSymbol("t") # unit: "year" #??? monthly
x = StateVariableTuple((C_ssl, C_msl, C_sl, C_ml, C_asom, C_ssom, C_7))
u = InputTuple((J_1*F_s, J_1*F_m, J_2*F_s, J_2*F_m, 0, 0, 0))
xi = f_T * f_W #environmental effects multiplier (DEFAG)
A = ([        -k_ssl,        0,            0,        0,            0,        0,        0],
     [           0,     -K_msl,            0,        0,            0,        0,        0],
     [           0,        0,         -k_sl,        0,            0,        0,        0],
     [           0,        0,            0,     -K_ml,            0,        0,        0],
     [alpha_asom_from_ssl*k_ssl, 0.45*K_msl, alpha_asom_from_sl*k_sl, 0.45*K_ml,         -k_asom, 0.42*K_ssom, 0.45*K_7],
     [alpha_ssom1*k_ssl,        0, alpha_ssom3*k_sl,        0, alpha_ssom5*k_asom,     -K_ssom,        0],
     [           0,        0,            0,        0,    0.004*k_asom, 0.03*K_ssom,     -K_7])
B = CompartmentalMatrix(xi * ImmutableMatrix(A))
#Original values without effects of temperature and soil moisture:
np1 = NumericParameterization(
    par_dict={
K_ssl: 0.076, K_msl: 0.28, K_sl: 0.094, K_ml: 0.35, K_asom: 0.14, K_ssom: 0.0038, K_7: 0.00013, xi: 1},
    func_dict=frozendict({})
)
nsv1 = NumericStartValueDict({
C_ssl: 100, C_msl: 200, C_sl: 00, C_ml: 0, C_asom: 0, C_ssom: 0, C_7: 0}) #faked values by Markus to make it run
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
        # fixme mm 01-20-2022: The parameterization is incompomplete
        # error: The following free symbols: {A_l, Ls, Tx} of the expression: {A_l, Tx, Ls} 
        # np1, 
        nsv1,
    #    ntimes
    },
    bgc_md2_computers()
)
