import numpy as np
from sympy import var, ImmutableMatrix, Min
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonStateVariableTuple,
    NumericParameterization,
    NumericStartValueDict,
    NumericSimulationTimes,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict = {
        'W_p': 'Mass of leaf proteins'
        ,'W_s': 'Mass of leaf structural components'
        ,'W_r': 'Mass of roots'
        ,'W_C': 'Substrate carbon'
        ,'W_N': 'Substrate nitrogen'
        ,'I_dens': '(I) photon flux density. Photosynthetically active radiation.' # "[\\mu mol\\, m^{-2}\\,s^{-1}]"
        ,'h_max': 'leaf max. thickness' # "[m]"
        ,'h_half': '$h_0.5$ leaf half thickness'
        ,'rho': 'leaf density'
        ,'h': ''
        ,'A': 'Area'
        ,'N_a': 'Leaf nitrogen concentration' # "gN*m^{-2}"
        ,'C_i': 'Intercellular CO_2 concentration' # "[\\mu l\\, l^{-1}]"
        ,'V_cmax': 'Maximum carboxylation velocity' # "[\\mu mol CO_2\\, m^{-2}\\,s^{-1}]"
        ,'J_max': 'Maximum rate of electron transport' # "[\\mu Eq\\, m^{-2}\\,s^{-1}]"
        ,'J': 'Rate of electron transport'
        ,'R_d': 'Dark respiration rate' # "[\\mu mol CO_2\\, m^{-2}\\,s^{-1}]"
        ,'A_1': 'RuBP saturated portion of the carbon dioxide response curve' # "[\\mu mol CO_2\\, m^{-2}\\,s^{-1}]"
        ,'A_2': 'RuBP limited portion of the carbon dioxide response curve' # "[\\mu mol CO_2\\, m^{-2}\\,s^{-1}]"
        ,'sigma_c': 'Photosynthetic rate per unit leaf'
        ,'W_g': 'Plant biomass'
        ,'kappa': 'growth rate coefficient'
        ,'sigma_r': 'specific root activity'  # "[g N (g root)^{-1} d^{-1}]"
        ,'f_C': 'Proportion of carbon'
        ,'f_N': 'Proportion of nitrogen'
        ,'Beta': 'Target whole plant nitrogen:carbon ratio'
        ,'f_cp': 'Proportion of carbon in leaf proteins'
        ,'f_cs': 'Proportion of carbon in leaf structural components'
        ,'f_cr': 'Proportion of carbon in roots'
        ,'f_np': 'Proportion of nitrogen in leaf proteins'
        ,'f_ns': 'Proportion of nitrogen in leaf structural components'
        ,'f_nr': 'Proportion of nitrogen in roots'
        ,'C': 'Substrate carbon concentration'
        ,'N': 'Substrate nitrogen concentration'
        ,'P': ''
        ,'Q': ''
        ,'lambda_p': ''
        ,'lambda_s': ''
        ,'lambda_r': ''
}

for name in sym_dict.keys():
    var(name)

h = h_max*I_dens/(h_half+I_dens)
A = W_s/(rho*h)
N_a = (f_np*W_p)/A
V_cmax = (35.76*N_a)+12.42 #k_3=35.76, k_4=12.42
J_max = (92.55*N_a)+13.85
J = (J_max*I_dens)/(I_dens+(2.1*J_max))
R_d = (0.775*N_a)-0.238
A_1 = (V_cmax*((C_i-31)/(C_i+827)))-R_d
A_2 = (J*((C_i-31)/((4.5*C_i)+(10.5*31))))-R_d
sigma_c = Min(A_1,A_2)
W_g = W_p + W_s + W_r
f_N = (sigma_r*W_r*f_C)/(sigma_c*A)
C = W_C/W_g
N = W_N/W_g
P = f_C*sigma_r*W_r/(f_N*sigma_c*A)
Q = f_N/(Beta*f_C)
lambda_p = P/(1+P+Q)
lambda_s = Q/(1+P+Q)
lambda_r = 1/(1+P+Q)

x = StateVariableTuple((W_N, W_C, W_p, W_s, W_r))
u = InputTuple((sigma_c*A,sigma_r*W_r,0,0,0))
A = CompartmentalMatrix(
[[-((f_cp*lambda_p*W_p)+(f_cs*lambda_s*W_s)+(f_cr*lambda_r*W_r))*((kappa*W_N)/W_g**2), 0, 0, 0, 0],
                              [0, -((f_np*lambda_p*W_p)+(f_ns*lambda_s*W_s)+(f_nr*lambda_r*W_r))*((kappa*W_C)/W_g**2), 0, 0, 0],
                              [0, (kappa*N*lambda_p*W_p)/W_g, 0, 0, 0],
                              [0, (kappa*N*lambda_s*W_s)/W_g, 0, 0, 0],
                              [0, (kappa*N*lambda_r*W_r)/W_g, 0, 0, 0]])
t = TimeSymbol("t") # unit: "day"

np1 = NumericParameterization(
    par_dict={
sigma_r: 0.01, I_dens: 1000, f_C: 0.5, Beta: 0.05, f_cp: 0.60, f_cr: 0.45, f_cs: 0.45, f_np: 0.16, f_nr: 0.03, f_ns: 0.005, h_max: 0.0003, h_half: 200, kappa: 300, rho: '5*10**5', C_i: 240},
    func_dict=frozendict({})
)
nsv1 = NumericStartValueDict({
N: 0.01, C: 0.15, W_p: 100, W_s: 100, W_r: 100, W_C: 45, W_N: 30
})
ntimes = NumericSimulationTimes(np.arange(0, 150, 25))
#ntimes = NumericSimulationTimes(np.arange(0, 1, 0.001)) #Fixme: There were 2 in the yaml file, not sure which one works best

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="",
            longName="", 
            version="1",
            entryAuthor="Verónika Ceballos-Núñez",
            entryAuthorOrcid="0000-0002-0046-1160",
            entryCreationDate="29/7/2015",
            doi="10.1093/oxfordjournals.aob.a088273",
            sym_dict=sym_dict
        ),
        A,  # the overall compartmental matrix
        u,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
        VegetationCarbonStateVariableTuple((W_N, W_C, W_p, W_s, W_r)),
        np1,
        nsv1,
        ntimes
    },
    bgc_md2_computers()
)
