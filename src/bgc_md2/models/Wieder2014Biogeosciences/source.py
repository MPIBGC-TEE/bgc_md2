import numpy as np
from sympy import var, ImmutableMatrix, exp
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
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet

sym_dict = {
        'LIT_m': 'Metabolic litter'
        ,'LIT_s': 'Structural litter'
        ,'MIC_r': 'Copiotrophic microbial biomass'
        ,'MIC_k': 'Oligotrophic microbial biomass'
        ,'SOM_p': 'Physically protected soil organic matter'
        ,'SOM_c': 'Chemically protected soil organic matter'
        ,'Vmax': '' 
        ,'Vmax_r1': ''
        ,'Vmax_r2': ''
        ,'Vmax_r3': ''
        ,'Vmax_r4': ''
        ,'Km': ''
        ,'Km_r1': ''
        ,'Km_r2': ''
        ,'Km_r3': ''
        ,'Km_r4': ''
        ,'F1': 'Flux from LIT_m to MIC_r'
        ,'F2': 'Flux from LIT_s to MIC_r'
        ,'F3': 'Flux from SOM_p to MIC_r'
        ,'F4': 'Flux from SOM_c to MIC_r'
        ,'F5': 'Flux from MIC_r to SOM_p'
        ,'F6': 'Flux from LIT_m to MIC_k'
        ,'F7': 'Flux from LIT_s to MIC_k'
        ,'F8': 'Flux from SOM_p to MIC_k'
        ,'F9': 'Flux from SOM_c to MIC_k'
        ,'F10': 'Flux from MIC_k to SOM_c'
        ,'S': 'Percent sand in soil' # percentage
}

for name in sym_dict.keys():
    var(name)
Vmax = exp(V_slope * T + V_int) * av
Vmax_r1 = Vmax * Vmod_r1 
Vmax_r2 = Vmax * Vmod_r2 
Vmax_r3 = Vmax * Vmod_r3 
Vmax_r4 = Vmax * Vmod_r4
Km = exp(K_slope * T + K_int) * ak 
Km_r1 = Km * Kmod_1 
Km_r2 = Km * Kmod_2 
Km_r3 = Km * Kmod_3 
Km_r4 = Km * Kmod_4 
F1 = MIC_r * Vmax_r1 * LIT_m /(Km_r1 + LIT_m)
F2 = MIC_r * Vmax_r2 * LIT_s /(Km_r2 + LIT_s)
F3 = MIC_r * Vmax_r3 * SOM_p /(Km_r3 + SOM_p)
F4 = MIC_r * Vmax_r4 * SOM_c /(Km_r4 + SOM_c)
F5 = MIC_r * tau_r
F6 = MIC_k * Vmax_k1 * LIT_m /(Km_k1 + LIT_m)
F7 = MIC_k * Vmax_k2 * LIT_s /(Km_k2 + LIT_s)
F8 = MIC_k * Vmax_k3 * SOM_p /(Km_k3 + SOM_p)
F9 = MIC_k * Vmax_k4 * SOM_c /(Km_k4 + SOM_c)
F10 = MIC_k * tau_k
t = TimeSymbol("t") # unit: "hour"
x = StateVariableTuple(())
u = InputTuple(())
B = CompartmentalMatrix(
)
np1 = NumericParameterization(
    par_dict={
},
    func_dict=frozendict({})
)

nsv1 = NumericStartValueDict({
})
ntimes = NumericSimulationTimes(np.arange())

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="MIMICS",
        longName="Microbial-Mineral Carbon Stabilization",
        version="1",
        entryAuthor="Carlos Sierra",
        entryAuthorOrcid="0000-0003-0009-4169",
        entryCreationDate="14/08/2018",
        doi="10.5194/bg-11-3899-2014",
        sym_dict=sym_dict
    ),
    B,  # the overall compartmental matrix
    u,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonStateVariableTuple(()),
    np1,
    nsv1,
    ntimes
})

