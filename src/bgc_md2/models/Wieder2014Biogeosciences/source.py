from sympy import var, exp
from bgc_md2.resolve.mvars import (
#    InFluxesBySymbol,
#    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

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
        ,'Vmax_k1': ''
        ,'Vmax_k2': ''
        ,'Vmax_k3': ''
        ,'Vmax_k4': ''
        ,'Km': ''
        ,'Km_r1': ''
        ,'Km_r2': ''
        ,'Km_r3': ''
        ,'Km_r4': ''
        ,'Km_k1': ''
        ,'Km_k2': ''
        ,'Km_k3': ''
        ,'Km_k4': ''
        ,'K_slope': ''
        ,'V_slope': ''
        ,'V_int': ''
        ,'K_int': ''
        ,'av': ''
        ,'ak': ''
        ,'T': ''
        ,'Vmod_r1': ''
        ,'Vmod_r2': ''
        ,'Vmod_r3': ''
        ,'Vmod_r4': ''
        ,'Kmod_1': ''
        ,'Kmod_2': ''
        ,'Kmod_3': ''
        ,'Kmod_4': ''
        ,'tau_r': ''
        ,'tau_k': ''
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
x = StateVariableTuple((LIT_m,LIT_s,MIC_r,MIC_k,SOM_p,SOM_c))

mvs = CMTVS(
    {
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
        TimeSymbol("t"), # unit: "hour"
        x,  # state vector of the complete system
    #    InFluxesBySymbol({})
    #    OutFluxesBySymbol({})
        InternalFluxesBySymbol({
            (LIT_m, MIC_r): F1,
            (LIT_s, MIC_r): F2,
            (SOM_p, MIC_r): F3,
            (SOM_c, MIC_r): F4,
            (MIC_r, SOM_p): F5,
            (LIT_m, MIC_k): F6,
            (LIT_s, MIC_k): F7,
            (SOM_p, MIC_k): F8,
            (SOM_c, MIC_k): F9,
            (MIC_k, SOM_c): F10
        })
    },
    bgc_md2_computers()
)    
