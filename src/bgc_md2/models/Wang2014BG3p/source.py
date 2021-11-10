from sympy import var
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict = {
        'C_l': 'litter carbon' # "g C m^{-2}"
        ,'C_s': 'soil organic matter' # "g C m^{-2}"
        ,'C_b': 'microbial biomass' # "g C m^{-2}"
        ,'alpha': 'fraction of carbon influx that directly enters the soil organic matter pool'
        ,'epsilon': 'microbial growth efficiency'
        ,'V_l': 'maximum rate of litter carbon assimilation per unit microbial biomass per year' # "year^{-1}"
        ,'K_l': 'half-saturation constant for litter carbon assimilation by microbial biomass' # "g C*m^{-2}"
        ,'V_s': 'maximum rate of soil carbon assimilation per unit microbial biomass per year' # "year^{-1}"
        ,'K_s': 'half-saturation constant for soil carbon assimilation by microbial biomass' # "g C*m^{-2}"
        ,'lamda_l': 'litter carbon decomposition rate'
        ,'lamda_s': 'soil organic matter decomposition rate'
        ,'mu_b': 'turnover rate of microbial biomass per year' # "year^{-1}"
        ,'F_NPP': 'carbon influx into soil' # "g C*m^{-2}*year^{-1}"
}

for name in sym_dict.keys():
    var(name)
lamda_l = (C_b*V_l)/(C_l+K_l)
lamda_s = (C_b*V_s)/(C_s+K_s)
t = TimeSymbol("t") # unit: "year"
x = StateVariableTuple((C_l, C_s, C_b))
u = InputTuple(((1-alpha)*F_NPP, alpha*F_NPP, 0))
B = CompartmentalMatrix(
[[       -lamda_l,               0,     0],
 [              0,        -lamda_s,  mu_b],
 [epsilon*lamda_l, epsilon*lamda_s, -mu_b]])

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="Three-pool microbial",
            longName="", 
            version="1",
            entryAuthor="Holger Metzler",
            entryAuthorOrcid="0000-0002-8239-1601",
            entryCreationDate="22/01/2018",
            doi="10.5194/bg-11-1817-2014",
            sym_dict=sym_dict
        ),
        B,  # the overall compartmental matrix
        u,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
    },
    bgc_md2_computers()
)
