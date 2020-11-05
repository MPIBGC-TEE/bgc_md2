from sympy import var, ImmutableMatrix
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
)
from ..BibInfo import BibInfo 
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet

sym_dict = {
        'C_s': 'soil organic matter' # "g C m^{-2}"
        ,'C_b': 'microbial biomass' # "g C m^{-2}"
        ,'epsilon': 'microbial growth efficiency'
        ,'V_s': 'maximum rate of soil carbon assimilation per unit microbial biomass per year' # "year^{-1}"
        ,'K_s': 'half-saturation constant for soil carbon assimilation by microbial biomass' # "g C m^{-2}"
        ,'lamda': 'soil organic matter decomposition rate'
        ,'mu_b': 'turnover rate of microbial biomass per year' # "year^{-1}"
        ,'F_NPP': 'carbon influx into soil' # "g C m^{-2} year^{-1}"
}

for name in sym_dict.keys():
    var(name)
lamda = (C_b*V_s)/(C_s+K_s)
t = TimeSymbol("t") # unit: "year"
x = StateVariableTuple((C_s, C_b))
u = InputTuple((F_NPP, 0))
B = CompartmentalMatrix([[       -lamda,  mu_b],
                         [epsilon*lamda, -mu_b]])

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="Two-pool microbial",
        longName="", 
        version="1",
        entryAuthor="Holger Metzler",
        entryAuthorOrcid="0000-0002-8239-1601",
        entryCreationDate="18/01/2018",
        doi="10.5194/bg-11-1817-2014",
        sym_dict=sym_dict
    ),
    B,  # the overall compartmental matrix
    u,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
})
