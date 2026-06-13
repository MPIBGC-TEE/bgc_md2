import numpy as np
from sympy import var, ImmutableMatrix, exp
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    NumericParameterization,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict = {
        'C_1': 'decomposable plant material pool (DPM)' # "t C*ha^{-1}"
        ,'C_2': 'resistant plant material pool (RPM)' # "t C*ha^{-1}"
        ,'C_3': 'microbial biomass pool (BIO)' # "t C*ha^{-1}"
        ,'C_4': 'humified organic matter pool (HUM)' # "t C*ha^{-1}"
        ,'C_5': 'inert organic matter pool (IOM)' # "t C*ha^{-1}"
        ,'k_1': 'decomposition rate of DPM' # "yr^{-1}"
        ,'k_2': 'decomposition rate of RPM' # "yr^{-1}"
        ,'k_3': 'decomposition rate of BIO' # "yr^{-1}"
        ,'k_4': 'decomposition rate of HUM' # "yr^{-1}"
        ,'pClay': 'percentage of clay in mineral soil'
        ,'DR': 'ratio of DPM to RPM'
        ,'x': '"CO$_2$ to (BIO+HUM) ratio"'
        ,'gamma': 'litter input partitioning coefficient'
        ,'J': 'mean annual carbon input' # "t C ha^{-1}yr^{-1}"
        ,'a': 'flux coefficient to BIO'
        ,'b': 'flux coefficient to HUM'
        ,'f_T': 'function of temperature'
        ,'f_W': 'function of soil moisture'
}

for name in sym_dict.keys():
    var(name)
x = 1.67 * (1.85 + 1.60 * exp(-0.0786 * pClay))
gamma = DR/(1+DR)
a = 0.46/(1+x)
b = 0.54/(1+x)
t = TimeSymbol("t") # unit: "year" # ??? parameters are yearly, but it is run on monthly basis
x = StateVariableTuple((C_1, C_2, C_3, C_4, C_5))
u = InputTuple((J * ImmutableMatrix(5, 1, [gamma, 1-gamma, 0, 0, 0])))
xi = f_T * f_W # environmental effects multiplier
A = ImmutableMatrix([[  -k_1,    0 ,          0,         0, 0],
                     [     0,  -k_2,          0,         0, 0],
                     [ a*k_1, a*k_2, -k_3+a*k_3,     a*k_4, 0],
                     [ b*k_1, b*k_2,      b*k_3,-k_4+b*k_4, 0],
                     [     0,     0,          0,         0, 0]])
B = CompartmentalMatrix(xi * A)
# Original values without effects of temperature and soil moisture
# doi = {10.1007/978-3-642-61094-3_17},
np1 = NumericParameterization(
    par_dict={
              k_1: 10
              ,k_2: 0.3
              ,k_3: 0.66
              ,k_4: 0.02
              ,pClay: 23.4
              ,DR: 1.44
              ,J: 1.7
              ,xi: 1
},
    func_dict=frozendict({})
)

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="RothC-26.3",
            longName="", 
            version="",
            entryAuthor="Holger Metzler",
            entryAuthorOrcid="0000-0002-8239-1601",
            entryCreationDate="10/03/2016",
            doi="10.1097/00010694-197705000-00005",
            further_references=BibInfo(doi="10.1007/978-3-642-61094-3_17"),
            sym_dict=sym_dict
        ),
        B,  # the overall compartmental matrix
        u,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
        np1,
    },
    bgc_md2_computers()
)
