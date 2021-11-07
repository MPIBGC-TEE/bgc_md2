from sympy import var, ImmutableMatrix
from frozendict import frozendict
from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS

sym_dict = {
        'S1': 'Substrate 1' # "mgC g^{-1}\\text{ soil}"
        ,'B1': 'Microbial biomass guild 1' # "mgC g^{-1}\\text{ soil}"
        ,'S2': 'Substrate 2' # "mgC g^{-1}\\text{ soil}"
        ,'B2': 'microbial biomass guild 2' # "mgC g^{-1}\\text{ soil}"
        ,'S3': 'Substrate 3' # "mgC g^{-1}\\text{ soil}"
        ,'B3': 'Microbial biomass guild 3' # "mgC g^{-1}\\text{ soil}"
        ,'k_s1': 'Decomposition rate of substrate 1' # "m^{3} d^{-1} g^{-1}"
        ,'k_s2': 'Decomposition rate of substrate 2' # "m^{3} d^{-1} g^{-1}"
        ,'k_s3': 'Decomposition rate of substrate 3' # "m^{3} d^{-1} g^{-1}"
        ,'k_b1': 'Microbial decay rate guild 1' # "d^{-1}"
        ,'k_b2': 'Microbial decay rate guild 2' # "d^{-1}"
        ,'k_b3': 'Microbial decay rate guild 3' # "d^{-1}"
        ,'r_1': 'respired C fraction 1' # 
        ,'r_2': 'respired C fraction 2' #
        ,'r_3': 'respired C fraction 3' # 
        ,'K_M1': 'Michaelis-Menten constant 1' # "g m^{-3}"
        ,'K_M2': 'Michaelis-Menten constant 2' # "g m^{-3}"
        ,'K_M3': 'Michaelis-Menten constant 3' # "g m^{-3}"
        ,'I_1': 'input to substrate 1' # "g m^{-3} d^{-1}"
        ,'I_2': 'input to substrate 2' # "g m^{-3} d^{-1}"
        ,'I_3': 'input to substrate 3' # "g m^{-3} d^{-1}"
        ,'T': 'transition operator'
        ,'N': 'decomposition operator'
}

for name in sym_dict.keys():
    var(name)
t = TimeSymbol("t") # unit: "day"
x = StateVariableTuple((S1, B1, S2, B2, S3, B3))
u = InputTuple((I_1, 0, I_2, 0, I_3, 0))
T = ImmutableMatrix([[   -1,   1,   0,   0,   0,   0],
           [1-r_1,  -1,   0,   0,   0,   0],
           [    0,   0,  -1,   1,   0,   0],
           [    0,   0,1-r_2,  -1,  0,   0],
           [    0,   0,   0,   0,  -1,   1],
           [    0,   0,   0,   0, 1-r_3, -1]])
N = ImmutableMatrix([[k_s1*B1/(K_M1+S1),              0,                0,        0,                0,    0],
           [                 0,            k_b1,               0,        0,                0,    0],
           [                 0,              0,k_s2*B2/(K_M2+S2),        0,                0,    0],
           [                 0,              0,                0,     k_b2,                0,    0],
           [                 0,              0,                0,        0,k_s2*B2/(K_M2+S2),    0],
           [                 0,              0,                0,        0,                0, k_b3]])
B = CompartmentalMatrix(T*N)

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="Three microbe-substrate pairs",
            longName="", 
            version="",
            entryAuthor="Carlos Sierra",
            entryAuthorOrcid="0000-0003-0009-4169",
            entryCreationDate="13/04/2016",
            doi="10.1890/12-0681.1",
            sym_dict=sym_dict
        ),
        B,  # the overall compartmental matrix
        u,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
    },
    bgc_md2_computers()
)
