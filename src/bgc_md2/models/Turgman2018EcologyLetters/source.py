from sympy import var, ImmutableMatrix
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
)
from ..BibInfo import BibInfo 
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet

sym_dict = {
        'C': 'Non-structural carbon (NSC), starch' # "kgC*m^{-2}" 
        ,'X': 'Xylem biomass ' # "kgC*m^{-2}" 
        ,'L': 'Leaf and fine roots biomass which are linearly related through a constant' # "kgC*m^{-2}" 

        ,'A_n': 'whole-plant net photosynthesis (classic photosynthetic model of CO2 demand for carbon-limited photosynthesis), function of the interstitial CO2 concentration, the daytime leaf respiration, the maximum rate of carboxylation, and empirical constants'
        ,'m_X': 'Turnover rate of the xylem' # "month^{-1}"
        ,'m_L': 'Turnover rate of the xylem' # "month^{-1}"
        ,'C_i': 'initial NSC that covers carbon for two crown and fine root flushes'
        ,'W_max': 'Maximum m sucrose loading rate and depends on tree size.'
        ,'k_c': 'Michaelis constant for the sucrose loading rate'
        ,'W': 'Sucrose loading rate from storage to the phloem'
        ,'L_opt': 'Optimal total leaf biomass'
        ,'U': 'Variable optimized in the system in the range [0,1]. Fraction of translocatable C invested in xylem reconstruction'
}

for name in sym_dict.keys():
    var(name)
W = (W_max * C)/((k_c*C_i)+C)
U = (X/(X+L_opt))*((W-((m_L+m_X)*L_opt))/W)
t = TimeSymbol("t") # unit: "month"
x = StateVariableTuple((C,X,L))
u = A_n
beta=(1,0,0)
Input = InputTuple(u*ImmutableMatrix(beta))
B = CompartmentalMatrix(
[[-W/C,0,0],
[W*U/C, -m_X, 0],
[W*(1-U)/C, 0, -m_L]])

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="",
        longName="", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="",
        doi="10.1111/ele.13136",
        sym_dict=sym_dict
    ),
    B,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(beta),
    VegetationCarbonStateVariableTuple((C,X,L)),
})
