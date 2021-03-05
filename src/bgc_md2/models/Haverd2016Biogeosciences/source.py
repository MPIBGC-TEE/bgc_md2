from sympy import var, ImmutableMatrix, exp
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
        'C_L': 'Leaf biomass' # "molC*m^{-2}" 
        ,'C_R': 'Fine roots biomass' # "molC*m^{-2}" 
        ,'C_stem': 'Trunk and coarse roots' # "molC*m^{-2}" 
        ,'F_Cgrowth': 'Growth, the flux of C to structural components, parameterized by a logistic curve. Function of soil water, C_L, and C_R relative to a prognostic carrying capacity C_max' # "molC*m^{-2}*day^{-1}" 
        ,'alpha_L': 'carbon allocation coefficient, Heaviside step function that depends on C_L and C_R'
        ,'alpha_R': 'carbon allocation coefficient, Heaviside step function that depends on C_L and C_R'
        ,'alpha_stem': 'carbon allocation coefficient, Heaviside step function that depends on C_L and C_R'
        ,'k_L': 'First-order rate constant' # "day^{-1}" 
        ,'k_R': 'First-order rate constant' # "day^{-1}" 
        ,'m_stem': 'Stem biomass turnover' # "molC*m^{-2}*day^{-1}" 
}

for name in sym_dict.keys():
    var(name)
t = TimeSymbol("t") # unit: "day"
x = StateVariableTuple((C_L,C_R,C_stem))
u = F_Cgrowth
beta = (alpha_L,alpha_R,alpha_stem)
alpha_L = 1 - (alpha_R + alpha_stem)
Input = InputTuple(u*ImmutableMatrix(beta))
B = CompartmentalMatrix(
                        [[-(k_L+(m_stem/C_stem)),0,0],
                         [0, -(k_R+(m_stem/C_stem)), 0],
                         [0, 0, -m_stem]])

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="HAVANA",
        longName="Hydrology and Vegetation-dynamics Algorithm for Northern Australia", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="",
        doi="10.5194/bg-13-761-2016",
        sym_dict=sym_dict
    ),
    B,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
#    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
#    VegetationCarbonInputPartitioningTuple(beta),
    VegetationCarbonStateVariableTuple((C_L,C_R,C_stem)),
})
