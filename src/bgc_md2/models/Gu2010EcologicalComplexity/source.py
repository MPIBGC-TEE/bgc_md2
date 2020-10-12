from sympy import var, ImmutableMatrix, diag
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
        'C_S': 'Carbon in stem'
        ,'C_R': 'Carbon in root'
        ,'C_L': 'Carbon in leaf'
        ,'NPP': 'Net Primary Production'
        ,'L': 'Scalar light availability'
        ,'W': 'Scalar water availability'
        ,'N_ava': 'effect of nitrogen availability on carbon allocation'
        ,'Omega': 'Sensitivity of allocation to changes in resources availability. If =0, partitioning is determined by constant allocation fractions.'
        ,'epsilon_S': 'Parameter relative to vegetation type'
        ,'epsilon_R': 'Parameter relative to vegetation type'
        ,'epsilon_L': 'Parameter relative to vegetation type'
        ,'a_S': 'Allocation fraction to stem'
        ,'a_R': 'Allocation fraction to root'
        ,'a_L': 'Allocation fraction to leaf'
        ,'gamma_S': 'Stem turnover rate'
        ,'gamma_R': 'Root turnover rate'
        ,'gamma_L': 'Stem turnover rate'
}

for name in sym_dict.keys():
    var(name)

epsilon_L = 1-epsilon_R-epsilon_S
a_S = (epsilon_S+(Omega*(1.5-L-(0.5*N_ava))))/(1+(Omega*(3-L-W-N_ava)))
a_R = (epsilon_R+(Omega*(1.5-W-(0.5*N_ava))))/(1+(Omega*(3-L-W-N_ava)))
a_L = epsilon_L/(1+(Omega*(3-L-W-N_ava))) #, "a_L = 1-a_S-a_R"

x = StateVariableTuple((C_S, C_R, C_L))
u = NPP
b = (a_S,a_R,a_L)
Input = InputTuple(u*ImmutableMatrix(b))
A = CompartmentalMatrix([[-gamma_S,          0,          0],
                               [        0,   -gamma_R,          0],
                               [        0,          0, -gamma_L]]
)
t = TimeSymbol("t")

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="CEVSA2  ",
        longName="", 
        version="2",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="",
        doi="10.1016/j.ecocom.2010.04.002",
        sym_dict=sym_dict
#        ,modApproach="process based"
#        ,partitioningScheme="dynamic"
#        ,claimedDynamicPart="yes"
#        ,spaceScale="forest"
#        ,timeResolution="daily"
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_S, C_R, C_L)),
})
