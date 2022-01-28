from sympy import var, ImmutableMatrix, diag
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
    NumericParameterization,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict = {
        'B_L': 'Carbon in foliage'
        ,'B_R': 'Carbon in roots'
        ,'B_S': 'Carbon in woody tissue'
        ,'Q_i': 'describes light availability'
        ,'G_i': 'describes water availability'
        ,'C_i': 'deviance of leaf biomass from $a_{0L}$'
        ,'A_CC': 'Leaf level photosynthetic rate'
        ,'R_g': 'Growth respiration' # PAR use efficiency
        ,'C': 'Canopy area of the plant'
        ,'C_delta': 'Net carbon gain'
    #CAU,'ON: Partitions in this work were the following: Leaves, Roots (including woody roots, not only fine roots?) and Stems
        ,'a_0R': 'fraction of carbon allocated to roots when resources are not limiting'
        ,'a_0S': 'fraction of carbon allocated to stems when resources are not limiting'
        ,'a_0L': 'fraction of carbon allocated to leaves when resources are not limiting'
        ,'a_L': ''
        ,'a_R': ''
        ,'a_S': ''
        ,'gamma_f': ''
        ,'gamma_r': ''
        ,'gamma_w': ''
}

for name in sym_dict.keys():
    var(name)
C_i = (B_L/(a_0L*(B_R+B_S+B_L)))
C_delta = A_CC*C-R_g
a_L = ((1-C_i)/(3+a_0R+a_0S-Q_i-G_i-C_i))
a_R = ((1+a_0R-G_i)/(3+a_0R+a_0S-Q_i-G_i-C_i))
a_S = ((1+a_0S-Q_i)/(3+a_0R+a_0S-Q_i-G_i-C_i))
t = TimeSymbol("t") # unit: "month"
x = StateVariableTuple((B_L, B_R, B_S))
u = C_delta
b = (a_L, a_R, a_S)
Input = InputTuple(u*ImmutableMatrix(b))
B = CompartmentalMatrix(
diag(-gamma_f, -gamma_r, -gamma_w))
np1 = NumericParameterization(
    par_dict={
a_0L: 0.3, a_0R: 0.5, a_0S: 0.2, Q_i: 1},
    #Q_i: 1 when light is highly available (See sup. material 1, pg 17)
    func_dict=frozendict({})
)

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="aDGVM",
            longName="", 
            version="",
    #       modApproach= "individuals based"
            entryAuthor="Verónika Ceballos-Núñez",
            entryAuthorOrcid="0000-0002-0046-1160",
            entryCreationDate="17/7/2015",
            doi="10.1111/j.1365-2486.2008.01838.x",
            sym_dict=sym_dict
        ),
        B,  # the overall compartmental matrix
        Input,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
        VegetationCarbonInputScalar(u),
        # vegetation carbon partitioning.
        VegetationCarbonInputPartitioningTuple(b),
        VegetationCarbonStateVariableTuple((B_L, B_R, B_S)),
        # fixme mm 01-20-2022 the parameterization is incomplete error: The
        # following free symbols: {gamma_f, gamma_r, gamma_w} of the
        # expression: {gamma_f, gamma_r, gamma_w} are not arguments.
        #np1,
    },
    bgc_md2_computers()
)
