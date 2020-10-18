from sympy import var, symbols, Symbol, ImmutableMatrix, diag, exp
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
from bgc_md2.helper import MVarSet

sym_dict = {
        'F': 'Foliage dry mass' #"kgC*m^{-2}" 
        ,'R': 'Fine roots dry mass' #"kgC*m^{-2}" 
        ,'W': 'Woody tissue dry mass' #"kgC*m^{-2}" 
        ,'k': 'Radiation extinction coefficient of canopy'
        ,'Phi_0': 'Incident PAR'
        ,'omega': 'Specific leaf area'
        ,'Phi': ' Annual photosynthetically active radiation (PAR) intercepted by the canopy' #"MJ*m*^{-2}*year^{-1}"
        ,'epsilon': 'Light utilization coefficient' #"kg*MJ^{-1}"
        ,'G': 'Rate of biomass production per unit ground area' #"kg*m^{-2}*year^{-1}" 
        ,'eta_f': 'Fraction of biomass production partitioned to leaves'
        ,'eta_r': 'Fraction of biomass production partitioned to roots'
        ,'eta_w': 'Fraction of biomass production partitioned to wood'
        ,'gamma_f': 'Senescence rate per unit foliage biomass'
        ,'gamma_r': 'Senescence rate per unit fine roots biomass'
}

for name in sym_dict.keys():
    var(name)

Phi = Phi_0*(1-exp(-k*omega*F))
G=Phi*epsilon
eta_w=1-eta_f-eta_r

x = StateVariableTuple((F, R, W))
u = G
b = (eta_f, eta_r, eta_w)
Input = InputTuple(u*ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-gamma_f, -gamma_r, 0)
)
t = TimeSymbol("t") # units: days, years for allocation

#Parameter set:  "Chosen based on the performance of Pinus radiata at Puruki, New Zeland":
## eta_ are not mentioned in the paper... probably the different values they gave to them were part of the hypothesis
np1 = NumericParameterization(
    par_dict={
        k: 0.5, 
        Phi_0: 2500, #"MJ*m*^{-2}*year^{-1}" 
        omega: 5, #"m*^2*kg^{-1}"
        gamma_r: 2, #"kg^{-1}"
        gamma_f: 0.5}, #"kg^{-1}"
    func_dict=frozendict({})
)

mvs=MVarSet({
    BibInfo(# Bibliographical Information
        name="",
        longName="", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="29/7/2015",
        doi="10.1093/treephys/12.2.119",
        sym_dict=sym_dict
    ),
    #
    # the following variables constitute the compartmental system:
    #
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((F,R,W)),
    np1
})
