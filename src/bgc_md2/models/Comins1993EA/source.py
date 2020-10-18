#import numpy as np
from sympy import var, symbols, Symbol, ImmutableMatrix, diag#, Rational
#from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
#    NumericParameterization,
#    NumericStartValueDict,
#    NumericSimulationTimes,
   )
from ..BibInfo import BibInfo 
from bgc_md2.helper import MVarSet

sym_dict={
        'F': 'Foliage carbon content per unit ground area at equilibrium'
        ,'R': 'Root carbon'
        ,'W': 'Carbon in woody tissue'
        ,'G': 'Net rate of plant carbon production'
        ,'eta_f': 'Allocation fraction to foliar biomass'
        ,'eta_r': 'Allocation fraction to roots biomass'
        ,'eta_w': 'Allocation fraction to wood (in stem, branches and large structurl roots) biomass'
        ,'gamma_f': 'Foliage senescence rate' #unit: "yr^{-1}" 
        ,'gamma_r': 'Roots senescence rate' #unit: "yr^{-1}" 
        ,'gamma_w': 'Wood senescence rate' #unit: "yr^{-1}" 
}
for name in sym_dict.keys():
    var(name)

x = StateVariableTuple((F, R, W ))
u = G
b = (eta_f, eta_w, eta_r)

Input = InputTuple(u * ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-gamma_f, -gamma_w, -gamma_r)
)
t = TimeSymbol("t") #'yr'

# Commented out the following lines because original publication only has 3 parameter values
#np1 = NumericParameterization(
#    par_dict={
#    G: , #"Mg*ha^{-1}*yr^{-1}"
#    eta_f: 'Rational(1,3)', 
#    eta_r: 'Rational(1,3)', 
#    eta_w: 'Rational(1,3)',
#},
#    func_dict=frozendict({})
#    # state_var_units=gram/kilometer**2,
#    # time_unit=day
#)
#nsv1 = NumericStartValueDict({
#    F: , #"Mg/ha"
#    W: , #"Mg/ha"
#    R: #"Mg/ha"
#})
#
#ntimes = NumericSimulationTimes(np.arange(, , ))

mvs=MVarSet({
    BibInfo(# Bibliographical Information
        name="G'DAY",
        longName="Generic Decomposition and Yield", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="27/1/2016",
        doi="10.2307/1942099",
        sym_dict=sym_dict
        
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((F, W, R)),
#    np1
})
