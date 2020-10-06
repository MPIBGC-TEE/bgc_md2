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
#        doi="",
#        modApproach: process based,
#        partitioningScheme: fixed,
#        claimedDynamicPart: "no",
#        spaceScale: global, 
#        #    unit: "1°",
#        timeResolution: monthly,
#         bibtex: "@article{Comins1993Ecological_Applications,
#                     author = {Comins, H. N. and McMurtrie, Ross E.},
#                     copyright = {Copyright {\\copyright} 1993 Ecological Society of America},
#                     journal = {Ecological Applications},
#                     language = {English},
#                     link = {http://www.jstor.org/stable/1942099},
#                     number = {4},
#                     pages = {666-681},
#                     publisher = {Ecological Society of America},
#                     title = {Long-Term Response of Nutrient-Limited Forests to CO$_2$ Enrichment; Equilibrium Behavior of Plant-Soil Models},
#                     volume = {3},
#                     year = {1993}
#                  }",
#        
#        abstract: "Established process-based models of forest biomass production in relation to atmospheric CO$_2$ concentration (McMurtrie 1991) and soil carbon/nutrient dynamics (Parton et al. 1987) are integrated to derive the \"Generic Decomposition and Yield\" model (G'DAY). The model is used to describe how photosynthesis and nutritional factors interact to determine the productivity of forests growing under nitrogen-limited conditions. A simulated instantaneous doubling of atmospheric CO$_2$ concentration leads to a growth response that is initially large (27% above productivity at current CO$_2$) but declines to <10% elevation within 5 yr. The decline occurs because increases in photosynthetic carbon gain at elevated CO$_2$ are not matched by increases in nutrient supply. Lower foliar N concentrations at elevated CO$_2$ have two countervailing effects on forest production: decreased rates of N cycling between vegetation and soils (with negative consequences for productivity), and reduced rates of N loss through gaseous emission, fire, and leaching. Theoretical analysis reveals that there is an enduring response to CO$_2$ enrichment, but that the magnitude of the long-term equilibrium response is extremely sensitive to the assumed rate of gaseous emission resulting from mineralization of nitrogen. Theory developed to analyze G'DAY is applicable to other published production-decomposition models describing the partitioning of soil carbon among compartments with widely differing decay-time constants.",
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
