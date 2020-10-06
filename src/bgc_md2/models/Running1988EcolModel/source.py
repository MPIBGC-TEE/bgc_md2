import numpy as np
from sympy import var, symbols, Symbol, ImmutableMatrix, diag
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
    NumericStartValueDict,
    NumericSimulationTimes,
   )
from ..BibInfo import BibInfo 
from bgc_md2.helper import MVarSet

sym_dict={
        'C_f': 'Carbon in foliage',
        'C_r': 'Carbon in roots',
        'C_w': 'Carbon in woody tissue',
        'eta_f': 'Fixed partitioning ratio of available carbon allocated to foliage',
        'eta_r': 'Fixed partitioning ratio of available carbon allocated to roots',
        'eta_w': 'Fixed partitioning ratio of available carbon allocated to wood',
        'gamma_f': 'Foliage turnover rate',
        'gamma_r': 'Roots turnover rate',
        'gamma_w': 'Wood turnover rate'
}
for name in sym_dict.keys():
    var(name)

t = TimeSymbol("t") #"days, years for allocation"
x = StateVariableTuple((C_f, C_r, C_w ))
u = u 
b = (eta_f, eta_w, eta_r)
Input = InputTuple(u * ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-gamma_f, -gamma_w, -gamma_r)
)

#        - "Original dataset of the publication":
#            values: {eta_f: 'Rational(25,100)', eta_r: 'Rational(40,100)', eta_w: 'Rational(35,100)', gamma_r: 'Rational(40,100)', gamma_f: 'Rational(33,100)', gamma_w: 'Rational(0,100)'}
#            doi: 10.1016/0304-3800(88)90112-3
#        - "Additional set 1":
#            values: {eta_f: 'Rational(20,100)', eta_r: 'Rational(55,100)', eta_w: 'Rational(25,100)', gamma_r: 'Rational(75,100)'}
#            doi: 10.1093/treephys/9.1-2.161 # Hunt1991TreePhysiol
#        - "Additional set 2":
#            values: {eta_f: 'Rational(48,100)', eta_r: 'Rational(37,100)',eta_w: 'Rational(15,100)', gamma_r: 'Rational(75,100)'}
#            doi: "10.1139/x91-151" # Korol1991CanJForRes

## The following are default values suggested by this entry's creator only to be able to run the model:
np1 = NumericParameterization(
    par_dict={u: 1400, eta_f: 0.48, eta_r: 0.44, eta_w: 0.49, gamma_r: 3.03, gamma_f: 23.32, gamma_w: 0.04},
    func_dict=frozendict({})
)

nsv1 = NumericStartValueDict({
    C_f: 200, 
    C_w: 5000, 
    C_r: 300
})

ntimes = NumericSimulationTimes(np.arange(0, 20000, 0.01))

mvs=MVarSet({
    BibInfo(# Bibliographical Information
        name="FOREST-BGC",
        longName="", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="17/7/2015",
        doi: 10.1016/0304-3800(88)90112-3,
#        modApproach: process based,
#        partitioningScheme: fixed,
#        claimedDynamicPart: "no",
#        spaceScale: forest, 
#        #    unit: "1°",
#        timeResolution: monthly,
#      - reviewer: Carlos Sierra
#        orcid: 0000-0003-0009-4169
#        date: 12/04/2016
#        desc: "Added a set of parameters that were missing from original publication."
#        type: shallow 
        sym_dict=sym_dict
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_f, C_w, C_r)),
    np1
})
