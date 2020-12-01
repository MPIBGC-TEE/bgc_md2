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
        'C_frl': 'Carbon in leaf/root litter',
        'C_s': 'Carbon in soil',
        'eta_f': 'Fixed partitioning ratio of available carbon allocated to foliage',
        'eta_r': 'Fixed partitioning ratio of available carbon allocated to roots',
        'eta_w': 'Fixed partitioning ratio of available carbon allocated to wood',
        'gamma_f': 'Foliage turnover rate',
        'gamma_r': 'Roots turnover rate',
        'gamma_w': 'Wood turnover rate',
        'l_s': 'Fraction of litter partitioned to the soil',
        'l_p': 'Fraction of litter available for plant',
        's_p': 'Fraction of soil carbon available for plant',
        'l_dr': 'Litter decomposition/respiration',
        's_dr': 'Soil decomposition/respiration',
        'u': 'scalar function of photosynthetic inputs'
}
for name in sym_dict.keys():
    var(name)

#Model based on Fig. 1 in page 127 (3 in the PDF) of doi= "10.1016/0304-3800(88)90112-3", no ODEs in publication
t = TimeSymbol("t") #"the model has a daily and a yearly component. Allocation occurs yearly"
x = StateVariableTuple((C_f, C_r, C_w,C_frl,C_s))
u 
b = (eta_f, eta_w, eta_r)
#Input = InputTuple(u * ImmutableMatrix(b))
Input = InputTuple((u * ImmutableMatrix(b),0,0))
A = CompartmentalMatrix(
#    diag(-gamma_f, -gamma_w, -gamma_r)
[[-gamma_f,    0   ,    0   ,   (l_p*eta_f)  ,(s_p*eta_f)],
 [    0   ,-gamma_r,    0   ,   (l_p*eta_r)  ,(s_p*eta_r)],
 [    0   ,    0   ,-gamma_w,   (l_p*eta_w)  ,(s_p*eta_w)],
 [ gamma_f, gamma_r, gamma_w, -(l_p+l_s+l_dr),      0    ],
 [    0   ,    0   ,    0   ,       l_s      ,-(s_p+s_dr)]
])

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
        doi= "10.1016/0304-3800(88)90112-3",
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
#    np1
})
