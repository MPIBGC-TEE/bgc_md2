import numpy as np
from sympy import var, ImmutableMatrix, diag
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonStateVariableTuple,
    NumericParameterization,
    NumericStartValueDict,
    NumericSimulationTimes,
)
from ..BibInfo import BibInfo 
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet

sym_dict = {
        'Y': 'young pool of soil carbon' # "kgCm^{-2}" 
        ,'O': 'old pool of soil carbon' # "kgCm^{-2}" 
        ,'k_1': 'decomposition rate of young pool' # "yr^{-1}"
        ,'k_2': 'decomposition rate of old pool' # "yr^{-1}"
        ,'i': 'mean annual carbon input' # "kgC m^{-2}yr^{-1}"
        ,'h': 'humification coefficient'
        ,'r': 'climatic and edaphic factors'
}

for name in sym_dict.keys():
    var(name)
xi = r #environmental effects multiplier
T = ImmutableMatrix([[-1,  0],
                     [ h, -1]]) #transition operator
N = diag(k_1, k_2) #decomposition operator
t = TimeSymbol("t") # unit: "year"
x = StateVariableTuple((Y, O))
u = InputTuple((i, 0))
B = CompartmentalMatrix(xi * T * N)
#np1 = NumericParameterization(
#    par_dict={
#},
#    func_dict=frozendict({})
#)
#
#nsv1 = NumericStartValueDict({
#})
#ntimes = NumericSimulationTimes(np.arange())

#model_run_data:
#    parameter_sets:
#        - "Bare fallow":
#            values:
#                k_1: 0.8
#                k_2: 0.00605
#                i: 0
#                h: 0.13
#                r: 1.32
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "+N +straw":
#            values:
#                k_1: 0.8
#                k_2: 0.00605
#                i: 0.285
#                h: 0.125
#                r: 1.00
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "-N +straw":
#            values: {k_1: 0.8, k_2: 0.00605, i: 0.248, h: 0.125, r: 1.22}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "-N -straw":
#            values: {k_1: 0.8, k_2: 0.00605, i: 0.057, h: 0.125, r: 1.17}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "+N -straw":
#            values: {k_1: 0.8, k_2: 0.00605, i: 0.091, h: 0.125, r: 1.07}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "Farmyard manure":
#            values: {k_1: 0.8, k_2: 0.00605, i: 0.272, h: 0.250, r: 1.10}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "Sewage sludge":
#            values: {k_1: 0.8, k_2: 0.00605, i: 0.296, h: 0.34, r: 0.97}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#
#    initial_values:
#        - "Bare fallow":
#            values: {Y: 0.3, O: 3.96}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "+N +straw":
#            values: {Y: 0.3, O: 4.11}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "-N +straw":
#            values: {Y: 0.3, O: 4.05}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "-N -straw":
#            values: {Y: 0.3, O: 3.99}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "+N -straw":
#            values: {Y: 0.3, O: 4.02}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "Farmyard manure":
#            values: {Y: 0.3, O: 3.99}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#        - "Sewage sludge":
#            values: {Y: 0.3, O: 4.14}
#            #doi: 10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2
#
#    run_times:
#        - RT1:
#            start: 0
#            end: 100
#            step_size: 0.1
#
#    possible_combinations:
#        - ["Bare fallow", "Bare fallow", RT1]
#        - ["+N +straw", "+N +straw", RT1]
#        - ["-N +straw", "-N +straw", RT1]
#        - ["-N -straw", "-N -straw", RT1]
#        - ["+N -straw", "+N -straw", RT1]
#        - ["Farmyard manure", "Farmyard manure", RT1]
#        - ["Sewage sludge", "Sewage sludge", RT1]
mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="ICBM",
        longName="Introductory Carbon Balance Model", 
        version="",
        entryAuthor="Holger Metzler",
        entryAuthorOrcid="0000-0002-8239-1601",
        entryCreationDate="09/03/2016",
        doi="10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2",
        sym_dict=sym_dict
    ),
    B,  # the overall compartmental matrix
    u,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
#    np1,
#    nsv1,
#    ntimes
})
