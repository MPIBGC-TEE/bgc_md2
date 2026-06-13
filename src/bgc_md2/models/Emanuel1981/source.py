import numpy as np
from sympy import var, ImmutableMatrix, Min
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
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict = {
        'x_1': 'Non-woody tree parts'
        ,'x_2': 'Woody tree parts '
        ,'x_3': 'Ground vegetation'
        ,'x_4': 'Detritus/Decomposers'
        ,'x_5': 'Active soil carbon'
        ,'I_1': 'Gross primary production of the non-woody pool' # "PgC yr^{-1}" 
        ,'I_3': 'Gross primary production of the ground vegetation pool' # "PgC yr^{-1}" 
        ,'F_1': 'Cycling rate of pool 1' # "yr^{-1}"
        ,'F_2': 'Cycling rate of pool 2' # "yr^{-1}"
        ,'F_3': 'Cycling rate of pool 3' # "yr^{-1}"
        ,'F_4': 'Cycling rate of pool 4' # "yr^{-1}"
        ,'F_5': 'Cycling rate of pool 5' # "yr^{-1}"
        ,'F_21': 'Transfer coefficient from non-woody to woody vegetation'
        ,'F_41': 'Transfer coefficient from non-woody to detritus'
        ,'F_42': 'Transfer coefficient from ground vegetaion to detritus'
        ,'F_52': 'Transfer coefficient from woody parts to soil'
        ,'F_43': 'Transfer coefficient from ground vegetation to deteritus'
        ,'F_53': 'Transfer coefficient from ground vegetation to soil'
        ,'F_54': 'Transfer coefficient from detritus to soil'
}

for name in sym_dict.keys():
    var(name)
t = TimeSymbol("t") # unit: "day"
x = StateVariableTuple((x_1, x_2, x_3, x_4, x_5))
u = InputTuple((I_1, 0, I_3, 0, 0))
B = CompartmentalMatrix([[-F_1,        0,       0,       0,        0],
                                [F_21,     -F_2,       0,       0,        0],
                                [   0,        0,    -F_3,       0,        0],
                                [F_41,     F_42,    F_43,    -F_4,        0],
                                [   0,     F_52,    F_53,    F_54,     -F_5]]
)
np1 = NumericParameterization(
    par_dict={I_1: 77, I_3: 36, F_1: 2.081, F_2: 0.0686, F_3: 0.5217, F_4: 0.5926, F_5: 9.813e-3, F_21: 0.8378, F_41:  0.5676, F_42: 0.0322, F_52: 4.425e-3, F_43: 0.1739, F_53: 0.0870, F_54: 0.0370
},
    func_dict=frozendict({})
)

nsv1 = NumericStartValueDict({x_1: 37, x_2: 452, x_3: 69, x_4: 81, x_5: 1121
})
ntimes = NumericSimulationTimes(np.arange(0,200,0.1))

mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="",
            longName="", 
            version="",
            entryAuthor="Carlos A. Sierra",
            entryAuthorOrcid="0000-0003-0009-4169",
            entryCreationDate="12/9/2016",
            doi="",
    #        bibtex= "@incollection{Emanuel1981,
    #                     author = {W. R. Emanuel and G. G. Killough and J. S. Olson},
    #                     booktitle = {Carbon Cycle Modelling},
    #                     editor = {Bert Bolin},
    #                     pages = {335--353},
    #                     publisher = {John Wiley and Sons},
    #                     series = {SCOPE 16},
    #                     title = {Modelling the circulation of carbon in the world's terrestrial ecosystems},
    #                     year = {1981}
    #                     abstract = {"A mathematical model for the circulation of carbon in the world terrestrial ecosystems is proposed. A five-compartment representation is developed which corresponds to the functional components studied by field ecologists. Rate coefficients for this linear dynamic model are calculated from estimates of the 1970 standingcrops and compartment exchanges of carbon. The model is analyzed in terms of response to a unit impulse, thereby displaying a transient time distribution. The response to a hypothetical pulse input through gross primary production is also simulated, illustrating the efficiency of the terrestrial carbon system in transferring carbon into longer storage components. Finally, the concept of CO$_2$ fertilization is examined by allowing gross primary production to increase in response to higher atmospheric concentrations. Although the standing crop of carbon in photosynthesizing compartments is induced to grow from a hypothetical preindustrial level to a specified 1970 level, the accompanying increase in other compartments is not as large as obtained in earlier model formulations which incorporate an input from the atmosphere directly to compartments containing carbon in woody material or soil."}
    #                }",
            sym_dict=sym_dict
        ),
        B,  # the overall compartmental matrix
        u,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
        VegetationCarbonStateVariableTuple((x_1,x_2)),
        np1,
        nsv1,
        ntimes
    },
    bgc_md2_computers()
)
