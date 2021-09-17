import numpy as np
from sympy import var, ImmutableMatrix, Min, sqrt, exp
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
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet
# questions 
# 1.) (for yuanyuan) What do the following variables mean exactly?
#     lig_leaf, lig_wood, clay, silt
#
# 2.) should the k_leaf ...k_passom be realy called k_ ?
#     d(because they are NOT the diagonal terms in the K matrix)
#
# 3.) I spent a lot of time transforming the coefficient names back
#     from numbers into intellegible names (e.g. f41 into f_leaf2metlit)
#     It is very important not to use numbers in the names for several reasons.
#     (remember the lessons learned part of the soilR talk)

sym_dict = {
        'C_leaf': 'Carbon in foliage',
        'C_root': 'Carbon in roots',
        'C_wood': 'Carbon in woody tissue',
        'C_metlit': 'Carbon in metabolic litter',
        'C_stlit': 'Carbon in structural litter',
        'CWD': 'Corse woody debris Carbon ?',
        'C_mic': 'Carbon ?',
        'C_slowsom': 'Carbon in slow SOM',
        'C_passsom': 'Carbon in passive SOM',
        'NPP': 'Photosynthetic rate (Carbon input) at time t',
        'beta_leaf': 'Fixed partitioning ratio (fraction) of available carbon allocated to foliage',
        'beta_root': 'Fixed partitioning ratio (fraction) of available carbon allocated to roots',
        'beta_wood': 'Fixed partitioning ratio (fraction) of available carbon allocated to wood',
        'k_leaf': 'Foliage cycling rate' ,
        'k_wood': 'Woody cycling rate' ,
        'k_root': 'Fine roots cycling rate' ,
        'k_metlit': 'Metabolic litter cycling rate' ,
        'k_stlit': 'Structural litter cycling rate' ,
        'k_CWD': 'cycling rate' ,
        'k_mic': 'Microbial SOM cycling rate' ,
        'k_slowsom': 'Slow SOM cycling rate' ,
        'k_passsom': 'Passive SOM cycling rate' ,
        'f_leaf2metlit': 'Transfer coefficient from Foliage to Metabilic Litter',
        #'f_leaf2stlit': 'Transfer coefficient from Foliage to Structural Litter',
        'f_wood2CWD': 'Transfer coefficient from Wood to CWD',
        'f_root2metlit': 'Transfer coefficient from Fine Roots to Metabolic Litter',
        #'f_root2stlit': 'Transfer coefficient from Fine Roots to Structural Litter',
        'f_metlit2mic': 'Transfer coefficient from Metabolic Litter to Fast SOM',
        'f_stlit2mic': 'Transfer coefficient from Structural Litter to Fast SOM',
        'f_stlit2slowsom': 'Transfer coefficient from Structural Litter to Slow SOM',
        'f_mic2slowsom': 'Transfer coefficient from Fast to Slow SOM',
        'f_mic2passsom': 'Transfer coefficient from Fast to Passive SOM',
        'f_CWD2slowsom': 'Transfer coefficient from CWD to Slow SOM',
        'f_CWD2passsom': 'Transfer coefficient from CWD to passive SOM',
        'f_slowsom2passsom': 'Transfer coefficient from Slow to Passive SOM',
        'lig_leaf': '?' ,
        'lig_wood': '?' ,
        'clay': '?',
        'silt': '?',
}

for name in sym_dict.keys():
    var(name)
    
beta_wood = 1.0 - (beta_leaf + beta_root)
f_leaf2stlit = 1.0 - f_leaf2metlit
f_root2stlit = 1.0 - f_root2metlit
f_stlit2mic = 0.45 * (1.0 - lig_leaf)
f_stlit2slowsom = 0.7 * lig_leaf
f_CWD2slowsom = 0.4 * (1.0 - lig_wood)
f_CWD2passsom = 0.7 * lig_wood
f_mic2slowsom = (0.85 - 0.68 * (clay+silt)) * (0.997 - 0.032 * clay)
f_mic2passsom = (0.85 - 0.68 * (clay+silt)) * (0.003 + 0.032 * clay)
f_slowsom2passsom = 0.45 * (0.003 + 0.009 * clay)

temp_leaf     = k_leaf
temp_wood     = k_wood
temp_root     = k_root
temp_metlit   = k_metlit
temp_stlit    = k_metlit / (5.75 * exp(-3.0 * lig_leaf))
temp_CWD      = k_metlit / 20.6
temp_mic      = k_mic
temp_slowsom  = k_slowsom
temp_passsom  = k_passsom

x = StateVariableTuple((C_leaf, C_root, C_wood, C_metlit, C_stlit, CWD, C_mic, C_slowsom, C_passsom))
K = ImmutableMatrix.diag([temp_leaf, temp_root, temp_wood, temp_metlit, temp_stlit, temp_CWD, temp_mic, temp_slowsom, temp_passsom] )
A = ImmutableMatrix(
        [
            [               -1,                 0,          0,                  0,                 0,             0,                    0,                 0,   0],
            [                0,                -1,          0,                  0,                 0,             0,                    0,                 0,   0],
            [                0,                 0,         -1,                  0,                 0,             0,                    0,                 0,   0],
            [ f_leaf2metlit, f_root2metlit,         0,                 -1,                 0,             0,                    0,                 0,   0],
            [  f_leaf2stlit,  f_root2stlit,         0,                  0,                -1,             0,                    0,                 0,   0],
            [                0,                 0, f_wood2CWD,                  0,                 0,            -1,                    0,                 0,   0],
            [                0,                 0,          0, f_metlit2mic, f_stlit2mic,             0,                   -1,                 0,   0],
            [                0,                 0,          0,                  0,   f_stlit2slowsom, f_CWD2slowsom,  f_mic2slowsom,                -1,   0],
            [                0,                 0,          0,                  0,                 0, f_CWD2passsom, f_mic2passsom, f_slowsom2passsom,  -1 ]
        ]
)   # tranfer

B = CompartmentalMatrix(A*K)
t = TimeSymbol("t")  # unit: "day"
u = NPP
b = ImmutableMatrix([beta_leaf, beta_root, beta_wood, 0, 0, 0, 0, 0, 0])
Input = InputTuple(u * b)

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="CABLE_yuanyuan",
        longName="Terrestrial Ecosystem Model", 
        version="",
        entryAuthor="Markus Müller",
        entryAuthorOrcid="0000-0003-0009-4169",
        entryCreationDate="08/24/2021",
        doi="",
        sym_dict=sym_dict
    ),
    B,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_leaf, C_root, C_wood)),
    #np1,
    #nsv1,
    #ntimes
})

