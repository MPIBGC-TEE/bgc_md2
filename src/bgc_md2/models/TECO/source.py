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

sym_dict = {
        'x_1': 'Carbon in foliage'
        ,'x_2': 'Carbon in roots'
        ,'x_3': 'Carbon in woody tissue'
        ,'x_4': 'Carbon in metabolic litter'
        ,'x_5': 'Carbon in structural litter'
        ,'x_6': 'Carbon in fast SOM'
        ,'x_7': 'Carbon in slow SOM'
        ,'x_8': 'Carbon in passive SOM'
        ,'T_k': 'Canopy temperature in Kelvin'
        ,'R': 'Universal gas constant'
        ,'E_p': 'Activation energy'
        ,'Arrhenius': 'Arrhenius multiplier'
        ,'C_a': 'Ambient CO2 concentration'
        ,'C_i': 'Leaf internal CO2 concentration'
        ,'Gamma': 'CO2 compensation point without dark respiration'
        ,'L': 'Leaf area index'
        ,'k': 'Light extinction coefficient'
        ,'I_0': 'Solar radiation at the top of the canopy'
        ,'I': 'Absorbed photosynthetically active radiation'
        ,'alpha_q': 'Quantum efficiency of photon capture'
        ,'J_m': 'Maximum electron transport rate'
        ,'J_e': 'Rate of light electron transport'
        ,'V_m': 'Maximum carboxylation rate'
        ,'O_x': 'Oxygen concentration in the air'
        ,'K_o': 'Michaelis-Menten constant for oxygenation'
        ,'K_c': 'Michaelis-Menten constant for carboxylation'
        ,'J_c': 'Rate of carboxilation with CO2 limitation'
        ,'R_d': 'Dark respiration'
        ,'A': 'Gross leaf CO2 assimilation rate'
        ,'D': 'Vapor pressure deficit'
        ,'D_0': 'Empirical coefficient'
        ,'g_l': 'Empirical coefficient'
        ,'G_s': 'Stomatal conductance'
        ,'A_n': 'Net photosynthesis rate at leaf level'
        ,'A_c': 'Canopy photosynthesis rate'
        ,'GPP': 'Photosynthetic rate (Carbon input) at time t'
        ,'T': 'Temperature'
        ,'Q_10': 'Temperature quotient that describes a change in decomposition rate for evey 10°C difference in temperature'
        ,'W': 'Volumetric soil moisture'
        ,'f_W': 'Function of W'
        ,'f_T': 'Function of T'
        ,'xi': 'Environmental scalar'
        ,'b_1': 'Fixed partitioning ratio (fraction) of available carbon allocated to foliage'
        ,'b_2': 'Fixed partitioning ratio (fraction) of available carbon allocated to roots'
        ,'b_3': 'Fixed partitioning ratio (fraction) of available carbon allocated to wood'
        ,'c_1': 'Foliage cycling rate' # "day^{-1}"
        ,'c_2': 'Woody cycling rate' # "day^{-1}"
        ,'c_3': 'Fine roots cycling rate' # "day^{-1}"
        ,'c_4': 'Metabolic litter cycling rate' # "day^{-1}"
        ,'c_5': 'Structural litter cycling rate' # "day^{-1}"
        ,'c_6': 'Fast SOM cycling rate' # "day^{-1}"
        ,'c_7': 'Slow SOM cycling rate' # "day^{-1}"
        ,'c_8': 'Passive SOM cycling rate' # "day^{-1}"
        ,'f_41': 'Transfer coefficient from Foliage to Metabilic Litter'
        ,'f_51': 'Transfer coefficient from Foliage to Structural Litter'
        ,'f_52': 'Transfer coefficient from Wood to Structural Litter'
        ,'f_43': 'Transfer coefficient from Fine Roots to Metabolic Litter'
        ,'f_53': 'Transfer coefficient from Fine Roots to Structural Litter'
        ,'f_64': 'Transfer coefficient from Metabolic Litter to Fast SOM'
        ,'f_65': 'Transfer coefficient from Structural Litter to Fast SOM'
        ,'f_75': 'Transfer coefficient from Structural Litter to Slow SOM'
        ,'f_76': 'Transfer coefficient from Fast to Slow SOM'
        ,'f_86': 'Transfer coefficient from Fast to Passive SOM'
        ,'f_67': 'Transfer coefficient from Slow to Fast SOM'
        ,'f_87': 'Transfer coefficient from Slow to Passive SOM'
        ,'f_68': 'Transfer coefficient from Passive to Fast SOM'
}

for name in sym_dict.keys():
    var(name)

R = 8.314
Arrhenius = exp((E_p * (T_k - 298))/(R * T_k * 298))
I = I_0 * exp(-k * L)
J_e = ((alpha_q * I * J_m)/(sqrt((J_m)^2 * (alpha_q)^2 * I^2))) * ((C_i - Gamma)/(4*(C_i + 2 * Gamma)))
J_c = (V_m * (C_i - Gamma))/(C_i + K_c *( 1 + (O_x/K_o) ))
A = Min(J_c, J_e) - R_d
g_l * A /((C_i - Gamma) *(1+(D/D_0)))
A_n = G_s*(C_a - C_i)
A_c = A_n *(1-exp(-k*L))/k
GPP = A_c*3600*12/1000000 # Scaling expression from TECO fortran code, line 667.' # gC*day^{-1} 
f_W = Min((0.5*W),1)
f_T = Q_10*((T-10)/10)
xi = f_W*f_T 
t = TimeSymbol("t") # unit: "day"
x = StateVariableTuple((x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8))
u = GPP
b = (b_1, b_2, b_3, 0, 0, 0, 0, 0)
Input = InputTuple(u * ImmutableMatrix(b))
B = CompartmentalMatrix(
[[  -c_1,     0,     0,     0,     0,     0,     0,     0],
[   0,    -c_2,     0,     0,     0,     0,     0,     0],
[   0,     0,    -c_3,     0,     0,     0,     0,     0],
[f_41,     0,  f_43,    -c_4,     0,     0,     0,     0],
[f_51,  f_52,  f_53,     0,    -c_5,     0,     0,     0],
[   0,     0,     0,  f_64,  f_65,    -c_6,  f_67,  f_68],
[   0,     0,     0,     0,  f_75,  f_76,    -c_7,     0],
[   0,     0,     0,     0,     0,  f_86,  f_87,    -c_8]])
np1 = NumericParameterization(
    par_dict={GPP: 3.370, b_1: 0.14, b_2: 0.26, b_3: 0.14, f_41: 0.9, f_51: 0.1, f_52: 1, f_43: 0.2, f_53: 0.8, f_64: 0.45, f_65: 0.275, f_75: 0.275, f_76: 0.296, f_86: 0.004, f_67: 0.42, f_87: 0.01, f_68: 0.45, c_1: 0.00258, c_2: 0.0000586, c_3: 0.002390, c_4: 0.0109, c_5: 0.00095, c_6: 0.0105, c_7: 0.0000995,c_8: 0.0000115
},
    func_dict=frozendict({})
)
# "Initial values as in Wang and Luo"
nsv1 = NumericStartValueDict({x_1: 250, x_2: 4145, x_3: 192, x_4: 93, x_5: 545, x_6: 146, x_7: 1585, x_8: 300
})
ntimes = NumericSimulationTimes(np.arange(0,200,1))

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="TECO",
        longName="Terrestrial Ecosystem Model", 
        version="",
        entryAuthor="Carlos A. Sierra",
        entryAuthorOrcid="0000-0003-0009-4169",
        entryCreationDate="12/4/2018",
        doi="",
#        bibtex= "@incollection{Luo2012TE,
#                 Address = {Berkeley},
#                 Author = {Yiqi Luo and Ensheng Weng and Yuanhe Yang},
#                 Booktitle = {Encyclopedia of Theoretical Ecology},
#                 Date-Added = {2015-05-05 15:20:40 +0000},
#                 Date-Modified = {2015-05-05 15:20:40 +0000},
#                 Editor = {Alan Hastings and Louis Gross},
#                 Pages = {219-229},
#                 Publisher = {University of California Press},
#                 Title = {Ecosystem Ecology},
#                 Year = {2012}}
#                 abstract = {"Ecosystem ecology is a subdiscipline of ecology that focuses on exchange of energy and materials between organisms and the environment. The materials that are commonly studied in ecosystem ecology include water, carbon, nitrogen, phosphorus, and other elements that organisms use as nutrients. The source of energy for most ecosystems is solar radiation. In this entry, material cy-cling and energy exchange are generally described before the carbon cycle is used as an example to illustrate our quantitative and theoretical understanding of ecosystem ecology."}",
#        sym_dict=sym_dict
    ),
    B,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((x_1, x_2, x_3)),
    np1,
    nsv1,
    ntimes
})

