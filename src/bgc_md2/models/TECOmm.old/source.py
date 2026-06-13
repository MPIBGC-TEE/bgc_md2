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
        'C_foliage': 'Carbon in foliage'
        ,'C_roots': 'Carbon in roots'
        ,'C_wood': 'Carbon in woody tissue'
        ,'C_metlit': 'Carbon in metabolic litter'
        ,'C_stlit': 'Carbon in structural litter'
        ,'C_fastsom': 'Carbon in fast SOM'
        ,'C_slowsom': 'Carbon in slow SOM'
        ,'C_passsom': 'Carbon in passive SOM'
        ,'T_k': 'Canopy temperature in Kelvin'
        ,'R': 'Universal gas constant'
        ,'E_p': 'Activation energy' # It may differ for each parameter
        ,'Arrhenius': 'Arrhenius multiplier'
        ,'C_a': 'Ambient CO2 concentration' # From external atmospheric CO2 concentration data
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
        ,'O_x': 'Oxygen concentration in the air' # value: 0.21
        ,'K_o': 'Michaelis-Menten constant for oxygenation'
        ,'K_c': 'Michaelis-Menten constant for carboxylation'
        ,'J_c': 'Rate of carboxilation with CO2 limitation'
        ,'R_d': 'Dark respiration' # Find expression in fortran code
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
        ,'b_foliage': 'Fixed partitioning ratio (fraction) of available carbon allocated to foliage'
        ,'b_roots': 'Fixed partitioning ratio (fraction) of available carbon allocated to roots'
        ,'b_wood': 'Fixed partitioning ratio (fraction) of available carbon allocated to wood'
        ,'cr_foliage': 'Foliage cycling rate' # "day^{-1}"
        ,'cr_wood': 'Woody cycling rate' # "day^{-1}"
        ,'cr_fineroot': 'Fine roots cycling rate' # "day^{-1}"
        ,'cr_metlit': 'Metabolic litter cycling rate' # "day^{-1}"
        ,'cr_stlit': 'Structural litter cycling rate' # "day^{-1}"
        ,'cr_fastsom': 'Fast SOM cycling rate' # "day^{-1}"
        ,'cr_slowsom': 'Slow SOM cycling rate' # "day^{-1}"
        ,'cr_passsom': 'Passive SOM cycling rate' # "day^{-1}"
        ,'f_foliage2metlit': 'Transfer coefficient from Foliage to Metabilic Litter'
        ,'f_foliage2stlit': 'Transfer coefficient from Foliage to Structural Litter'
        ,'f_wood2stlit': 'Transfer coefficient from Wood to Structural Litter'
        ,'f_fineroots2metlit': 'Transfer coefficient from Fine Roots to Metabolic Litter'
        ,'f_fineroots2stlit': 'Transfer coefficient from Fine Roots to Structural Litter'
        ,'f_metlit2fastsom': 'Transfer coefficient from Metabolic Litter to Fast SOM'
        ,'f_stlit2fastsom': 'Transfer coefficient from Structural Litter to Fast SOM'
        ,'f_stlit2slowsom': 'Transfer coefficient from Structural Litter to Slow SOM'
        ,'f_fastsom2slowsom': 'Transfer coefficient from Fast to Slow SOM'
        ,'f_fastsom2passsom': 'Transfer coefficient from Fast to Passive SOM'
        ,'f_slowsom2fastsom': 'Transfer coefficient from Slow to Fast SOM'
        ,'f_slowsom2passsom': 'Transfer coefficient from Slow to Passive SOM'
        ,'f_passsom2fastsom': 'Transfer coefficient from Passive to Fast SOM'
}

for name in sym_dict.keys():
    var(name)

R = 8.314
Arrhenius = exp((E_p * (T_k - 298))/(R * T_k * 298))
I = I_0 * exp(-k * L)
J_e = 1#((alpha_q * I * J_m)/(sqrt((J_m)^2 * (alpha_q)^2 * I^2))) * ((C_i - Gamma)/(4*(C_i + 2 * Gamma))) #Fixme: had to set J_e to 1 because problem with sqrt AttributeError: 'Not' object has no attribute '_eval_power'
J_c = (V_m * (C_i - Gamma))/(C_i + K_c *( 1 + (O_x/K_o) ))
A = Min(J_c, J_e) - R_d
g_l * A / ((C_i - Gamma) * (1 + (D / D_0)))
A_n = G_s*(C_a - C_i)
A_c = A_n * (1 - exp(- k * L)) / k
GPP = A_c*3600*12/1000000 # Scaling expression from TECO fortran code, line 667.' # gC*day^{-1} 
f_W = Min((0.5*W), 1)
f_T = Q_10*((T-10)/10)
xi = f_W*f_T
t = TimeSymbol("t")  # unit: "day"
#x = StateVariableTuple((C_foliage, C_roots, C_wood, C_metlit, C_stlit, C_fastsom, C_slowsom, C_passsom))
x = StateVariableTuple((C_foliage, C_wood, C_roots, C_metlit, C_stlit, C_fastsom, C_slowsom, C_passsom))
u = GPP
b = (b_foliage, b_roots, b_wood, 0, 0, 0, 0, 0)
Input = InputTuple(u * ImmutableMatrix(b))
B = CompartmentalMatrix([
[  -cr_foliage,     0,     0,     0,     0,     0,     0,     0],
[   0,    -cr_wood,     0,     0,     0,     0,     0,     0],
[   0,     0,    -cr_fineroot,     0,     0,     0,     0,     0],
[f_foliage2metlit,     0,  f_fineroots2metlit,    -cr_metlit,     0,     0,     0,     0],
[f_foliage2stlit,  f_wood2stlit,  f_fineroots2stlit,     0,    -cr_stlit,     0,     0,     0],
[   0,     0,     0,  f_metlit2fastsom,  f_stlit2fastsom,    -cr_fastsom,  f_slowsom2fastsom,  f_passsom2fastsom],
[   0,     0,     0,     0,  f_stlit2slowsom,  f_fastsom2slowsom,    -cr_slowsom,     0],
[   0,     0,     0,     0,     0,  f_fastsom2passsom,  f_slowsom2passsom,    -cr_passsom]])
np1 = NumericParameterization(
    par_dict={GPP: 3.370, b_foliage: 0.14, b_roots: 0.26, b_wood: 0.14, f_foliage2metlit: 0.9, f_foliage2stlit: 0.1, f_wood2stlit: 1, f_fineroots2metlit: 0.2, f_fineroots2stlit: 0.8, f_metlit2fastsom: 0.45, f_stlit2fastsom: 0.275, f_stlit2slowsom: 0.275, f_fastsom2slowsom: 0.296, f_fastsom2passsom: 0.004, f_slowsom2fastsom: 0.42, f_slowsom2passsom: 0.01, f_passsom2fastsom: 0.45, cr_foliage: 0.00258, cr_wood: 0.0000586, cr_fineroot: 0.002390, cr_metlit: 0.0109, cr_stlit: 0.00095, cr_fastsom: 0.0105, cr_slowsom: 0.0000995,cr_passsom: 0.0000115
},
    func_dict=frozendict({})
)
# "Initial values as in Wang and Luo"
nsv1 = NumericStartValueDict({C_foliage: 250, C_roots: 4145, C_wood: 192, C_metlit: 93, C_stlit: 545, C_fastsom: 146, C_slowsom: 1585, C_passsom: 300
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
        sym_dict=sym_dict
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
    ),
    B,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_foliage, C_roots, C_wood)),
    np1,
    nsv1,
    ntimes
})

