#import numpy as np
from sympy import var, symbols, Symbol, ImmutableMatrix, diag, exp
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
        'C_f': 'Foliage carbon content per unit ground area at equilibrium'
        ,'C_r': 'Fine root carbon'
        ,'C_w': 'Carbon in woody tissue'
        ,'N_f': 'Nitrogen in foliage'
        ,'N_pore': 'Soil pore water nutrient pool'
        ,'F_i': 'Herbivore functional response'#unit: "gN*m^{-2}*day^{-1}"
        ,'N_r': 'Nitrogen in roots'
        ,'v_f': 'N:C ratio in foliage'
        ,'v_r': 'N:C ratio in fine roots'
        ,'N_w': 'Nitrogen in wood'
        ,'v_w': 'N:C ratio in wood'
        ,'v_m': 'N:C ratio for reproductive propagules'
        ,'G_0': 'Maximum possible primary production, assuming all light is captured and photosynthesizing material (foliage) is operating optimally.'
        ,'b_f': 'Converts carbon per square meter to LAI'
        ,'k_f': 'Foliage light-extinction (Beer-Lambert law) coefficient'
        ,'v_0': 'Half-saturation constant for the effect of foliar nitrogen concentration on primary production'
        ,'G': 'Net carbon production or growth per unit time'#unit: "gC*m^{-2}*day^{-1}"
        ,'g_N': 'Maximum possible nutrient uptake rate'
        ,'k_N': 'Half-saturation constant for uptake of soil porewater N'
        ,'k_r': 'Coefficient analogous to k$_f$'
        ,'b_r': 'Coefficient of fine root length per unit C'
        ,'U': 'Nutrient uptake rate of plant available nutrient. Saturated response of uptake to soil porewater concentration is assumed'#unit: "gN*m^{-2}*day^{-1}"
        ,'s_f': 'Allocation ratio of wood to foliage'
        ,'s_r': 'Allocation ratio of wood to fine roots'
        ,'eta_f': 'Allocation fraction to foliar biomass'
        ,'eta_r': 'Allocation fraction to roots biomass'
        ,'eta_w': 'Allocation fraction to wood (in stem, branches and large structurl roots) biomass'
        ,'eta_m': 'Allocation fraction to reproduction'
        ,'eta_d': 'Allocation fraction to plant defense'
        ,'gamma_f': 'Foliage senescence rate'#unit: "day^{-1}" 
        ,'gamma_r': 'Roots senescence rate'#unit: "day^{-1}" 
        ,'gamma_w': 'Wood senescence rate'#unit: "day^{-1}" 
}

for name in sym_dict.keys():
    var(name)

v_f = N_f/C_f
v_r = N_r/C_r
v_w = N_w/C_w
G = G_0*(1-(exp(-k_f*b_f*C_f)))*(v_f/(v_0+v_f))
U = (g_N*N_pore/(k_N+ N_pore))*(1-exp(-k_r*b_r*C_r)) # See page 4
eta_w = (s_f*eta_f)+(s_r*eta_r)
eta_d = 1 - (eta_f + eta_r + eta_w + eta_m)
u = G
x = StateVariableTuple((C_f, C_r, C_w, N_f))
b = (eta_f, eta_r, eta_w, (U/G)-eta_r*v_r-eta_w*v_w-eta_m*v_m)
Input = InputTuple(u * ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-gamma_f-(F_i/N_f), -gamma_r, -gamma_w, -gamma_f-(F_i/N_f))
)
t = TimeSymbol("t") #'yr'
# Commented out the following lines because original publication only has a few parameter values. See table 1 in page 7
#np1 = NumericParameterization(
#    par_dict={
#    v_r: 0.002,
#    v_w: 0.00005,
#    v_m: 0.005,
#    gamma_f: 0.005,
#    gamma_r: 0.01,
#    gamma_w: 0.001,
#    G_0: 30,
#    v_0: 0.02,
#    g_N: 15,
#    k_N: 5,
#    b_f: 0.004,
#    k_f: 0.2,
#    b_r: 0.001,
#    k_r: 0.15,
#    s_f: 0.5,
#    s_r: 0.5
#},
#    func_dict=frozendict({})
#    # state_var_units=gram/kilometer**2,
#    # time_unit=day
#)
#nsv1 = NumericStartValueDict({
#    C_f: , #"g*m^{-2}"
#    C_r: , #"g*m^{-2}"
#    C_w: , #"g*m^{-2}"
#    N_f: #"g*m^{-2}"
#})
#
#ntimes = NumericSimulationTimes(np.arange(, , ))

mvs=MVarSet({
    BibInfo(# Bibliographical Information
        name="",
#        basedOn="G'DAY",
        longName="", 
        version="",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="15/3/2016",
        doi="10.1007/s12080-011-0135-z",
        sym_dict=sym_dict
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_f, C_r, C_w, N_f)),
#    np1
})

