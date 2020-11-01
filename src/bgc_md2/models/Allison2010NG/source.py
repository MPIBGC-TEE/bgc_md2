import numpy as np
from sympy import var, ImmutableMatrix, exp
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
#    NumericParameterization,
#    NumericStartValueDict,
#    NumericSimulationTimes,
)
from ..BibInfo import BibInfo 
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet

sym_dict = {
        'S': 'soil organic carbon pool' # "mgC cm^{-3}"
        ,'D': 'dissolved organic carbon pool' # "mgC cm^{-3}"
        ,'B': 'microbial biomass pool' # "mgC cm^{-3}"
        ,'E': 'enzyme pool' # "mgC cm^{-3}"
        ,'R': 'ideal gas constant'
        ,'r_B': '"microbial biomass turnover rate"' # "h^{-1}"
        ,'r_E': 'enzyme production rate' # "h^{-1}"
        ,'r_L': 'enzyme loss rate' # "h^{-1}"
        ,'a_BS': 'fraction of dead microbial biomass transferred to soil organic matter' # 
        ,'V_Umax': 'reference for maximum rate of dissolved organic carbon uptake' # "h^{-1}"
        ,'V_max': 'reference for maximum rate of soil organic carbon decomposition' # "h^{-1}"
        ,'E_aU': 'activation energy to convert substrate into product' # "kJ mol^{-1}"
        ,'E_a': 'activation energy for soil organic carbon decomposition' # "kJ mol^{-1}"
        ,'epsilon_0': 'base carbon uptake efficiency' #
        ,'epsilon_s': 'carbon uptake efficieny slope' # "°C^{-1}"
        ,'K_U0': 'base half saturation constant for carbon uptkae' # "mg C cm^{-3}"
        ,'K_Us': 'half saturation constant slope for carbon uptake' # "mg C cm^{-3} °C^{-1}"
        ,'K_0': 'base half saturation constant for soil organic carbon decomposition' # "mg C cm^{-3}"
        ,'K_s': 'half saturation constant slope for soil organic carbon decomposition' # "mg C cm^{-3} °C^{-1}"
        ,'T': 'temperature' # "°C"
        ,'V_U': 'maximum dissolved organic carbon uptake rate'
        ,'V': 'maximum decomosition rate of soil orgacic carbon'
        ,'E_C': 'carbon uptake efficiency'
        ,'K_U': 'half saturation constant for carbon uptake'
        ,'K': 'half saturation constant for soil organic carbon decomposition'
        ,'I_S': 'soil organic carbon input rate' # "mg C cm^{-3} h^{-1}"
        ,'I_D': 'dissolved organic carbon input rate' # "mg C cm^{-3} h^{-1}"
}

for name in sym_dict.keys():
    var(name)
R = 0.008314 # "kJ mol^{-1} K^{-1}"
V_U = V_Umax * exp(-E_aU/(R*(T+273))) # 
V = V_max * exp(-E_a/(R*(T+273))) #
E_C = epsilon_0 + epsilon_s * T #
K_U = K_U0 + K_Us * T # "mg C cm^{-3}"
K = K_0 + K_s * T # "mg C cm^{-3}"
t = TimeSymbol("t") # unit: "hour"
x = StateVariableTuple((S, D, B, E))
u = InputTuple((I_S, I_D, 0, 0))
T_M = ImmutableMatrix([[-1,   0,     a_BS*r_B/(r_B+r_E),  0],
             [ 1,  -1, (1-a_BS)*r_B/(r_B+r_E),  1],
             [ 1, E_C,                     -1,  0],
             [ 0,   0,          r_E/(r_B+r_E), -1]])
N = ImmutableMatrix([[V*E/(K+S),             0,       0,   0],
           [         0, V_U*B/(K_U+D),       0,   0],
           [         0,             0, r_B+r_E,   0],
           [         0,             0,       0, r_L]])
B = CompartmentalMatrix(T_M * N)
#np1 = NumericParameterization(
#    par_dict={
#},
#    func_dict=frozendict({})
#)
#
#nsv1 = NumericStartValueDict({
#})
#ntimes = NumericSimulationTimes(np.arange())
# still don't know how to bring parameter sets here
# different model descriptions in Li and in original paper
# different units

# we would need to be able to split the yaml file into two versions of the model

#model-run-data:
#    parameter_sets:
#        - "Li":
#            values:
#                - I_S: 0.00014
#                - I_D: 0.00001
#                - V_max: 1
#                - V_Umax: 0.01
#                - 
#                
#            desc: "unit used here ist $mgC g^{-1}\\text{ soil}$ instead of $mgC cm^{-3}$"
#            bibtex: "@article{Li2014Biogeochemistry,
#                         author = {Li, Jianwei and Wang, Gangsheng and Allison, Steven D and Mayes, Melanie A and Luo, Yiqi},
#                         title = {Soil carbon sensitivity to temperature and carbon use efficiency compared across microbial-ecosystem models of varying complexity},
#                         journal = {Biogeochemistry},
#                         volume = {119},
#                         number = {1-3},
#                         pages = {67--84},
#                         year = {2014},
#                         publisher = {Springer}
#                        }"
##        - "Spinup":
##            values: 
##                - I_S: 0.0005
##                - I_D: 0.0005
##                - r_B: 0.0002
##                - r_E: 0.00005
##                - r_L: 0.001
##                
##                
##            desc: spinup data for enzyme model runs from supplementary material
#
#    initial_values:
##        - "Spinup":
##            values: {S: 100, D: 0.5, B: 0.5, E: 0.01}
##            desc: spinup data for enzyme model runs from supplementary material
#
#    run_times:
##        - "Spinup":
##            start: 0
##            end: 24000000
##            step_size: 
##           interval: 240000
#    possible_combinations:
#        - [,,]


mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="AWB",
        longName="", 
        version="",
        entryAuthor="Holger Metzler",
        entryAuthorOrcid="0000-0002-8239-1601",
        entryCreationDate="17/03/2016",
        doi="10.1038/ngeo846",
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
