from sympy import var, symbols, Symbol, exp, log
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonStateVariableTuple,
)
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict={
        'x_1': 'Leaves' # Pg C
        ,'x_2': 'Roots' # Pg C
        ,'x_3': 'Wood' # Pg C
        ,'x_4': 'Litter 1' # Pg C
        ,'x_5': 'Litter 2' # Pg C
        ,'x_6': 'Litter 3' # Pg C
        ,'x_7': 'Soil 1' # Pg C
        ,'x_8': 'Soil 2' # Pg C
        ,'x_9': 'Soil 3' # Pg C
        ,'b_11': 'Cycling rate for leaf pool' # "yr^{-1}"
        ,'b_22': 'Cycling rate for root pool' # "yr^{-1}"
        ,'b_33': 'Cycling rate for wood pool' # "yr^{-1}"
        ,'b_44': 'Cycling rate for litter pool 1' # "yr^{-1}"
        ,'b_55': 'Cycling rate for litter pool 2' # "yr^{-1}"
        ,'b_66': 'Cycling rate for litter pool 3' # "yr^{-1}"
        ,'b_77': 'Cycling rate for soil pool 1' # "yr^{-1}"
        ,'b_88': 'Cycling rate for soil pool 2' # "yr^{-1}"
        ,'b_99': 'Cycling rate for soil pool 3' # "yr^{-1}"
        ,'sigma': 'Sensitivity of global temperatures to atmospheric CO2  $x_a$'
        ,'T_s0': 'Mean land surface temperature in 1850'
        ,'x_a': 'Atmospheric CO2 content'
        ,'b_41': 'Transfer rate from Leaves to Litter 1' # "yr^{-1}"
        ,'b_51': 'Transfer rate from Leaves to Litter 2' # "yr^{-1}"
        ,'b_42': 'Transfer rate from Roots to Litter 1' # "yr^{-1}"
        ,'b_52': 'Transfer rate from Roots to Litter 2' # "yr^{-1}"
        ,'b_63': 'Transfer rate from Wood to Litter 3' # "yr^{-1}"
        ,'b_74': 'Transfer rate from Litter 1 to Soil 1' # "yr^{-1}"
        ,'b_75': 'Transfer rate from Litter 2 to Soil 1' # "yr^{-1}"
        ,'b_85': 'Transfer rate from Litter 2 to Soil 2' # "yr^{-1}"
        ,'b_76': 'Transfer rate from Litter 3 to Soil 1' # "yr^{-1}"
        ,'b_86': 'Transfer rate from Litter 3 to Soil 2' # "yr^{-1}"
        ,'b_87': 'Transfer rate from Soil 1 to Soil 2' # "yr^{-1}"
        ,'b_97': 'Transfer rate from Soil 1 to Soil 3' # "yr^{-1}"
        ,'b_78': 'Transfer rate from Soil 2 to Soil 1' # "yr^{-1}"
        ,'b_98': 'Transfer rate from Soil 2 to Soil 3' # "yr^{-1}"
        ,'b_79': 'Transfer rate from Soil 3 to Soil 1' # "yr^{-1}"
        ,'b_89': 'Transfer rate from Soil 3 to Soil 2' # "yr^{-1}"
        ,'T_s': 'Surface temperature'
        ,'xi_b': 'Scaling of decomposition rates at 20 degrees Celsius'
        ,'f_i': 'proportion of carbon input going to different carbon pools'
        ,'alpha': 'proportion of gross primary production that remains after respiration'
        ,'rho': 'ratio of intercellular CO$_2$ to $x_a$' # in the text it says x ???
        ,'Gamma': ''
        ,'beta': 'sensitivity of $s$ to $x_a$'
        ,'s_0': ''
        ,'s_i': 'general input function'
        ,'s_1': 'input to pool 1'
        ,'s_2': 'input to pool 2'
        ,'s_3': 'input to pool 3'
        ,'t': 'time symbol'
}
for name in sym_dict.keys():
    var(name)

x_a = 1715*exp(0.0305*t)/(1715+exp(0.0305*t)-1)+284 # from code
T_s = T_s0 + sigma/log(2)*log(x_a/285)
Gamma = 42.7 + 1.68*(T_s-25) + 0.012*(T_s-25)**2
beta = 3*rho*x_a*Gamma/((rho*x_a-Gamma)*(rho*x_a+2*Gamma))
s_i = f_i*alpha*s_0*(1+2.5*beta*log(x_a/285)) # 2.5 from code (betak)
s_1 = s_i
s_2 = s_i
s_3 = s_i
            
xi = xi_b**(0.1*T_s-1.5) #environmental effects multiplier, not for all pools. Paper says -2, code says -1.5, parameter vals say -1.5
t = TimeSymbol("t")
s = InputTuple((s_1, s_2, s_3, 0, 0, 0, 0, 0, 0))
B = CompartmentalMatrix(
[[-b_11,     0,     0,        0,        0,        0,        0,        0,        0], 
                                [    0, -b_22,     0,        0,        0,        0,        0,        0,        0],
                                [    0,     0, -b_33,        0,        0,        0,        0,        0,        0],
                                [ b_41,  b_42,     0, -b_44*xi,        0,        0,        0,        0,        0],
                                [ b_51,  b_52,     0,        0, -b_55*xi,        0,        0,        0,        0],
                                [    0,     0,  b_63,        0,        0, -b_66*xi,        0,        0,        0],
                                [    0,     0,     0,  b_74*xi,  b_75*xi,  b_76*xi, -b_77*xi,  b_78*xi,  b_79*xi], 
                                [    0,     0,     0,        0,  b_85*xi,  b_86*xi,  b_87*xi, -b_88*xi,  b_89*xi],
                                [    0,     0,     0,        0,        0,        0,  b_97*xi,  b_98*xi, -b_99*xi]])
C = StateVariableTuple((x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9))

mvs=CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="CASA",
            longName="Carnegie-Ames-Stanford approach", 
            version="2016",
            entryAuthor="Holger Metzler",
            entryAuthorOrcid="0000-0002-8239-1601",
            entryCreationDate="18/01/2018",
            doi="10.1007/s00285-016-0990-8",
            sym_dict=sym_dict
            
        ),
        # the following variables constitute the compartmental system:
        s,  # the overall input
        B,  # the overall compartmental matrix
        C, #state vector of the complete system
        VegetationCarbonStateVariableTuple((x_1, x_2, x_3))
    },
    bgc_md2_computers()
)
