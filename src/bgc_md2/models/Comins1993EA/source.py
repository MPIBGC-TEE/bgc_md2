#import numpy as np
from sympy import var, symbols, Symbol, ImmutableMatrix, exp#, Rational
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
#    NumericStartValueDict,
#    NumericSimulationTimes,
   )
from ..BibInfo import BibInfo 
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict={
        'C_f': 'Foliage carbon content'
        ,'C_r': 'Coarse and fine root carbon'
        ,'C_w': 'Carbon in woody stem, branches, and large structural roots'
        ,'C_m': 'Carbon in surface metabolic litter'
        ,'C_n': 'Carbon in soil metabolic litter'
        ,'C_u': 'Carbon in surface structural litter'
        ,'C_v': 'Carbon in soil structural litter'
        ,'C_a': 'Carbon in active soil organic matter'
        ,'C_s': 'Carbon in slow soil organic matter'
        ,'C_p': 'Carbon in passive soil organic matter'
        ,'G': 'Net rate of plant carbon production'
        ,'eta_f': 'Allocation fraction to foliar biomass'
        ,'eta_r': 'Allocation fraction to roots biomass'
        ,'eta_w': 'Allocation fraction to wood (in stem, branches and large structurl roots) biomass'
        ,'gamma_f': 'Foliage senescence rate' #unit: "yr^{-1}" 
        ,'gamma_r': 'Roots senescence rate' #unit: "yr^{-1}" 
        ,'gamma_w': 'Wood senescence rate' #unit: "yr^{-1}" 
        ,'p_uf': 'Partition coefficient'
        ,'p_mf': 'Partition coefficient'
        ,'p_vr': 'Partition coefficient'
        ,'p_nr': 'Partition coefficient'
        ,'p_su': 'Partition coefficient'
        ,'p_sv': 'Partition coefficient'
        ,'p_sa': 'Partition coefficient'
        ,'p_ap': 'Partition coefficient'
        ,'p_as': 'Partition coefficient'
        ,'p_an': 'Partition coefficient'
        ,'p_av': 'Partition coefficient'
        ,'p_am': 'Partition coefficient'
        ,'p_au': 'Partition coefficient'
        ,'p_pa': 'Partition coefficient'
        ,'p_ps': 'Partition coefficient'
        ,'d_1': 'Decay rate'
        ,'d_2': 'Decay rate'
        ,'d_3': 'Decay rate'
        ,'d_4': 'Decay rate'
        ,'d_5': 'Decay rate'
        ,'d_6': 'Decay rate'
        ,'d_7': 'Decay rate'
        ,'L_fl': 'Lignin to biomass ratio in leaf litter'
        ,'L_rl': 'Lignin to biomass ratio in root litter'
        ,'omega': 'Carbon content of biomass'
        ,'lambda_f': 'Ratio of litter N:C to live leaf'
        ,'upsilon_f': 'N:C ratio in foliage'
        ,'lambda_r': 'Ratio of litter N:C to live root'
        ,'upsilon_r': 'N:C ratio in roots'
        ,'T': 'Soil texture'
        ,'T_soil': 'Average soil temperature'
        ,'AT_soil': 'Soil-temperature activity factor'
}
for name in sym_dict.keys():
    var(name)
eta_w = 1-eta_f-eta_r #Added by Vero
#SEE page 7 for equations of G, upsilon...
p_mf = 0.85 - 0.018*L_fl/(omega*lambda_f*upsilon_f)
p_uf = 1 - p_mf
p_nr = 0.85 - 0.018*L_rl/(omega*lambda_r*upsilon_r)
p_vr = 1 - p_nr 
p_au = 0.55*(1-L_fl)
p_su = 0.7*L_fl
p_av = 0.45*(1-L_fl)
p_sv = 0.7*L_fl
p_sa = 0.996-(0.85-0.68*T)
d_1 = 0.076*exp(-3*L_fl)*AT_soil
d_2 = 0.28*AT_soil 
d_3 = 0.094*exp(-3*L_rl)*AT_soil
d_4 = 0.35*AT_soil
d_5 = 0.14*(1-0.75*T)*AT_soil
d_6 = 0.0038*AT_soil
d_7 = 0.00013*AT_soil
AT_soil = 0.0326 + 0.00351*(T_soil)**(1.652) - (T_soil/41.748)**7.19
x = StateVariableTuple((C_f, C_w, C_r, C_u, C_m, C_v, C_n, C_a, C_s, C_p))
u = G
b = ImmutableMatrix((eta_f, eta_w, eta_r))
t = TimeSymbol("t") #'yr'
Input = InputTuple(tuple(u*b)+(0,0,0,0,0,0,0))
A = CompartmentalMatrix(
[[  -gamma_f  ,    0   ,      0     ,    0   ,    0   ,    0   ,    0   ,    0   ,    0   ,    0   ]
,[      0     ,-gamma_w,      0     ,    0   ,    0   ,    0   ,    0   ,    0   ,    0   ,    0   ]
,[      0     ,    0   ,  -gamma_r  ,    0   ,    0   ,    0   ,    0   ,    0   ,    0   ,    0   ]
,[gamma_f*p_uf, gamma_w,      0     ,  -d_1  ,    0   ,    0   ,    0   ,    0   ,    0   ,    0   ]
,[gamma_f*p_mf,    0   ,      0     ,    0   ,  -d_2  ,    0   ,    0   ,    0   ,    0   ,    0   ]
,[      0     ,    0   ,gamma_r*p_vr,    0   ,    0   ,  -d_3  ,    0   ,    0   ,    0   ,    0   ]
,[      0     ,    0   ,gamma_r*p_nr,    0   ,    0   ,    0   ,  -d_4  ,    0   ,    0   ,    0   ]
,[      0     ,    0   ,      0     ,d_1*p_au,d_2*p_am,d_3*p_av,d_4*p_an,  -d_5  ,d_6*p_as,d_7*p_ap]
,[      0     ,    0   ,      0     ,d_1*p_su,    0   ,d_3*p_sv,    0   ,d_5*p_sa,  -d_6  ,    0   ]
,[      0     ,    0   ,      0     ,    0   ,    0   ,    0   ,    0   ,d_5*p_pa,d_6*p_ps,  -d_7  ]
])
# Commented out the following lines because original publication only has 3 parameter values
np1 = NumericParameterization(
    par_dict={
#    G: , #"Mg*ha^{-1}*yr^{-1}"
    eta_f: 0.3 
    ,eta_r: 0.3 
    ,gamma_f: 0.5 #"yr^{-1}"
    ,gamma_r: 1.5 #"yr^{-1}"
    ,gamma_w: 0.01 #"yr^{-1}"
    ,lambda_f: 1
    ,lambda_r: 1
    ,omega: 0.45
    ,L_fl: 0.2
    ,L_rl: 0.16
    ,T: 0.5
    ,p_am: 0.45
    ,p_an: 0.45
    ,p_pa: 0.004
    ,p_as: 0.42
    ,p_ps: 0.03
    ,p_ap: 0.45
},
    func_dict=frozendict({})
    # state_var_units=gram/kilometer**2,
    # time_unit=day
)
#nsv1 = NumericStartValueDict({
#    F: , #"Mg/ha"
#    W: , #"Mg/ha"
#    R: #"Mg/ha"
#})
#
#ntimes = NumericSimulationTimes(np.arange(, , ))

mvs=CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="G'DAY",
            longName="Generic Decomposition and Yield", 
            version="1",
            entryAuthor="Verónika Ceballos-Núñez",
            entryAuthorOrcid="0000-0002-0046-1160",
            entryCreationDate="27/1/2016",
            doi="10.2307/1942099",
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
    },
    bgc_md2_computers()
)
