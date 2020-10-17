from sympy import var, ImmutableMatrix, diag, exp, Min, integrate
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
)
from ..BibInfo import BibInfo 
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet

sym_dict = {
        'C_il': 'Carbon in leaves of plant functional type (PFT) i'
        ,'C_is': 'Carbon in transport tissue (mainly stems) of PFT$_i$'
        ,'C_ir': 'Carbon in fine roots of PFT$_i$'
        ,'Q_p': 'Flux density of photosynthetically active radiation absorbed by the leaf # Einstein*m^{-2}*s^{-1}'
        ,'alpha_3': 'Intrinsic quantum efficiency for CO_2 uptake in C_3 plants' # mol CO_2*Einstein^{-1}
        ,'alpha_4': ''# C4
        ,'O_2': 'Atmospheric [O_2] (value: 0.209)' # mol*mol^{-1}
        ,'tau': 'Ratio of kinetic parameters describing the partitioning of enzyme activity to carboxylase or oxygenase function'
        ,'Gamma': 'Gamma^* is the compensation point for gross photosynthesis' # mol*mol^{-1}
        ,'C_i': '[CO_2] in the intercellular air spaces of the leaf' # mol*mol^{-1}
        ,'J_e': 'Light-limited rate of photoynthesis'
        ,'J_e4': 'Rubisco-limited rate of photosynthesis (C4 plants)'
        ,'V_m': 'Maximum capacity of Rubisco to perform the carboxylase fuction' # mol CO_2*m^{-2}*s^{-1}
        ,'K_c': 'Michaelis-Menten coefficient for CO$_2$' # mol*mol^{-1}
        ,'K_o': 'Michaelis-Menten coefficient for O$_2$' # mol*mol^{-1}
        ,'J_c': 'Rubisco-limited rate of photosynthesis'
        ,'k':''  # C4
        ,'J_c4': 'CO$_2$-limited rate of photosynthesis at low [CO$_2$] (C4 plants)'
        ,'T': 'Rate of triose phosphate utilization'
        ,'J_p': ''
        ,'J_s': 'Triose phosphate-limited rate of photosynthesis'
        ,'J_i': 'Light-limited rate of photosynthesis (C4 plants)'
        ,'A_g': 'Gross photosynthesis rate per unit of area' # mol CO_2*m^{-2}*s^{-2}
        ,'gamma': 'Leaf respiration cost of Rubisco acivity'
        ,'Beta_stem': 'Maintenance respiration coefficient defined at 15°C'
        ,'Beta_root': 'Maintenance respiration coefficient defined at 15°C'
        ,'lambda_sapwood': 'Sapwood fraction of the total stem biomass (estimated from an assumed sap velocity and the maximum rate of transpiration experienced during the previous year)'
        ,'E_0': 'Temperature sensitivity factor'
        ,'T_0': 'Set to absolute zero (-273.16 °C)'
        ,'T_stem': 'Stem temperature' # °C
        ,'T_soil': 'Temperature of the soil in the rooting zone' # °C
        ,'fT_stem': 'f(T) is the Arrenhius temperature function'
        ,'fT_soil': 'f(T) is the Arrenhius temperature function'
        ,'R_leaf': 'Leaf maintenance respiration' # mol CO_2*m^{-2}*s^{-1}
        ,'R_stem': 'Stem maintenance respiration'
        ,'R_root': 'Root maintenance respiration'
        ,'A_n': 'Net leaf assimilation rate' # mol CO_2*m^{-2}*s^{-1}
        ,'GPP-i': 'Gross primary productivity'
        ,'eta': 'Fraction of carbon lost in the construction of net plant material because of growth respiration (value 0.33)'
        ,'NPP_i': 'Net Primary Production for PFT$_i$'
        ,'a_il': 'Fraction of annual NPP allocated to leaves for PFT$_i$'
        ,'a_is': 'Fraction of annual NPP allocated to stem for PFT$_i$'
        ,'a_ir': 'Fraction of annual NPP allocated to roots for PFT$_i$'
        ,'tau_il': 'Residence time of carbon in leaves for PFT$_i$'
        ,'tau_is': 'Residence time of carbon in stem for PFT$_i$'
        ,'tau_ir': 'Residence time of carbon in roots for PFT$_i$'
        ,'t': 'TimeSymbol("t")' #Fixme: Needed to include it here so that the equation for GPP would work
}

for name in sym_dict.keys():
    var(name)

Gamma = (O_2/(2*tau))
J_e = (alpha_3*Q_p*((C_i-Gamma)/(C_i+(2*Gamma))))
J_e4 = V_m
J_c = ((V_m*(C_i-Gamma))/(C_i + (K_c*(1+(O_2/K_o)))))
J_c4 = k*C_i # the compensation point is taken to be 0 for C4 plants
T = V_m/8.2
J_s = ((3*T*(1-(Gamma/C_i)))+(J_p*Gamma/C_i))
J_i = alpha_4*Q_p
A_g = Min(J_e,J_c,J_s)
# A_g = Min(J_i,J_e,J_c) #C4
fT_stem = exp(E_0*((1/(15-T_0))-(1/(T_stem - T_0))))
R_leaf = gamma*V_m
fT_soil = exp(E_0*((1/(15-T_0))-(1/(T_soil - T_0))))
R_stem = Beta_stem * lambda_sapwood * C_is * fT_stem
R_root = Beta_root * C_ir * fT_soil
A_n = A_g - R_leaf
GPP_i = integrate(A_g, t)
NPP_i = ((1-eta)*(integrate(A_g-R_leaf-R_stem-R_root,t)))

x = StateVariableTuple((C_il, C_is, C_ir))
u = NPP_i
b = (a_il, a_is, a_ir)
Input = InputTuple(u*ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-1/tau_il, -1/tau_is, -1/tau_ir)
)
t = TimeSymbol("t")

#model_run_data:
#    parameter_sets:
#        - "Tropical evergreen trees":
#            values: {a_il: 0.25,a_is: 0.5,a_ir: 0.25}

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="IBIS",
        longName="Integrated Biosphere Simulator", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="26/1/2016",
        doi="10.1029/96GB02692 ",
        #further_references=BibInfo(doi=""),
        sym_dict=sym_dict
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_il, C_is, C_ir)),
})
