from sympy import var, ImmutableMatrix, diag, exp, Piecewise
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

sym_dict={
        'age': 'Age of the stand' #changed symbol in comparison to yaml to separate from TimeSymbol
        ,'t_1': 'Age of the stand at which $\\epsilon_{0}$ begins to decline'
        ,'t_2': 'Age of the stand at which $\\epsilon_{0}$ reaches a minimum'
        ,'C_f': 'Foliar carbon mass' #unit: "kgC*m^{-2}"
        ,'C_r': 'Root carbon' #unit: "kgC*m^{-2}"
        ,'C_w': 'Carbon in woody tissue' #unit: "kgC*m^{-2}"
        ,'C_sw': 'Sapwood carbon content' #unit: "kgC*m^{-2}"
        ,'N_f': 'Nitrogen content of foliage' #unit: "kgC*m^{-2}"
        ,'N_r': 'Nitrogen content of fine roots'
        ,'n_f': 'Foliar N:C ratio'
        ,'n_crit': 'Foliar N:C ratio below which production is N-limited'
        ,'T_a': 'Mean air temperature'
        ,'Q_10': ''
        ,'Q_010': ''
        ,'R_c': 'Total construction respiration'
        ,'R_0': 'Respiration rate per unit nitrogen content corresponding to a temperature of 0°C' #unit: "kgC*kgN^{-1}*year^{-1}"
        ,'R_mf': 'Annual maintenance respiration rate of foliage (dark period only)' 
        ,'R_mr': 'Annual maintenance respiration rate of fine roots (dark period only)'
        ,'R_msw': 'Annual maintenance respiration rate of sapwood (dark period only)'
        ,'R_m': 'Total maintenance respiration'
        ,'I_0': 'Incident PAR' #unit: "GJ*m^{-2}"
        ,'sigma': 'Leaf area per unit carbon' #unit: "m^{2}*kgC^{-1}"
        ,'k': 'Light extinction coefficient' #unit: "kgC*m^{-2}"
        ,'APAR': 'Absorbed photosynthetically active radiation'
        ,'E_nf': 'Function that represents the dependence of NPP on foliar N:C ratio (n_f)'
        ,'epsilon_young': 'Maximum gross PAR utilization efficiency of young stands' #unit: "gC*MJ^{-1}"
        ,'epsilon_old': 'Maximum gross PAR utilization efficiency of old stands' #unit: "gC*MJ^{-1}"
        ,'epsilon_0': 'Maximum gross PAR utilization efficiency' #unit: "gC*MJ^{-1}"
        ,'GPP': 'Gross primary production'
        ,'NPP': 'Annual net primary production' #unit: "kgC*m^{-2}*year^{-1}"
        ,'a_f': 'Allocation fraction to foliar biomass'
        ,'a_r': 'Allocation fraction to roots biomass'
        ,'a_w': 'Allocation fraction to wood (in stem, branches and large structurl roots) biomass'
        ,'gamma_f': 'Foliage senescence rate' #unit: "yr^{-1}" 
        ,'gamma_r': 'Roots senescence rate' #unit: "yr^{-1}" 
        ,'gamma_w': 'Wood senescence rate' #unit: "yr^{-1}" 
}

for name in sym_dict.keys():
    var(name)

C_sw = 1.11*C_w**0.77
R_mf = 0.5*R_0*N_f*Q_10**(T_a/10)
R_mr = R_0*N_r*Q_10**(T_a/10)
R_msw = 0.00876*C_sw*Q_010**(T_a/10)
R_m = R_mf + R_mr + R_msw
APAR = I_0*(1-exp(-k*sigma*C_f))
E_nf = Piecewise((((((1.84*n_f)-0.01)/(0.017+n_f))/(((1.84*n_crit)-0.01)/(0.017+n_crit))),n_f<n_crit),(1,n_f>n_crit))
epsilon_0 = Piecewise((epsilon_young,age<=t_1),(Piecewise(((epsilon_young - ((epsilon_young-epsilon_old)*((age-t_1)/(t_2-t_1)))),t_1<age),(Piecewise(((epsilon_young - ((epsilon_young-epsilon_old)*((age-t_1)/(t_2-t_1)))),age<t_2),(epsilon_old,age>=t_2)),True)),True))
GPP = epsilon_0*E_nf*APAR
NPP = GPP -(R_c+R_m)
a_w = 1-a_f-a_r
x = StateVariableTuple((C_f, C_r, C_w))
u = NPP
b = (a_f, a_r, a_w)
Input = InputTuple(u*ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-gamma_f, -gamma_r, -gamma_w)
)
t = TimeSymbol("t")

# Parameter sets not working because some of the symbols have 2 values.
#model_run_data:
#    parameter_sets:
#        - "Original dataset of the publication":
#            values: {R_c: 0.25*NPP,Q_10: 2,Q_010: 1.94,n_crit: 0.034,I_0: 1.164,k: 0.5,R_0: 27,T_a: 3.8,epsilon_0: (1.05,1.25),sigma: 7.6,a_f: (0.16,0.19),a_r: (0.42,0.58),a_w: (0.23,0.42)}
#            doi: 10.1016/S0304-3800(00)00345-8

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="",
        longName="", 
        version="",
#        basedOn="G'DAY"
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="16/3/2016",
        doi="10.1016/S0304-3800(00)00345-8",
#        ,modApproach="process based"
#        ,# Age-dependent
#        ,partitioningScheme="fixed"
#        ,#partitioningScheme="semi_dynamic"
#        ,claimedDynamicPart="no"
#        ,spaceScale="global"
#        ,#    unit= "1°"
#        ,timeResolution="yearly"
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_f, C_r, C_w)),
})
