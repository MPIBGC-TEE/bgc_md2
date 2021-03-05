from sympy import var, ImmutableMatrix, diag, exp, Piecewise
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
        ,'C_wl': 'Carbon content of woody litter' #unit: "kgC*m^{-2}"
        ,'C_a': 'Carbon in active SOM' # See Fig. 1
        ,'C_s': 'Carbon content of slow SOM' #unit: "kgC*m^{-2}"
        ,'C_p': 'Carbon content of passive SOM' #unit: "kgC*m^{-2}"
        ,'C_u': 'Carbon in surface structural litter' # See Fig. 1
        ,'C_v': 'Carbon in soil structural litter' # See Fig. 1
        ,'C_m': 'Carbon in surface metabolic litter' # See Fig. 1
        ,'C_n': 'Carbon in soil metabolic litter' # See Fig. 1
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
        ,'buf': 'Fraction of C released from foliage entering surface structural litter'
        ,'bvr': 'Fraction of C released from foliage entering soil structural litter'
        ,'bmf': 'Fraction of C released from foliage entering surface metabolic litter'
        ,'bmr': 'Fraction of C released from fine roots entering surface metabolic litter'
        ,'bam': 'Fraction of C released from surface metabolic litter entering active SOM'
        ,'ban': 'Fraction of C released from soil metabolic litter entering active SOM'
        ,'bau': 'Fraction of C released from surface structural litter entering active SOM'
        ,'bav': 'Fraction of C released from soil structural litter entering active SOM'
        ,'baw': 'Fraction of C released from woody tissue entering active SOM'
        ,'bas': 'Fraction of C released from slow SOM entering active SOM'
        ,'bsa': 'Fraction of C released from active SOM entering slow SOM'
        ,'bsu': 'Fraction of C released from surface structural litter entering slow SOM'
        ,'bsv': 'Fraction of C released from soil structural litter entering slow SOM'
        ,'bsw': 'Fraction of C released from woody tissue entering slow SOM'
        ,'bps': 'Fraction of C released from slow SOM entering passive SOM'
        ,'bpa': 'Fraction of C released from active SOM entering passive SOM'
        ,'bap': 'Fraction of C released from passive SOM entering active SOM'
        ,'s_f': 'Foliage senescence rate' #unit: "yr^{-1}" 
        ,'s_r': 'Roots senescence rate' #unit: "yr^{-1}" 
        ,'s_w': 'Wood senescence rate' #unit: "yr^{-1}" 
        ,'d_a': 'Intrinsic decomposition rate of active SOM'
        ,'d_s': 'Intrinsic decomposition rate of slow SOM'
        ,'d_p': 'Intrinsic decomposition rate of passive SOM'
        ,'d_wl': 'Intrinsic decomposition rate of woody litter'
        ,'d_u': 'Intrinsic decomposition rate of surface structural litter'
        ,'d_v': 'Intrinsic decomposition rate of soil structural litter'
        ,'d_m': 'Intrinsic decomposition rate of surface metabolic litter'
        ,'d_n': 'Intrinsic decomposition rate of soil metabolic litter'
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
x = StateVariableTuple((C_f, C_r, C_w, C_wl, C_u, C_m, C_v, C_n, C_a, C_s, C_p))
u = NPP
b = ImmutableMatrix((a_f, a_r, a_w))
Input = InputTuple(tuple(u*b)+(0,0,0,0,0,0,0,0))
A = CompartmentalMatrix(
[[   -s_f,      0,   0,       0,      0,      0,      0,      0,      0,      0,      0]
,[      0,   -s_r,   0,       0,      0,      0,      0,      0,      0,      0,      0]
,[      0,      0,-s_w,       0,      0,      0,      0,      0,      0,      0,      0]
,[      0,      0, s_w,   -d_wl,      0,      0,      0,      0,      0,      0,      0]
,[buf*s_f,      0,   0,       0,   -d_u,      0,      0,      0,      0,      0,      0]
,[bmf*s_f,bmr*s_r,   0,       0,      0,   -d_m,      0,      0,      0,      0,      0]
,[      0,bvr*s_r,   0,       0,      0,      0,   -d_v,      0,      0,      0,      0]
,[      0,      0,   0,       0,      0,      0,      0,   -d_n,      0,      0,      0]
,[      0,      0,   0,baw*d_wl,bau*d_u,bam*d_m,bav*d_v,ban*d_n,   -d_a,bas*d_s,bap*d_p]
,[      0,      0,   0,bsw*d_wl,bsu*d_u,      0,bsv*d_v,      0,bsa*d_a,   -d_s,      0]
,[      0,      0,   0,       0,      0,      0,      0,      0,bpa*d_a,bps*d_s,   -d_p]
])
t = TimeSymbol("t")

np1 = NumericParameterization(
    par_dict={
    a_f: 0.16
#    ,buf: 
#    ,bvr: 
#    ,bmf: 
#    ,bmr: 
#    ,bam: 
#    ,ban: 
#    ,bau: 
#    ,bav: 
    ,baw: 0.413 
    ,bas: 0.42
    ,bsa: 0.35
#    ,bsu: 
#    ,bsv: 
    ,bsw: 0.175
    ,bps: 0.032
    ,bpa: 0.004
    ,bap: 0.45
    ,s_f: 0.12 #"year^{-1}"
    ,s_r: 1.0 #"year^{-1}"
    ,s_w: 0.0069 #"year^{-1}"
},
    func_dict=frozendict({})
    # state_var_units=kgC*m^{-2}
    # time_unit=year
)


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
