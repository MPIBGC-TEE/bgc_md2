from sympy import var, ImmutableMatrix, Piecewise
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CarbonCompartmentalMatrix,
    CarbonInputTuple,
    TimeSymbol,
    StateVariableTuple,
    CarbonStateVariableTuple,
    # VegetationCarbonInputScalar,
    # VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
    NumericParameterization,
#    NumericStartValueDict,
#    NumericSimulationTimes,
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    NitrogenInFluxesBySymbol,
    NitrogenOutFluxesBySymbol,
    NitrogenInternalFluxesBySymbol,
    CarbonInFluxesBySymbol,
    CarbonOutFluxesBySymbol,
    CarbonInternalFluxesBySymbol,
)

# the next lines are unusual but possible
# since we use the conversion capabilities 
# explicitly, while they usually are applied
# by the framework
from bgc_md2.resolve.computers import (
    carbon_in_fluxes_by_symbol_1,
    carbon_out_fluxes_by_symbol_1,
    carbon_internal_fluxes_by_symbol_1
)

from ..BibInfo import BibInfo 
from ...  import helper as h
from bgc_md2.resolve.MVarSet import MVarSet

sym_dict = {
        'C_leaf': 'Carbon in foliage'
        ,'C_wood': 'Carbon in wood'
        ,'C_root': 'Carbon in roots'
        ,'C_labile': 'Labile carbon'
        ,'C_bud': 'Bud carbon'
        ,'C_labileRa': 'Maintenance respiration pool'
        ,'C_litter': 'Carbon in litter'
        ,'C_soil': 'Carbon in soil'
        ,'C_cwd': 'Carbon in coarse woody debris (cwd)'
        ,'N_leaf': 'Nitrogen in foliage'
        ,'N_wood': 'Nitrogen in wood'
        ,'N_root': 'Nitrogen in roots'
        ,'N_labile': 'Labile nitrogen'
        ,'N_bud': 'Bud nitrogen'
        ,'N_litter': 'Nitrogen in litter'
        ,'N_soil': 'Nitrogen in soil'
        ,'N_cwd': 'Nitrogen in coarse woody debris'
        ,'N_NH4': '"Mineral N (NH$_4^+$)"'
        ,'N_NO3': '"Mineral N (NO$_3^-$)"'
        ,'GPP': 'Photosynthesis; based on ACM model (see article for description)' # "gC*day^{-1}"
        ,'a_budC2leaf': 'Allocation from bud C pool to leaf C' # "gC*m^{-2}*day^{-1}"
        ,'a_woodC': 'Allocation from labile C to wood C' # "gC*m^{-2}*day^{-1}"
        ,'a_rootC': 'Allocation from labile C to root C' # "gC*m^{-2}*day^{-1}"
        ,'a_budC2Ramain': 'Allocation of bud C pool to maintenance respiration pool when maintain respiration pool reaches zero; represents forgoing future leaf C to prevent carbon starvation.' # "gC*m^{-2}*day^{-1}"
        ,'a_budC': 'Allocation of labile C to bud C; a fraction of the potential maximum leaf C' # "gC*m^{-2}*day^{-1}"
        ,'a_Ramain': 'Allocation of labile C to future maintenance respiration; helps prevent carbon starvation during periods of negative NPP' # "gC*m^{-2}*day^{-1}"
        ,'a_budN2leaf': 'Allocation from bud N pool to leaf C (???); bud N is set in previous year' # "gN*m^{-2}*day^{-1}"
        ,'a_budN2Ramain': 'When bud C is used for maintenance respiration (a$_budC2Ramain$ > 0), bud N is returned to the labile N pool' # "gN*m^{-2}*day^{-1}"
        ,'a_budN': 'Allocation of labile N to bud N; in seasonal environments it occurs in year prior to being displayed as leaf N' # "gN*m^{-2}*day^{-1}"
        ,'a_woodN': 'Allocation from labile N to wood N' # "gN*m^{-2}*day^{-1}"
        ,'a_rootN': 'Allocation from labile N to root N (???)' # "gN*m^{-2}*day^{-1}"
        ,'a_labileRamain': 'Allocation of labile C to respiration of living tissues' # "gC*m^{-2}*day^{-1}"
        ,'Ndep_NH4': '"Input of N deposition to NH$_4^+$"' # "gN*m^{-2}*day^{-1}"
        ,'Ndep_NO3': '"Input of N deposition to NO$_3^-$"' # "gN*m^{-2}*day^{-1}"
        ,'U_NH4': '"Uptake of NH$_4^+$ from mineral soil NH$_4^+$"' # "gN*m^{-2}*day^{-1}"
        ,'U_NH4_immob': '"Immobilization of NH$_4^+$ to soil N associated with the turnover of litter C and N"' # "gN*m^{-2}*day^{-1}"#
        ,'U_NO3': '"Uptake of NO$_3^-$ from mineral soil NO$_3^-$"' # "gN*m^{-2}*day^{-1}"
        ,'U_NO3_immob': '"Immobilization of NO$_3^-$ to soil N associated with the turnover of litter C and N"' # "gN*m^{-2}*day^{-1}"#
        ,'U_Nfix': '"Fixation of N from N$_2$; function of Ra$_excess$ flux, temperature, N demand, and C cost"' # "gN*m^{-2}*day^{-1}"
        ,'DOY': 'Day Of Year' # "day"
        ,'DOY_senesc': 'Day Of Year that growth ends and leaf fall begins' # "day"
        ,'tau_leaf': 'Turnover of leaf (C and N)' # "day^{-1}"
        ,'tau_wood': 'Turnover of wood (C and N)' # "day^{-1}"
        ,'tau_root': 'Turnover of root (C and N)' # "day^{-1}"
        ,'tau_excessC': 'Turnover of labile C when pool exceeds the maximum size of the labile C pool' # "day^{-1}"
        ,'tau_litter': 'Litter turnover rate' # "day^{-1}"
        ,'tau_cwd': 'Coarse woody debris turnover rate' # "day^{-1}"
        ,'tau_soil': 'Soil turnover rate' # "day^{-1}"
        ,'t_leafC': 'Turnover of leaf C to litter C; constant over year in humid tropics; seasonal otherwise' # "gC*m^{-2}*day^{-1}"
        ,'t_woodC': 'Turnover of wood C to CWDC pool; occurs throughout year' # "gC*m^{-2}*day^{-1}"
        ,'t_rootC': 'Turnover of root C to litter C; occurs throughout year' # "gC*m^{-2}*day^{-1}"
        ,'t_CWDC': 'Turnover of coarse woody debris into litter C pool' # "gC*m^{-2}*day^{-1}"
        ,'Pot_litterC_soilC': 'Potential turnover of litter C with flux to soil'
        ,'t_litterC_soilC': 'Turnover of litter C pool to soil C pool' # "gC*m^{-2}*day^{-1}"
        ,'Pot_litterC_atm': 'Potential turnover of litter C with flux to soil'
        ,'t_litterC_atm': 'Turnover of litter C pool released as heterotrophic respiration' # "gC*m^{-2}*day^{-1}"
        ,'t_soilC_atm': 'Turnover of litter C pool released as heterotrophic respiration' # "gC*m^{-2}*day^{-1}"
        ,'t_retransN': 'Reabsorption of N from leaves to labile N' # "gN*m^{-2}*day^{-1}"
        ,'t_leafN': 'Turnover of leaf N to litter N; constant over year in humid tropics; seasonal otherwise' # "gN*m^{-2}*day^{-1}"
        ,'t_woodN': 'Turnover of wood N to CWDN pool; occurs throughout year' # "gN*m^{-2}*day^{-1}"
        ,'t_rootN': 'Turnover of root N to litter N; occurs throughout year' # "gN*m^{-2}*day^{-1}"
        ,'t_litterN': 'Turnover of litter N to soil N' # "gN*m^{-2}*day^{-1}"
        ,'t_CWDN': 'Turnover of coarse woody debris N to litter N pool' # "gN*m^{-2}*day^{-1}"
        ,'Ra_growth': 'Growth respiration that occurs when tissue is allocated; a constant fraction of carbon allocated to tissue' # "gC*m^{-2}*day^{-1}"
        ,'Ra_excess': 'Respiration that occurs when labile C exceeds a maximum labile C store; used for N fixation' # "gC*m^{-2}*day^{-1}"
        ,'Ra_main': 'Respiration of living tissues; a function of N content and temperature' # "gC*m^{-2}*day^{-1}"
        ,'m_resp_frac': 'Proportion of litter C turnover respired' #unitless
        ,'Q_h': 'Soil respiration Q_10' #unitless
        ,'T_a': 'Daily air temperature' # ◦C
        ,'fx_N': 'Ratio of actual to potential immobilizations (process whereby mineral N is incorporated into organic, soil N  by microbial action), fx_N = (NH_4immob + NH_3immob)/total_immob' #See page 14, equations 66-70
        ,'DON_leach_prop': 'Proportion of soil N turnover lost through DON leaching'
        ,'L_DON': 'Production and leaching of dissolved organic N'
        ,'L_NO3': '"Leaching of NO$_3^-$"'
        ,'t_soilN': '"Mineralization of soil N to NH$_4^+$"' # see eq 72, pg 14 # "gN*m^{-2}*day^{-1}"
        ,'nitr': 'Nitrification ratio' # see equation 74 (pg 14)
        ,'Retrans_frac': 'Proportion of foliar N retranslocated to labile plant N pool'
        ,'leach_rate': 'NO$_3^-$ leaching rate'
        ,'nitr_rate': 'Nitrification rate'
}

for name in sym_dict.keys():
    var(name)

t_leafC = Piecewise((C_leaf*tau_leaf, DOY>DOY_senesc), (0, DOY<=DOY_senesc)) # "gC*m^{-2}*day^{-1}"
t_woodC = C_wood*tau_wood # "gC*m^{-2}*day^{-1}"
t_rootC = C_root*tau_root # "gC*m^{-2}*day^{-1}"
g_T = Q_h**((T_a-20)/10)
t_CWDC = C_cwd*tau_cwd*g_T # See equation 61, page 13
Pot_litterC_soilC = C_litter * tau_litter * g_T * (1 - m_resp_frac)
Pot_litterC_atm = C_litter * tau_litter * g_T * m_resp_frac
t_litterC_soilC = Pot_litterC_soilC * fx_N
t_litterC_atm = Pot_litterC_atm * fx_N
t_leafN = Piecewise((N_leaf*tau_leaf*(1-Retrans_frac), DOY>DOY_senesc), (0, DOY<=DOY_senesc))
t_retransN = Piecewise((N_leaf*tau_leaf*Retrans_frac, DOY>DOY_senesc), (0, DOY<=DOY_senesc))
t_CWDN = N_cwd*tau_cwd*g_T
t_litterN = N_litter*tau_litter*g_T # "gN*m^{-2}*day^{-1}"
t_rootN = N_root*tau_root # "gN*m^{-2}*day^{-1}"
t_woodN = N_wood*tau_wood # "gN*m^{-2}*day^{-1}"
t_soilN = N_soil*tau_soil*g_T*(1-DON_leach_prop)
nitr = N_NH4*nitr_rate*g_T #In equation <nitr_ratio>, in table 6 <nitr_rate>
L_NO3 = N_NO3*leach_rate
L_DON = N_soil*tau_soil*g_T*DON_leach_prop

xc = CarbonStateVariableTuple((C_labile, C_bud, C_leaf, C_wood, C_root, C_labileRa, C_litter, C_soil, C_cwd))
x = StateVariableTuple((C_labile, C_bud, C_leaf, C_wood, C_root, C_labileRa, C_litter, C_soil, C_cwd, N_leaf, N_wood, N_root, N_labile, N_bud, N_litter, N_soil, N_cwd, N_NH4, N_NO3))
u = GPP
b = ImmutableMatrix((1,0,0,0,0,0))
c_in_t= CarbonInputTuple(tuple(u*b)+(0,0,0))
#Input = InputTuple(tuple(u*b)+(0,0,0,0,0,0,0,0,0,0,0,0,0))
A_c = CarbonCompartmentalMatrix(
[[-(a_budC+a_rootC+a_woodC+a_labileRamain+Ra_growth+Ra_excess)/C_labile,0,0,0,0,0,0,0,0],
[                            a_budC/C_labile                           ,-(a_budC2leaf+a_budC2Ramain)/C_bud,0,0,0,0,0,0,0],
[                                   0                                  ,         a_budC2leaf/C_bud        ,-t_leafC/C_leaf,0,0,0,0,0,0],
[                            a_woodC/C_labile                          ,                  0               ,    0    ,-tau_wood,0,0,0,0,0],
[                            a_rootC/C_labile                          ,                  0               ,    0    ,    0    ,-tau_root,0,0,0,0],
[                            a_labileRamain/C_labile                   ,         a_budC2Ramain/C_bud      ,    0    ,    0    ,    0    ,-Ra_main/C_labileRa,0,0,0],
[                                   0                                  ,                  0               , t_leafC/C_leaf,    0    , tau_root,          0        ,-(t_litterC_soilC+t_litterC_atm)/C_litter,0,t_CWDC/C_cwd],
[                                   0                                  ,                  0               ,    0    ,    0    ,    0    ,          0        ,          t_litterC_soilC/C_litter       ,-t_soilC_atm/C_soil,0],
[                                   0                                  ,                  0               ,    0    , tau_wood,    0    ,          0        ,                      0                  ,0,-t_CWDC/C_cwd]
])
t = TimeSymbol("t")

np1 = NumericParameterization(
    par_dict={
    tau_leaf: 0.0027 
    ,tau_wood: 5e-05
    ,tau_root: 0.002
    ,tau_excessC: 0.05
    ,tau_litter: 0.029
    ,tau_cwd: 0.001
    ,tau_soil: 1e-04
    ,Q_h: 1.4
    ,m_resp_frac: 0.5
    ,DON_leach_prop: 0.0015
    ,leach_rate: 0.00001 # day^{-1}
    ,nitr_rate: 0.0001 # day^{-1}
},
    func_dict=frozendict({})
    # state_var_units= gC*m^{-2}
    # time_unit=day
)
in_fl_c = carbon_in_fluxes_by_symbol_1(c_in_t, xc)
internal_fl_c = carbon_internal_fluxes_by_symbol_1(A_c, xc)
out_fl_c = carbon_out_fluxes_by_symbol_1(A_c, xc)

in_fl_n = NitrogenInFluxesBySymbol({N_labile: U_Nfix, N_NH4: Ndep_NH4, N_NO3: Ndep_NO3})
out_fl_n = NitrogenOutFluxesBySymbol({N_soil: L_DON, N_NH4: U_NH4, N_NO3: U_NO3+L_NO3})
internal_fl_n = NitrogenInternalFluxesBySymbol({
    (N_bud, N_leaf): a_budN2leaf,
    (N_bud, N_labile): a_budN2Ramain,
    (N_labile, N_bud): a_budN,
    (N_leaf, N_litter): t_leafN,
    (N_leaf, N_labile): t_retransN,
    (N_labile, N_wood): a_woodN,
    (N_wood, N_cwd): t_woodN,
    (N_root, N_labile): a_rootN,
    (N_root, N_litter): t_rootN,
    (N_cwd, N_litter): t_CWDN,
    (N_litter, N_soil): t_litterN,
    (N_NH4, N_soil): U_NH4_immob,
    (N_NO3, N_soil): U_NO3_immob,
    (N_soil, N_NH4): t_soilN,
    (N_NH4, N_NO3): nitr
})
mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="ACONITE",
        longName="A new, simple model of ecosystem C–N cycling and interactions", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="29/3/2016",
        doi="10.5194/gmd-7-2015-2014",
        sym_dict=sym_dict
    ),
    InFluxesBySymbol(h.combine(in_fl_c, in_fl_n)),
    OutFluxesBySymbol(h.combine(out_fl_c, out_fl_n)),
    InternalFluxesBySymbol(h.combine(internal_fl_c, internal_fl_n)), 
    in_fl_n,
    out_fl_n,
    internal_fl_n,
    # I,
    # O,
    # Int,
    #Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    xc, # state vector of the carbon sub system
    
    A_c,  # the carbon compartmental matrix
    # in_fl_c, #alternatively to the CarbonCompartmentalMatrix
    # out_fl_c,# 
    internal_fl_c,
    VegetationCarbonStateVariableTuple((C_labile, C_bud, C_leaf, C_wood, C_root, C_labileRa)),
    # the following can be computed automatically
    # VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    # VegetationCarbonInputPartitioningTuple(b),
})
