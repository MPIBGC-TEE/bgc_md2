from sympy import var, ImmutableMatrix
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
        'C_leaf': 'Carbon in foliage'
        ,'C_wood': 'Carbon in wood'
        ,'C_root': 'Carbon in roots'
        ,'C_labile': 'Labile carbon'
        ,'C_bud': 'Bud carbon'
        ,'C_labileRa': 'Maintenance respiration pool'
#        ,'N_leaf': ''
#        ,'N_wood': ''
#        ,'N_root': ''
#        ,'N_labile': ''
#        ,'N_bud': ''
        ,'GPP': 'Photosynthesis; based on ACM model (see article for description)' # "gC*day^{-1}"
        ,'a_budC2leaf': 'Allocation from bud C pool to leaf C' # "gC*m^{-2}*day^{-1}"
        ,'a_woodC': 'Allocation from labile C to wood C' # "gC*m^{-2}*day^{-1}"
        ,'a_rootC': 'Allocation from labile C to root C' # "gC*m^{-2}*day^{-1}"
        ,'a_budC2Ramain': 'Allocation of bud C pool to maintenance respiration pool when maintain respiration pool reaches zero; represents forgoing future leaf C to prevent carbon starvation.' # "gC*m^{-2}*day^{-1}"
        ,'a_budC': 'Allocation of labile C to bud C; a fraction of the potential maximum leaf C' # "gC*m^{-2}*day^{-1}"
        ,'a_Ramain': 'Allocation of labile C to future maintenance respiration; helps prevent carbon starvation during periods of negative NPP' # "gC*m^{-2}*day^{-1}"
#        ,'a_budN2leaf': 'Allocation from bud N pool to leaf C (???); bud N is set in previous year' # "gN*m^{-2}*day^{-1}"
#        ,'a_budN2Ramain': 'When bud C is used for maintenance respiration (a$_budC2Ramain$ > 0), bud N is returned to the labile N pool' # "gN*m^{-2}*day^{-1}"
#        ,'a_budN': 'Allocation of labile N to bud N; in seasonal environments it occurs in year prior to being displayed as leaf N' # "gN*m^{-2}*day^{-1}"
#        ,'a_woodN': 'Allocation from labile N to wood N' # "gN*m^{-2}*day^{-1}"
#        ,'a_rootN': 'Allocation from labile N to root N (???)' # "gN*m^{-2}*day^{-1}"
        ,'a_labileRamain': 'Allocation of labile C to respiration of living tissues' # "gC*m^{-2}*day^{-1}"
#        ,'U_NH4': '"Uptake of NH$_4^+$ from mineral soil NH$_4^+$"' # "gN*m^{-2}*day^{-1}"
#        ,'  doi: 10.1007/BF00015315
#        ,'U_NO3': '"Uptake of NO$_3^-$ from mineral soil NO$_3^-$"' # "gN*m^{-2}*day^{-1}"
#        ,'  doi: 10.1007/BF00015315
#        ,'U_Nfix': '"Fixation of N from N$_2$; function of Ra$_excess$ flux, temperature, N demand, and C cost"' # "gN*m^{-2}*day^{-1}"
        ,'tau_leaf': 'Turnover of leaf (C and N) ' # "day^{-1}"
        ,'tau_wood': 'Turnover of wood (C and N) ' # "day^{-1}"
        ,'tau_root': 'Turnover of root (C and N)' # "day^{-1}"
        ,'t_leafC': 'Turnover of leaf C to litter C; constant over year in humid tropics; seasonal otherwise'
        ,'t_woodC': 'Turnover of wood C to CWDC pool; occurs throughout year'
        ,'t_rootC': 'Turnover of root C to litter C; occurs throughout year'
#        ,'t_retransN': 'Reabsorption of N from leaves to labile N' # "gN*m^{-2}*day^{-1}"
#        ,'t_leafN': 'Turnover of leaf N to litter N; constant over year in humid tropics; seasonal otherwise' # "gN*m^{-2}*day^{-1}"
#        ,'t_woodN': 'Turnover of wood N to CWDN pool; occurs throughout year'
#        ,'t_rootN': 'Turnover of root N to litter N; occurs throughout year'
        ,'Ra_growth': 'Growth respiration that occurs when tissue is allocated; a constant fraction of carbon allocated to tissue' # "gC*m^{-2}*day^{-1}"
        ,'Ra_excess': 'Respiration that occurs when labile C exceeds a maximum labile C store; used for N fixation' # "gC*m^{-2}*day^{-1}"
        ,'Ra_main': 'Respiration of living tissues; a function of N content and temperature' # "gC*m^{-2}*day^{-1}"
}

for name in sym_dict.keys():
    var(name)

t_leafC = C_leaf*tau_leaf # if Day Of the Year (DOY) > DOY_senesc. t_leafC = 0, otherwise. # "gC*m^{-2}*day^{-1}"
t_woodC = C_wood*tau_wood # "gC*m^{-2}*day^{-1}"
t_rootC = C_root*tau_root # "gC*m^{-2}*day^{-1}"
#t_rootN = N_root*tau_root # "gN*m^{-2}*day^{-1}"
#t_woodN = N_wood*tau_wood # "gN*m^{-2}*day^{-1}"

x = StateVariableTuple((C_labile, C_bud, C_leaf, C_wood, C_root, C_labileRa))
#x = StateVariableTuple((C_leaf, C_wood, C_root, C_labile, C_bud, C_labileRa, N_leaf, N_wood, N_root, N_labile, N_bud))
u = GPP
#            exprs: "u = Matrix(11,1,[, , , GPP, , + , a_budN2leaf, a_woodN, a_rootN, U_NH4+U_NO3+U_Nfix+t_retransN+a_budN2Ramain, a_budN2leaf])"
b = (1,0,0,0,0,0)
Input = InputTuple(u*ImmutableMatrix(b))
A = CompartmentalMatrix(
[[-(a_budC+a_rootC+a_woodC+a_labileRamain+Ra_growth+Ra_excess)/C_labile,0,0,0,0,0],
                               [a_budC/C_labile,-(a_budC2leaf+a_budC2Ramain)/C_bud,0,0,0,0],
                               [0, a_budC2leaf/C_bud,-tau_leaf,0,0,0],
                               [a_woodC/C_labile,0,0,-tau_wood,0,0],
                               [a_rootC/C_labile,0,0,0,-tau_root,0],
                               [a_labileRamain/C_labile, a_budC2Ramain/C_bud, 0, 0, 0, -Ra_main/C_labileRa]])
t = TimeSymbol("t")

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
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_labile, C_bud, C_leaf, C_wood, C_root, C_labileRa)),
})
