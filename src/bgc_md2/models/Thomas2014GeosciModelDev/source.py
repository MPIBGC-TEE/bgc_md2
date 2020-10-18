from sympy import var, ImmutableMatrix, diag
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
        - C_leaf:
            desc: Carbon in foliage
        - C_wood: 
            desc: Carbon in wood 
        - C_root:
            desc: Carbon in roots 
        - C_labile:
            desc: Labile carbon 
        - C_bud:
            desc: Bud carbon 
        - C_labileRa:
            desc: Maintenance respiration pool
#        - N_leaf:
#        - N_wood:
#        - N_root:
#        - N_labile:
#        - N_bud:
    
    - photosynthetic_parameters:
        - GPP:
            desc: Photosynthesis; based on ACM model (see article for description) 
            unit: "gC*day^{-1}"
            exprs: 
            type: variable
            key: GPP
    - allocation_fluxes:
        - a_budC2leaf:
            desc: Allocation from bud C pool to leaf C
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
        - a_woodC:
            desc: Allocation from labile C to wood C 
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
        - a_rootC:
            desc: Allocation from labile C to root C
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
        - a_budC2Ramain:
            desc: Allocation of bud C pool to maintenance respiration pool when maintain respiration pool reaches zero; represents forgoing future leaf C to prevent carbon starvation.
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
        - a_budC:
            desc: Allocation of labile C to bud C; a fraction of the potential maximum leaf C
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
        - a_Ramain:
            desc: Allocation of labile C to future maintenance respiration; helps prevent carbon starvation during periods of negative NPP 
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
#        - a_budN2leaf:
#            desc: Allocation from bud N pool to leaf C (???); bud N is set in previous year
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
#        - a_budN2Ramain:
#            desc: When bud C is used for maintenance respiration (a$_budC2Ramain$ > 0), bud N is returned to the labile N pool 
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
#        - a_budN:
#            desc: Allocation of labile N to bud N; in seasonal environments it occurs in year prior to being displayed as leaf N 
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
#        - a_woodN:
#            desc: Allocation from labile N to wood N
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
#        - a_rootN:
#            desc: Allocation from labile N to root N (???)
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
        - a_labileRamain:
            desc: Allocation of labile C to respiration of living tissues # Added by Vero
            unit: "gC*m^{-2}*day^{-1}" # Added by Vero
            type: variable
#    - nitrogen_uptake_and_fixation:
#        - U_NH4:
#            desc: "Uptake of NH$_4^+$ from mineral soil NH$_4^+$"
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
#            doi: 10.1007/BF00015315
#        - U_NO3:
#            desc: "Uptake of NO$_3^-$ from mineral soil NO$_3^-$"
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
#            doi: 10.1007/BF00015315
#        - U_Nfix:
#            desc: "Fixation of N from N$_2$; function of Ra$_excess$ flux, temperature, N demand, and C cost"
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
    - turnover_fluxes:
        - tau_leaf:
            desc: Turnover of leaf (C and N) 
            unit: "day^{-1}"
            type: parameter
        - tau_wood:
            desc: Turnover of wood (C and N) 
            unit: "day^{-1}"
            type: parameter
        - tau_root:
            desc: Turnover of root (C and N)
            unit: "day^{-1}"
            type: parameter
        - t_leafC:
            desc: Turnover of leaf C to litter C; constant over year in humid tropics; seasonal otherwise
            exprs: "t_leafC = C_leaf*tau_leaf" # if Day Of the Year (DOY) > DOY_senesc. t_leafC = 0, otherwise. 
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
        - t_woodC:
            desc: Turnover of wood C to CWDC pool; occurs throughout year
            exprs: "t_woodC = C_wood*tau_wood"
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
        - t_rootC:
            desc: Turnover of root C to litter C; occurs throughout year
            exprs: "t_rootC = C_root*tau_root"
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
#        - t_retransN:
#            desc: Reabsorption of N from leaves to labile N
#            unit: "gN*m^{-2}*day^{-1}"
#        - t_leafN:
#            desc: Turnover of leaf N to litter N; constant over year in humid tropics; seasonal otherwise
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
#        - t_woodN:
#            desc: Turnover of wood N to CWDN pool; occurs throughout year
#            exprs: "t_woodN = N_wood*tau_wood"
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
#        - t_rootN:
#            desc: Turnover of root N to litter N; occurs throughout year
#            exprs: "t_rootN = N_root*tau_root"
#            unit: "gN*m^{-2}*day^{-1}"
#            type: variable
    - respiration_fluxes:
        - Ra_growth:
            desc: Growth respiration that occurs when tissue is allocated; a constant fraction of carbon allocated to tissue
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
        - Ra_excess:
            desc: Respiration that occurs when labile C exceeds a maximum labile C store; used for N fixation
            unit: "gC*m^{-2}*day^{-1}"
            type: variable
        - Ra_main:
            desc: Respiration of living tissues; a function of N content and temperature
            unit: "gC*m^{-2}*day^{-1}"
