from sympy import  Symbol, Function
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import module_computers
from bgc_md2.models.BibInfo import BibInfo
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

# Make a small dictionary for the variables we will use
sym_dict={
    "r_vl_2_vw": "internal flux rate from leaf to wood", 
    "r_vw_2_vl": "internal flux rate from wood to leaf", 
    'C_soil_fast': '',
    'C_soil_slow': '',
    'C_soil_passive': '',
    'C_leaf': '',
    'C_root': '',
    'C_wood': '',
    'C_leaf_litter': '',
    'C_root_litter': '',
    'C_wood_litter': '',
    'r_C_leaf_2_C_leaf_litter': '',
    'r_C_root_2_C_root_litter': '',
    'r_C_wood_2_C_wood_litter': '',
    'r_C_leaf_litter_rh': '',
    'r_C_root_litter_rh': '',
    'r_C_wood_litter_rh': '',
    'r_C_soil_fast_rh': '',
    'r_C_soil_slow_rh': '',
    'r_C_soil_passive_rh': '',
    'r_C_leaf_litter_2_C_soil_fast': '',
    'r_C_leaf_litter_2_C_soil_slow': '',
    'r_C_leaf_litter_2_C_soil_passive': '',
    'r_C_wood_litter_2_C_soil_fast': '',
    'r_C_wood_litter_2_C_soil_slow': '',
    'r_C_wood_litter_2_C_soil_passive': '',
    'r_C_root_litter_2_C_soil_fast': '',
    'r_C_root_litter_2_C_soil_slow': '',
    'r_C_root_litter_2_C_soil_passive': '',
    'env_modifier': '',
    'tas': 'air temperature',
    'mrso': '',
    'T_0': '',
    'E': '',
    'KM': '',
    'beta_leaf': '',
    'beta_wood': '',
}
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    'xi': 'a scalar function of temperature and moisture and thereby ultimately of time',
    'NPP': '',
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")
beta_root = 1.0- (beta_leaf+beta_wood)
mvs = CMTVS(
    {
        t,
        StateVariableTuple((
            C_leaf,
	        C_wood,
	        C_root,
	        C_leaf_litter,
	        C_wood_litter,
	        C_root_litter,
	        C_soil_fast,
	        C_soil_slow,
	        C_soil_passive,
        )),
        InFluxesBySymbol(
            {
                C_leaf: NPP(t) * beta_leaf, 
                C_root: NPP(t) * beta_root, 
                C_wood: NPP(t) * beta_wood
            }
        ),
        OutFluxesBySymbol(
            {
                C_leaf_litter: r_C_leaf_litter_rh*C_leaf_litter*xi(t)*env_modifier,
                C_wood_litter: r_C_wood_litter_rh*C_wood_litter*xi(t)*env_modifier,
                C_root_litter: r_C_root_litter_rh*C_root_litter*xi(t)*env_modifier,
                C_soil_fast: r_C_soil_fast_rh*C_soil_fast*xi(t)*env_modifier,
                C_soil_slow: r_C_soil_slow_rh*C_soil_slow*xi(t)*env_modifier,
                C_soil_passive: r_C_soil_passive_rh*C_soil_passive*xi(t)*env_modifier,
            }
        ),
        InternalFluxesBySymbol(
            {
                (C_leaf, C_leaf_litter): r_C_leaf_2_C_leaf_litter*C_leaf, 
                (C_wood, C_wood_litter): r_C_wood_2_C_wood_litter*C_wood, 
                (C_root, C_root_litter): r_C_root_2_C_root_litter*C_root, 
                (C_leaf_litter, C_soil_fast)    : r_C_leaf_litter_2_C_soil_fast * C_leaf_litter*xi(t)*env_modifier,
                (C_leaf_litter, C_soil_slow)    : r_C_leaf_litter_2_C_soil_slow * C_leaf_litter*xi(t)*env_modifier,
                (C_leaf_litter, C_soil_passive) : r_C_leaf_litter_2_C_soil_passive * C_leaf_litter*xi(t)*env_modifier,
                (C_wood_litter, C_soil_fast)    : r_C_wood_litter_2_C_soil_fast * C_wood_litter*xi(t)*env_modifier,
                (C_wood_litter, C_soil_slow)    : r_C_wood_litter_2_C_soil_slow * C_wood_litter*xi(t)*env_modifier,
                (C_wood_litter, C_soil_passive) : r_C_wood_litter_2_C_soil_passive * C_wood_litter*xi(t)*env_modifier,
                (C_root_litter, C_soil_fast)    : r_C_root_litter_2_C_soil_fast * C_root_litter*xi(t)*env_modifier,
                (C_root_litter, C_soil_slow)    : r_C_root_litter_2_C_soil_slow * C_root_litter*xi(t)*env_modifier,
                (C_root_litter, C_soil_passive) : r_C_root_litter_2_C_soil_passive * C_root_litter*xi(t)*env_modifier,
            }
        ),
        BibInfo(# Bibliographical Information
            name="Visit",
            longName="",
            version="1",
            entryAuthor="Konstiantyn Viatkin",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        ),
    },


    computers=module_computers(bgc_c)
)
