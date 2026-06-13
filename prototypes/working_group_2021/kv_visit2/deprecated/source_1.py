from sympy import  Symbol, Function, exp, diff, Piecewise
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import module_computers
from bgc_md2.models.BibInfo import BibInfo
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    Temperature,
    StateVariableTuple,
    VegetationCarbonStateVariableTuple,
    SoilCarbonStateVariableTuple,
    CarbonStateVariableTuple,
    LuoXiBySymbol
)
import bgc_md2.resolve.computers as bgc_c

# Make a small dictionary for the variables we will use
sym_dict={
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
    'beta_leaf': '',
    'beta_wood': '',
    'mrso': 'soil moisture',
    'TAS': 'air temperature',
}
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    'xi_a': 'a scalar function of temperature and moisture and thereby ultimately of time',
    'NPP': 'net primary production',
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")
xi=xi_a(TAS,mrso)
beta_root = 1.0 - (beta_leaf + beta_wood)
svt= (
            C_leaf,
	        C_wood,
	        C_root,
	        C_leaf_litter,
	        C_wood_litter,
	        C_root_litter,
	        C_soil_fast,
	        C_soil_slow,
	        C_soil_passive,
)
mvs = CMTVS(
    {
        t,
        Temperature(TAS),
        StateVariableTuple(svt),
        CarbonStateVariableTuple(svt), # the same since we here have only carbon pools
        VegetationCarbonStateVariableTuple((
            C_leaf,
	        C_wood,
	        C_root,
        )),
        SoilCarbonStateVariableTuple((
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
                C_leaf_litter: r_C_leaf_litter_rh*C_leaf_litter*xi,
                C_wood_litter: r_C_wood_litter_rh*C_wood_litter*xi,
                C_root_litter: r_C_root_litter_rh*C_root_litter*xi,
                C_soil_fast: r_C_soil_fast_rh*C_soil_fast*xi,
                C_soil_slow: r_C_soil_slow_rh*C_soil_slow*xi,
                C_soil_passive: r_C_soil_passive_rh*C_soil_passive*xi,
            }
        ),
        #LuoXiBySymbol(
        #    {
        #        C_leaf_litter: xi,
        #        C_wood_litter: xi,
        #        C_root_litter: xi,
        #        C_soil_fast: xi,
        #        C_soil_slow: xi,
        #        C_soil_passive: r_C_soil_passive_rh*C_soil_passive*xi,
        #    }
        #),
        InternalFluxesBySymbol(
            {
                (C_leaf, C_leaf_litter): r_C_leaf_2_C_leaf_litter*C_leaf, 
                (C_wood, C_wood_litter): r_C_wood_2_C_wood_litter*C_wood, 
                (C_root, C_root_litter): r_C_root_2_C_root_litter*C_root, 
                (C_leaf_litter, C_soil_fast)    : r_C_leaf_litter_2_C_soil_fast * C_leaf_litter*xi,
                (C_leaf_litter, C_soil_slow)    : r_C_leaf_litter_2_C_soil_slow * C_leaf_litter*xi,
                (C_leaf_litter, C_soil_passive) : r_C_leaf_litter_2_C_soil_passive * C_leaf_litter*xi,
                (C_wood_litter, C_soil_fast)    : r_C_wood_litter_2_C_soil_fast * C_wood_litter*xi,
                (C_wood_litter, C_soil_slow)    : r_C_wood_litter_2_C_soil_slow * C_wood_litter*xi,
                (C_wood_litter, C_soil_passive) : r_C_wood_litter_2_C_soil_passive * C_wood_litter*xi,
                (C_root_litter, C_soil_fast)    : r_C_root_litter_2_C_soil_fast * C_root_litter*xi,
                (C_root_litter, C_soil_slow)    : r_C_root_litter_2_C_soil_slow * C_root_litter*xi,
                (C_root_litter, C_soil_passive) : r_C_root_litter_2_C_soil_passive * C_root_litter*xi,
            }
        ),
        BibInfo(# Bibliographical Information
            name="Visit",
            longName="",
            version="1",
            entryAuthor="Kostiantyn Viatkin",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        ),
    },


    computers=module_computers(bgc_c)
)
