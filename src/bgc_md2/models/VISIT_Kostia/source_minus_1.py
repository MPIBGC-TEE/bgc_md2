from sympy import symbols, Function, exp, var
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers
from bgc_md2.models.BibInfo import BibInfo 
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)

# from sympy.vector import CoordSysND,express
# fixme mm:
# add this boilerplatecode automatically
# from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
sym_dict= {
    'C_soil_fast': '',
    'C_soil_slow': '',
    'C_soil_passive': '',
    'C_leaf': '',
    'C_root': '',
    'C_wood': '',
    'C_leaf_litter': '',
    'C_root_litter': '',
    'C_wood_litter': '',
    'k_leaf_2_leaf_litter': '',
    'k_root_2_root_litter': '',
    'k_wood_2_wood_litter': '',
    'k_leaf_litter_rh': '',
    'k_root_litter_rh': '',
    'k_wood_litter_rh': '',
    'k_soil_fast_rh': '',
    'k_soil_slow_rh': '',
    'k_soil_passive_rh': '',
    'k_leaf_litter_2_soil_fast': '',
    'k_leaf_litter_2_soil_slow': '',
    'k_leaf_litter_2_soil_passive': '',
    'k_wood_litter_2_soil_fast': '',
    'k_wood_litter_2_soil_slow': '',
    'k_wood_litter_2_soil_passive': '',
    'k_root_litter_2_soil_fast': '',
    'k_root_litter_2_soil_slow': '',
    'k_root_litter_2_soil_passive': '',
    'tsl': '',
    'mrso': '',
    't': '',
    'T0': '',
    'E': '',
    'KM': '',
    'NPP': '',
    'beta_leaf': '',
    'beta_wood': '',
    'xi': '',
}
for name in sym_dict.keys():
    var(name)
#xi = exp(E * (1 / (10 - T0) - 1 / (tsl(t) - T0))) * mrso(t) / (KM + mrso(t))
# the keys of the internal flux dictionary are tuples (source_pool,target_pool)

# srm:SmoothReservoirModel
# srm=SmoothReservoirModel.from_state_variable_indexed_fluxes(
#     in_fluxes
#    ,out_fluxes
#    ,internal_fluxes
# )

beta_root = 1.0- (beta_leaf+beta_wood)
# specialVars = {
mvs = CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="Visit",
            longName="Vegetation Intergrative Simulator for Trace Gases", 
            version="1",
            entryAuthor="Konstantin Viatkin",
            entryAuthorOrcid="",
            entryCreationDate="11/24/2021",
            doi="",
            further_references=BibInfo(doi="10.1016/j.agrformet.2009.05.002"),
            #  Also from PDF in Reflex experiment
            sym_dict=sym_dict
            
        ),
        InFluxesBySymbol(
            {
                C_leaf: NPP * beta_leaf, 
                C_root: NPP * beta_root, 
                C_wood: NPP * beta_wood
            }
        ),
        OutFluxesBySymbol(
            {
                C_leaf_litter: k_leaf_litter_rh*C_leaf_litter*xi,
                C_wood_litter: k_wood_litter_rh*C_wood_litter*xi,
                C_root_litter: k_root_litter_rh*C_root_litter*xi,
                C_soil_fast: k_soil_fast_rh*C_soil_fast*xi,
                C_soil_slow: k_soil_slow_rh*C_soil_slow*xi,
                C_soil_passive: k_soil_passive_rh*C_soil_passive*xi,
            }
        ),
        InternalFluxesBySymbol(
            {
                (C_leaf, C_leaf_litter): k_leaf_2_leaf_litter*C_leaf, 
                (C_wood, C_wood_litter): k_wood_2_wood_litter*C_wood, 
                (C_root, C_root_litter): k_root_2_root_litter*C_root, 
                (C_leaf_litter, C_soil_fast)    : k_leaf_litter_2_soil_fast * C_leaf_litter*xi,
                (C_leaf_litter, C_soil_slow)    : k_leaf_litter_2_soil_slow * C_leaf_litter*xi,
                (C_leaf_litter, C_soil_passive) : k_leaf_litter_2_soil_passive * C_leaf_litter*xi,
                (C_wood_litter, C_soil_fast)    : k_wood_litter_2_soil_fast * C_wood_litter*xi,
                (C_wood_litter, C_soil_slow)    : k_wood_litter_2_soil_slow * C_wood_litter*xi,
                (C_wood_litter, C_soil_passive) : k_wood_litter_2_soil_passive * C_wood_litter*xi,
                (C_root_litter, C_soil_fast)    : k_root_litter_2_soil_fast * C_root_litter*xi,
                (C_root_litter, C_soil_slow)    : k_root_litter_2_soil_slow * C_root_litter*xi,
                (C_root_litter, C_soil_passive) : k_root_litter_2_soil_passive * C_root_litter*xi,
            }
        ),
        TimeSymbol('t'),
        StateVariableTuple(
            (
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
        )
        #srm
    },
    bgc_md2_computers()

)
