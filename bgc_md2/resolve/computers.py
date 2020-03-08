from sympy import Symbol
from typing import Tuple
from .mvars import InFluxesBySymbol,OutFluxesBySymbol,InternalFluxesBySymbol,TimeSymbol,StateVariableTuple
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

def smooth_reservoir_model_from_fluxes(
        in_fluxes:InFluxesBySymbol
        ,out_fluxes:OutFluxesBySymbol
        ,internal_fluxes:InternalFluxesBySymbol
        ,time_symbol:TimeSymbol
        ,state_variable_tuple:StateVariableTuple
    )->SmoothReservoirModel:
    pass

