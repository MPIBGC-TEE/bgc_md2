from sympy import Symbol,Matrix
import numpy as np
from typing import Tuple
from .mvars import (
    InFluxesBySymbol
	,OutFluxesBySymbol
	,InternalFluxesBySymbol
	,TimeSymbol
	,StateVariableTuple
	,CompartmentalMatrix
)
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

def smooth_reservoir_model_from_fluxes(
        in_fluxes:InFluxesBySymbol
        ,out_fluxes:OutFluxesBySymbol
        ,internal_fluxes:InternalFluxesBySymbol
        ,time_symbol:TimeSymbol
        ,state_variable_tuple:StateVariableTuple
    )->SmoothReservoirModel:
    n=len(state_variable_tuple)
    return SmoothReservoirModel.from_state_variable_indexed_fluxes(
        state_vector=list(state_variable_tuple)
        ,time_symbol=time_symbol
        ,input_fluxes=in_fluxes
        ,output_fluxes=out_fluxes
        ,internal_fluxes=internal_fluxes
   )

def compartmental_matrix_from_smooth_reservoir_model(
        smr:SmoothReservoirModel
   )->CompartmentalMatrix:
   return CompartmentalMatrix(smr.compartmental_matrix)

