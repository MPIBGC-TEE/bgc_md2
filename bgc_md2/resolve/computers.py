from sympy import Symbol,Matrix
import numpy as np
from typing import Tuple
from .mvars import (
    InFluxesBySymbol,
	OutFluxesBySymbol,
	InternalFluxesBySymbol,
	TimeSymbol,
	StateVariableTuple,
	CompartmentalMatrix,  
    InputTuple
)
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel


def smooth_reservoir_model_from_fluxes(
        in_fluxes:              InFluxesBySymbol,
        out_fluxes:             OutFluxesBySymbol,
        internal_fluxes:        InternalFluxesBySymbol,
        time_symbol:            TimeSymbol,
        state_variable_tuple:   StateVariableTuple
    ) -> SmoothReservoirModel:
    return SmoothReservoirModel.from_state_variable_indexed_fluxes(
        state_vector=list(state_variable_tuple),
        time_symbol=time_symbol,
        input_fluxes=in_fluxes,
        output_fluxes=out_fluxes,
        internal_fluxes=internal_fluxes
    )

def smooth_reservoir_model_from_input_tuple_and_matrix(
        u:              InputTuple,
        B:              CompartmentalMatrix,
        time_symbol:    TimeSymbol
        #state_variable_tuple: StateVariableTuple
    ) -> SmoothReservoirModel:
    return SmoothReservoirModel.from_B_u(
        state_vector=list(StateVariableTuple),
        time_symbol=time_symbol,
        B=B,
        u=u
    )

def compartmental_matrix_from_smooth_reservoir_model(
        smr: SmoothReservoirModel
    ) -> CompartmentalMatrix:
    return CompartmentalMatrix(smr.compartmental_matrix)

