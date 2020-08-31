from sympy import Symbol, ImmutableMatrix
import numpy as np
from typing import Tuple
from sympy.physics.units import Quantity
from .mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
    CompartmentalMatrix,
    InputTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonInputTuple,
    VegetationCarbonCompartmentalMatrix,
    NumericSimulationTimes,
    NumericParameterization,
    NumericStartValueArray,
    NumericStartValueDict,
    QuantityParameterization,
    QuantityParameterizedModel,
    QuantitySimulationTimes,
    NumericParameterizedSmoothReservoirModel,
    QuantityStartValueDict,
    StateVarUnitTuple,
)
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun


def smooth_reservoir_model_from_fluxes(
    in_fluxes: InFluxesBySymbol,
    out_fluxes: OutFluxesBySymbol,
    internal_fluxes: InternalFluxesBySymbol,
    time_symbol: TimeSymbol,
    state_variable_tuple: StateVariableTuple,
) -> SmoothReservoirModel:
    return SmoothReservoirModel.from_state_variable_indexed_fluxes(
        state_vector=list(state_variable_tuple),
        time_symbol=time_symbol,
        input_fluxes=in_fluxes,
        output_fluxes=out_fluxes,
        internal_fluxes=internal_fluxes,
    )


def smooth_reservoir_model_from_input_tuple_and_matrix(
    u: InputTuple,
    B: CompartmentalMatrix,
    time_symbol: TimeSymbol,
    state_variable_tuple: StateVariableTuple,
) -> SmoothReservoirModel:
    return SmoothReservoirModel.from_B_u(
        state_vector=ImmutableMatrix(state_variable_tuple),
        time_symbol=time_symbol,
        B=B,
        u=ImmutableMatrix(u),
    )


def compartmental_matrix_from_smooth_reservoir_model(
    smr: SmoothReservoirModel,
) -> CompartmentalMatrix:
    return CompartmentalMatrix(smr.compartmental_matrix)


def vegetation_carbon_input_tuple_from_vegetation_carbon_input_partinioning_tuple_and_vegetation_carbon_input_scalar(
    u: VegetationCarbonInputScalar, b: VegetationCarbonInputPartitioningTuple
) -> VegetationCarbonInputTuple:
    return VegetationCarbonInputTuple(ImmutableMatrix(b) * u)


# def vegetation_carbon_compartmental_matrix_from_compartmental_matrix_and_vegetation_carbon_state_variable_tuple(
#       B:       CompartmentalMatrix,
#       svt:     StateVariableTuple,
#       vcsvt:   VegetationCarbonStateVariableTuple
#    ) ->
#    return CompartmentalMatrix(smr.compartmental_matrix)


def numeric_model_run_1(
    npsrm: NumericParameterizedSmoothReservoirModel,
    start_values_num: NumericStartValueArray,
    times_num: NumericSimulationTimes,
) -> SmoothModelRun:
    return SmoothModelRun(
        npsrm.srm,
        npsrm.parameterization.par_dict,
        start_values_num,
        times_num,
        npsrm.parameterization.func_dict,
    )


# def numeric_model_run_2(
#        srm:        SmoothReservoirModel,
#        para_num:   NumericParameterization,
#        start_values_num: NumericStartValueDict,
#        times_num: NumericSimulationTimes
#    ) -> SmoothModelRun:
#    return SmoothModelRun(
#        srm,
#        para_num.par_dict,
#        start_values_num,
#        times_num,
#        para_num.func_dict
#    )
#
def numeric_parameterized_smooth_reservoir_model_1(
    srm: SmoothReservoirModel, para_num: NumericParameterization,
) -> NumericParameterizedSmoothReservoirModel:
    return NumericParameterizedSmoothReservoirModel(srm, para_num)


def numeric_start_value_aray_1(
    nsvd: NumericStartValueDict,
    svt: StateVariableTuple
) -> NumericStartValueArray:
    tup = tuple(nsvd[k] for k in svt)
    return NumericStartValueArray(tup)

def quantity_parameterization_1(
        np: NumericParameterization,
        state_var_units: StateVarUnitTuple,
        time_unit: Quantity
    ) -> QuantityParameterization:
    return QuantityParameterization(
        np.par_dict, 
        np.func_dict,
        state_var_units,
        time_unit
    )
