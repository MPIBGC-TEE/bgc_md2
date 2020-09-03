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
    NumericParameterizedSmoothReservoirModel,
    NumericSolutionArray,
    QuantityParameterization,
    QuantitySimulationTimes,
    QuantityParameterizedSmoothReservoirModel,
    QuantityStartValueDict,
    QuantityStartValueArray,
    QuantityModelRun,
    QuantitySolutionArray,
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


def numeric_parameterized_smooth_reservoir_model_1(
    srm: SmoothReservoirModel, para_num: NumericParameterization,
) -> NumericParameterizedSmoothReservoirModel:
    return NumericParameterizedSmoothReservoirModel(srm, para_num)

def numeric_start_value_array_1(
    nsvd: NumericStartValueDict,
    svt: StateVariableTuple
) -> NumericStartValueArray:
    tup = tuple(nsvd[k] for k in svt)
    return NumericStartValueArray(tup)

def numeric_start_value_array_2(
    smr: SmoothModelRun
) -> NumericStartValueArray:
    tup = tuple(nsvd[k] for k in svt)
    return NumericStartValueArray(smr.start_values)

def numeric_start_value_dict(
    nsva: NumericStartValueArray,
    svt: StateVariableTuple
) -> NumericStartValueDict:
    return NumericStartValueDict({sv:nsva[i]  for i,sv in enumerate(svt)})

def numeric_solution_array_1(
    smr: SmoothModelRun
    )->NumericSolutionArray:
    return NumericSolutionArray(smr.solve())

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

def quantity_parameterized_smooth_reservoir_model_1(
    srm: SmoothReservoirModel,
    para_q: QuantityParameterization
) -> QuantityParameterizedSmoothReservoirModel:
    return QuantityParameterizedSmoothReservoirModel(srm, para_q)

def quantity_start_value_array_1(
    qsvd: QuantityStartValueDict,
    svt: StateVariableTuple
) -> QuantityStartValueArray:
    tup = tuple(qsvd[k] for k in svt)
    return QuantityStartValueArray(tup)

def quantity_model_run_1(
    qpsrm: QuantityParameterizedSmoothReservoirModel,
    start_values_q: QuantityStartValueArray,
    times_q: QuantitySimulationTimes,
) -> QuantityModelRun:
    return QuantityModelRun(
        qpsrm.srm,
        qpsrm.parameterization.par_dict,
        start_values_num,
        times_q,
        qpsrm.parameterization.func_dict,
    )
def quantity_solution_array_1(
    qmr: QuantityModelRun
    )->QuantitySolutionArray:
    return QuantitySolutionArray(qmr.solve())
