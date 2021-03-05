"""
Module for defining variable types for the description of compartmental models and model runs
"""
from typing import Tuple
from frozendict import frozendict
import numpy as np
from sympy import (
    Symbol,
    # symbols,
    Function,
    prod,
    sin,
    cos,
    pi,
    lambdify,
    simplify,
    factor,
    ImmutableMatrix,
    Expr,
)
from sympy.core import Basic, Dict, Integer, Tuple
from sympy.core.cache import cacheit
from sympy.core.sympify import converter as sympify_converter, _sympify
from sympy.matrices.dense import DenseMatrix
from sympy.matrices.expressions import MatrixExpr
from sympy.matrices.matrices import MatrixBase

from sympy.physics.units import Quantity
from sympy.physics.units.systems import SI
from sympy.physics.units import time
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun


from bgc_md2.described_quantities import to_number
 
from testinfrastructure.helpers import pp, pe

# fixme: mm 03-12-2020
# At the moment the classes are just defined to provide
# the vocabulary for the computer signatures and model
# descriptions.
# Ultimately the classes should have sensible constructors
# that check their arguments thoroughly to expose
# inadequate imput early in the process
# part of this should be a dimensional analysis
# which requires that the model authors provide dimensions
# for there variables


class InFluxesBySymbol(frozendict):
    pass


class OutFluxesBySymbol(frozendict):
    pass


class InternalFluxesBySymbol(frozendict):
    pass


class CarbonInFluxesBySymbol(frozendict):
    pass


class CarbonOutFluxesBySymbol(frozendict):
    pass


class VegetationCarbonInFluxesBySymbol(frozendict):
    pass


class VegetationCarbonOutFluxesBySymbol(frozendict):
    pass


class VegetationCarbonInternalFluxesBySymbol(frozendict):
    pass


class CarbonInternalFluxesBySymbol(frozendict):
    pass


class NitrogenInFluxesBySymbol(frozendict):
    pass


class NitrogenOutFluxesBySymbol(frozendict):
    pass


class NitrogenInternalFluxesBySymbol(frozendict):
    pass

class TimeSymbol(Symbol):
    # should become a quantity with dimension time
    pass

class MatrixLike(ImmutableMatrix):
    @classmethod
    def _new(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], cls):
            return args[0]
        if kwargs.get('copy', True) is False:
            if len(args) != 3:
                raise TypeError("'copy=False' requires a matrix be initialized as rows,cols,[list]")
            rows, cols, flat_list = args
        else:
            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
            flat_list = list(flat_list) # create a shallow copy

        obj = Basic.__new__(cls,
            Integer(rows),
            Integer(cols),
            Tuple(*flat_list))
        obj._rows = rows
        obj._cols = cols
        obj._mat = flat_list
        return obj

class ColumnVectorLike(MatrixLike):
    # fixme add some check that there is only one column...
    pass

class StateVariableTuple(ColumnVectorLike):
    # fixme mm 02/14/2021
    # Check that we only allow column vectors
    """Defines the ordering of the state variables.
    This defines the coordinate system and thereby
    the interpretation of all tuple or matrix related variables.
    Although two models with identical fluxes between there pools are identical on a physical
    level the coordinate based variables like the InputTuple or the CompartmentalMatrix depend 
    on the ordering and have to be consistent."""
    pass

class VegetationCarbonStateVariableTuple(ColumnVectorLike):
    # fixme mm 02/14/2021
    # Check that we only allow column vectors
    """Defines the ordering of the vegetation carbon state variables.  The
    order defines the order of the axes in a vegetation related coordinate
    system and thereby the interpretation of all tuple or matrix related carbon
    vegetation variables.  Although two models with identical vegetation carbon
    fluxes between their pools are identical on a physical level the coordinate
    based variables like the VegetationInputTuple or the
    VegetationInputTupleCompartmentalMatrix depend on the ordering and have to
    be consistent."""
    pass

class CarbonStateVariableTuple(ColumnVectorLike):
    # fixme mm 02/14/2021
    # Check that we only allow column vectors
    """Defines the ordering of the carbon state variables.  The order defines
    the order of the axes in a carbon related coordinate system and thereby the
    interpretation of all tuple or matrix related carbon variables.  Although
    two models with identical carbon fluxes between their pools are identical
    on a physical level the coordinate based variables like the
    VegetationInputTuple or the VegetationInputTupleCompartmentalMatrix depend
    on the ordering and have to be consistent."""
    pass


class NitrogenStateVariableTuple(ColumnVectorLike):
    # fixme mm 02/14/2021
    # Check that we only allow column vectors
    """Defines the ordering of the nitrogen state variables.  The order defines
    the order of the axes in a nitrogen related coordinate system and thereby the
    interpretation of all tuple or matrix related nitrogen variables.  Although
    two models with identical nitrogen fluxes between their pools are identical
    on a physical level the coordinate based variables like the
    VegetationInputTuple or the VegetationInputTupleCompartmentalMatrix depend
    on the ordering and have to be consistent."""
    pass


class InputTuple(ColumnVectorLike):
    # fixme mm 02/14/2021
    # Check that we only allow column vectors
    """The coordinates of the input flux vector
    The interpretation requires the statevector since the nth entry is interpreted as input to 
    the pool denoted by the nth state variable"""
    pass
 

class CarbonInputTuple(ColumnVectorLike):
    # fixme mm 02/14/2021
    # Check that we only allow column vectors
    """The coordinates of the carbon input flux vector
    The interpretation requires the carbon statevector since the nth entry is interpreted as input to 
    the pool denoted by the nth carbon state variable"""
    pass

class NitrogenInputTuple(ColumnVectorLike):
    # fixme mm 02/14/2021
    # Check that we only allow column vectors
    """The coordinates of the nitrogen input flux vector
    The interpretation requires the nitrogen statevector since the nth entry is interpreted as input to 
    the pool denoted by the nth nitrogen state variable"""
    pass

# vegetation specific variables
class VegetationCarbonInputTuple(ColumnVectorLike):
    # fixme mm 02/14/2021
    # Check that we only allow column vectors
    """The coordinates of the vegetation part of the the Carbon input flux vector
    Note that this will affect 
    """
    pass


class VegetationCarbonInputScalar(Expr):
    pass


#class CompartmentalMatrix(ImmutableMatrix, MatrixExpr):
class CompartmentalMatrix(MatrixLike):
    """Create a CompartmentalMatrix (which is immutable).

    Examples
    ========

    >>> from sympy import eye
    >>> from bgc_md2.resolve.mvars import CompartmentalMatrix
    >>> CompartmentalMatrix(eye(3))
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])
    >>> _[0, 0] = 42
    Traceback (most recent call last):
    ...
    TypeError: Cannot set values of CompartmentalMatrix
    """
    pass


class CarbonCompartmentalMatrix(MatrixLike):
    """Create a CompartmentalMatrix (which is immutable).

    Examples
    ========

    >>> from sympy import eye
    >>> from bgc_md2.resolve.mvars import CarbonCompartmentalMatrix
    >>> CarbonCompartmentalMatrix(eye(3))
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])
    >>> _[0, 0] = 42
    Traceback (most recent call last):
    ...
    TypeError: Cannot set values of CompartmentalMatrix
    """
    pass

class NitrogenCompartmentalMatrix(MatrixLike):
    """Create a CompartmentalMatrix (which is immutable).

    Examples
    ========

    >>> from sympy import eye
    >>> from bgc_md2.resolve.mvars import NitrogenCompartmentalMatrix
    >>> NitrogenCompartmentalMatrix(eye(3))
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])
    >>> _[0, 0] = 42
    Traceback (most recent call last):
    ...
    TypeError: Cannot set values of CompartmentalMatrix
    """
    pass

class VegetationCarbonInputPartitioningTuple(ColumnVectorLike):
    # fixme mm 02/14/2021
    # Check that we only allow column vectors
    """The partitioning components vector
    The order will affect the VegetationCarbonInput"""
    pass


class VegetationCarbonCompartmentalMatrix(MatrixLike):  # cycling matrix
    """Create a CompartmentalMatrix (which is immutable) for the vegetation part only.
    Note that for the same physical model this matrix can take different shapes, depending
    on the order of the vegetation state variables (as defined via the VegetationStateVariableTuple)
    
    It is more likely that you will not create this matrix itself since it can easily be computed
    from the VegetationCarbonStateVariables and the CompartmentalMatrix for the whole system. (or the carbon cycle)

    Examples
    ========

    >>> from sympy import eye
    >>> from bgc_md2.resolve.mvars import VegetationCarbonCompartmentalMatrix
    >>> VegetationCarbonCompartmentalMatrix(eye(3))
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])
    >>> _[0, 0] = 42
    Traceback (most recent call last):
    ...
    TypeError: Cannot set values of ImmutableDenseMatrix
    """
    pass


class NumericStartValueDict(frozendict):
    pass


class StateVarUnitTuple(tuple):
    pass


# extending ndarray is special
# https://numpy.org/doc/stable/user/basics.subclassing.html
class NumericSimulationTimes(np.ndarray):
    def __new__(cls, input_array):
        # Input array is an already formed ndarray instance
        # We cast to be our class type
        obj = np.asarray(input_array).view(cls)
        obj.flags.writeable = False
        return obj

    def __hash__(self):
        return hash(tuple(self))


class NumericStartValueArray(np.ndarray):
    def __new__(cls, input_array):
        # Input array is an already formed ndarray instance
        # We cast to be our class type
        obj = np.asarray(input_array).view(cls)
        obj.flags.writeable = False
        return obj

    def __hash__(self):
        return hash(tuple(self))


class NumericParameterization:
    # Note:
    # A parameterization implicitly refers to a unique specific
    # symbolic model:
    #
    # 1.)   Obviously the keys of the par_dict and func_dict can
    #       only be substituted in a model with these symbols present
    # 2.)   Even if a symbolic model had the same free symbols
    #       these symbols aquire meaning from their use in the symbolic
    #       expressions.
    #
    #
    # An instance would naturally contain a referece to the model (self.model=).
    # This reference is omitted on purpose since the model might be given
    # only implicitly by variables defined in a model describing source.py
    def __init__(self, par_dict, func_dict):
        self.par_dict = frozendict(par_dict)
        self.func_dict = frozendict(func_dict)


class NumericParameterizedSmoothReservoirModel:
    def __init__(self, srm, parameterization):
        self.srm = srm
        self.parameterization = parameterization


class QuantityStartValueArray(np.ndarray):
    def __new__(cls, input_array):
        # Input array is an already formed ndarray instance
        # We cast to be our class type
        obj = np.asarray(input_array).view(cls)
        obj.flags.writeable = False
        return obj

    def __hash__(self):
        return hash(tuple(self))


class QuantityStartValueDict(frozendict):
    pass


class QuantitySimulationTimes(np.ndarray):
    def __new__(cls, input_array):
        # Input array is an already formed ndarray instance
        # We cast to be our class type
        obj = np.asarray(input_array).view(cls)
        obj.flags.writeable = False
        return obj

    def __hash__(self):
        return hash(tuple(self))


class NumericSolutionArray(np.ndarray):
    # just a wrapper class to distinguish the array as a solution
    def __new__(cls, input_array):
        # Input array is an already formed ndarray instance
        # We cast to be our class type
        obj = np.asarray(input_array).view(cls)
        obj.flags.writeable = False
        return obj

    def __hash__(self):
        return hash(tuple(self))

class QuantitySolutionArray(np.ndarray):
    # just a wrapper class to distinguish the array as a solution (with units) via the type
    def __new__(cls, input_array):
        # Input array is an already formed ndarray instance
        # We cast to be our class type
        obj = np.asarray(input_array).view(cls)
        obj.flags.writeable = False
        return obj

    def __hash__(self):
        return hash(tuple(self))



class QuantityParameterization(NumericParameterization):
    # If a numeric parameterization has some physical meanimg.
    # the purely numerical parameter values have implicit units.
    # Also purely numerical functions have to be interpreted in the
    # context of the units of their arguments and return value.
    #
    # Instead of seperately attaching units to
    # ALL parameters, function arguments and  return values and expressions
    # and letting sympy take care of the conversions,
    # we require the user to choose the units of
    # - start values and
    # - times
    # in a simulation using the paramterization.
    # Thus the user has to provide a numerical parameterization for the
    # WHOLE compartmental model to be consistent with these units.
    # Possibly parameters have to be converted manually by the user to work
    # for the units specified.
    # This  keeps the numerical computations several orders of magnitude faster
    # than fully automatic unit derivation in sympy and can be seen as a kind
    # of cached unit computation.

    def __init__(self, par_dict, func_dict, state_var_units, time_unit):
        super().__init__(par_dict, func_dict)
        self.state_var_units = state_var_units
        self.time_unit = time_unit


        

class QuantityParameterizedSmoothReservoirModel:
    def __init__(self, srm, parameterization):
        self.srm = srm
        self.parameterization = parameterization


class QuantityModelRun:
    def __init__(
        self, qpm, start_values_quant, times_quant,
    ):
        self.qpm = qpm
        p = qpm.parameterization
        if len(start_values_quant) != len(p.state_var_units):
            raise Exception('size inconsistency')
        times_num = np.array([to_number(tv, p.time_unit) for tv in times_quant])
        start_values_num = np.array(
            [to_number(sv, p.state_var_units[i]) for i, sv in enumerate(start_values_quant)]
        )
        self.smr = SmoothModelRun(
            qpm.srm, p.par_dict, start_values_num, times_num, p.func_dict
        )

    def solve(self):
        # compute the solution with respect to the state_var_units given in
        # the model parameterization
        sol_num = self.smr.solve()

        # the result is correct since it comes with unit
        # and can be converted in any other unit.
        sol_quant = sol_num * self.qpm.parameterization.state_var_units

        return sol_quant
