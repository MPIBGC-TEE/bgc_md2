from typing import List
from sympy import Symbol, symbols
from sympy.matrices import ImmutableMatrix
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

##################################################
# Mvars...
class CompartmentalMatrix(ImmutableMatrix):
    # The constructor could perform some checks
    # but in general it is hard to check symbolicly if a
    # matrix is compartmental
    pass


class InputTuple(ImmutableMatrix):
    # we could have the constructor check for a one
    # dimensional and positive input
    pass


class StateTuple(ImmutableMatrix):
    # we could have the constructor check for a one
    # dimensional and purely symbolic input
    pass


class TimeSymbol(Symbol):
    pass


##################################################
# Computers ...
def reservoirModel(
    sv: StateTuple, t: TimeSymbol, A: CompartmentalMatrix, I: InputTuple
) -> SmoothReservoirModel:
    return SmoothReservoirModel.from_B_u(sv, t, A, I)


# ...many more functions comprising our domain knowledge


##################################################
# start inspections
from inspect import signature

# This is how we  build the computability graph from
# the analysis of the computers and the user code
# first we inspect the implemented computers

# The following would have to be done for every function in a module reserved for the computers.
sig = signature(reservoirModel)
[val.annotation.__name__ for key, val in sig.parameters.items()]
# get the types of the parameters
# ['StateTuple', 'TimeSymbol', 'CompartmentalMatrix', 'InputTuple']
# and the type of the return value
sig.return_annotation.__name__
# 'SmoothReservoirModel'


def input_mvars(computer):
    params = signature(computer).parameters.values()
    return {param.annotation for param in params}


def output_mvar(computer):
    return signature(computer).return_annotation


# allComputers is the union of all functions in the computers
# module
# The set of all Mvartypes is the union of the signatures of all
# computers
# we can compute the computability graph from those two sets


# to get the set of user defined variables import the user code
# normaly from a model specific file like:
# from .models.markus1.source import special_vars
# here for demonstration just:

##################################################
# user Code (normaly stored in a model specific directory (module)
# has to implement a variable with a specified name
# here "special_vars"

a, b, t = symbols("a,b,t")
sv = StateTuple([a, b])
cm = CompartmentalMatrix([[1, a], [a, a]])
i = InputTuple([3, 2])
special_vars = [sv, cm, i]

user_defined_vars = special_vars
# the set of Mvars for a model is then gives as
[type(v).__name__ for v in user_defined_vars]


# testcode
# should be constructed from a chosen path through the graph
a, b, t = symbols("a,b,t")
sv = StateTuple([a, b])
cm = CompartmentalMatrix([[1, a], [a, a]])
i = InputTuple([3, 2])
mod = reservoirModel(sv, t, cm, i)
