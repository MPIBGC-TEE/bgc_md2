from frozendict import frozendict
from sympy import Symbol,Quantity,ImmutableMatrix
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

class TimeSymbol(Symbol):
    # should become a quantity with dimension time
    pass

class NPP(Quantity):
    # should become a quantity with dimension time
    pass

class StateVariableTuple(tuple):
    pass

class InputTuple(tuple):
    pass

class VegetationCarbonInputPartitioningTuple(tuple):
    pass

class CompartmentalMatrix(ImmutableMatrix):
    pass

