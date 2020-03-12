from frozendict import frozendict
from sympy import Symbol,ImmutableMatrix
# fixme: mm 03-12-2020
# At the moment the classes are just defined to provide
# the vocabulary for the computer signatures and model 
# descriptions.
# Ultimately the classes should have sensible constructors
# that check their arguments thoroughly to expose
# inadequate imput early in the process 

class InFluxesBySymbol(frozendict):
    pass

class OutFluxesBySymbol(frozendict):
    pass

class InternalFluxesBySymbol(frozendict):
    pass

class TimeSymbol(Symbol):
    pass

class StateVariableTuple(tuple):
    pass

class CompartmentalMatrix(ImmutableMatrix):
    pass
