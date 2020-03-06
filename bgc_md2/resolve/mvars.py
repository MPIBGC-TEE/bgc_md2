from frozendict import frozendict
from sympy import Symbol

# The classes should have sensible constructors
# that check their arguments thuroughly to expose
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

