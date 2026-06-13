from typing import NamedTuple
from .Obs import Observations
from ..general_helpers import EP

class NT(NamedTuple):
    a: float
    b: float
    c: float
    d: float

class EstimatedParameters(NT, EP):
    pass


def param2res(ep: EstimatedParameters) -> Observations:
    return Observations(x= ep.a+ep.c+ep.d, y=ep.a-ep.c-ep.d)
