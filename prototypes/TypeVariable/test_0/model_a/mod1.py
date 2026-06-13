from typing import NamedTuple
from .Obs import Obs

class EP(NamedTuple):
    a: float
    b: float

def param2res(ep:EP) -> Obs:
    return Obs(x= ep.a+ep.b, y=ep.a-ep.b)

