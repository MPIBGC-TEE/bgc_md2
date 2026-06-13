from typing import NamedTuple
from .Obs import Obs

class EP(NamedTuple):
    a: float
    c: float
    d: float


def param2res(ep:EP) -> Obs:
    return Obs(x= ep.a+ep.c+ep.d, y=ep.a-ep.c-ep.d)
