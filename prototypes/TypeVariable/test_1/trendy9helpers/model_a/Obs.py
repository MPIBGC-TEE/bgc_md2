from typing import NamedTuple
from ..general_helpers import Obs 

class NT(NamedTuple):
    x: float
    y: float

class Observations(NT, Obs):
    pass
