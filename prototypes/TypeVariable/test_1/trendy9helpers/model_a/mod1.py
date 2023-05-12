from typing import NamedTuple
from pathlib import Path
from .Obs import Observations
from .. general_helpers import EP


class NT(NamedTuple):
    a: float
    b: float


class EstimatedParameters(NT, EP):
    pass


def param2res(ep:EstimatedParameters) -> Observations:
    return Observations(x=ep.a + ep.b, y=ep.a-ep.b)

def show_file_content():
    pp = Path(__file__).parent
    dp=pp.joinpath("data.txt")
    with dp.open("r") as f:
        print("\n".join(f.readlines()))
