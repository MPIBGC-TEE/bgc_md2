from typing import TypeVar, Callable, NamedTuple, Generic
class Obs(NamedTuple):
    x: float
    y: float

T=TypeVar("T")
class EP(Generic[T]):
    def __init__(self,t:T):
        self.t=t
