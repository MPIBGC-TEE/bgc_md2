import unittest
from bgc_md2.resolve.mvars import (
    TimeSymbol,
    StateVariableTuple,
    InputTuple,
    CompartmentalMatrix,
)


class TestInputTuple:
    def test_shape(self):
        I = InputTuple((1, 1))
        print(I, type(I))
