import unittest
from sympy import ImmutableMatrix
from bgc_md2.resolve.mvars import (
    TimeSymbol,
    StateVariableTuple,
    InputTuple,
    CompartmentalMatrix,
)


class TestInputTuple(unittest.TestCase):
    def test_shape(self):
        I = InputTuple((1, 1))
        print(I, type(I))

class TestCompartmentalMatrix(unittest.TestCase):
    def test_creation(self):
        A = ImmutableMatrix(1,1,[1])
        B = ImmutableMatrix(1,1,[1])
        C = CompartmentalMatrix(A*B)
        print(C, type(C))
        self.assertEqual(type(C),CompartmentalMatrix)

class TestInputTuple(unittest.TestCase):
    def test_creation(self):
        I1 = ImmutableMatrix(2,1,[1,3])
        I2 = ImmutableMatrix(2,1,[1,3])
        I = InputTuple(I1+I2)
        print(I, type(I))
        self.assertEqual(type(I),InputTuple)
