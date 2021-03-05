import unittest
from sympy import Matrix, ImmutableMatrix, symbols, Symbol
from bgc_md2.resolve.mvars import (
    TimeSymbol,
    StateVariableTuple,
    InputTuple,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputTuple,
    CompartmentalMatrix,
    # CompartmentalMatrixStructure,
)


class TestStateVariableTuple(unittest.TestCase):
    def test_creation(self):
        a, b, c = symbols('a b c')
        X = StateVariableTuple((a, b))
        self.assertEqual(type(X), StateVariableTuple)
        print(X, type(X))
        Y = X.subs({b: c})
        print(Y, type(Y))
        self.assertEqual(type(Y), StateVariableTuple)

    def test_immutability(self):
        a, b, c = symbols('a b c')
        T = StateVariableTuple((a, b, c))
        with self.assertRaises(TypeError):
            T[1] = 1


class TestVegetationCarbonInputScalar(unittest.TestCase):
    def test_creation(self):
        a, b, c = symbols('a b c')

        u = VegetationCarbonInputScalar(a)
        self.assertEqual(type(u), VegetationCarbonInputScalar)

        v = u.subs({a: b})
        self.assertEqual(type(v), VegetationCarbonInputScalar)

        u = VegetationCarbonInputScalar(1)
        self.assertEqual(type(u), VegetationCarbonInputScalar)


class TestInputTuple(unittest.TestCase):
    def test_creation(self):
        a, b, c = symbols('a b c')
        In = InputTuple((a, b))
        self.assertEqual(type(In), InputTuple)

        J = In.subs({b: c})
        self.assertEqual(type(J), InputTuple)

        In = InputTuple((1, 1))
        self.assertEqual(type(In), InputTuple)

        M = Matrix((1, 1))
        In = InputTuple(M)
        self.assertEqual(type(In), InputTuple)

        I1 = ImmutableMatrix(2, 1, [1, 3])
        I2 = ImmutableMatrix(2, 1, [1, 3])
        In = InputTuple(I1+I2)
        self.assertEqual(type(In), InputTuple)
        self.assertEqual(In[0, 0], 2)
        self.assertEqual(In[1, 0], 6)

    def test_immutability(self):
        T = InputTuple([1, 0, 1, 0])
        with self.assertRaises(TypeError):
            T[1] = 1


class TestVegetationCarbonInputTuple(unittest.TestCase):
    def test_creation(self):
        a, b, c = symbols('a b c')
        In = VegetationCarbonInputTuple((a, b))
        self.assertEqual(type(In), VegetationCarbonInputTuple)

        J = In.subs({b: c})
        self.assertEqual(type(J), VegetationCarbonInputTuple)

        In = VegetationCarbonInputTuple((1, 1))
        self.assertEqual(type(In), VegetationCarbonInputTuple)

        M = Matrix((1, 1))
        In = VegetationCarbonInputTuple(M)
        self.assertEqual(type(In), VegetationCarbonInputTuple)

        I1 = ImmutableMatrix(2, 1, [1, 3])
        I2 = ImmutableMatrix(2, 1, [1, 3])
        In = VegetationCarbonInputTuple(I1+I2)
        self.assertEqual(type(In), VegetationCarbonInputTuple)
        self.assertEqual(In[0, 0], 2)
        self.assertEqual(In[1, 0], 6)

    def test_immutability(self):
        T = VegetationCarbonInputTuple([1, 0, 1, 0])
        with self.assertRaises(TypeError):
            T[1] = 1


class TestVegetationCarbonInputPartitioningTuple(unittest.TestCase):
    def test_creation(self):
        a, b, c = symbols('a b c')
        In = VegetationCarbonInputPartitioningTuple((a, b))
        self.assertEqual(type(In), VegetationCarbonInputPartitioningTuple)

        J = In.subs({b: c})
        self.assertEqual(type(J), VegetationCarbonInputPartitioningTuple)

        In = VegetationCarbonInputPartitioningTuple((1, 1))
        self.assertEqual(type(In), VegetationCarbonInputPartitioningTuple)

        M = Matrix((1, 1))
        In = VegetationCarbonInputPartitioningTuple(M)
        self.assertEqual(type(In), VegetationCarbonInputPartitioningTuple)

        I1 = ImmutableMatrix(2, 1, [1, 3])
        I2 = ImmutableMatrix(2, 1, [1, 3])
        In = VegetationCarbonInputPartitioningTuple(I1+I2)
        self.assertEqual(type(In), VegetationCarbonInputPartitioningTuple)
        self.assertEqual(In[0, 0], 2)
        self.assertEqual(In[1, 0], 6)

    def test_immutability(self):
        T = VegetationCarbonInputPartitioningTuple([1, 0, 1, 0])
        with self.assertRaises(TypeError):
            T[1] = 1


class TestCompartmentalMatrix(unittest.TestCase):
    # fixme mm 2/11/2021:
    # The tests up to now only check technical properties Although it is
    # difficult to decide if a particular (symbolic) matrix is compartmental we
    # at least can sometimes infer that it is not.  We should warn the user
    # about an attempt to build a noncompartmental instance of
    # CompartmentalMatrix and possibly even refuse to do so.  The tests should
    # include such cases but at the time do not.
    def test_creation_from_list(self):
        C = CompartmentalMatrix(2, 2, [1, 2, 3, 4])
        print(C, type(C))
        self.assertEqual(type(C), CompartmentalMatrix)

    def test_creation_from_index_function(self):
        C = CompartmentalMatrix(2, 2, lambda i, j: i + j)
        print(C, type(C))
        self.assertEqual(type(C), CompartmentalMatrix)

    def test_creation_from_matrix(self):
        A = ImmutableMatrix(1, 1, [1])
        B = ImmutableMatrix(1, 1, [1])
        D = A*B
        C = CompartmentalMatrix(D)
        print(C, type(C))
        self.assertEqual(type(C), CompartmentalMatrix)

    def test_creation_from_expression(self):
        A = ImmutableMatrix(1, 1, [1])
        B = ImmutableMatrix(1, 1, [1])
        C = CompartmentalMatrix(A*B)
        print(C, type(C))
        self.assertEqual(type(C), CompartmentalMatrix)

    def test_immutability(self):
        A = CompartmentalMatrix(2, 2, [1, 0, 1, 0])
        with self.assertRaises(TypeError):
            A[1, 1] = 1

