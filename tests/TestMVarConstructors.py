import unittest
from sympy import Matrix,ImmutableMatrix
from bgc_md2.resolve.mvars import (
    TimeSymbol,
    StateVariableTuple,
    InputTuple,
    CompartmentalMatrix,
    CompartmentalMatrixStructure,
)


class TestInputTuple(unittest.TestCase):
    def test_shape(self):
        I = InputTuple((1, 1))
        print(I, type(I))


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
        A = CompartmentalMatrixStructure(2, 2, [1, 0, 1, 0])
        with self.assertRaises(TypeError):
            A[1, 1] = 1


class TestCompartmentalMatrixStructure(unittest.TestCase):
    # Fixme mm 2/11/2021
    # This class is merely defined as a RESULT of a cumputation
    # and the only argument is likely to be an instance of
    # CompartmentalMatrix.
    # Nearly nothing can be computed FROM it.
    # In the computability graph it is a dead end.
    # In this light it seems a bit odd to provide
    # the whole plethora of sympy ways to construct such a thing.
    # This might be a recurring theme for other computations 
    # that are merely projections, since they loose information.

    def test_creation_from_list(self):
        C = CompartmentalMatrixStructure(2, 2, [1, 0, 1, 0])
        print(C, type(C))
        self.assertEqual(type(C), CompartmentalMatrixStructure)
        
        with self.assertRaises(TypeError):
            # Adjacency Matrices contain only zeros and ones
            C = CompartmentalMatrixStructure(2, 2, [2, 0, 1, 0])

    def test_creation_from_index_function(self):
        C = CompartmentalMatrixStructure(2, 2, lambda i, j: i * j)
        print(C, type(C))
        self.assertEqual(type(C), CompartmentalMatrixStructure)
        with self.assertRaises(TypeError):
            # Adjacency Matrices contain only zeros and ones
            C = CompartmentalMatrixStructure(2, 2, lambda i, j: i + j)

    def test_creation_from_matrix(self):
        A = Matrix(1, 1, [1])
        B = Matrix(1, 1, [1])
        D = A*B
        C = CompartmentalMatrixStructure(D)
        print(C, type(C))
        self.assertEqual(type(C), CompartmentalMatrixStructure)
        with self.assertRaises(TypeError):
            # Adjacency Matrices contain only zeros and ones
            A = Matrix(1, 1, [1])
            B = Matrix(1, 1, [2])
            D = A*B
            C = CompartmentalMatrixStructure(D)

    def test_creation_from_expression(self):
        A = Matrix(1, 1, [1])
        B = Matrix(1, 1, [1])
        C = CompartmentalMatrixStructure(A*B)
        print(C, type(C))
        self.assertEqual(type(C), CompartmentalMatrixStructure)

    def test_immutability(self):
        A = CompartmentalMatrixStructure(2, 2, [1, 0, 1, 0])
        with self.assertRaises(TypeError):
            A[1, 1] = 1


class TestInputTuple(unittest.TestCase):
    def test_creation(self):
        I1 = ImmutableMatrix(2,1,[1,3])
        I2 = ImmutableMatrix(2,1,[1,3])
        I = InputTuple(I1+I2)
        print(I, type(I))
        self.assertEqual(type(I),InputTuple)
