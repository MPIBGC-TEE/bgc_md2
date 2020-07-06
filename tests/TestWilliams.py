import unittest
from bgc_md2.models.helpers import (
    provided_mvars
    ,computable_mvars
    ,path_dict_to_single_mvar
    ,get_single_mvar_value
)
from bgc_md2.resolve.mvars import (
    TimeSymbol
    ,StateVariableTuple
    ,CompartmentalMatrix
    ,InputPartitioningTuple
)

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

class TestWilliams(unittest.TestCase):
    def setUp(self):
        self.mn='Williams2005GCB'

    def test_provided_mvars(self):
        mvs=provided_mvars(self.mn)
        res=frozenset([
            CompartmentalMatrix
            ,TimeSymbol
            ,StateVariableTuple
            ,InputPartitioningTuple
        ])
        self.assertSetEqual(mvs,res)

    def test_computable_mvars(self):
        res=frozenset([
           # InFluxesBySymbol
           # ,OutFluxesBySymbol
           # ,InternalFluxesBySymbol
            TimeSymbol
            ,StateVariableTuple
            ,CompartmentalMatrix
            ,InputPartitioningTuple
        ])
        mvs=computable_mvars(self.mn)
        
        self.assertSetEqual(mvs,res)
        list_str="\n".join(["<li> "+str(var.__name__)+" </li>" for var in mvs])
        print(list_str)
        for var in mvs:
            print(get_single_mvar_value(var,self.mn))

    def test_get_mvar_value(self):
        # first get a provided value
        t=get_single_mvar_value(TimeSymbol,self.mn)
        self.assertEqual(t,TimeSymbol('t'))
        # now get a variable that is not provided directly but computable in one step

        #res=get_single_mvar_value(SmoothReservoirModel,self.mn)
        # now get a variable that is not provided directly but computable in two steps
        A=get_single_mvar_value(CompartmentalMatrix,self.mn)
        print(A)
        #self.assertTrue(False)

