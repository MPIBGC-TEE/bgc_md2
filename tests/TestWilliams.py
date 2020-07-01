
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
)

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

class TestWilliams(unittest.TestCase):
    def setUp(self):
        self.mn='Williams2005GCB'

    def test_provided_mvars(self):
        mvts=provided_mvars(self.mn)
        res=frozenset([
            CompartmentalMatrix
            ,TimeSymbol
            ,StateVariableTuple
        ])
        self.assertSetEqual(mvts,res)

    def test_computable_mvars(self):
        res=frozenset([
           # InFluxesBySymbol
           # ,OutFluxesBySymbol
           # ,InternalFluxesBySymbol
            TimeSymbol
            ,StateVariableTuple
            ,CompartmentalMatrix
        ])
        mvts=computable_mvars(self.mn)
        self.assertSetEqual(mvts,res)

    def test_get_mvar_value(self):
        # first get a provided value
        res=get_single_mvar_value(TimeSymbol,self.mn)
        self.assertEqual(res,TimeSymbol('t'))
        # now get a variable that is not provided directly but computable in one step

        #res=get_single_mvar_value(SmoothReservoirModel,self.mn)
        # now get a variable that is not provided directly but computable in two steps
        res=get_single_mvar_value(CompartmentalMatrix,self.mn)
        print(res)
        #self.assertTrue(False)

