import unittest
from sympy import symbols
from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
    CompartmentalMatrix,
)
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

class TestMVarSet(unittest.TestCase):
    def setUp(self):
        I_vl, I_vw = symbols("I_vl I_vw")
        vl, vw = symbols("vl vw")
        k_vl, k_vw = symbols("k_vl k_vw")
        self.provided_values=frozenset(
            [
                InFluxesBySymbol({vl: I_vl, vw: I_vw}),
                OutFluxesBySymbol({vl: k_vl * vl, vw: k_vw * vw}),
                InternalFluxesBySymbol({(vl, vw): k_vl * vl, (vw, vl): k_vw * vw}),
                TimeSymbol("t"),
                StateVariableTuple((vl, vw))
            ]
        )
        self.mvs=MVarSet(self.provided_values)
    
    def test_from_model_name(self):
        # alternative constructor
        mvs= MVarSet.from_model_name("testVectorFree")
        res = frozenset(
            [
                InFluxesBySymbol,
                OutFluxesBySymbol,
                InternalFluxesBySymbol,
                TimeSymbol,
                StateVariableTuple,
            ]
        )
        self.assertSetEqual(mvs.provided_mvar_types, res)

    def test_provided_mvar_values(self):
        mvs = self.mvs
        pmvs = mvs.provided_mvar_values
        self.assertSetEqual(pmvs, self.provided_values)
        
    def test_provided_mvar_types(self):
        mvs = self.mvs
        pmvts = mvs.provided_mvar_types
        self.assertSetEqual(
            pmvts,
            frozenset(type(v) for v in self.provided_values)
        )
        

    # @unittest.skip
    def test_computable_mvar_types(self):
        ''' This test also depends on the content of bgc_md2.resolve.computers
            since more variables become computable if we add computers...
        '''
        mvs = self.mvs
        res = frozenset(
            [
                InFluxesBySymbol,
                OutFluxesBySymbol,
                InternalFluxesBySymbol,
                TimeSymbol,
                StateVariableTuple,
                SmoothReservoirModel,
                CompartmentalMatrix,
            ]
        )
        cmvs = mvs.computable_mvar_types()
        print("###############################################")
        print(cmvs)
        self.assertSetEqual(cmvs, res)

    
    def test_paths_to_single_mvar(self):
        mvs = self.mvs
        pd = mvs.path_dict_to_single_mvar(SmoothReservoirModel)
        startNode = frozenset(
            [
                InFluxesBySymbol,
                OutFluxesBySymbol,
                InternalFluxesBySymbol,
                TimeSymbol,
                StateVariableTuple,
            ]
        )
        paths = pd[startNode]
        res = paths[0]
        ref = [startNode, frozenset({SmoothReservoirModel})]
        self.assertEqual(res, ref)
    
    def test_get_mvar_value(self):
        mvs = self.mvs
        res = mvs._get_single_mvar_value(TimeSymbol)
        self.assertEqual(res, TimeSymbol("t"))
        
        ## now get a variable that is not provided directly but computable in one step
        res = mvs._get_single_mvar_value(SmoothReservoirModel)
        ## now get a variable that is not provided directly but computable in two steps
        res = mvs._get_single_mvar_value(CompartmentalMatrix)
        #print(res)
        ## self.assertTrue(False)

