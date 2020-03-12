# The purpose of this Schema is to work backwards from the minimal requirement that the 
# an Instance of CompartmentalModel can be created.
# So a model consists at minimum constructor call.
# and possibly some variable definitions to populate the namespace in which the constructor is called. 


import unittest
from testinfrastructure.helpers import pe
from testinfrastructure.InDirTest import InDirTest

from sympy import Symbol,Number
from typing import List
#from CompartmentalSystems import smooth_reservoir_model 
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol
    ,OutFluxesBySymbol
    ,InternalFluxesBySymbol
    ,TimeSymbol
    ,StateVariableTuple
    ,CompartmentalMatrix
)
import bgc_md2.resolve.computers as bgc_computers
from bgc_md2.models.helpers import (
    provided_mvars
    ,computable_mvars
    ,path_dict_to_single_mvar
    ,get_single_mvar_value
)
import sys
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
sys.setrecursionlimit(2500)

class TestModelFiles(unittest.TestCase):
    # The aim is a proof of concept implementation for the retrieval of the information that is neccessary to build the 
    # compartmental Matrix 
    # Here we execute a python script in a special sandbox environment
    
    def test_provided_mvars(self):
        mvts=provided_mvars('testVectorFree')
        res=frozenset([
            InFluxesBySymbol
            ,OutFluxesBySymbol
            ,InternalFluxesBySymbol
            ,TimeSymbol
            ,StateVariableTuple
        ])
        self.assertSetEqual(mvts,res)
        
    def test_computable_mvars(self):
        res=frozenset([
            InFluxesBySymbol
            ,OutFluxesBySymbol
            ,InternalFluxesBySymbol
            ,TimeSymbol
            ,StateVariableTuple
            ,SmoothReservoirModel
            ,CompartmentalMatrix
        ])
        mvts=computable_mvars('testVectorFree')
        self.assertSetEqual(mvts,res)

    def test_paths_to_single_mvar(self):
        pd=path_dict_to_single_mvar(SmoothReservoirModel,'testVectorFree')
        startNode=frozenset([
            InFluxesBySymbol
            ,OutFluxesBySymbol
            ,InternalFluxesBySymbol
            ,TimeSymbol
            ,StateVariableTuple
        ])
        paths=pd[startNode]
        res=paths[0]
        ref=[startNode,frozenset({SmoothReservoirModel})]
        self.assertEqual(res,ref)

    def test_get_mvar_value(self):
        # first get a provided value
        res=get_single_mvar_value(TimeSymbol,'testVectorFree') 
        self.assertEqual(res,TimeSymbol('t'))
        # now get a variable that is not provided directly but computable in one step 
        
        res=get_single_mvar_value(SmoothReservoirModel,'testVectorFree') 
        # now get a variable that is not provided directly but computable in two steps
        res=get_single_mvar_value(CompartmentalMatrix,'testVectorFree') 
        print(res)
        #self.assertTrue(False)

    @unittest.skip
    def test_many(self):
        #https://docs.python.org/3/library/unittest.html#distinguishing-test-iterations-using-subtests
        #for md in ['miniCable','Sujan']:
        for md in ['Sujan']:
            with self.subTest(md=md):
                #assertions would be possible too
                srm=get3(var_name="smooth_reservoir_model",allMvars=myMvars,allComputers=myComputers,model_id=md)
                


