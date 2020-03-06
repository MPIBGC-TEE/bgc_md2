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
)
import bgc_md2.resolve.computers as bgc_computers
from bgc_md2.models.helpers import (
        available_mvars
        ,computable_mvars
        )
import sys
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
sys.setrecursionlimit(2500)

class TestModelFiles(unittest.TestCase):
    # The aim is a proof of concept implementation for the retrieval of the information that is neccessary to build the 
    # compartmental Matrix 
    # Here we execute a python script in a special sandbox environment
    
    def test_available_mvars(self):
        mvts=available_mvars('testVectorFree')
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
        ])
        mvts=computable_mvars('testVectorFree')
        self.assertSetEqual(mvts,res)

    def test_get_mvar(self):
        self.assertTrue(False)

    @unittest.skip
    def test_many(self):
        #https://docs.python.org/3/library/unittest.html#distinguishing-test-iterations-using-subtests
        #for md in ['miniCable','Sujan']:
        for md in ['Sujan']:
            with self.subTest(md=md):
                #assertions would be possible too
                srm=get3(var_name="smooth_reservoir_model",allMvars=myMvars,allComputers=myComputers,model_id=md)
                


