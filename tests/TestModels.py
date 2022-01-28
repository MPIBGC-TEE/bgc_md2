from time import time as now
from unittest import TestCase, skip
from ComputabilityGraphs.rec_graph_helpers import fast_graph
from ComputabilityGraphs.helpers import all_mvars
from ComputabilityGraphs.TypeSet import TypeSet

import bgc_md2.helper as h 
from bgc_md2.resolve.mvars import *

class TestModels(TestCase):
    #@skip
    def test_all_computable_mvars_for_all_models(self):
        model_names = h.list_models(
                #explicit_exclude_models=frozenset({'Hilbert1991AnnBot', 'Thomas2014GeosciModelDev'})
                explicit_exclude_models=frozenset()
        )
        #print(model_names)
        for mn in model_names:
            # https://docs.python.org/3/library/unittest.html#distinguishing-test-iterations-using-subtests
            with self.subTest(mn=mn):
                mvs = h.CMTVS_from_model_name(mn)
                mvars = mvs.computable_mvar_types()
                list_str = "\n".join(["<li> " + str(var.__name__) + " </li>" for var in mvars])
                print(list_str)
                for var in mvars:
                    print("########################################")
                    print(str(var.__name__))
                    #print(mvs._get_single_value(var))
                    print(mvs._get_single_value_by_depgraph(var))
