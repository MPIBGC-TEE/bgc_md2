import unittest
import matplotlib.pyplot as plt

import ipywidgets as widgets
import bgc_md2.helper as h 

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun

from testinfrastructure.InDirTest import InDirTest

from bgc_md2.resolve.helpers import (
    bgc_md2_computers,
    bgc_md2_computer_aliases,
    bgc_md2_mvar_aliases,
)
from bgc_md2.resolve.graph_helpers import sparse_powerset_graph
from bgc_md2.resolve.non_graph_helpers import (
    all_mvars
)
from bgc_md2.resolve.mvars import (
    TimeSymbol,
    StateVariableTuple,
    CompartmentalMatrix,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonInputTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonStateVariableTuple,
    InputTuple,
    NumericParameterization,
    NumericStartValueDict,
    NumericSimulationTimes,
    NumericParameterizedSmoothReservoirModel,
)
from bgc_md2.models.BibInfo import BibInfo 

from bgc_md2.resolve.computers import numeric_parameterized_smooth_reservoir_model_1
from bgc_md2.resolve.graph_plotting import (
    AGraphComputerSetMultiDiGraph,
    AGraphComputerMultiDiGraph,
    draw_update_sequence,
    draw_ComputerSetMultiDiGraph_matplotlib,
    # ,draw_Graph_with_computers_svg
)
from bgc_md2.resolve.MVarSet import MVarSet
from testinfrastructure.helpers import pp

class TestModels(InDirTest):
    def test_list_target_models(self):
        #li = h.list_target_models(
        #    explicit_exclude_models=frozenset({'CARDAMOM'})
        #)
        #print(li)

        # li = h.list_target_models(
        #     target_classes=frozenset(
        #        {CompartmentalMatrix, StateVariableTuple}
        #     ),
        #     explicit_exclude_models=frozenset({'CARDAMOM'})
        # )
        #print(li)
        li = h.list_target_models(
            target_classes=frozenset(
                {
                    CompartmentalMatrix,
                    StateVariableTuple,
                    VegetationCarbonInputTuple
                }
            ),
            # explicit_exclude_models=frozenset({'CARDAMOM'})
        )
        print(li)

    def test_list_models(self):
        li = h.list_models(explicit_exclude_models=frozenset({'CARDAMOM'}))
        print(li)

    def test_graph(self):
        # for debugging we draw the sparse_powerset_graph of the actually present computers and mvars
        spsg = sparse_powerset_graph(bgc_md2_computers())
        f = plt.figure()
        ax = f.add_subplot(1, 1, 1)
        draw_ComputerSetMultiDiGraph_matplotlib(
                ax,
                spsg,
                bgc_md2_mvar_aliases(),
                bgc_md2_computer_aliases(),
                targetNode=frozenset({SmoothModelRun})
        )
        f.savefig("spgs.pdf")
        fig = plt.figure()
        draw_update_sequence(bgc_md2_computers(), max_it=8, fig=fig)
        fig.savefig("c1.pdf")


    def test_all_computable_mvars_for_all_models(self):
        # https://docs.python.org/3/library/unittest.html#distinguishing-test-iterations-using-subtests
        model_names = h.list_models(
                #explicit_exclude_models=frozenset({'CARDAMOM'})
        )
        for mn in model_names:
            with self.subTest(mn=mn):
                mvs = MVarSet.from_model_name(mn)
                mvars = mvs.computable_mvar_types()
                list_str = "\n".join(["<li> " + str(var.__name__) + " </li>" for var in mvars])
                print(list_str)
                for var in mvars:
                    print("########################################")
                    print(str(var.__name__))
                    print(mvs._get_single_mvar_value(var))
# to do:
# 1.) create a jupyter notebook that implements something like the last three
# tests and recognizes an interactively changed model file. (Supports working
# on the notebook and the model file in parallel
#
# 2.) Create a table of models with columns:name,compartmental matrix
