import unittest
import matplotlib.pyplot as plt
from string import ascii_lowercase, ascii_uppercase

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun

from testinfrastructure.InDirTest import InDirTest

from bgc_md2.models.helpers import (
    provided_mvars,
    computable_mvars,
    path_dict_to_single_mvar,
    get_single_mvar_value,
    bgc_md2_computers,
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
    InputTuple,
    NumericParameterization,
    NumericStartValueDict,
    NumericSimulationTimes,
    NumericParameterizedSmoothReservoirModel,
    QuantityModelRun,
    QuantityParameterizedModel,
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



def computer_color_func(allComputers):
    colordict = TABLEAU_COLORS
    color_names = [n for n in colordict.keys()]
    color_dict = {c.__name__: color_names[i] for i, c in enumerate(allComputers)}

    def cf(c):
        col = colordict[color_dict[c.__name__]]
        return col

    return cf


class TestWilliams(InDirTest):
    def setUp(self):
        self.mn = "Williams2005GCB"
        self.ref_provided_mvars = frozenset(
            [
                CompartmentalMatrix,
                TimeSymbol,
                StateVariableTuple,
                VegetationCarbonInputPartitioningTuple,
                VegetationCarbonInputScalar,
                InputTuple,
                NumericStartValueDict,
                NumericSimulationTimes,
                NumericParameterization,
                BibInfo,
                # QuantityModelRun,
                # QuantityParameterizedModel
            ]
        )

    def test_provided_mvars(self):
        mvs = provided_mvars(self.mn)
        self.assertSetEqual(mvs, self.ref_provided_mvars)

    def test_computable_mvars(self):
        # fixme: the next lines want to move to resolve/computers 
        self.bgc_computers = bgc_md2_computers()
        self.computer_aliases = {
            v.__name__: ascii_lowercase[i] for i, v in enumerate(self.bgc_computers)
        }
        allVars = all_mvars(self.bgc_computers)
        self.mvar_aliases = {
            name: ascii_uppercase[i]
            for i, name in enumerate(sorted(map(lambda v: v.__name__, allVars)))
        }
        # end_fixme
        spsg=sparse_powerset_graph(bgc_md2_computers())
        f = plt.figure()
        ax = f.add_subplot(1,1,1)
        draw_ComputerSetMultiDiGraph_matplotlib(
                ax,
                spsg, 
                self.mvar_aliases, 
                self.computer_aliases,
                targetNode=frozenset({SmoothModelRun})
        )
        f.savefig("spgs.pdf")
        fig = plt.figure()
        draw_update_sequence(bgc_md2_computers(), max_it=8, fig=fig)
        fig.savefig("c1.pdf")
        res = frozenset(
            [
                VegetationCarbonInputTuple,
                SmoothReservoirModel,
                SmoothModelRun,
                NumericParameterizedSmoothReservoirModel,
            ]
        ).union(self.ref_provided_mvars)
        mvs = computable_mvars(self.mn)
        # self.assertSetEqual(mvs,res)
        list_str = "\n".join(["<li> " + str(var.__name__) + " </li>" for var in mvs])
        print(list_str)
        # for var in mvs:
        #    print("########################################")
        #    print(str(var.__name__))
        #    print(get_single_mvar_value(var,self.mn))
        mr = get_single_mvar_value(SmoothModelRun, "Williams2005GCB")
