import unittest
import matplotlib.pyplot as plt
from string import ascii_lowercase, ascii_uppercase

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
    InputTuple,
    NumericParameterization,
    NumericStartValueDict,
    NumericSimulationTimes,
    NumericParameterizedSmoothReservoirModel,
    QuantityParameterization,
    QuantityStartValueDict,
    QuantitySimulationTimes,
    QuantityParameterizedSmoothReservoirModel,
    QuantityModelRun,
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
from testinfrastructure.helpers import pp
from bgc_md2.resolve.MVarSet import MVarSet


class TestWilliams(InDirTest):
    def setUp(self):
        self.mn = "Williams2005GCB"
        self.mvs = MVarSet.from_model_name(self.mn)
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
                QuantityStartValueDict,
                QuantitySimulationTimes,
                QuantityParameterization,
                BibInfo,
                # QuantityModelRun,
                # QuantityParameterizedModel
            ]
        )

    def test_provided_mvars(self):
        mvs = self.mvs 
        self.assertSetEqual(mvs.provided_mvar_types, self.ref_provided_mvars)

