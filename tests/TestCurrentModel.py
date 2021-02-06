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
#from bgc_md2.models.helpers import (
#    provided_mvars,
#    computable_mvars,
#    path_dict_to_single_mvar,
#    get_single_mvar_value,
#    bgc_md2_computers,
#    bgc_md2_computer_aliases,
#    bgc_md2_mvar_aliases,
#)
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
from testinfrastructure.helpers import pp
from bgc_md2.resolve.MVarSet import MVarSet


class TestCurrentModel(InDirTest):
#class TestLuo(InDirTest):
    def setUp(self):
        self.mn = "CARDAMOM"
#        self.mn = "Murty2000EcolModell"
#        self.mn = "ElMasri2013AgricForMeteorol" #Paper shows results for soil carbon, and fig. 1 has litter and C pools, but the equations on table A2 (although very detailed) don't include them
#        self.mn = "Scheiter2009GlobalChangeBiol" #No soil compartments
#        self.mn = "Turgman2018EcologyLetters" #No soil compartments
#        self.mn = "Haverd2016Biogeosciences" #No soil compartments
#        self.mn = "Foley1996GBC" #No equations for litter and soil, but the figure has those compartments
#        self.mn = "Gu2010EcologicalComplexity" #No equations for litter and soil, but the model description (CEVSA) metions them
#        self.mn = "King1993TreePhysiol" #No soil compartments
#        self.mn = "DeAngelis2012TheorEcol" #No soil compartments (model based on G’Day, but removed litter and soil compartments)
#        self.mn = "Potter1993GlobalBiogeochemicalCycles" #No equations for litter and soil, but the model description (CEVSA) metions them
        self.mvs = MVarSet.from_model_name(self.mn)
        self.ref_provided_mvars = frozenset(
            [
                CompartmentalMatrix,
                TimeSymbol,
                StateVariableTuple,
                VegetationCarbonInputPartitioningTuple,
                VegetationCarbonInputScalar,
                VegetationCarbonStateVariableTuple,
                InputTuple,
#                # NumericStartValueDict,
#                # NumericSimulationTimes,
#                 NumericParameterization,
#                # QuantityStartValueDict,
#                # QuantitySimulationTimes,
#                # QuantityParameterization,
                BibInfo,
#                # QuantityModelRun,
#                # QuantityParameterizedModel
            ]
        )

#    def test_provided_mvars(self):
#        mvs = self.mvs 
#        self.assertSetEqual(mvs.provided_mvar_types, self.ref_provided_mvars)

#    @unittest.skip
    def test_computable_mvars(self):
        spsg=sparse_powerset_graph(bgc_md2_computers())
        f = plt.figure()
        ax = f.add_subplot(1,1,1)
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
        mvs=self.mvs
        mvars = mvs.computable_mvar_types()
        list_str = "\n".join(["<li> " + str(var.__name__) + " </li>" for var in mvars])
        print(list_str)
        for var in mvars:
            print("########################################")
            print(str(var.__name__))
            print(mvs._get_single_mvar_value(var))
