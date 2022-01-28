import matplotlib.pyplot as plt
from networkx.algorithms.operators.binary import difference
import bgc_md2.helper as bh 
from unittest import TestCase
from testinfrastructure.InDirTest import InDirTest
import ComputabilityGraphs.fast_graph_helpers as fgh
import ComputabilityGraphs.rec_graph_helpers as rgh
import ComputabilityGraphs.dep_graph_helpers as dgh
import ComputabilityGraphs.helpers as h
import ComputabilityGraphs.graph_helpers as gh
import ComputabilityGraphs.graph_plotting as gp
from ComputabilityGraphs.helpers import all_mvars
from ComputabilityGraphs.TypeSet import TypeSet

from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import *

class TestSpecialModels(InDirTest):
    def test_bgc_graph_computation(self):
        #possible_types = all_mvars(bgc_md2_computers())
        #print(TypeSet(possible_types))
        possible_types={
            #NumericSolutionArray,
	    #TimeSymbol,
	    #VegetationCarbonInputTuple,
	    #CompartmentalMatrix,
	    #QuantityStartValueDict,
	    #NumericParameterizedSmoothReservoirModel,
	    #CarbonInputTuple,
	    #NitrogenCompartmentalMatrix,
	    #SmoothReservoirModel,
	    #NitrogenInFluxesBySymbol,
	    #NumericStartValueArray,
	    #CarbonStateVariableTuple,
	    #VegetationCarbonCompartmentalMatrix,
	    #VegetationCarbonInternalFluxesBySymbol,
	    #StateVarUnitTuple,
	    #InternalFluxesBySymbol,
	    #VegetationCarbonInFluxesBySymbol,
	    #QuantityParameterizedSmoothReservoirModel,
	    #InputTuple,
	    #VegetationCarbonStateVariableTuple,
	    #QuantitySolutionArray,
	    #CarbonOutFluxesBySymbol,
	    #NitrogenInputTuple,
	    #QuantitySimulationTimes,
	    #NitrogenInternalFluxesBySymbol,
	    #InFluxesBySymbol,
	    #VegetationCarbonInputScalar,
	    #QuantityStartValueArray,
	    #CarbonCompartmentalMatrix,
	    #NitrogenOutFluxesBySymbol,
	    #NumericParameterization,
	    VegetationCarbonInputPartitioningTuple,
	    #CarbonInternalFluxesBySymbol,
	    #NumericSimulationTimes,
	    #Quantity,
	    #CarbonInFluxesBySymbol,
	    #NitrogenStateVariableTuple,
	    #VegetationCarbonOutFluxesBySymbol,
	    #NumericStartValueDict,
	    #QuantityModelRun,
	    #OutFluxesBySymbol,
	    #StateVariableTuple,
	    #QuantityParameterization,
	    #SmoothModelRun,
	    #_empty
        }

        for root in possible_types:
            with self.subTest(root=root):
                print(root)
                g = rgh.fast_graph(
                    cs=bh.bgc_md2_computers(),
                    root_type=root,
                    given=frozenset()
                )
        

    def test_all_computable_mvars_for_all_models(self):
        # https://docs.python.org/3/library/unittest.html#distinguishing-test-iterations-using-subtests
        model_names = bh.list_models(
                #explicit_exclude_models=frozenset({'CARDAMOM'})
        )
        print(model_names)
        #for mn in model_names:
        for mn in [
            #'ACGCA',
            #'ACGCASoilModel',
            #'ACGCAWoodProductModel',
            #'Allison2010NG',
            #'Andren1997EA',
            #'Arora2005GCB-1',
            #'CARDAMOM',
            #'Castanho2013Biogeosciences',
            'Comins1993EA',
            #'DeAngelis2012TheorEcol',
            #'ElMasri2013AgricForMeteorol',
            #'Emanuel1981',
            #'Foley1996GBC',
            #'Fontaine2005Ecologyletters',
            #'Fontaine2005Ecologyletters_1',
            #'Fontaine2005Ecologyletters_2',
            #'Fontaine2005Ecologyletters_3_1',
            #'Fontaine2005Ecologyletters_3_2',
            #'Fontaine2005Ecologyletters_4_1',
            #'Fontaine2005Ecologyletters_4_2',
            #'Gu2010EcologicalComplexity',
            #'Haverd2016Biogeosciences',
            #'Henin1945AA',
            # 'Hilbert1991AnnBot',
            #'Jenkinson1977SoilScience',
            #'King1993TreePhysiol',
            #'Luo2012TE',
            #'Murty2000EcolModell',
            'Parton1987SoilSciSocAmJ',
            #'Pavlick2013Biogeosciences',
            #'Potter1993GlobalBiogeochemicalCycles',
            #'Rasmussen2016JMB',
            #'Running1988EcolModel',
            'Scheiter2009GlobalChangeBiol',
            #'TECO',
            #'TECOmm',
            #'Thomas2014GeosciModelDev',
            #'Turgman2018EcologyLetters',
            #'Wang2010Biogeosciences',
            #'Wang2013EcologicalApplications',
            #'Wang2014BG2p',
            #'Wang2014BG3p',
            #'Wieder2014Biogeosciences',
            #'Williams2005GCB',
            #'Zelenev2000MicrobialEcology',
            #'cable_all',
            #'sixPairsModel',
            #'testVectorFree'
        ]:

            with self.subTest(mn=mn):
                mvs = bh.CMTVS_from_model_name(mn)
                mvars = mvs.computable_mvar_types()
                list_str = "\n".join(["<li> " + str(var.__name__) + " </li>" for var in mvars])
                print(list_str)
                for var in mvars:
                    print("########################################")
                    print(str(var.__name__))
                    print(mvs._get_single_value(var))

    

    def test_special_mvar_for_special_model(self):
        mn = 'Thomas2014GeosciModelDev'
        root =  VegetationCarbonInputPartitioningTuple
        mvs = bh.CMTVS_from_model_name(mn)
        fg=rgh.fast_graph(root_type=root,cs=mvs.computers, given=mvs.provided_mvar_values)
        fig=plt.figure(figsize=(40,40));
        ax=fig.subplots(1,1)
        fg.draw_matplotlib(ax)
        fig.savefig("fast_graph.pdf")
        plt.close(fig)
        
        fig=plt.figure(figsize=(40,40));
        ax=fig.subplots(1,1)
        gp.draw_ComputerSetMultiDiGraph_matplotlib(
            ax,
            fgh.project_to_multiDiGraph(fg)
        )
        fig.savefig("graph.pdf")
        plt.close(fig)
        print(mvs._get_single_value(root))
        
        g=dgh.shortest_computable_dep_graph(root_type=root,cs=mvs.computers, given=mvs.provided_mvar_types)
        fig=plt.figure()
        ax=fig.subplots(1,1)
        g.draw_matplotlib(ax)
        fig.savefig("dep_graph.pdf")
        plt.close(fig)
        print(mvs._get_single_value_by_depgraph(root))
