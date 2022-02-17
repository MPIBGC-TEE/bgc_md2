from time import time as now
from unittest import TestCase, skip
import matplotlib.pyplot as plt
from testinfrastructure.InDirTest import InDirTest
from ComputabilityGraphs.rec_graph_helpers import fast_graph
from ComputabilityGraphs.dep_graph_helpers import duplicated_computer_dict, all_dep_graphs
from ComputabilityGraphs.helpers import all_mvars
from ComputabilityGraphs.TypeSet import TypeSet
from ComputabilityGraphs.graph_plotting import (
    draw_update_sequence,
    draw_ComputerSetDiGraph_matplotlib,
    draw_ComputerSetMultiDiGraph_matplotlib,
)

import bgc_md2.helper as h 
from bgc_md2.resolve.mvars import *

class TestDependencyGraphComputation(InDirTest):
    def test_duplicated_computers(self):
        dcd = duplicated_computer_dict(h.bgc_md2_computers())
        for key,val in dcd.items():
            print(key,len(val))
            for f in val:
                print(f.__name__) 
        # fixme mm 1/17/2022
        # check if the duplicated computers are redundant
        # If the arguments of computer c1 can always be computed from the arguments 
        # of computer c2 then computer c2 is redundant. 
        # we could try to implement this check and run it here to inform 
        # the author of the computermodule that they could get rid of a function
        # (This has several benefits:
        # - The graph computations become much faster since 
        #   the number of computercombies decreases.
        # - Contradictions are avoided since there is only one way to compute 
        #   a desired result from the same set of given values.

        # To implement the check we could do the following.
        # - for every argument of c1 compute all dep_graphs 
        #   (the arguments themselfes might be computable in different ways)
        #   compute all possible unions of all the leafs of these graphs (thus
        #   provide all the combinations of uncomputable (given) values
        #   sufficient to compute the arguments of c1 and apply it) the result
        #   is a set of sets of types.
        # - do the same for c2
        # - do the same for c3 ...
        # - assert the sets are not equal 
        #   this has to be done for every combination of computers with the same result
        #   for 3 computers computing A we would have to have 3-1+1+3 tests
        #   (leafsetset(c1)!=leafsetset(c2); leafsetset(c1)!=leafsetset(c3);
        #   leafsetset(c3)!=leafsetset(c3) 
        #   for 4 computers with the same result we would have (4-1)+(3-1)+(2-1) tests
       
    def test_bgc_dep_graph_computation(self):
        cs=h.bgc_md2_computers()
        possible_types = all_mvars(cs).difference({VegetationCarbonInputPartitioningTuple})
        #print(TypeSet(possible_types))
        for t in possible_types:
            print(t.__name__)

        for root in possible_types:
            before_dg=now()
            with self.subTest(root=root):
                gs= list(
                    all_dep_graphs(
                        root_type=root,
                        cs=cs,
                        given=frozenset()
                    )
                )
                fig=plt.figure()
                axs = fig.subplots(len(gs),1)
                for i,ax in enumerate(axs):
                    g = gs[i]
                    g.draw_matplotlib(ax)
                fig.savefig("depgraphs_" + root.__name__ + '.pdf')


            
            after_dg=now()
            time=(after_dg - before_dg)
            #if time >1:
            if True:
                print("depgraph time",root, time,len(gs))


    def test_bgc_fast_graph_computation(self):
        cs=h.bgc_md2_computers()
        possible_types = all_mvars(cs).difference({VegetationCarbonInputPartitioningTuple})
        #print(TypeSet(possible_types))
        for t in possible_types:
            print(t.__name__)

        for root in possible_types:
            before_dg=now()
            with self.subTest(root=root):
                before_fg=now()
                fg = fast_graph(
                    cs=cs,
                    root_type=root,
                    given=frozenset()
                )
                after_fg=now()
                time=(after_fg - before_fg)
                if time >1:
                    print("fast_graph time", root, time)
                
                fig=plt.figure(figsize=(20,20))
                ax1 = fig.add_subplot(1, 1, 1)
                fg.draw_matplotlib(
                    ax1
                
                )
                fig.savefig(root.__name__+'.pdf')
            after_dg=now()
            time=(after_dg - before_dg)
            #if time >1:
            if True:
                print("fast_graph time",root, time)
                


