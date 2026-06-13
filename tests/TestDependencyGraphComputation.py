from time import time as now
from unittest import TestCase, skip
import matplotlib.pyplot as plt
from testinfrastructure.InDirTest import InDirTest
from ComputabilityGraphs.rec_graph_helpers import fast_graph
from ComputabilityGraphs.dep_graph_helpers import (
    duplicated_computer_dict,
    all_dep_graphs,
)
from ComputabilityGraphs.helpers import all_mvars
# from ComputabilityGraphs.TypeSet import TypeSet
from ComputabilityGraphs.or_graph_helpers import t_tree
# from ComputabilityGraphs.graph_plotting import (
#     draw_update_sequence,
#     draw_ComputerSetDiGraph_matplotlib,
#     draw_ComputerSetMultiDiGraph_matplotlib,
# )

import bgc_md2.helper as h
from bgc_md2.resolve.mvars import *


class TestDependencyGraphComputation(InDirTest):
    
    @skip
    def test_duplicated_computers(self):
        dcd = duplicated_computer_dict(h.bgc_md2_computers())
        for key, val in dcd.items():
            print(key, len(val))
            for f in val:
                print(f.__name__)

    def test_bgc_or_graph_computation(self):
        cs = h.bgc_md2_computers()
        possible_types = all_mvars(cs).difference(
            {VegetationCarbonInputPartitioningTuple}
        )
        # print(TypeSet(possible_types))
        for t in possible_types:
            print(t.__name__)

        for root in possible_types:
            before_og = now()
            with self.subTest(root=root):
                at = frozenset({})
                tt = t_tree(
                    root_type=root,
		            available_computers =cs,
		            avoid_types=at
                    #given=frozenset({})
                )
                fig = plt.figure()
                ax = fig.subplots(1, 1)
                og = tt.to_networkx_graph(avoid_types=frozenset({}))
                og.draw_matplotlib(ax)
                fig.savefig("or_graphs_" + root.__name__ + ".pdf")
                plt.close(fig)

            after_dg = now()
            time = after_dg - before_og
            # if time >1:
            if True:
                print("or_graph time", root, time)

    @skip
    def test_bgc_dep_graph_computation(self):
        cs = h.bgc_md2_computers()
        possible_types = all_mvars(cs).difference(
            {VegetationCarbonInputPartitioningTuple}
        )
        # print(TypeSet(possible_types))
        for t in possible_types:
            print(t.__name__)

        for root in possible_types:
            before_dg = now()
            with self.subTest(root=root):
                gs = list(all_dep_graphs(root_type=root, cs=cs, given=frozenset()))
                fig = plt.figure()
                axs = fig.subplots(len(gs), 1)
                for i, ax in enumerate(axs):
                    g = gs[i]
                    g.draw_matplotlib(ax)
                fig.savefig("depgraphs_" + root.__name__ + ".pdf")
                plt.close(fig)

            after_dg = now()
            time = after_dg - before_dg
            # if time >1:
            if True:
                print("depgraph time", root, time, len(gs))

    @skip
    def test_bgc_fast_graph_computation(self):
        cs = h.bgc_md2_computers()
        possible_types = all_mvars(cs).difference(
            {VegetationCarbonInputPartitioningTuple}
        )
        # print(TypeSet(possible_types))
        for t in possible_types:
            print(t.__name__)

        for root in possible_types:
            before_dg = now()
            with self.subTest(root=root):
                before_fg = now()
                fg = fast_graph(cs=cs, root_type=root, given=frozenset())
                after_fg = now()
                time = after_fg - before_fg
                if time > 1:
                    print("fast_graph time", root, time)

                # fig=plt.figure(figsize=(20,20))
                # ax1 = fig.add_subplot(1, 1, 1)
                # fg.draw_matplotlib(
                #    ax1
                #
                # )
                # fig.savefig(root.__name__+'.pdf')
            after_dg = now()
            time = after_dg - before_dg
            # if time >1:
            if True:
                print("fast_graph time", root, time)
