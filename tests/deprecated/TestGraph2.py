#!/usr/bin/env python3
import matplotlib.pyplot as plt
from frozendict import frozendict
from copy import copy, deepcopy
import networkx as nx
from testinfrastructure.InDirTest import InDirTest
from testinfrastructure.helpers import pp, pe
from unittest import skip
from copy import copy, deepcopy

from bgc_md2.resolve.graph_helpers import (
    arg_set_graph,
    initial_sparse_powerset_graph,
    minimal_startnodes_for_single_var,
    minimal_startnodes_for_node,
    update_step,
    toDiGraph,
    equivalent_singlegraphs,
    equivalent_multigraphs,
    node_2_string,
    nodes_2_string,
    edge_2_string,
    product_graph,
    product_graph_2,
    sparse_powerset_graph,
    update_generator,
    # ,draw_multigraph_plotly
    # ,draw_Graph_svg
)
from bgc_md2.resolve.graph_plotting import (
    draw_update_sequence,
    draw_sequence,
    draw_ComputerSetDiGraph_matplotlib,
    draw_ComputerSetMultiDiGraph_matplotlib,
)
from bgc_md2.resolve.non_graph_helpers import arg_set_set, all_mvars


class A:
    pass


class B:
    pass


class C:
    pass


class D:
    pass


class E:
    pass


class F:
    pass


class G:
    pass


class H:
    pass


class I:
    pass


class J:
    pass


class K:
    pass


def a_from_i(i: I) -> A:
    return A()


def d_from_c_i(c: C, i: I) -> D:
    return F()


def f_from_j_g(j: J, g: G) -> F:
    return G()


def h_from_d_e_f(d: D, e: E, f: F) -> H:
    return H()


def i_from_a_b_k_j(a: A, b: B, k: K, j: J) -> I:
    return I()


# for easier debugging in ipython
computers = frozenset({a_from_i, d_from_c_i, f_from_j_g, h_from_d_e_f, i_from_a_b_k_j})


class TestGraph2(InDirTest):
    def setUp(self):
        # fixme mm 02-03-2020:
        # assemble this set from a file of class definitions
        self.computers = computers

    def test_arg_set_graph(self):
        asg_I = arg_set_graph(I, self.computers)
        asg_B = arg_set_graph(B, self.computers)
        # For compatibility arg_set_graph returns a multigraph
        # although we do not have more than one edge between a pair
        # of nodes.

        ref_I = nx.MultiDiGraph()
        ref_I.add_edge(
            frozenset({A, B, K, J}),
            frozenset({I}),
            computers=frozenset({i_from_a_b_k_j}),
        )

        ref_B = nx.MultiDiGraph()
        ref_B.add_node(frozenset({B}))
        # picture for manual check
        fig = plt.figure(figsize=(20, 20))
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)
        draw_ComputerSetMultiDiGraph_matplotlib(ax1, ref_I)
        draw_ComputerSetMultiDiGraph_matplotlib(ax2, asg_I)
        draw_ComputerSetMultiDiGraph_matplotlib(ax3, ref_B)
        draw_ComputerSetMultiDiGraph_matplotlib(ax4, asg_B)
        fig.savefig("arg_set_graph.pdf")

        self.assertTrue(equivalent_multigraphs(asg_I, ref_I))
        self.assertTrue(equivalent_multigraphs(asg_B, ref_B))

    def test_product_graph(self):

        asg_C = arg_set_graph(C, computers)
        asg_E = arg_set_graph(E, computers)
        asg_G = arg_set_graph(G, computers)
        asg_I = arg_set_graph(I, computers)
        asg_J = arg_set_graph(J, computers)
        # cart_CE=nx.cartesian_product(asg_C,asg_E)
        # pg2_CE=product_graph_2(asg_C,asg_E)
        # pg_CE=product_graph(asg_C,asg_E)
        pg_CI = product_graph(asg_C, asg_I)
        pg_CEGIJ = product_graph(asg_C, asg_E, asg_G, asg_I, asg_J)
        fig1 = plt.figure(figsize=(20, 100))
        draw_sequence(
            fig1,
            tups=(
                (asg_C, "arg_set_graph(C)"),
                (asg_E, "arg_set_graph(E)"),
                (asg_G, "arg_set_graph(G)"),
                (asg_I, "arg_set_graph(I)"),
                (asg_J, "arg_set_graph(J)"),
                # (cart_CE,"cart_CE"),
                # (pg2_CE,"pg2_CE"),
                # (pg_CE,"pg_CE")
                (pg_CI, "pg_CI"),
                (pg_CEGIJ, "pg_CEGIJ"),
            ),
        )

        fig1.savefig("AB_Z.pdf")
        plt.close(fig1)
        # fig2=plt.figure(figsize=(20,20))
        # nr=1
        # axs=fig.subplots(nr,1,sharex=True,sharey=True)
        # nx.draw(cart_CE)

    def test_minimal_startnodes_for_single_var(self):
        spsg = sparse_powerset_graph(self.computers)
        ## After the graph has been computed we can use it
        ## to infer computability of all Mvars
        fig1 = plt.figure(figsize=(20, 20))
        axs = fig1.subplots(1, 1)
        draw_ComputerSetMultiDiGraph_matplotlib(axs, spsg)
        fig1.savefig("spsg.pdf")
        plt.close(fig1)

        fig2 = plt.figure(figsize=(20, 100))
        draw_update_sequence(self.computers, max_it=8, fig=fig2)
        fig2.savefig("c1.pdf")
        plt.close(fig2)
        res = minimal_startnodes_for_single_var(spsg, H)
        print(nodes_2_string(res))
        path = nx.single_source_shortest_path(
            spsg, source=frozenset({A, B, C, E, G, J, K})
        )
        print("path", nodes_2_string(path))
        # self.assertSetEqual(
        #    res,
        #    frozenset({
        #        frozenset({E,F}),
        #        frozenset({C,D}),
        #        frozenset({H,C,G})
        #    })
        # )


#
#
#    def test_minimal_startnodes_for_node(self):
#        spsg=sparse_powerset_graph(self.computers)
#        targetVars=frozenset({A,B})
#        res=minimal_startnodes_for_node(
#                 spsg
#                ,targetVars
#        )
#        #print('###############################')
#        #print('miminal startsets for: ', node_2_string(targetVars),nodes_2_string(res))
#        #print('###############################')
#        self.assertSetEqual(
#            res,
#            frozenset({
#                frozenset({E,F,I}),
#                frozenset({C,D,I}),
#                frozenset({H,C,G,I})
#            })
#        )
#
