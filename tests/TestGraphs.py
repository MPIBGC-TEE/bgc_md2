#!/usr/bin/env python3
import matplotlib.pyplot as plt
from frozendict import frozendict
from copy import copy,deepcopy
import networkx as nx
from testinfrastructure.InDirTest import InDirTest
from testinfrastructure.helpers import pp,pe
from unittest import skip
from copy import copy,deepcopy

from bgc_md2.resolve.graph_helpers import (
    arg_set_graph
    ,initial_sparse_powerset_graph
    ,minimal_startnodes_for_single_var
    ,minimal_startnodes_for_node
    ,update_step
    ,toDiGraph
    ,equivalent_singlegraphs
    ,equivalent_multigraphs
    ,node_2_string
    ,nodes_2_string
    ,edge_2_string
    ,product_graph
    ,sparse_powerset_graph
    ,update_generator
    #,draw_multigraph_plotly
    #,draw_Graph_svg
)
from bgc_md2.resolve.graph_plotting import (
    draw_update_sequence
    ,draw_ComputerSetDiGraph_matplotlib
    ,draw_ComputerSetMultiDiGraph_matplotlib
)
from bgc_md2.resolve.non_graph_helpers import  arg_set_set,all_mvars

class A_minus_1:
    pass

class A_minus_2:
    pass



class A:
    """A variable we assume to be given"""
    pass

class B:
    pass

class A_minus_2:
    pass

class A_minus_1:
    pass

class A0:
    pass

class A1:
    pass

class A2:
    pass

class A3:
    pass

class B_minus_2:
    pass


class B_minus_1:
    pass

class B0:
    pass

class B1:
    pass

class B2:
    pass

class B3:
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

class X:
    pass

class Y:
    pass


class Z:
    pass

def a_from_b_c(b: B, c: C) -> A:
    return A()

def a_from_b_d(b: B, d: D) -> A:
    return A()

def a_from_x(x:X)->A:
    return A()

def b_from_x(x:X)->B:
    return B()

def a_from_y(y:Y)->A:
    return A()

def b_from_y(y:Y)->B:
    return B()

def a_from_z(z:Z)->A:
    return A()

def b_from_z(z:Z)->B:
    return B()

def c_from_z(z:Z)->C:
    return C()

def a_from_i(i: I) -> A:
    """Computes a from i"""
    return A()


def b_from_c_d(c: C, d: D) -> B:
    return B()

def b_from_d_e(d: D,e: E) -> B:
    return B()

def b_from_e_f(e: E, f: F) -> B:
    return B()


def c_from_b(b: B) -> C:
    """Computes c from b"""
    return C()

def d_from_a(a: A) -> D:
    return D()


def d_from_b(b: B) -> D:
    """Computes d from b"""
    return D()


def d_from_g_h(g: G, h: H) -> D:
    """Computes d from g and h"""
    return D()


def e_from_b(b: B) -> E:
    """Computes e from b"""
    return E()


def f_from_b(b: B) -> F:
    """Computes f from b"""
    return F()

def a_minus_1_from_a_minus_2(x:A_minus_2)->A_minus_1:
    return A_minus_1()

def a0_from_a_minus_1(x:A_minus_1)->A0:
    return A0()

def a1_from_a0(a0:A0) -> A1:
    return A1()

def a2_from_a1(a1:A1) -> A2:
    return A2()

def a3_from_a2(a2:A2) -> A3:
    return A3()

def a0_from_b0(x:B0)->A0:
    return A0()

def b_minus_1_from_b_minus_2(x:B_minus_2)->B_minus_1:
    return B_minus_1()

def b0_from_b_minus_1(x:B_minus_1)->B0:
    return B0()

def b1_from_b0(b0:B0) -> B1:
    return B1()

def b2_from_b1(b1:B1) -> B2:
    return B2()

def b3_from_b2(b2:B2) -> B3:
    return B3()

#for easier debugging in ipython
computers = frozenset({
            a_from_i,
            b_from_c_d,
            b_from_e_f,
            c_from_b,
            d_from_b,
            d_from_g_h,
            e_from_b,
            f_from_b,
})

class TestGraphs(InDirTest):
    def setUp(self):
        # fixme mm 02-03-2020:
        # assemble this set from a file of class definitions
        self.mvars = {
            A,
            B,
            C,
            D,
            E,
            F,
            G,
            H,
            I,
        }
        # fixme mm 02-03-2020:
        # assemble this set from a file of annotated functions in
        # special file
        self.computers = computers

    def test_toDiGraph(self):
        g_multi=nx.MultiDiGraph()
        src =  frozenset({X,Y})
        dest = frozenset({A,B})
        computers_1=frozenset({ a_from_x,b_from_y})

        computers_2=frozenset({ a_from_y,b_from_x})
        g_multi.add_edge(
            src
            ,dest
            ,computers=computers_1
        )
        g_multi.add_edge(
            src
            ,dest
            ,computers=computers_2
        )
        g_single=toDiGraph(g_multi)
        self.assertSetEqual(
            g_single.get_edge_data(src,dest)['computers'],
            frozenset({
                computers_1,
                computers_2
            })
        )

    def test_equivalent_multigraphs(self):
        g1=nx.MultiDiGraph()
        g1.add_edge(
            frozenset({X,Y})
            ,frozenset({A,B})
            ,computers=frozenset({ a_from_x,b_from_y})
        )
        g1.add_edge(
            frozenset({X,Y})
            ,frozenset({A,B})
            ,computers=frozenset({ a_from_y,b_from_x})
        )

        g2=deepcopy(g1)
        # now we exchange the edge data
        g2.add_edge(
            frozenset({X,Y})
            ,frozenset({A,B})
            ,computers=frozenset({ a_from_y,b_from_x})
        )
        g2.add_edge(
            frozenset({X,Y})
            ,frozenset({A,B})
            ,computers=frozenset({ a_from_x,b_from_y})
        )

        self.assertTrue(equivalent_multigraphs(g1,g2))

    def test_arg_set_graph(self):
        asg=arg_set_graph(D,self.computers)
        # For compatibility arg_set_graph returns a multigraph
        # although we do not have more than one edge between a pair
        # of nodes.

        ref=nx.MultiDiGraph()
        ref.add_edge(
            frozenset({B})
            ,frozenset({D})
            ,computers=frozenset({d_from_b})
        )
        ref.add_edge(
            frozenset({G, H})
            ,frozenset({D})
            ,computers=frozenset({d_from_g_h})
        )

        #picture for manual check
        fig=plt.figure(figsize=(20,20))
        ax1=fig.add_subplot(1,2,1)
        ax2=fig.add_subplot(1,2,2)
        draw_ComputerSetMultiDiGraph_matplotlib(ax1,ref)
        draw_ComputerSetMultiDiGraph_matplotlib(ax2,asg)
        fig.savefig("arg_set_graph.pdf")

        self.assertTrue(equivalent_multigraphs(
            asg , ref
        ))

    def test_product_graph(self):
        computers=frozenset({a_from_y,a_from_z,b_from_y,b_from_z})

        asg_A=arg_set_graph(A,computers)
        asg_B=arg_set_graph(B,computers)
        pg_A_B=product_graph(asg_A,asg_B)

        fig1=plt.figure(figsize=(20,100))

        ax1=fig1.add_subplot(411,frame_on=True,title="arg_set_graph(A)")
        ax2=fig1.add_subplot(412,frame_on=True,title="arg_set_graph(B)")
        ax3=fig1.add_subplot(413,frame_on=True,title="product_graph(A,B)")
        #ax4=fig1.add_subplot(414,frame_on=True,title="ref")
        draw_ComputerSetMultiDiGraph_matplotlib(ax1,    asg_A)
        draw_ComputerSetMultiDiGraph_matplotlib(ax2,    asg_B)
        draw_ComputerSetMultiDiGraph_matplotlib(ax3,   pg_A_B)
        fig1.savefig('AB_Z.pdf')

        computers=frozenset({a_from_y,b_from_y,b_from_z})
        asg_A=arg_set_graph(A,computers)
        asg_B=arg_set_graph(B,computers)
        pg_A_B=product_graph(asg_A,asg_B)

        fig1=plt.figure(figsize=(10,30))
        ax1=fig1.add_subplot(411,frame_on=True,title="arg_set_graph(A)")
        ax2=fig1.add_subplot(412,frame_on=True,title="arg_set_graph(B)")
        ax3=fig1.add_subplot(413,frame_on=True,title="product_graph(A,B)")
        #ax4=fig1.add_subplot(414,frame_on=True,title="ref")
        draw_ComputerSetMultiDiGraph_matplotlib(ax1,    asg_A)
        draw_ComputerSetMultiDiGraph_matplotlib(ax2,    asg_B)
        draw_ComputerSetMultiDiGraph_matplotlib(ax3,   pg_A_B)
        fig1.savefig('A_y_B_Y_B_Z.pdf')

        #    {a(z)}          {b(z)}          {a(z),b(z)}

        # {z}  ->  {a} x {z}  ->   {b}   = {z} -> {a,b}

        computers=frozenset({ a_from_z,b_from_z})
        asg_A=arg_set_graph(A,computers)
        asg_B=arg_set_graph(B,computers)
        pg_A_B=product_graph(asg_A,asg_B)

        fig1=plt.figure(figsize=(10,30))
        ax1=fig1.add_subplot(411,frame_on=True,title="arg_set_graph(A)")
        ax2=fig1.add_subplot(412,frame_on=True,title="arg_set_graph(B)")
        ax3=fig1.add_subplot(413,frame_on=True,title="product_graph(A,B)")
        #ax4=fig1.add_subplot(414,frame_on=True,title="ref")
        draw_ComputerSetMultiDiGraph_matplotlib(ax1,     asg_A)
        draw_ComputerSetMultiDiGraph_matplotlib(ax2,     asg_B)
        draw_ComputerSetMultiDiGraph_matplotlib(ax3,    pg_A_B)
        fig1.savefig('A_Z_B_Y.pdf')

        computers=frozenset({ a_from_z,b_from_z,c_from_b})
        asg_A=arg_set_graph(A,computers)
        asg_B=arg_set_graph(B,computers)
        asg_C=arg_set_graph(C,computers)
        asg_C=arg_set_graph(C,computers)
        prod=product_graph(*[asg_A,asg_B,asg_C])

        fig1=plt.figure(figsize=(10,30))
        ax1=fig1.add_subplot(511,frame_on=True,title="arg_set_graph(A)")
        ax2=fig1.add_subplot(512,frame_on=True,title="arg_set_graph(B)")
        ax3=fig1.add_subplot(513,frame_on=True,title="arg_set_graph(C)")
        ax4=fig1.add_subplot(514,frame_on=True,title="product_graph(A,B,C)")
        #ax4=fig1.add_subplot(414,frame_on=True,title="ref")
        draw_ComputerSetMultiDiGraph_matplotlib(ax1,    asg_A)
        draw_ComputerSetMultiDiGraph_matplotlib(ax2,    asg_B)
        draw_ComputerSetMultiDiGraph_matplotlib(ax3,    asg_C)
        draw_ComputerSetMultiDiGraph_matplotlib(ax4,    prod)
        fig1.savefig('ABC.pdf')






    def test_Markus_graph_creation(self):
        # Now we build the directed Graph we can use to compute connectivity
        # the Nodes are sets of Mvars (elemenst of the powerset of all Mvars)
        # and a connection between two sets indicates computability of the target set from
        # the source set.
        # The complete graph would contain all elements of the powerset of allMvars and all
        # possible connections, which is prohibitively expensive.
        # Instead we will compute a subgraph where we start with one element sets as targets
        # and infer the predecessors of those sets and then the predecessors of the predecessors and so on until we do not find new nodes=start_sets.
        #spsg=sparse_powerset_graph(self.mvars,self.computers)

        ################# linear A1->A2->A3
        spsg=sparse_powerset_graph(
            frozenset({a2_from_a1,a3_from_a2})
        )
        self.assertSetEqual(
            set(spsg.nodes()),
            {
                frozenset({A1}),
                frozenset({A2}),
                frozenset({A3})
            }
        )
        self.assertSetEqual(
            set(spsg.edges()),
            {
                (frozenset({A1}),frozenset({A2})),
                (frozenset({A2}),frozenset({A3}))
            }
        )
        ################## cross
        #       B-2->B-1->B0->B1->B2
        #                  ||
        #                  \/
        #       A-2->A-1->A0->A1->A2
        spsg=sparse_powerset_graph(
            frozenset({
                b_minus_1_from_b_minus_2,
                b0_from_b_minus_1,
                b1_from_b0,
                b2_from_b1,
                b3_from_b2,
                #
                a0_from_b0,
                #
                a_minus_1_from_a_minus_2,
                a0_from_a_minus_1,
                a1_from_a0,
                a2_from_a1,
                a3_from_a2,
             })
        )
        self.assertSetEqual(
            set(spsg.nodes()),
            {
                frozenset({B_minus_1}),
                frozenset({B_minus_2}),
                frozenset({B0}),
                frozenset({B1}),
                frozenset({B2}),
                frozenset({B3}),
                frozenset({A_minus_1}),
                frozenset({A_minus_2}),
                frozenset({A0}),
                frozenset({A1}),
                frozenset({A2}),
                frozenset({A3}),
            }
        )
        self.assertSetEqual(
            set(spsg.edges()),
            {
                (frozenset({B_minus_2}),frozenset({B_minus_1})),
                (frozenset({B_minus_1}),frozenset({B0})),
                (frozenset({B0}),frozenset({B1})),
                (frozenset({B1}),frozenset({B2})),
                (frozenset({B2}),frozenset({B3})),
                #
                (frozenset({B0}),frozenset({A0})),
                #
                (frozenset({A_minus_2}),frozenset({A_minus_1})),
                (frozenset({A_minus_1}),frozenset({A0})),
                (frozenset({A0}),frozenset({A1})),
                (frozenset({A1}),frozenset({A2})),
                (frozenset({A2}),frozenset({A3})),
            }
        )
    def test_minimal_startnodes_for_single_var(self):
        spsg=sparse_powerset_graph(self.computers)
        fig1=plt.figure()
        axs=fig1.subplots(1,1)
        draw_ComputerSetMultiDiGraph_matplotlib(axs,spsg)
        fig1.savefig('spsg.pdf')
        plt.close(fig1)
        ## After the graph has been computed we can use it
        ## to infer computability of all Mvars
        res=minimal_startnodes_for_single_var(spsg,B)
        ##print(nodes_2_string(res))
        self.assertSetEqual(
            res,
            frozenset({
                frozenset({E,F}),
                frozenset({C,D}),
                frozenset({H,C,G})
            })
        )


    def test_minimal_startnodes_for_node(self):
        spsg=sparse_powerset_graph(self.computers)
        targetVars=frozenset({A,B})
        res=minimal_startnodes_for_node(
                 spsg
                ,targetVars
        )
        #print('###############################')
        #print('miminal startsets for: ', node_2_string(targetVars),nodes_2_string(res))
        #print('###############################')
        self.assertSetEqual(
            res,
            frozenset({
                frozenset({E,F,I}),
                frozenset({C,D,I}),
                frozenset({H,C,G,I})
            })
        )

