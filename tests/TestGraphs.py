#!/usr/bin/env python3
import matplotlib.pyplot as plt
from bgc_md2.resolve.graph_helpers import ( 
    direct_prerequisites,
    #,minimal_startnodes_for_single_var
    #,minimal_startnodes_for_node
    #,update_step
    #,remove_supersets
    #,node_2_string
    #,nodes_2_string
    #,edge_2_string
    #,cartesian_union
    #,create_multigraph
    #,draw_multigraph_graphviz
    #,draw_multigraph_matplotlib
    #,draw_multigraph_plotly
    sparse_powerset_graph,
    #,draw_Graph_png
    #,powerlist
)
from copy import copy,deepcopy
import networkx as nx
from unittest import TestCase,skip
from copy import copy



class A:
    """A variable we assume to be given"""


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


def a_from_i(i: I) -> A:
    """Computes a from i"""
    return i**2


def b_from_c_d(c: C, d: D) -> B:
    """Computes b from c and d"""
    return 3


def b_from_e_f(e: E, f: F) -> B:
    """Computes b from e and f"""
    return 3


def c_from_b(b: B) -> C:
    """Computes c from b"""
    return 2*b


def d_from_b(b: B) -> D:
    """Computes d from b"""
    return 5


def d_from_g_h(g: G, h: H) -> D:
    """Computes d from g and h"""
    return 5


def e_from_b(b: B) -> E:
    """Computes e from b"""
    return 4


def f_from_b(b: B) -> F:
    """Computes f from b"""
    return 4


class TestGraphs(TestCase):
    # we produce a small set of Mvars with a loop (b could be something like a CompartmentalSystem that can be computed # in different ways and can also be queried about its constituents.

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
        self.computers = {
            a_from_i,
            b_from_c_d,
            b_from_e_f,
            c_from_b,
            d_from_b,
            d_from_g_h,
            e_from_b,
            f_from_b,
        }
    
    def test_direct_prerequisites(self):
        g = sparse_powerset_graph(self.mvars, self.computers)
        
        self.assertSetEqual(
            direct_prerequisites(g, A), {(frozenset({I}), a_from_i)})
        
        self.assertSetEqual(
            direct_prerequisites(g, B),
            {
                (frozenset({C, D}), b_from_c_d),
                (frozenset({E, F}), b_from_e_f),
            })

        # pns=direct_predecessor_nodes(
        #         frozenset({ self.Mvars['c'],self.Mvars['d']}) ,self.Mvars,self.Computers) # this function should also return the computers it used 
        # self.assertSetEqual(
        #      pns
        #     ,frozenset({
        #          frozenset({self.Mvars['b']})
        #         ,frozenset({self.Mvars['h'],self.Mvars['c'],self.Mvars['g']})
        #         })
        # )

    @skip(" until new Version of MVars and Computers  are derived from files")
    def test_update_step(self):
        new_nodes=frozenset([frozenset({v}) for v in self.Mvars])
        G=nx.DiGraph()
        G.add_nodes_from(new_nodes)
        draw_Graph_png(G,"original_Graph")
        G_new,new_nodes=update_step(G,new_nodes,self.Mvars,self.Computers)
        draw_Graph_png(G_new,"updated_Graph_1")
        
    @skip(" until new Version of MVars and Computers  are derived from files")
    def test_minimal_startnodes_for_single_var(self):
        G=sparse_powerset_graph(self.Mvars,self.Computers)
        res=minimal_startnodes_for_single_var(G,self.Mvars['b'])
        print(nodes_2_string(res))
        
        
    @skip(" until new Version of MVars and Computers  are derived from files")
    def test_minimal_startnodes_for_node(self):
        G=sparse_powerset_graph(self.Mvars,self.Computers)
        targetVars=frozenset({self.Mvars['a'],self.Mvars['b']})
        res=minimal_startnodes_for_node(
                 G
                ,targetVars
        )
        print('###############################')
        print('miminal startsets for: ', node_2_string(targetVars),nodes_2_string(res))
        print('###############################')
        
    @skip(" until new Version of MVars and Computers  are derived from files")
    def test_draw_multigraph_graphviz(self):
        draw_multigraph_graphviz(self.Mvars,self.Computers)
    
    @skip(" until new Version of MVars and Computers  are derived from files")
    def test_draw_multigraph_matplotlib(self):
        draw_multigraph_matplotlib(self.Mvars,self.Computers)
    
    @skip(" until new Version of MVars and Computers  are derived from files")
    def test_draw_multigraph_plotly(self):
        draw_multigraph_plotly(self.Mvars,self.Computers)
    
    def test_creation(self):
        # only for visualization draw the connections of the mvars via computers
        # note that this is not a graph we can query for connectivity
        
        # Now we build the directed Graph we can use to compute connectivity
        # the Nodes are sets of Mvars (elemenst of the powerset of all Mvars)
        # and a connection between two sets indicates computability of the target set from 
        # the source set.
        # The complete graph would contain all elements of the powerset of allMvars and all
        # possible connections, which is prohibitively expensive.
        # Instead we will compute a subgraph where we start with one element sets as targets
        # and infer the predecessors of those sets and then the predecessors of the predecessors and so on until we do not find new nodes=start_sets.
         
        g = sparse_powerset_graph(self.mvars, self.computers)
        # After the graph has been computed we can use it
        # to infer computability of all Mvars

        self.assertEqual({
            frozenset({E}),
            frozenset({D}),
            frozenset({C}),
            frozenset({B}),
            frozenset({A}),
            frozenset({G}),
            frozenset({I}),
            frozenset({F}),
            frozenset({H}),
            frozenset({H, G}),
            frozenset({C, D}),
            frozenset({E, F}),
        }, set(g.nodes()))

        self.assertEqual({
            (frozenset({B}), frozenset({D})),
            (frozenset({B}), frozenset({E})),
            (frozenset({B}), frozenset({F})),
            (frozenset({B}), frozenset({C})),
            (frozenset({I}), frozenset({A})),
            (frozenset({H, G}), frozenset({D})),
            (frozenset({C, D}), frozenset({B})),
            (frozenset({E, F}), frozenset({B})),
        }, g.edges())
