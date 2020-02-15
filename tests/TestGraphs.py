#!/usr/bin/env python3
import matplotlib.pyplot as plt
from bgc_md2.resolve.graph_helpers import ( 
    # Thomas's functions
    direct_prerequisites_Thomas
    ,graph_Thomas
    #
    #
    # Markus's functions
    ,sparse_powerset_graph
    ,direct_predecessor_nodes
    ,minimal_startnodes_for_single_var
    ,minimal_startnodes_for_node
    ,update_step
    #,remove_supersets
    ,node_2_string
    ,nodes_2_string
    ,edge_2_string
    #,cartesian_union
    #,create_multigraph
    ,draw_multigraph_graphviz
    ,draw_multigraph_matplotlib
    #,draw_multigraph_plotly
    ,draw_Graph_svg
    #,powerlist
)
from copy import copy,deepcopy
import networkx as nx
from unittest import TestCase,skip
from copy import copy

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


def a_from_i(i: I) -> A:
    """Computes a from i"""
    return A()


def b_from_c_d(c: C, d: D) -> B:
    """Computes b from c and d"""
    return B()


def b_from_e_f(e: E, f: F) -> B:
    """Computes b from e and f"""
    return B()


def c_from_b(b: B) -> C:
    """Computes c from b"""
    return C()


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
        self.computers = frozenset({
            a_from_i,
            b_from_c_d,
            b_from_e_f,
            c_from_b,
            d_from_b,
            d_from_g_h,
            e_from_b,
            f_from_b,
        })
    
        # pns=direct_predecessor_nodes(
        #         frozenset({ self.mvars['c'],self.mvars['d']}) ,self.mvars,self.computers) # this function should also return the computers it used 
        # self.assertSetEqual(
        #      pns
        #     ,frozenset({
        #          frozenset({self.mvars['b']})
        #         ,frozenset({self.mvars['h'],self.mvars['c'],self.mvars['g']})
        #         })
        # )
    

    def test_direct_predecessor_nodes(self):
        self.assertSetEqual(
            direct_predecessor_nodes(
                frozenset({A})
                ,self.computers
            ) 
            ,frozenset({
                 frozenset({I})
            })
        )

        self.assertSetEqual(
             direct_predecessor_nodes(
                 frozenset({B})
                 ,self.computers
            )
            ,frozenset({
                 frozenset({C,D})
                ,frozenset({F,E})
            })
        )

        self.assertSetEqual(
            direct_predecessor_nodes(
                frozenset({C,D})
                ,self.computers
            ) 
            ,frozenset({
                 frozenset({B})
                ,frozenset({H,C,G})
            })
        )

    def test_update_step(self):
        new_nodes=frozenset([frozenset({v}) for v in self.mvars])
        G=nx.DiGraph()
        G.add_nodes_from(new_nodes)
        draw_Graph_svg(G,"test_update_step_original_Graph")
        G,new_nodes=update_step(G,new_nodes,self.computers)
        draw_Graph_svg(G,"test_update_step_updated_Graph_1")
        G,new_nodes=update_step(G,new_nodes,self.computers)
        draw_Graph_svg(G,"test_update_step_updated_Graph_2")
        G,new_nodes=update_step(G,new_nodes,self.computers)
        draw_Graph_svg(G,"test_update_step_updated_Graph_3")
        G,new_nodes=update_step(G,new_nodes,self.computers)
        draw_Graph_svg(G,"test_update_step_updated_Graph_4")
        G,new_nodes=update_step(G,new_nodes,self.computers)
        draw_Graph_svg(G,"test_update_step_updated_Graph_5")
        G,new_nodes=update_step(G,new_nodes,self.computers)
        draw_Graph_svg(G,"test_update_step_updated_Graph_6")
        
    def test_draw_multigraph_graphviz(self):
        draw_multigraph_graphviz(self.mvars,self.computers)
    
    def test_draw_multigraph_matplotlib(self):
        draw_multigraph_matplotlib(self.mvars,self.computers)
    
    @skip("very immature and nearly manual, but maybe neccessary to make the connections clickable") 
    def test_draw_multigraph_plotly(self):
        draw_multigraph_plotly(self.mvars,self.computers)
        
   


    def test_Thomas_graph_creation(self):
        g = graph_Thomas(self.mvars, self.computers)

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
    

    def test_direct_prerequisites_Thomas(self):
        g = graph_Thomas(self.mvars, self.computers)
        
        self.assertSetEqual(
            direct_prerequisites_Thomas(g, A), {(frozenset({I}), a_from_i)}
        )
        
        self.assertSetEqual(
            direct_prerequisites_Thomas(g, B),
            {
                (frozenset({C, D}), b_from_c_d),
                (frozenset({E, F}), b_from_e_f)
            }
        )

    
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
            frozenset({A1,A2,A3}),
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
                B_minus_2
                ,B_minus_1
                ,B0
                ,B1
                ,B2
                ,A_minus_2
                ,A_minus_1
                ,A0
                ,A1
                ,A2
                })
            ,frozenset({
                 b_minus_1_from_b_minus_2
                ,b0_from_b_minus_1
                ,b1_from_b0
                ,b2_from_b1
                ,b3_from_b2
                #
                ,a0_from_b0
                #
                ,a_minus_1_from_a_minus_2
                ,a0_from_a_minus_1
                ,a1_from_a0
                ,a2_from_a1
                ,a3_from_a2
             })
        )
        self.assertSetEqual(
            set(spsg.nodes()),
            {
                frozenset({B_minus_1})
                ,frozenset({B_minus_2})
                ,frozenset({B0})
                ,frozenset({B1})
                ,frozenset({B2})
                ,frozenset({A_minus_1})
                ,frozenset({A_minus_2})
                ,frozenset({A0})
                ,frozenset({A1})
                ,frozenset({A2})
            }
        )
        self.assertSetEqual(
            set(spsg.edges()),
            {
                (frozenset({B_minus_2}),frozenset({B_minus_1}))
                ,(frozenset({B_minus_1}),frozenset({B0}))
                ,(frozenset({B0}),frozenset({B1}))
                ,(frozenset({B1}),frozenset({B2}))
                #
                ,(frozenset({B0}),frozenset({A0})) 
                #
                ,(frozenset({A_minus_2}),frozenset({A_minus_1}))
                ,(frozenset({A_minus_1}),frozenset({A0}))
                ,(frozenset({A0}),frozenset({A1}))
                ,(frozenset({A1}),frozenset({A2}))
            }
        )
    def test_minimal_startnodes_for_single_var(self):
        spsg=sparse_powerset_graph(self.mvars,self.computers)
        # After the graph has been computed we can use it
        # to infer computability of all Mvars
        res=minimal_startnodes_for_single_var(spsg,B)
        #print(nodes_2_string(res))
        self.assertSetEqual(
            res,
            frozenset({ 
                frozenset({E,F}),
                frozenset({C,D}),
                frozenset({H,C,G})
            })
        )
        
        
    def test_minimal_startnodes_for_node(self):
        spsg=sparse_powerset_graph(self.mvars,self.computers)
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
        
