#!/usr/bin/env python3
import matplotlib.pyplot as plt
from frozendict import frozendict
from copy import copy,deepcopy
import networkx as nx
from testinfrastructure.InDirTest import InDirTest
from unittest import skip
from copy import copy

from bgc_md2.resolve.graph_helpers import ( 
    # Thomas's functions
    direct_prerequisites_Thomas
    ,graph_Thomas
    #
    #
    # Markus's functions
    ,immutable_edge
    ,sparse_powerset_graph
    ,direct_predecessor_nodes
    ,direct_predecessor_graph
    ,arg_set_graph
    ,minimal_startnodes_for_single_var
    ,minimal_startnodes_for_node
    ,update_step
    #,remove_supersets
    ,node_2_string
    ,nodes_2_string
    ,edge_2_string
    ,product_graph
    ,product_edge_triple
    #,cartesian_union_graph
    #,create_multigraph
    ,draw_multigraph_graphviz
    ,draw_multigraph_matplotlib
    #,draw_multigraph_plotly
    ,draw_SetMultiDiGraph
    ,draw_SetDiGraph
    ,draw_Graph_svg
    ,draw_Graph_with_computers_svg
    #,powerlist
)

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

def a_from_y(y:Y)->A:
    return A()

def b_from_y(y:Y)->B:
    return B()

def a_from_z(z:Z)->A:
    return A()

def b_from_z(z:Z)->B:
    return B()

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
    # we produce a small set of Mvars with a loop (b could be something like a CompartmentalSystem that can be computed # in different ways and can also be queried about its constituents.

    # we should build a graph comparison
    def assertGraphEqual(self,g1,g2):
        self.assertSetEqual(set(g1.nodes),set(g2.nodes))
        self.assertSetEqual(
            set(immutable_edge(e) for e in g1.edges(data=True))
            ,set(immutable_edge(e) for e in g2.edges(data=True))
        )

        
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

    def test_assertGraphEqual(self):
        # this should be a test for a derived testClass
        g1=nx.DiGraph()
        g1.add_nodes_from([A,B,Z])
        g1.add_edge(Z,A,computer=frozenset({a_from_z}))
        g1.add_edge(Z,B,computer=frozenset({b_from_z}))

        g2=nx.DiGraph()
        g2.add_nodes_from([A,B,Z])
        g2.add_edge(Z,A,computer=frozenset({a_from_z}))
        g2.add_edge(Z,B,computer=frozenset({b_from_z}))
        
        self.assertGraphEqual(g1,g2)
        g2_1=deepcopy(g2)
        self.assertGraphEqual(g1,g2_1)
        g2_1.add_node(D)
        with self.assertRaises(Exception):
            self.assertGraphEqual(g1,g2_1)

    def test_arg_set_graph(self):
        asg=arg_set_graph(D,self.computers)
        
        ref=nx.DiGraph()
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


        self.assertGraphEqual(
            asg , ref
        )
        draw_Graph_with_computers_svg(ref,'ref')
        draw_Graph_with_computers_svg(asg,'asg')
    
    
    def test_product_edge_triple(self):
        self.assertEqual(
            product_edge_triple(
                 ({I},{A},{'computers':{a_from_i}})
                ,({B},{D},{'computers':{d_from_b}})
            )
            ,({I,B},{A,D},{a_from_i,d_from_b})
        )
        
        
        # take an example where we could loose a costructor
        # by just writing edges with the same start and destination
        # so that only the edge that is written last remains.
        # computers:a(z),b(z) 
        #
        #    {a(z)}          {b(z)}          {a(z),b(z)}
        # {z}  ->  {a} x {z}  ->   {b}   = {z} -> {a,b}

        self.assertEqual(
            product_edge_triple(
                 ({Z},{A},{'computers':{a_from_z}})
                ,({Z},{B},{'computers':{b_from_z}})
            )
            ,({Z},{A,B},{a_from_z,b_from_z})
        )

    def test_product_graph(self):
        computers=frozenset({a_from_i,d_from_b})
        
        asg_A=arg_set_graph(A,computers)
        asg_D=arg_set_graph(D,computers)
        
        pg_A_D=product_graph(asg_A,asg_D)

        ref_pg_A_D=nx.MultiDiGraph()
        ref_pg_A_D.add_edge(
            frozenset({B,I})
            ,frozenset({D,A})
            ,computers=frozenset({a_from_i,d_from_b})
        )
        self.assertGraphEqual(pg_A_D,ref_pg_A_D)

        # Now take an example where we could loose a costructor
        # by just writing edges with the same start and destination
        # so that only the edge that is written last remains.
        # computers:a(z),b(z) 
        computers=frozenset({a_from_z,b_from_z})
        #
        #    {a(z)}          {b(z)}          {a(z),b(z)}
        # {z}  ->  {a} x {z}  ->   {b}   = {z} -> {a,b}
        asg_A=arg_set_graph(A,computers)
        asg_B=arg_set_graph(B,computers)
        pg_A_B=product_graph(asg_A,asg_B)
        ref_pg_A_B=nx.MultiDiGraph()
        ref_pg_A_B.add_edge(
            frozenset({Z})
            ,frozenset({A,B})
            ,computers=frozenset({a_from_z,b_from_z})
        ) 
        self.assertGraphEqual(pg_A_B,ref_pg_A_B)
        fig1=plt.figure(figsize=(5,20))
        ax1=fig1.add_subplot(411,frame_on=True,title="arg_set_graph(A)")
        ax2=fig1.add_subplot(412,frame_on=True,title="arg_set_graph(B)")
        ax3=fig1.add_subplot(413,frame_on=True,title="product_graph(A,B)")
        #ax4=fig1.add_subplot(414,frame_on=True,title="ref")
        draw_SetDiGraph(     asg_A,ax=ax1)
        draw_SetDiGraph(     asg_B,ax=ax2)
        draw_SetMultiDiGraph(    pg_A_B,ax=ax3)
        #draw_SetMultiDiGraph(ref_pg_A_B,ax=ax4)
        fig1.savefig('AB_Z.pdf')
        
        # Now take another example where we could loose a costructor
        # by just writing edges with the same start and destination
        # so that only the edge that is written last remains.
        # computers: a(y),a(z)  ,b(y),b(z) 
        # the interesting edge is {z,y} -> {a,b}
        # which can be computed by the computer combinations
        # {a(z),b(y)} and {a(y),
        computers=frozenset({a_from_y,a_from_z,b_from_y,b_from_z})

        asg_A=arg_set_graph(A,computers)
        asg_B=arg_set_graph(B,computers)
        pg_A_B=product_graph(asg_A,asg_B)
        ref_pg_A_B=nx.MultiDiGraph()
        
        ref_pg_A_B.add_edge(
            frozenset({Z})
            ,frozenset({A,B})
            ,computers=frozenset({a_from_z,b_from_z})
        ) 
        ref_pg_A_B.add_edge(
            frozenset({Y})
            ,frozenset({A,B})
            ,computers=frozenset({a_from_y,b_from_y})
        ) 
        ref_pg_A_B.add_edge(
            frozenset({Y,Z})
            ,frozenset({A,B})
            ,computers=frozenset({a_from_y,b_from_z})
        ) 
		# same nodes different computers (multigraph)
        ref_pg_A_B.add_edge(
            frozenset({Y,Z})
            ,frozenset({A,B})
            ,computers=frozenset({a_from_z,b_from_y})
        ) 
        self.assertGraphEqual(pg_A_B,ref_pg_A_B)
        fig2=plt.figure(figsize=(5,20))
        ax1=fig2.add_subplot(4,1,1,title="arg_set_graph(A)") 
        ax2=fig2.add_subplot(4,1,2,title="arg_set_graph(B)") 
        ax3=fig2.add_subplot(4,1,3,title="product_graph(A,B)")
        ax4=fig2.add_subplot(4,1,4,title="ref")
        draw_SetDiGraph(     asg_A,ax=ax1)
        draw_SetDiGraph(     asg_B,ax=ax2)
        draw_SetMultiDiGraph(    pg_A_B,ax=ax3)
        draw_SetMultiDiGraph(ref_pg_A_B,ax=ax4)
        fig2.savefig('AB_YZ.pdf')
    
    def test_direct_predecessor_graph(self):
        # for a node N with a singel mvar M={mvar} this is the same graph as 
        # returned by arg_set_graph(mvar)
        dpg=direct_predecessor_graph(
            frozenset({D})
            ,self.computers
        ) 
        #draw_Graph_with_computers_svg(dpg,'dpg')

        #self.assertSetEqual(
        #    set(asg.nodes())
        #    ,{
        #        frozenset({D})
        #        ,frozenset({B})
        #        ,frozenset({G, H})
        #    }
        # )
       
        #def immutable_edge(edge):
        #    s,d,dat=edge
        #    return (s,d,frozendict(dat))

        #edge_set=set(immutable_edge(e) for e in dpg.edges(data=True))
        #print(edge_set)
        #self.assertSetEqual(
        #    edge_set
        #    ,{
        #        
        #        (
        #            frozenset({B})
        #            ,frozenset({D})
        #            ,frozendict({'computers':frozenset({d_from_b})})
        #        )
        #        ,(
        #            frozenset({G, H})
        #            ,frozenset({D})
        #            ,frozendict({'computers':frozenset({d_from_g_h})})
        #        )
        #    }
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
        
