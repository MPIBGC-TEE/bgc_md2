#!/usr/bin/env python3
import matplotlib.pyplot as plt
from bgc_md2.resolve.graph_helpers import ( 
    # direct_predecessor_nodes
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
    sparse_powerset_Graph
    #,draw_Graph_png
    #,powerlist
)
#from bgc_md.resolve.MVar import MVar
#from bgc_md.resolve.Computer import Computer
#from bgc_md.resolve.functions import srm_from_B_u_tens
#from bgc_md.resolve.IndexedSet import IndexedSet
from copy import copy,deepcopy
import networkx as nx
from unittest import TestCase,skip
from copy import copy




class TestGraphs(TestCase):
    # we produce a small set of Mvars with a loop (b could be something like a CompartmentalSystem that can be computed # in different ways and can also be queried about its constituents.

    def setUp(self):
        # fixme mm 02-03-2020: 
        # assemble this set from a file of class definitions
        self.Mvars=IndexedSet({
            MVar('a',description= """ a varible we assume to be given """)
            ,MVar('b')
            ,MVar('c')
            ,MVar('d')
            ,MVar('e')
            ,MVar('f') 
            ,MVar('g')
            ,MVar('h')
            ,MVar('i')
        })
        # fixme mm 02-03-2020: 
        # assemble this set from a file of annotated functions in 
        # special file
        self.Computers=IndexedSet({
            Computer(
                'a(i)'
                ,func=lambda e: e**2
                ,description="""computes f from e"""
            )
            ,Computer(
                'b(c,d)'
                ,func=lambda c,d: 3
                ,description="""computes b from c and d"""
            )
            ,Computer(
                'b(e,f)'
                ,func=lambda c,d: 3
                ,description="""computes b from e and f"""
            )
            ,Computer(
                'c(b)'
                ,func=lambda b:2*b
                ,description="""computes c from b"""
            )
            ,Computer(
                'd(b)'
                ,func=lambda b:5 # we make it consistent but atually the result is not important for the test
                ,description="""computes d from b """
            )
            ,Computer(
                'd(g,h)'
                ,func=lambda g,h:5
                ,description="""computes d from g  and h """
            )
            ,Computer(
                'e(b)'
                ,func=lambda b: 4
                ,description="""computes e from b"""
            )
            ,Computer(
                'f(b)'
                ,func=lambda b: 4
                ,description="""computes f from b"""
            )
        })
    
    @skip(" until new Version of MVars and Computers  are derived from files")
    def test_direct_predecessor_nodes(self):
        pns=direct_predecessor_nodes(frozenset({self.Mvars['a']}),self.Mvars,self.Computers) # this function should also return the computers it used 
        self.assertSetEqual(
             pns
            ,frozenset({
                 frozenset({self.Mvars['i']})
                })
        )
        #print(nodes_2_string(pns))
        
        pns=direct_predecessor_nodes(frozenset({self.Mvars['b']}),self.Mvars,self.Computers) # this function should also return the computers it used 
        self.assertSetEqual(
             pns
            ,frozenset({
                 frozenset({self.Mvars['c'],self.Mvars['d']})
                ,frozenset({self.Mvars['f'],self.Mvars['e']})
                })
        )
        pns=direct_predecessor_nodes(
                frozenset({ self.Mvars['c'],self.Mvars['d']}) ,self.Mvars,self.Computers) # this function should also return the computers it used 
        self.assertSetEqual(
             pns
            ,frozenset({
                 frozenset({self.Mvars['b']})
                ,frozenset({self.Mvars['h'],self.Mvars['c'],self.Mvars['g']})
                })
        )

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
        G=sparse_powerset_Graph(self.Mvars,self.Computers)
        res=minimal_startnodes_for_single_var(G,self.Mvars['b'])
        print(nodes_2_string(res))
        
        
    @skip(" until new Version of MVars and Computers  are derived from files")
    def test_minimal_startnodes_for_node(self):
        G=sparse_powerset_Graph(self.Mvars,self.Computers)
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
         
        G=sparse_powerset_Graph(self.Mvars,self.Computers)
        # After the graph has been computed we can use it
        # to infer computability of all Mvars
