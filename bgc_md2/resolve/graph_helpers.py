from functools import lru_cache,reduce
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS,BASE_COLORS,TABLEAU_COLORS
from pygraphviz.agraph import AGraph
from typing import List,Set,Tuple,Callable
from copy import deepcopy

from .non_graph_helpers import  (
    computable_mvars
	,directly_computable_mvars
	,input_mvars
	,output_mvar
	,arg_set
	,arg_set_set
	,all_mvars
	,applicable_computers
	,all_computers_for_mvar
    ,pretty_name
)
def powerlist(S):
    # We do not want to rely on the set union operation (which necessiates the creation of sets
    # in the first place which is a O(n^2) operation)
    # So we avoid the manyfold occurence of a sublit in the  resultlist 'manually' by not creating
    # it in the first place
    # we start with an empty list
    initial=[[]]
    # and gradually enhance it
    return reduce(lambda acc,el:acc+[ subset +[el] for subset in acc],S,initial)

def node_2_string(node):
    return '{'+",".join([pretty_name(v) for v in node])+'}'

def nodes_2_string(node):
    return '[ '+",".join([node_2_string(n) for n in node])+' ]'

def edge_2_string(e):
    return "("+node_2_string(e[0])+','+node_2_string(e[1])+')'

def cartesian_product(l:List[Set])->Set[Tuple]:
    left_tupels=frozenset([tuple(el) for el in l[0]])
    if len(l)==1:
        return left_tupels
    else:
        right_tupels=cartesian_product(l[1:])
        return frozenset([lt+rt for lt in left_tupels for  rt in right_tupels ])

def cartesian_union(l:List[Set])->Set[Set]:
    #pe('l',locals())
    return frozenset([frozenset(t) for t in cartesian_product(l)])

def remove_supersets_once(sets):
    key_func=lambda s:len(s)
    sets=sorted(sets,key=key_func)
    #print('Startnodes:')
    #print([node_2_string(val)  for val in sets])
    #print('##############################:')

    minimal_sets=[]
    for n in sets:
        if not(any([m.issubset(n) for m in minimal_sets])):
            minimal_sets.append(n)

    return frozenset(minimal_sets)

def remove_supersets(sets):
    new_nodes=remove_supersets_once(sets)
    
    if new_nodes==sets:
        return(new_nodes)
    else:
        return remove_supersets(new_nodes)
    
def direct_predecessor_nodes(
        node        :Set[type],
        allComputers:Set[Callable]
    )->Set[Set[type]]:
    # assume that we want to compute a set of MVars (a node in our graph) from other sets of Mvars
    # let s_a be the set of nodes from which we can reach the set {a} (where a is a MVar} 
    # and s_b the set of nodes from which we can reach the node {b} (where b is an Mvar
    # to reach the node set {a,b} we can combine any startnode from s_a with any startnode from s_b
    # in fact we can reach {a,b} from all nodes in the set {s: s=n_a v n_b for na in sa v {a}  for nb in sb v {b} }
    # we call this set the 'cartesian_union(A,B) with  A=s_a v {a}  and B=s_b v {b} 
    # This can be generalized to an arbitrary list of sets. We build the cartesian product of the sets A,B,... and then
    # transform every tupel of the product to a set (thereby removing the duplicates and order)
    res=cartesian_union(
        [ {frozenset({v})}.union(arg_set_set(v,allComputers)) for v in node]
    )
    #pe('node',locals())

    # note that the cartesian product contains the original node
    # we remove all nodes that are just supersets of it
    # and afterwards the node itself
    #return res.difference(frozenset({node}))
    return remove_supersets(res).difference(frozenset({node}))

# fixme: Markus 2-14 2020
# I do not know what we can use this graph for
# I renamed it since it is not the "sparse_powerset_graph" (my = Markus's) tests. 
def graph_Thomas(mvars, computers):
    g = nx.DiGraph()
    g.add_nodes_from(frozenset({v}) for v in mvars)
    for computer in computers:
        _input = frozenset(input_mvars(computer))
        _output = frozenset({output_mvar(computer)})
        g.add_node(_input)
        g.add_edge(_input, _output, computer=computer)
    return g


# fixme: Markus 2-14 2020
# added the suffix _Thomas to distinguish this function from the rest, since its
# functionality is not used at the moment.
# The information it provides is presently retrieved from the computers directly
# by the function 'arg_set_set'
# The function could possibly be used to reimplement parts of the creation of the 
# original (Markus's) "sparse_powerset_graph"
def direct_prerequisites_Thomas(graph, mvar):
    node = frozenset({mvar})
    return set(
        (pre, data['computer'])
        for pre, _, data in graph.in_edges(node, data=True))


def draw_multigraph_graphviz(allMvars,allComputers):
    # build initial multigraph
    # for visualization draw the directed Multigraph with the MVars as nodes
    # unfortunately it is not useful for connectivity computations
    # since all 'edges' defined by a computer c are misleading in the sense that 
    # we need the union of all the source variables of c to go to the target Mvar of c 
    # while the arrows suggest ways from any of the arguments...
    # for visualization it would helpful to draw all arcs belonging to the same computer
    # in the same color.
    # Since we do not compute anything from this graph we actually do not need a 
    # graph library
    # but can visualize it immediately with graphviz
    # We use a uniqud .draw('Multigraph.svg',prog="circo") # draw using circo
    colordict=TABLEAU_COLORS
    color_names=[n for n in colordict.keys()]
    computer_colors={c.__name__:color_names[i] for i,c in enumerate(allComputers)}
    A=AGraph(directed=True)
    A.node_attr['style']='filled'
    A.node_attr['shape']='circle'
    A.node_attr['fixedsize']='false'
    A.node_attr['fontcolor']='#FFFFFF'
    A.node_attr['color']='black'
    A.add_nodes_from([pretty_name(v) for v in allMvars])
    cols=['blue','red','green','orange'] 
    for v in allMvars:
        for c in all_computers_for_mvar(v,allComputers):
            ans=input_mvars(c)
            edges=[(pretty_name(an),pretty_name(v))  for an in ans]
            for e in edges:
                A.add_edge(e)
                Ae=A.get_edge(e[0],e[1])
                Ae.attr['color']=colordict[computer_colors[c.__name__]]
                Ae.attr['fontcolor']=colordict[computer_colors[c.__name__]]
                Ae.attr['label']=c.__name__
    #print(A.string()) # print to screen
    A.draw('Multigraph.svg',prog="circo") # draw using circo




def draw_multigraph_matplotlib(allMvars,allComputers):
    # only for visualization draw the connections of the mvars via computers
    # note that this is not a graph we can query for connectivity
    colordict=TABLEAU_COLORS
    color_names=[n for n in colordict.keys()]
    computer_colors={c.__name__:color_names[i] for i,c in enumerate(allComputers)}
    G=create_multigraph(allMvars,allComputers)
    # possibly create new Graph with text nodes
    g=G
    # create a colorlist for the edges using the computer attribute 
    edgelist=[e for e in g.edges]
    computers=[g[e[0]][e[1]]['computer'] for e in edgelist]
    edge_color=[computer_colors[c.__name__] for c in computers]
    fig=plt.figure(figsize=(5,5))
    axes=fig.add_subplot(1,1,1)
    nx.draw_networkx(
            g
            ,edgelist=edgelist
            ,edge_color=edge_color
            ,ax=axes)
    fig.savefig("Multigraph_matplotlib.pdf")

def create_multigraph(allMvars,allComputers):
    G=nx.DiGraph()
    for v in allMvars:
        for c in all_computers_for_mvar(v,allComputers):
            ans=input_mvars(c)
            for an in ans:
                G.add_edge(
                     pretty_name(an)
                    ,pretty_name(v)
                    ,computer=c
                )
    return G



def update_step(G,extendable_nodes,allComputers):
    
    # update the Graph by looking at the new_nodes and adding all their possible predecessors as new nodes 
    # The nodes representing one-element sets have been already treated by the first step 
    # (Thein predecessors found by their computers)
    # now we have to find the predecessors of the sets with more than one element
    # we do this by looking at the tensor product of the predecessor-sets  of the elements
    G=deepcopy(G)
    present_nodes=frozenset(G.nodes)
    present_edges=frozenset(G.edges)
    new_global_nodes=frozenset({})
    for n in extendable_nodes:
        #print('\nextendable node:'+node_2_string(n))
        # it should actually be possible to infer the nodes that can 
        # be computed from other nodes from the graph G alone
        # Up to now we only use the computability of mvars
        # Actually the graph  can also provide information
        # about the computablitiy of sets (by using the union of the sources of single vars)
        pns=direct_predecessor_nodes(n,allComputers) # this function should also return the computers it used 
        for pn in pns:
            G.add_node(pn) 
            e=(pn,n)
            #if not(e in present_edges):
            if not(e in G.edges()):
                #print('new_edge:'+edge_2_string(e))
                G.add_edge(pn,n) 
        #print('direct_predecessor_nodes:'+nodes_2_string(pns))
        arn=present_nodes.intersection(pns)
        #print('allready known:'+nodes_2_string(arn))
        new_local_nodes= pns.difference(arn)
        new_global_nodes=new_global_nodes.union(new_local_nodes)
        #print('new_local_node:'+nodes_2_string(new_global_nodes))

    #print('new_global_nodes:'+nodes_2_string(new_global_nodes))
    #extendable_nodes=new_nodes
    #new_minimal_nodes=new_minimal_nodes#.difference(present_nodes)
    return (G,new_global_nodes)


def updated_Graph(G,extendable_nodes,allMvars,allComputers,counter=0):
    G_new,extendable_nodes_new=update_step(G,extendable_nodes,allComputers)
    if len(extendable_nodes)==0: 
        return (G,extendable_nodes)
    else:
        #draw_Graph_svg(G,"updated_Graph"+str(counter))
        return updated_Graph(G_new,extendable_nodes_new,allMvars,allComputers,counter+1)


def sparse_powerset_graph(allMvars,allComputers):
    new_nodes=frozenset([frozenset({v}) for v in allMvars])
    G=nx.DiGraph()
    G.add_nodes_from(new_nodes)
    G_final,new_nodes=updated_Graph(G,new_nodes,allMvars,allComputers)
    # new_nodes is now empty 
    return G_final


def minimal_startnodes_for_single_var(
         spg:nx.Graph
        ,targetVar:type
    ):
    ''' spg is a sparse powerset Graph, which means that it only contains all one element sets as targets.'''
    # We first create a graph with the direction of the edges reversed
    targetSet=frozenset({targetVar})
    GR=spg.reverse()
    res=[p for p in nx.all_pairs_shortest_path(GR) if p[0]==targetSet]
    possible_startnodes=frozenset([n for n in res[0][1].keys()])
    print("possible_startnodes for",node_2_string(targetSet)," including supersets:",[node_2_string(n) for n in possible_startnodes])
    minimal_startnodes=remove_supersets(possible_startnodes)
    minimal_startnodes=[n for n in filter(lambda n: not(n.issuperset(targetSet)),minimal_startnodes)]
    return frozenset(minimal_startnodes)
    

def minimal_startnodes_for_node(
         spg:nx.Graph
        ,targetNode:Set[type]
    ):
    single_sets=[
            minimal_startnodes_for_single_var(spg,var) for var in targetNode
    ]
    possible_startnodes=cartesian_union(single_sets)
    minimal_startnodes=remove_supersets(possible_startnodes)
    return frozenset(minimal_startnodes)
    
    def filter_func(n):
        # remove every set that contains one of the variables we are looking for ...
        return not(any([ (v in n) for v in targetNode]))

    minimal_startnodes=[n for n in filter(filter_func,minimal_startnodes)]




def draw_Graph_svg(G,file_name_trunk):
    # the next line is the standard translation 
    # We could do this using the attributes of the edges to
    # make much niceer pictures representing the different computers in different
    # colors or annotate them....
    A=nx.nx_agraph.to_agraph(G)
    A=AGraph(directed=True)
    A.node_attr['style']='filled'
    A.node_attr['shape']='rectangle'
    A.node_attr['fixedsize']='false'
    A.node_attr['fontcolor']='black'
    
    for node in G.nodes:
        A.add_node(node_2_string(node))
    for edge in G.edges:
        A.add_edge(node_2_string(edge[0]),node_2_string(edge[1]))
    #print(A.string()) # print to screen
    #A.draw(file_name_trunk+'.png',prog="circo") # draw to png using circo
    A.draw(file_name_trunk+'.svg',prog='circo') 

