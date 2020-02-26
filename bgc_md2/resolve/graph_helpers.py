
from functools import lru_cache,reduce
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS,BASE_COLORS,TABLEAU_COLORS
from pygraphviz.agraph import AGraph
from typing import List,Set,Tuple,Callable
from copy import deepcopy
from frozendict import frozendict
from testinfrastructure.helpers import pp,pe

from .non_graph_helpers import  (
    #computable_mvars
	#,directly_computable_mvars
	input_mvars
	#,output_mvar
	,arg_set
	#,arg_set_set
	,all_mvars
	#,applicable_computers
	,all_computers_for_mvar
    ,pretty_name
)
def compset_2_string(compset):
    return '{'+",".join([pretty_name(c) for c in compset])+'}'

def node_2_string(node):
    return '{'+",".join([pretty_name(v) for v in node])+'}'

def nodes_2_string(node):
    return '[ '+",".join([node_2_string(n) for n in node])+' ]'

def edge_2_string(e):
    return "("+node_2_string(e[0])+','+node_2_string(e[1])+')'

def immutable_edge(edge):
    s,d,dat=edge
    return (s,d,frozendict(dat))

def graphs_equal(
        g1_multi:nx.MultiDiGraph
        ,g2_multi:nx.MultiDiGraph
    )->bool:
    # since we deal with a multigraph
    # it is possible to get several edges between two nodes.
    # The method get_edge_data returns a dictionary
    # with the numbers of these edges as keys.
    # But we want to consider two graphs equivalent if the resulting
    # SET is equal, in other words: 
    # If a graph has two edges EACH from : A->B  
    # we do not care which of the edges has which computerset
    # Therefore we compare the set of computersets belonging


    g1_single=toDiGraph(g1_multi) 
    g2_single=toDiGraph(g2_multi)
    return all(
            [g1_single.get_edge_data(*e)
                ==g2_single.get_edge_data(*e) 
                for e in g1_single.edges()
            ]
            + [g1_single.get_edge_data(*e)
                ==g2_single.get_edge_data(*e) 
                for e in g2_single.edges()]
            ) & (g1_single.nodes()==g2_single.nodes())

@lru_cache(maxsize=None) 
def arg_set_graph( 
        mvar :type,
        allComputers:Set[Callable]
    )->nx.MultiGraph:
    # return the subgraph of arg_name_sets for all computers that
    # return this mvar
    target=frozenset({mvar})
    g = nx.MultiDiGraph()
    g.add_node(target)
    for c in all_computers_for_mvar(mvar,allComputers):
        g.add_edge(arg_set(c), target,computers=frozenset({c}))
    return g

def product_graph(
        *graphs:Tuple[nx.MultiDiGraph]
    )->nx.MultiDiGraph:
    return reduce(lambda u,v:product_graph_2(u,v),graphs)

def product_graph_2(
        g1:nx.MultiDiGraph
        ,g2:nx.MultiDiGraph
    )->nx.MultiDiGraph:
    cp=nx.cartesian_product(g1,g2)
    prod=nx.MultiDiGraph()


    for edge in cp.edges(data=True):
        s_tup,d_tup,data=edge
        prod.add_edge(
            reduce(lambda acc,n :acc.union(n), s_tup)
            ,reduce(lambda acc,n :acc.union(n), d_tup)
            ,computers=data['computers'])
        
    return prod

def initial_sparse_powerset_graph(computers:Set[Callable])->nx.MultiDiGraph:
    spsg=nx.MultiDiGraph()
    allMvars = all_mvars(computers) 
    for v in allMvars:
        spsg.add_edges_from(arg_set_graph(v,computers).edges(data=True))
    return spsg

def update_step(
        spsg:nx.MultiDiGraph,
        computers:Set[Callable]
    )->nx.MultiDiGraph:
    new=deepcopy(spsg)
    start_nodes= [n for n in new.nodes() if len(new.in_edges(n))==0]
    #print(nodes_2_string(start_nodes)) 
    for node in start_nodes:
        #print(tuple(v for v in node))
        pg=product_graph(*[arg_set_graph(v,computers) for v in node])
        #new.add_edges_from(pg.edges(data=True))
        new=nx.compose(new,pg)
    
    return new

def sparse_powerset_graph(
        computers:Set[Callable]
    )->nx.MultiDiGraph:
    old=initial_sparse_powerset_graph(computers)
    new=update_step(old,computers)
    while not(graphs_equal(old,new)):
        new=update_step(old,computers)
        old=new
        print(graphs_equal(old,new))
    return new

def update_generator(
        computers:Set[Callable]
        ,max_it:int
    )->List[nx.MultiDiGraph]:

    if max_it<0:
        raise IndexError("update sequence indices have to be larger than 0")
    
    val=initial_sparse_powerset_graph(computers)
    yield val
    old=deepcopy(val)
    val=update_step(val,computers)

    print(graphs_equal(old,val))

    counter=1
    while max_it>counter and not(graphs_equal(old,val)):
        yield val
        old=deepcopy(val)
        val=update_step(val,computers)
        counter +=1
        print(counter,graphs_equal(old,val))



def draw_update_sequence(computers,max_it,file_name=None):
    lg=[ g for g in update_generator(computers,max_it=max_it)]
    nr=len(lg)
    print("nr=",nr)
    pos = nx.spring_layout(lg[-1], iterations=20)
    #pos = nx.circular_layout(lg[-1] )
    #pos = nx.kamada_kawai_layout (lg[-1])
    #pos = nx.planar_layout (lg[-1])
    #pos = nx.random_layout (lg[-1])
    #pos = nx.shell_layout (lg[-1])
    #pos = nx.spectral_layout (lg[-1])
    #pos = nx.spiral_layout (lg[-1])
    fig=plt.figure(figsize=(20,20*nr))   
    axs=fig.subplots(nr,1,sharex=True,sharey=True)
    for i in range(nr):
        draw_SetMultiDiGraph(lg[i],axs[i],pos=pos)
    
    fig.savefig(file_name)

def toDiGraph(
        g_multi:nx.MultiDiGraph
    )->nx.DiGraph:
    def edgeDict_to_set(ed):
        target='computers'
        comp_set_set=frozenset(
            [v[target] for v in ed.values() if target in v.keys()]
        )
        return comp_set_set

    g_single=nx.DiGraph()
    for e in g_multi.edges():
        s,t=e
        edgeDict=g_multi.get_edge_data(s,t)
        comp_set_set=edgeDict_to_set(edgeDict)
        if (g_single.has_edge(s,t)):
            comp_set_set=comp_set_set.union( 
                g_single.get_edge_data(s,t)['computers']
            )
        g_single.add_edge(s,t,computers=comp_set_set)
    return g_single



def minimal_startnodes_for_single_var(
         spg:nx.Graph
        ,targetVar:type
    ):
    ''' spg is a sparse powerset Graph, which means that it only contains all one element sets as targets.'''
    # We first create a graph with the direction of the edges reversed
    targetSet=frozenset({targetVar})
    res=nx.shortest_path(spg,target=targetSet) 
    possible_startnodes=frozenset(res.keys())
    minimal_startnodes=[n for n in filter(lambda n: not(n.issuperset(targetSet)),possible_startnodes)]
    return frozenset(minimal_startnodes)
    
def target_subgraph(
         spg:nx.Graph
        ,targetNode:Set[type]
    ):
    paths=nx.shortest_path(spg,target=targetNode).values()
    connected_nodes=reduce(lambda acc, p : acc.union(p),paths,frozenset({}))
    return spg.subgraph(connected_nodes)



def minimal_target_subgraph_for_single_var(
         spg:nx.Graph
        ,targetVar:type
    ):
    # all minimal starting points
    start_nodes=minimal_startnodes_for_single_var(spg,targetVar)
    path_list_list=[
        nx.shortest_path(spg,target=targetVar, source=sn)
        for sn in start_nodes
    ]
    connected_nodes=reduce(lambda acc, pl : acc.union(pl),path_list_list,frozenset({}))
    return spg.subgraph(connected_nodes)


def minimal_startnodes_for_node(
         spg:nx.Graph
        ,targetNode:Set[type]
    )->Set[Set]:
    if spg.has_node(targetNode):
        # The targetNode is already part of the spg. 
        # (because it has been added in one of the update 
        # steps as an argument set of a computer) 
        # in which case we can simply return the startnodes
        # of the paths leading to it.
        path_dict=nx.shortest_path(spg,target=targetNode)
        possible_startnodes=frozenset(path_dict.keys())
    else:
        # Although we do not find the node itself
        # we can find the nodes for the single element sets 
        # of the mvars in the node, since we have built the graph
        # starting wiht them. E.g if node {A,B,C} is not part of 
        # the graph we know at least that the nodes {A}, {B} amd {C}
        # are in the graph. 
        # For each of them we can compute the subgraph of 
        # spg that leads to it. If we compute the product of these
        # subgraphs it will contain the desired node.

        # fixme: mm 02-26-2020
        # We could make this more efficient by looking for all the
        # disjoint unions of the targetNode Mvars and compute
        # the product of the graphs leading to the subsets.
        # If the subsets are not one element sets we need fewer
        # multiplications.
        prod_g=product_graph(*[target_subgraph(spg,frozenset({v})) for v in targetNode])
        prod_path_dict=nx.shortest_path(prod_g,target=targetNode)
        possible_startnodes=frozenset(prod_path_dict.keys())
        
    
    def filter_func(n):
        # remove every set that contains one of the variables we are looking for ...
        return not(any([ (v in n) for v in targetNode]))

    minimal_startnodes=frozenset([n for n in filter(filter_func,possible_startnodes)])
    return minimal_startnodes




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


def draw_SetMultiDiGraph(spsg,ax,pos=None,**kwargs):
    if pos is None:
        pos=nx.spring_layout(spsg)
    
    nx.draw(
        spsg
        ,labels={n:node_2_string(n) for n in spsg.nodes()}
        ,ax=ax
        ,node_size=1000
        ,node_shape='s'
        ,pos=pos
        ,**kwargs
    )
    # at the moment it is not possible to draw
    # more than one edge (egde_lables) between nodes
    # directly (no edgelabels for MultiDiGraphs)
    # therefore we draw only one line for all computersets
    # and assemble the label from the different edges 
    def edgeDict_to_string(ed):
        target='computers'
        comp_sets=[v[target] for v in ed.values() if target in v.keys()]
        comp_set_strings=[compset_2_string(cs) for cs in comp_sets]
        res="\n".join(comp_set_strings)
        #print(res)
        return res
    
    edge_labels={e:edgeDict_to_string(spsg.get_edge_data(*e)) for e in spsg.edges()}  

    nx.draw_networkx_edge_labels(
        spsg 
        ,ax=ax
        ,edge_labels=edge_labels
        ,pos=pos
    )

def draw_Graph_with_computers_svg(spsg,file_name_trunk):
    # the next line is the standard translation 
    # We could do this using the attributes of the edges to
    # make much niceer pictures representing the different computers in different
    # colors or annotate them....
    allComputers=reduce(
        lambda acc,edge:acc.union(edge[2]['computers'],acc)
        ,spsg.edges(data=True)
        ,set()
    )
    colordict=TABLEAU_COLORS
    color_names=[n for n in colordict.keys()]
    computer_colors={c.__name__:color_names[i] for i,c in enumerate(allComputers)}
    A=nx.nx_agraph.to_agraph(spsg)
    A=AGraph(directed=True)
    A.node_attr['style']='filled'
    A.node_attr['shape']='rectangle'
    A.node_attr['fixedsize']='false'
    A.node_attr['fontcolor']='black'
    
    for node in spsg.nodes:
        A.add_node(node_2_string(node))
    edges=spsg.edges(data=True)
    for edge in edges:
        s,t,data_dict=edge
        computer_set=data_dict['computers']
        print('computer_set')
        print(computer_set)
        ss,st=tuple(map(node_2_string,(s,t)))
        A.add_edge(ss,st)
        Ae=A.get_edge(ss,st)
        #Ae.attr['color']=colordict[computer_colors[c.__name__]]
        #Ae.attr['fontcolor']=colordict[computer_colors[c.__name__]]
        Ae.attr['label']="\n".join([c.__name__ for c in computer_set]) 
    #print(A.string()) # print to screen
    #A.draw(file_name_trunk+'.png',prog="circo") # draw to png using circo
    A.draw(file_name_trunk+'.svg',prog='circo') 

