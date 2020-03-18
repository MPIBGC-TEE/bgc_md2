from pygraphviz.agraph import AGraph
from typing import Callable
from functools import reduce
import networkx as nx

from bgc_md2.resolve.graph_helpers import ( 
    update_generator
)
from .non_graph_helpers import  (
    #,
    pretty_name
)

def compset_2_string(compset):
    return '{'+",".join([pretty_name(c) for c in compset])+'}'

def node_2_string(node):
    return '{'+",".join([pretty_name(v) for v in node])+'}'

def nodes_2_string(node):
    return '[ '+",".join([node_2_string(n) for n in node])+' ]'

def edge_2_string(e):
    return "("+node_2_string(e[0])+','+node_2_string(e[1])+')'



def draw_update_sequence(computers,max_it,fig):
    lg=[ g for g in update_generator(computers,max_it=max_it)]
    nr=len(lg)
    fig.set_size_inches(20,20*nr)
    pos = nx.spring_layout(lg[-1] )
    # layout alternatives
    #pos = nx.spring_layout(lg[-1], iterations=20)
    #pos = nx.circular_layout(lg[-1] )
    #pos = nx.kamada_kawai_layout (lg[-1])
    #pos = nx.planar_layout (lg[-1])
    #pos = nx.random_layout (lg[-1])
    #pos = nx.shell_layout (lg[-1])
    #pos = nx.spectral_layout (lg[-1])
    #pos = nx.spiral_layout (lg[-1])
    axs=fig.subplots(nr,1,sharex=True,sharey=True)
    for i in range(nr):
        draw_ComputerSetMultiDiGraph_matplotlib(lg[i],axs[i],pos=pos)
    
def draw_ComputerSetDiGraph_matplotlib(
        spsg:nx.DiGraph,
        ax,
        pos=None,
        **kwargs):
    if pos is None:
        pos=nx.spring_layout(spsg)
        #pos = nx.circular_layout(spsg)
    
    nx.draw(
        spsg
        ,labels={n:node_2_string(n) for n in spsg.nodes()}
        ,ax=ax
        ,node_size=2000
        ,node_shape='s'
        ,pos=pos
        ,**kwargs
    )
    for e in spsg.edges():
        print(spsg.get_edge_data(*e))
        
    edge_labels= {e:compset_2_string(spsg.get_edge_data(*e)['computers']) for e in spsg.edges()}
    nx.draw_networkx_edge_labels(
        spsg 
        ,ax=ax
        ,edge_labels=edge_labels
        ,pos=pos
    )

def draw_ComputerSetMultiDiGraph_matplotlib(spsg,ax,pos=None,**kwargs):
    if pos is None:
        pos=nx.spring_layout(spsg)
        #pos = nx.circular_layout(spsg)
    
    nx.draw(
        spsg
        ,labels={n:node_2_string(n) for n in spsg.nodes()}
        ,ax=ax
        ,node_size=2000
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

def AGraphComputerSetMultiDiGraph(
        spsg:nx.MultiDiGraph
        ,cf:Callable
    )->AGraph:
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
        ss,st=tuple(map(node_2_string,(s,t)))
        A.add_edge(ss,st)
        Ae=A.get_edge(ss,st)
        Ae.attr['label']="\n".join(
                [c.__name__ for c in computer_set]
        ) 
    return A
    
def AGraphComputerMultiDiGraph(
        spsg:nx.MultiDiGraph
        ,cf:Callable
    )->AGraph:
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
        for c in computer_set:
            ss,st=tuple(map(node_2_string,(s,t)))
            A.add_edge(ss,st)
            Ae=A.get_edge(ss,st)
            Ae.attr['color']=cf(c)
            Ae.attr['fontcolor']=cf(c)
            Ae.attr['label']= c.__name__ 
    
    return A
