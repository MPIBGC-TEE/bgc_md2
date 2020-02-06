from .inspect_type_hints import input_mvars, output_mvar
from functools import lru_cache,reduce
import networkx as nx


def powerlist(S):
    # We do not want to rely on the set union operation (which necessiates the creation of sets
    # in the first place which is a O(n^2) operation)
    # So we avoid the manyfold occurence of a sublit in the  resultlist 'manually' by not creating
    # it in the first place
    # we start with an empty list
    initial=[[]]
    # and gradually enhance it
    return reduce(lambda acc,el:acc+[ subset +[el] for subset in acc],S,initial)

def sparse_powerset_graph(mvars, computers):
    g = nx.DiGraph()
    g.add_nodes_from(frozenset({v}) for v in mvars)
    for computer in computers:
        _input = frozenset(input_mvars(computer))
        _output = frozenset({output_mvar(computer)})
        g.add_node(_input)
        g.add_edge(_input, _output, computer=computer)
    return g


def draw_multigraph_graphviz(allMvars,allComputers):
    # build initial multigraph
    # for visualization draw the directed Multigraph with the MVars as nodes
    # unfortunately it is not useful for connectivity computations
    # since all 'edges' defined by a computer c are misleading in the sense that 
    # we need the union of all the source variables of c to go to the target Mvar of c 
    # while the arrows suggest ways from any of the arguments...
    # for visualization it would helpful to draw all arcs belonging to the same computer
    # in the same color.
    # Since we do not compute anything from this graph we actually do not need a graphlibrary
    # but can visualize it immediately with graphviz
    # We use a unique color for every computer
    #colordict=CSS4_COLORS
    colordict=TABLEAU_COLORS
    color_names=[n for n in colordict.keys()]
    computer_colors={c.name:color_names[i] for i,c in enumerate(allComputers)}
    A=AGraph(directed=True)
    A.node_attr['style']='filled'
    A.node_attr['shape']='circle'
    A.node_attr['fixedsize']='false'
    A.node_attr['fontcolor']='#FFFFFF'
    A.node_attr['color']='black'
    A.add_nodes_from([v.name for v in allMvars])
    cols=['blue','red','green','orange'] 
    for v in allMvars:
        for c in v.computers(allComputers):
            ans=c.arg_name_set
            edges=[(an,v.name)  for an in ans]
            for e in edges:
                A.add_edge(e)
                Ae=A.get_edge(e[0],e[1])
                Ae.attr['color']=colordict[computer_colors[c.name]]
                Ae.attr['fontcolor']=colordict[computer_colors[c.name]]
                Ae.attr['label']=c.name
    #print(A.string()) # print to screen
    A.draw('Multigraph.svg',prog="circo") # draw using circo
