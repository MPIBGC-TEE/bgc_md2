# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import ipycytoscape
import ipywidgets as widgets
import networkx as nx

G = nx.complete_graph(5)
directed = ipycytoscape.CytoscapeWidget()
directed.graph.add_graph_from_networkx(G, directed=True)
directed


#
# # Custom networkx Node Objects
#
# The most common choices for Nodes in networkx are numbers or strings as shown above. A node can also be any hashable object (except None) which work as well.
#

# +
class Node:
    def __init__(self, name):
        self.name = name
        
    def __str__(self):
        return "Node: " + str(self.name)

n1 = Node("node 1")
n2 = Node("node 2")
        
G = nx.Graph()

G.add_node(n1)
G.add_node(n2)

G.add_edge(n1, n2)

w = ipycytoscape.CytoscapeWidget()
w.graph.add_graph_from_networkx(G)
w


# -

#
# # Custom networkx Node Objects that inherit from ipycytoscape.Node
#
# While custom networkx Node objects work, they do not allow as much control over formatting as you may need. The easiest way to achieve customization with custom Node objects is to subclass ipycytoscape.Node as show below.
#

# +
class CustomNode(ipycytoscape.Node):
    def __init__(self, name, classes=''):
        super().__init__()
        self.data['id'] = name
        self.classes = classes

n1 = CustomNode("node 1", classes='class1')
n2 = CustomNode("node 2", classes='class2')
        
G = nx.Graph()

G.add_node(n1)
G.add_node(n2)

G.add_edge(n1, n2)

custom_inherited = ipycytoscape.CytoscapeWidget()
custom_inherited.graph.add_graph_from_networkx(G)
custom_inherited.set_style([
                        {
                            'selector': 'node.class1',
                            'css': {
                                'background-color': 'red'
                            }
                        },
                        {
                            'selector': 'node.class2',
                            'css': {
                                'background-color': 'green'
                            }
                        }])
custom_inherited
# -
# Fully customized using a dictionary (that we would have to build from our graph)

data = {
    'nodes': [
        { 'data': { 'id': 'desktop', 'name': 'C_leaf', 'href': 'http://cytoscape.org' } },
        { 'data': { 'id': 'a', 'name': 'C_root', 'href': 'http://cytoscape.org' } },
        { 'data': { 'id': 'b', 'name': 'C_wood', 'href': 'http://cytoscape.org' } },
        { 'data': { 'id': 'c', 'name': 'C_fast_soil', 'href': 'http://cytoscape.org' } },
        { 'data': { 'id': 'js', 'name': 'C_slow_soil', 'href': 'http://js.cytoscape.org' } }
    ],
    'edges': [
        {'data': { 'source': 'desktop', 'target': 'js' ,'name': 'my edge'}},
        {'data': { 'source': 'a', 'target': 'b' }},
        {'data': { 'source': 'a', 'target': 'c' }},
        {'data': { 'source': 'b', 'target': 'c' }},
        {'data': { 'source': 'js', 'target': 'b' }}
    ]
}




cytoscapeobj = ipycytoscape.CytoscapeWidget()
cytoscapeobj.graph.add_graph_from_json(data)


cytoscapeobj.set_style([
    {
                        'selector': 'node',
                        'css': {
                            'content': 'data(name)',
                            'text-valign': 'center',
                            'color': 'white',
                            'text-outline-width': 2,
                            'text-outline-color': 'green',
                            'background-color': 'green'
                        }
    },
    {
                        'selector': 'edge',
                        'css': {
                            'content': 'data(name)',
                            'text-valign': 'center',
                            'target-arrow-shape': 'triangle',
                            'target-arrow-color': 'black',
                            'source-arrow-color': 'black',
                            'line-color': '#333',
                            'width': 1.5,
                            'curve-style': 'bezier'
                        }
    },
    {
                        'selector': ':selected',
                        'css': {
                            'background-color': 'black',
                            'line-color': 'black',
                            'target-arrow-color': 'black',
                            'source-arrow-color': 'black',
                            'text-outline-color': 'black'
    }
    }
])
cytoscapeobj

# To do:
#  - How to plot subgraphs **INTO** the existing graph?
#  - How to use latex in edge or node descriptions.
#  





