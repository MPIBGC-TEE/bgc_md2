# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# +
from bokeh.plotting import figure, show
import numpy as np
import networkx as nx
from functools import reduce
from typing import List,Tuple

def veclist_to_coordlists(veclist):
    xs=[v[0] for v in veclist]
    ys=[v[1] for v in veclist]
    return xs,ys
    #return list(zip(*veclist))

def norm(vec):
    return np.sqrt(np.dot(vec,vec))

def with_arrow_head(
        xs: List[float],
        ys: List[float], 
        width: float = 0.5,
        length: float =0.5,
    )->Tuple[List[float],List[float]]:
    
    last_point=np.array([
        xs[-1],
        ys[-1]
    ])
    second_last_point=np.array([
        xs[-2],
        ys[-2]
    ])
    last_vector=last_point-second_last_point
    last_vector_n=last_vector/norm(last_vector)
    print('last_vector',last_vector)
    left_head_vector=np.array([
        -last_vector_n[1],
         last_vector_n[0]
    ])    
    print('left_head_vector',left_head_vector)
    right_head_vector = -left_head_vector
    left_head_point=last_point -length*last_vector_n +width*left_head_vector
    
    right_head_point = last_point - length*last_vector_n  +  width*right_head_vector
    #we go first from the tip to left head
    left_stroke=[left_head_point,last_point]
    right_stroke=[right_head_point,last_point]

    
    left_xs,left_ys = veclist_to_coordlists(left_stroke)
    right_xs,right_ys = veclist_to_coordlists(right_stroke)
    #xss= [
    #    xs, 
    #    left_xs, 
    #    right_xs
    #]
    #yss=[
    #    ys, 
    #    left_ys, 
    #    right_ys
    #]
    xss=xs+left_xs+right_xs
    yss=ys+left_ys+right_ys
    return xss,yss

# create a new plot with a title and axis labels
p = figure(
    title="Simple line example", 
    x_axis_label='x', 
    y_axis_label='y'
)
# prepare some data
xs_1 = [ 3, 4, 5, 6, 7]
ys_1 = [ 1, 1, 1, 1, 1]
xs_2 = [ 3, 4, 5, 6, 7]
ys_2 = [ 2, 2, 2, 2, 2]
xss_1,yss_1=with_arrow_head(xs_1,ys_1,0.2,0.5)
xss_2,yss_2=with_arrow_head(xs_2,ys_2)
p.multi_line(
    [xss_1,xss_2],
    [yss_1,yss_2],
    line_color=["blue","red"], 
    line_width=4,
)
show(p)

# +
import math

from bokeh.io import output_file, show, output_notebook
from bokeh.models import Ellipse, GraphRenderer, StaticLayoutProvider
from bokeh.palettes import Spectral8
from bokeh.plotting import figure

N = 5#<8

node_indices = list(range(N))

plot = figure(title="Graph Layout Demonstration", x_range=(-1.1,1.1), y_range=(-1.1,1.1),
              tools="", toolbar_location=None)

graph = GraphRenderer()

graph.node_renderer.data_source.add(node_indices, 'index')
graph.node_renderer.data_source.add(Spectral8[:N], 'color')
graph.node_renderer.glyph = Ellipse(height=0.1, width=0.2, fill_color="color")

graph.edge_renderer.data_source.data = dict(
    start=[0]*N,
    end=node_indices)

### start of layout code
circ = [i*2*math.pi/N for i in node_indices]
x = [math.cos(i) for i in circ]
y = [math.sin(i) for i in circ]
graph_layout = dict(zip(node_indices, zip(x, y)))
graph.layout_provider = StaticLayoutProvider(graph_layout=graph_layout)

### Draw quadratic bezier paths
def bezier(start, end, control, steps):
    return [(1-s)**2*start + 2*(1-s)*s*control + s**2*end for s in steps]

def double_bezier(start_point,end_point,control_vector,steps):
    start_x,start_y = start_point
    end_x,end_y = end_point
    control_x,control_y = control_vector
    xs = bezier(start_x,end_x,control_x,steps)
    ys = bezier(start_y,end_y,control_y,steps)
    xss,yss=with_arrow_head(xs,ys,0.1,0.1)
    return (
        xss,
        yss
    )

n_steps=20
steps = [i/float(n_steps) for i in range(n_steps)]


# +
def add_line(acc,node_index):
    xss,yss=acc
    axs,ays= double_bezier(
        start_point=graph_layout[0],
        end_point=graph_layout[node_index],
        control_vector=(0,0),
        steps=steps
    )
    return( 
        xss+[axs],
        yss+[ays]
    )

xss,yss=reduce(add_line,node_indices,([],[]))
graph.edge_renderer.data_source.data['xs'] = xss
graph.edge_renderer.data_source.data['ys'] = yss

plot.renderers.append(graph)
output_notebook()
#output_file("graph.html")

show(plot)
# -



#xss[2]
#
## +
#import networkx as nx
#
#from bokeh.io import output_file, show
#from bokeh.models import (BoxSelectTool, Circle, EdgesAndLinkedNodes, HoverTool,
#                          MultiLine, NodesAndLinkedEdges, Plot, Range1d, TapTool)
#from bokeh.palettes import Spectral4
#from bokeh.plotting import from_networkx
#
#G=nx.karate_club_graph()
#
#plot = Plot(width=400, height=400,
#            x_range=Range1d(-1.1,1.1), y_range=Range1d(-1.1,1.1))
#plot.title.text = "Graph Interaction Demonstration"
#
#plot.add_tools(HoverTool(tooltips=None), TapTool(), BoxSelectTool())
#
#graph_renderer = from_networkx(G, nx.circular_layout, scale=1, center=(0,0))
#
#graph_renderer.node_renderer.glyph = Circle(size=15, fill_color=Spectral4[0])
#graph_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
#graph_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])
#
#graph_renderer.edge_renderer.glyph = MultiLine(line_color="#CCCCCC", line_alpha=0.8, line_width=5)
#graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
#graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)
#
#graph_renderer.selection_policy = NodesAndLinkedEdges()
#graph_renderer.inspection_policy = EdgesAndLinkedNodes()
#
#plot.renderers.append(graph_renderer)
#
#output_file("interactive_graphs.html")
#show(plot)


# -


