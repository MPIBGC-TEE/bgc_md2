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
# prepare some data
x = [ 3, 4, 5]
y = [ 3, 4, 5]

last_point=np.array([
    x[-1],
    y[-1]
])
second_last_point=np.array([
    x[-2],
    y[-2]
])
last_vector=last_point-second_last_point
left_head_vector=np.array([
    -last_vector[0],
     last_vector[1]
])    
right_head_vector=np.array([
     last_vector[0],
    -last_vector[1]
])    
width = 0.5
left_head_point=second_last_point+width*left_head_vector
right_head_point=second_last_point+width*right_head_vector
#we go first from the tip to left head
left_stroke=[left_head_point,last_point]
right_stroke=[right_head_point,last_point]

def veclist_to_coordlists(veclist):
    xs=[v[0] for v in veclist]
    ys=[v[1] for v in veclist]
    return xs,ys
    
left_xs,left_ys = veclist_to_coordlists(left_stroke)
right_xs,right_ys = veclist_to_coordlists(right_stroke)

# create a new plot with a title and axis labels
p = figure(
    title="Simple line example", 
    x_axis_label='x', 
    y_axis_label='y'
)
p.multi_line(
    [
        x, 
        left_xs, 
        right_xs
    ],
    [
        y, 
        left_ys, 
        right_ys
    ],
    color="blue", 
    line_width=4,
)
show(p)
np.dot(right_head_vector,last_vector)

# +
import math

from bokeh.io import output_file, show
from bokeh.models import Ellipse, GraphRenderer, StaticLayoutProvider
from bokeh.palettes import Spectral8
from bokeh.plotting import figure

N = 8
node_indices = list(range(N))

plot = figure(title="Graph Layout Demonstration", x_range=(-1.1,1.1), y_range=(-1.1,1.1),
              tools="", toolbar_location=None)

graph = GraphRenderer()

graph.node_renderer.data_source.add(node_indices, 'index')
graph.node_renderer.data_source.add(Spectral8, 'color')
graph.node_renderer.glyph = Ellipse(height=0.1, width=0.2, fill_color="color")

graph.edge_renderer.data_source.data = dict(
    start=[0]*N,
    end=node_indices)

### start of layout code
circ = [i*2*math.pi/8 for i in node_indices]
x = [math.cos(i) for i in circ]
y = [math.sin(i) for i in circ]
graph_layout = dict(zip(node_indices, zip(x, y)))
graph.layout_provider = StaticLayoutProvider(graph_layout=graph_layout)

### Draw quadratic bezier paths
def bezier(start, end, control, steps):
    return [(1-s)**2*start + 2*(1-s)*s*control + s**2*end for s in steps]

xs, ys = [], []
sx, sy = graph_layout[0]
steps = [i/100. for i in range(100)]
for node_index in node_indices:
    ex, ey = graph_layout[node_index]
    xs.append(bezier(sx, ex, 0, steps))
    ys.append(bezier(sy, ey, 0, steps))
graph.edge_renderer.data_source.data['xs'] = xs
graph.edge_renderer.data_source.data['ys'] = ys

plot.renderers.append(graph)

output_file("graph.html")
show(plot)
# -



sx,sy,ex,ey

bezier(1.0,2.0,0,[1/10*i for i in range(10)])


