import math
from functools import reduce
from bokeh.io import output_file, show, output_notebook
from bokeh.models import Ellipse, GraphRenderer, StaticLayoutProvider
from bokeh.palettes import Spectral8
from bokeh.plotting import figure
from helpers import with_arrow_head
import networkx as nx

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
#output_notebook()
#output_file("graph.html")

show(plot)
