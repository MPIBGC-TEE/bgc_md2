import networkx as nx
import numpy as np
import inspect
from bokeh.io import output_file, show
from bokeh.plotting import figure#, from_networkx
from bokeh.models import (
    Ellipse,
    GraphRenderer,
    StaticLayoutProvider,
    BoxSelectTool,
    Circle,
    EdgesAndLinkedNodes,
    EdgesOnly,
    NodesOnly,
    HoverTool,
    MultiLine,
    NodesAndLinkedEdges,
    Plot,
    Range1d,
    TapTool,
    ColumnDataSource,
    Label,
    LabelSet,
    Range1d
)
from bokeh.palettes import Spectral4,Cividis256
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

def bezier(start, end, control, steps):
    return [(1-s)**2*start + 2*(1-s)*s*control + s**2*end for s in steps]

def double_bezier(start_point,end_point,control_vector,steps):
    start_x,start_y = start_point
    end_x,end_y = end_point
    control_x,control_y = control_vector
    xs = bezier(start_x,end_x,control_x,steps)
    ys = bezier(start_y,end_y,control_y,steps)
    xss,yss=with_arrow_head(xs,ys,0.5,0.5)
    return (
        xss,
        yss
    )



# +
def from_networkx(graph, layout_function, **kwargs):
        '''
        Generate a ``GraphRenderer`` from a ``networkx.Graph`` object and networkx
        layout function. Any keyword arguments will be passed to the
        layout function.

        Only two dimensional layouts are supported.

        Args:
            graph (networkx.Graph) : a networkx graph to render
            layout_function (function or dict) : a networkx layout function or mapping of node keys to positions.
            The position is a two element sequence containing the x and y coordinate.

        Returns:
            instance (GraphRenderer)

        .. note::
            Node and edge attributes may be lists or tuples. However, a given
            attribute must either have *all* lists or tuple values, or *all*
            scalar values, for nodes or edges it is defined on.

        .. warning::
            Node attributes labeled 'index' and edge attributes labeled 'start' or 'end' are ignored.
            If you want to convert these attributes, please re-label them to other names.

        Raises:
            ValueError

        '''

        # Handles nx 1.x vs 2.x data structure change
        # Convert node attributes
        node_dict = dict()
        node_attr_keys = [attr_key for node in list(graph.nodes(data=True))
                          for attr_key in node[1].keys()]
        node_attr_keys = list(set(node_attr_keys))

        for attr_key in node_attr_keys:
            values = [node_attr[attr_key] if attr_key in node_attr.keys() else None
                      for _, node_attr in graph.nodes(data=True)]

            values = _handle_sublists(values)

            node_dict[attr_key] = values

        if 'index' in node_attr_keys:
            from warnings import warn
            warn("Converting node attributes labeled 'index' are skipped. "
                 "If you want to convert these attributes, please re-label with other names.")

        node_dict['index'] = list(graph.nodes())
        N=graph.number_of_nodes()
        node_dict['color'] = Cividis256[:N]
        n_steps=100

        steps = [i/float(n_steps) for i in range(n_steps)]

        graph_renderer = GraphRenderer()
        if callable(layout_function):
            graph_layout = layout_function(graph, **kwargs)
        else:
            graph_layout = layout_function

            node_keys = graph_renderer.node_renderer.data_source.data['index']
            if set(node_keys) != set(layout_function.keys()):
                from warnings import warn
                warn("Node keys in 'layout_function' don't match node keys in the graph. "
                     "These nodes may not be displayed correctly.")

        # Convert edge attributes
        edge_dict = dict()
        edge_attr_keys = [attr_key for edge in graph.edges(data=True)
                          for attr_key in edge[2].keys()]
        edge_attr_keys = list(set(edge_attr_keys))

        for attr_key in edge_attr_keys:
            values = [edge_attr[attr_key] if attr_key in edge_attr.keys() else None
                      for _, _, edge_attr in graph.edges(data=True)]

            values = _handle_sublists(values)

            edge_dict[attr_key] = values

        if 'start' in edge_attr_keys or 'end' in edge_attr_keys:
            from warnings import warn
            warn("Converting edge attributes labeled 'start' or 'end' are skipped. "
                 "If you want to convert these attributes, please re-label them with other names.")

        edge_dict['start'] = [x[0] for x in graph.edges()]
        edge_dict['end'] = [x[1] for x in graph.edges()]
        xs_ys=[ 
            double_bezier(
                start_point=graph_layout[sp],
                end_point=graph_layout[ep],
                control_vector=(1,1),
                steps=steps
            )
            for (sp,ep) in graph.edges()
        ]    
        edge_dict['xs'] = [x for x,y in xs_ys]
        edge_dict['ys'] = [y for x,y in xs_ys]


        #graph_renderer = GraphRenderer()
        #graph_renderer.node_renderer.data_source.add(Cividis256[:N], 'color')
        #from IPython import embed;embed()
        #graph_renderer.node_renderer.glyph = Ellipse(height=0.1, width=0.2, fill_color="color")
        graph_renderer.node_renderer.data_source.data = node_dict
        graph_renderer.edge_renderer.data_source.data = edge_dict

        #def add_line(acc,edge):
        #    xss,yss=acc
        #    sp,ep=edge
        #    axs,ays= double_bezier(
        #        start_point=graph_layout[sp],
        #        end_point=graph_layout[ep],
        #        control_vector=(0,0),
        #        steps=steps
        #    )
        #    return( 
        #        xss+[axs],
        #        yss+[ays]
        #    )
        #
        #xss,yss=reduce(add_line,graph.edges,([],[]))
        #graph_renderer.edge_renderer.data_source.data['xs'] = xss
        #graph_renderer.edge_renderer.data_source.data['ys'] = yss

        #graph_renderer.edge_renderer.data_source.data['start'] = [x[0] for x in graph.edges()]
        #graph_renderer.edge_renderer.data_source.data['end'] = [x[1] for x in graph.edges()]
        #graph_renderer.node_renderer.data_source.data['xs'] = [v[0] for v in graph_layout.values()]
        #graph_renderer.node_renderer.data_source.data['ys'] = [v[1] for v in graph_layout.values()]

        graph_renderer.layout_provider = StaticLayoutProvider(graph_layout=graph_layout)
        #graph_renderer.node_renderer.glyph = Circle(size=15, fill_color=Spectral4[0])
        graph_renderer.node_renderer.glyph = Circle(size=30, fill_color='color')
        graph_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
        graph_renderer.node_renderer.hover_glyph = Circle(size=25, fill_color=Spectral4[1])
        
        graph_renderer.edge_renderer.glyph = MultiLine(line_color="#CCCCCC", line_alpha=0.8, line_width=5)
        graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
        graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)
        
        graph_renderer.inspection_policy = EdgesAndLinkedNodes()
        #graph_renderer.selection_policy = NodesAndLinkedEdges()
        #graph_renderer.inspection_policy = EdgesOnly()
        #graph_renderer.inspection_policy = NodesOnly()
        #graph_renderer.inspection_policy = NodesAndLinkedEdges()

        return graph_renderer

G=nx.DiGraph()
G.add_edges_from([("A",n) for n in ("B","C","D","E","F")])
plot = figure(title="Networkx Integration Demonstration")
#, x_range=(-1.1,1.1), y_range=(-1.1,1.1),
#              tools="", toolbar_location=None)
node_hover_tool = HoverTool(tooltips=[("index", "$index"), ("start","@start"),("end","@end")])
plot.add_tools(node_hover_tool, TapTool(), BoxSelectTool())
layout_function=nx.spring_layout
graph_layout = layout_function(G)
names= list(graph_layout.keys())
source = ColumnDataSource(
        data=dict(
            xs=[graph_layout[n][0] for n in names],
            ys=[graph_layout[n][1] for n in names],
            names=names 
        )
)
labels = LabelSet(
        x='xs', 
        y='ys', 
        text='names',
        x_offset=-15, 
        y_offset=-15, 
        source=source, 
        border_line_color='white',
        background_fill_color='white',
        render_mode='canvas'
)

graph = from_networkx(G, graph_layout, scale=2, center=(0,0))


plot.renderers.append(graph)
plot.add_layout(labels)
output_file("networkx_graph.html")
show(plot)
