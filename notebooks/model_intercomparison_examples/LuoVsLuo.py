# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from IPython.display import HTML

display(HTML("<style>.container { width:100% !important; }</style>"))

# #%load_ext autoreload
# %autoreload 2
import bgc_md2.helper as h

model_inspection = h.MvarSetInspectionBox()

from bgc_md2.resolve.mvars import CompartmentalMatrix, StateVariableTuple, VegetationCarbonInputPartitioningTuple,VegetationCarbonInputTuple
from bgc_md2.resolve.MVarSet import MVarSet

# This time we are only interested in Vegetation models and Carbon input partitioning. Therefore we look for models for which the variable
# `VegetationCarbonInputPartitioningTuple` is defined or computable.

# +
#li = h.list_target_models(
#    target_classes=frozenset(
#        {
#            CompartmentalMatrix,
#            StateVariableTuple,
#            VegetationCarbonInputPartitioningTuple,
#            VegetationCarbonInputTuple
#            
#        }
#    ),
#    # explicit_exclude_models=frozenset({'CARDAMOM'})
#)
#li    
# -


# From these we chose two models to investigate more thoroughly.
#

mvs_Luo =  MVarSet.from_model_name('Luo2012TE')
mvs_TECO =  MVarSet.from_model_name('TECO')

mvs_Luo.get_InputTuple()

mvs_TECO.get_InputTuple()

B_Luo = mvs_Luo.get_CompartmentalMatrix();B_Luo

B_TECO = mvs_TECO.get_CompartmentalMatrix();B_TECO


# The matrices look structurally similar. Lets check if this is really the case.

# +
from sympy import Matrix
def MatrixStructure(M):
    return Matrix(M.rows,M.cols,lambda i, j : 1 if M[i,j]!=0 else 0)

S_Luo = MatrixStructure(B_Luo)
S_TECO = MatrixStructure(B_TECO)
for x in (S_Luo,S_TECO,S_Luo == S_TECO):
    display(x)

# -

# So the two matrices look very similar.
# However this could still be  misleading since  the state variables behind this description could be different.  
# Lets look at the state variables and their order as represented by the `StatevariableTuple`

mvs_Luo.get_StateVariableTuple()

mvs_TECO.get_StateVariableTuple()

# Ok not very informative. Lets investigate the additional information that the translator of the model provided about the meaning of these symbols

bib_Luo=mvs_Luo.get_BibInfo();bib_Luo.sym_dict


bib_TECO=mvs_TECO.get_BibInfo();bib_TECO.sym_dict


# We can see that although the matrix structures are identical the meaning of the state variables differs between the two models!

[bib_Luo.sym_dict[str(sym)] for sym in mvs_Luo.get_StateVariableTuple()]


[bib_TECO.sym_dict[str(sym)] for sym in mvs_TECO.get_StateVariableTuple()]

# So the positions of **root**  and **wood** carbon pools is exchanged in the state vectors. 
# To check more easily what is going on, we choose more descriptive names for the state variables that are not related to a specific indexing scheme.
# The names are chosen to reflect the meaning as recorded in bib_TECO.sym_dict (as printed above).


rename_dict = {
        "x_1":    "C_foliage",
        "x_2":    "C_roots",
        "x_3":    "C_wood",
        "x_4":    "C_metlit",
        "x_5":    "C_stlit",
        "x_6":    "C_fastsom",
        "x_7":    "C_slowsom",
        "x_8":    "C_passsom",
        "b_1":    "b_foliage",
        "b_2":    "b_roots",
        "b_3":    "b_wood",
        "c_1":    "cr_foliage",
        "c_2":    "cr_wood",
        "c_3":    "cr_fineroot",
        "c_4":    "cr_metlit",
        "c_5":    "cr_stlit",
        "c_6":    "cr_fastsom",
        "c_7":    "cr_slowsom",
        "c_8":    "cr_passsom",
        "f_41":   "f_foliage2metlit",
        "f_51":   "f_foliage2stlit",
        "f_52":   "f_wood2stlit",
        "f_43":   "f_fineroots2metlit",
        "f_53":   "f_fineroots2stlit",
        "f_64":   "f_metlit2fastsom",
        "f_65":   "f_stlit2fastsom",
        "f_75":   "f_stlit2slowsom",
        "f_76":   "f_fastsom2slowsom",
        "f_86":   "f_fastsom2passsom",
        "f_67":   "f_slowsom2fastsom",
        "f_87":   "f_slowsom2passsom",
        "f_68":   "f_passsom2fastsom",
}
from sympy import var, Symbol
for old, new in rename_dict.items():
    var(old)
    var(new)
subs_dict = { 
    Symbol(old): Symbol(new)
    for old, new in rename_dict.items()
}
# We want to compute the fluxes with renamed symbols.
# We first check what we have to substitute.
# +

mvs_TECO.provided_mvar_types


# -


mvs_TECO.get_StateVariableTuple().subs(subs_dict)

mvs_subs=MVarSet({var.subs(subs_dict) for var in {mvs_TECO.get_CompartmentalMatrix(),mvs_TECO.get_StateVariableTuple(),mvs_TECO.get_InputTuple()} })

for v in mvs_subs.computable_mvar_types():
    display(mvs_subs._get_single_mvar_value(v))


# This description already that the CompartmentalMatrix and Statevector are probably not consistent.
# We can make this even more obvious by computing the outflux from the `C_root` pool.

 mvs_subs.get_OutFluxesBySymbol()[C_roots]

# and the internal flux from `C_roots` to `C_stlit`.

mvs_subs.get_InternalFluxesBySymbol()[(C_roots,C_stlit)]

# We can probably repair this by exchanging the positions of `C_roots` and `C_woods` as has been done in the following version of the model


mvs_mm =  MVarSet.from_model_name('TECOmm')

for key,fl in mvs_mm.get_InternalFluxesBySymbol().items():
    print(key);display(fl) 

in_fluxes, internal_fluxes, out_fluxes = mvs_mm.get_InFluxesBySymbol(),mvs_mm.get_InternalFluxesBySymbol(),mvs_mm.get_OutFluxesBySymbol()


in_flux_targets, out_flux_sources = [[str(k) for k in d.keys()] for d in (in_fluxes, out_fluxes)] 

internal_connections = [(str(s),str(t)) for s,t in internal_fluxes.keys()]                                                                

internal_connections

import networkx as nx
import matplotlib.pyplot as plt


# +
virtual_in_flux_sources=["virtual_in_" + str(t) for t in in_flux_targets]
virtual_in_flux_sources
GVI=nx.DiGraph()
for n in virtual_in_flux_sources:
    GVI.add_node(n,virtual=True)
for n in in_flux_targets:
    GVI.add_nodes_from(in_flux_targets)
for i in range(len(in_flux_targets)):
    GVI.add_edge(virtual_in_flux_sources[i],in_flux_targets[i])

GVI.nodes,GVI.edges

# +
virtual_out_flux_targets=["virtual_out_" + str(t) for t in out_flux_sources]
GVO=nx.DiGraph()
for n in virtual_out_flux_targets:
    GVO.add_node(n,virtual=True)
for n in out_flux_sources:
    GVO.add_nodes_from(out_flux_sources)
for i in range(len(out_flux_sources)):
    GVO.add_edge(out_flux_sources[i], virtual_out_flux_targets[i])

GVO.nodes,GVO.edges
# -

GINT=nx.DiGraph()
for c in internal_connections:
    GINT.add_edge(c[0],c[1])
GINT.edges
GINT.nodes


set(GVI.nodes).intersection(GINT.nodes)
G1=nx.compose(GVI,GINT)
G=nx.compose(G1,GVO)




# +
#import CompartmentalSystems.helpers_reservoir as hr
#G, GVI, GINT, GVO = hr.nxgraphs(mvs_mm.get_StateVariableTuple,in_fluxes,internal_fluxes,out_fluxes)
# -

# # matplotlib plotting

# +
ax = plt.axes()
#pos=nx.spiral_layout(G)
#pos=nx.circular_layout(G)
#pos=nx.spring_layout(G)
#pos=nx.kamada_kawai_layout(G)

#pos=nx.planar_layout(G)
#pos=nx.shell_layout(G)
#pos=nx.spectral_layout(G)
pos=nx.spring_layout(G)
virtual_node_options={
    'node_size': 10,
}
real_node_options={
    'node_color': 'black',
    'node_size': 100,
}
ax=plt.axes()
nx.draw_networkx_nodes(ax=ax,pos=pos,G=G,nodelist=virtual_in_flux_sources,**virtual_node_options,node_color='blue')
nx.draw_networkx_nodes(ax=ax,pos=pos,G=G,nodelist=in_flux_targets,**real_node_options)

nx.draw_networkx_nodes(ax=ax,pos=pos,G=G,nodelist=GINT.nodes,**real_node_options)
nx.draw_networkx_nodes(ax=ax,pos=pos,G=G,nodelist=virtual_out_flux_targets,**virtual_node_options,node_color='red')

nx.draw_networkx_edges(ax=ax,pos=pos,G=G)
# -

list(G.nodes());G.nodes['virtual_in_C_foliage']

# +
 
from bokeh.io import output_file, show
from bokeh.models import (BoxSelectTool, Circle, Text, EdgesAndLinkedNodes, HoverTool,
                          MultiLine, ArrowHead, NodesAndLinkedEdges, Plot, Range1d, TapTool,)
from bokeh.palettes import Spectral4
from bokeh.plotting import from_networkx

#G = nx.karate_club_graph()


plot = Plot(plot_width=400, plot_height=400,
            x_range=Range1d(-1.1,1.1), y_range=Range1d(-1.1,1.1))
plot.title.text = "Graph Interaction Demonstration"

plot.add_tools(HoverTool(tooltips=None), TapTool(), BoxSelectTool())

graph_renderer = from_networkx(G, nx.circular_layout, scale=1, center=(0,0))

graph_renderer.node_renderer.glyph = Circle(size=15, fill_color=Spectral4[0])
#graph_renderer.node_renderer.glyph = Text()
graph_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
graph_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])

graph_renderer.edge_renderer.glyph = MultiLine(line_color="#CCCCCC", line_alpha=0.8, line_width=5)
graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)

graph_renderer.selection_policy = NodesAndLinkedEdges()
graph_renderer.inspection_policy = EdgesAndLinkedNodes()


plot.renderers.append(graph_renderer)

output_file("networkx_graph.html")
show(plot)


# +
VC, RC = "black", "red"

G.nodes.data()
# -

node_attrs = {}
for n,y in G.nodes(data=True):
    #print(n)
    #print(y)
    if 'virtual' in y.keys():
        node_attrs[n]=VC
    #edge_color = SAME_CLUB_COLOR if G.nodes[start_node]["club"] == G.nodes[end_node]["club"] else DIFFERENT_CLUB_COLOR
    #edge_attrs[(start_node, end_node)] = edge_color
node_attrs
nx.set_node_attributes(G,node_attrs,'node_color')

# +
plot = Plot(plot_width=400, plot_height=400,
            x_range=Range1d(-1.1,1.1), y_range=Range1d(-1.1,1.1))
plot.title.text = "Graph Interaction Demonstration"

plot.add_tools(HoverTool(tooltips=None), TapTool(), BoxSelectTool())

graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0,0))

graph_renderer.node_renderer.glyph = Circle(size=15, fill_color="node_color")
#graph_renderer.node_renderer.glyph = Text()
graph_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
graph_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])

graph_renderer.edge_renderer.glyph = MultiLine(line_color="#CCCCCC", line_alpha=0.8, line_width=5)


graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5) 
graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)

graph_renderer.selection_policy = NodesAndLinkedEdges()
graph_renderer.inspection_policy = EdgesAndLinkedNodes()


plot.renderers.append(graph_renderer)

output_file("networkx_graph.html")
show(plot)
# -

graph_renderer.edge_renderer

# To do:
#     implement an edge_renderer.glyph that can draw arrows instead of lines (directed graphs)


