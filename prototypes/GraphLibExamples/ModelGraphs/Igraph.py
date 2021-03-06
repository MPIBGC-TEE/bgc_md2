# ---
# jupyter:
#   jupytext:
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
from bgc_md2.resolve.mvars import CompartmentalMatrix, StateVariableTuple, VegetationCarbonInputPartitioningTuple,VegetationCarbonInputTuple
from bgc_md2.resolve.MVarSet import MVarSet
display(HTML("<style>.container { width:100% !important; }</style>"))
mvs_mm =  MVarSet.from_model_name('TECOmm')

import igraph
g=igraph.Graph(directed=True)

in_fluxes, internal_fluxes, out_fluxes = mvs_mm.get_InFluxesBySymbol(),mvs_mm.get_InternalFluxesBySymbol(),mvs_mm.get_OutFluxesBySymbol()

in_flux_targets, out_flux_sources = [[str(k) for k in d.keys()] for d in (in_fluxes, out_fluxes)]

internal_connections = [(str(s),str(t)) for s,t in internal_fluxes.keys()]


# +
G=igraph.Graph(directed=True)
for i,s in enumerate(mvs_mm.get_StateVariableTuple()):
    G.add_vertex(str(s))# the name attribute is set automatically
    
virtual_in_flux_sources=["virtual_in_" + str(t) for t in in_flux_targets ]
for vs in virtual_in_flux_sources:
    G.add_vertex(vs,virtual=True)
    
for i in range(len(in_flux_targets)):
    # edges can be defined via node names
    # G.add_edge("virtual_in_C_foliage","C_foliage")
    G.add_edge(virtual_in_flux_sources[i],in_flux_targets[i],type='in')

for s,t in internal_connections:
    G.add_edge(s,t,type='internal')

virtual_out_flux_targets=["virtual_out_" + str(t) for t in out_flux_sources]
for vs in virtual_out_flux_targets:
    G.add_vertex(vs,virtual=True)

for i in range(len(out_flux_sources)):
    G.add_edge(out_flux_sources[i], virtual_out_flux_targets[i],type='out')


# +
vertex_size = [1 if v['virtual'] else 50 for v in G.vs]
labels = [v['name'] if not v['virtual'] else '' for v in G.vs]
edge_color_dict = {'in':'blue','internal':'black','out':'red'}
edge_colors = [edge_color_dict[e['type']] for e in G.es]
layout = G.layout('sugiyama')

#layout = G.layout('grid')
#layout = G.layout('kk')
igraph.plot(
    G,
    layout=layout,
    vertex_size=vertex_size,
    vertex_label=labels,
    edge_color=edge_colors    
)

# -
mvs_mm.get_InFluxesBySymbol()

layout.




