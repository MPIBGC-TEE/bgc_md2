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

in_fluxes, internal_fluxes, out_fluxes = mvs_mm.get_InFluxesBySymbol(),mvs_mm.get_InternalFluxesBySymbol(),mvs_mm.get_OutFluxesBySymbol()

in_flux_targets, out_flux_sources = [[str(k) for k in d.keys()] for d in (in_fluxes, out_fluxes)]

internal_connections = [(str(s),str(t)) for s,t in internal_fluxes.keys()]


import CompartmentalSystems.helpers_reservoir as hr
Gnx, GVI, GINT, GVO = hr.nxgraphs(mvs_mm.get_StateVariableTuple(),in_fluxes,internal_fluxes,out_fluxes)
virtual_in_flux_sources=[n for n in GVI.nodes if GVI.in_degree(n)==0]
virtual_out_flux_targets=[n for n in GVO.nodes if GVO.out_degree(n)==0]



Gnx.edges.data

import igraph as ig
import networkx
G = ig.Graph.from_networkx(Gnx))

Gnx.edges[ ('virtual_in_C_wood', 'C_wood')]['type']


# +
vertex_size = [1 if v['virtual'] else 50 for v in G.vs]
labels = [v['_nx_name'] if not v['virtual'] else '' for v in G.vs]

edge_color_dict = {'in':'blue','internal':'black','out':'red'}
edge_colors = [edge_color_dict[e['type']] for e in G.es]
layout = G.layout('sugiyama')

#layout = G.layout('grid')
#layout = G.layout('kk')

ig.plot(
    G,
    layout=layout,
    vertex_size=vertex_size,
    vertex_label=labels,
    edge_color=edge_colors    
)
# -
# We have a function for this now
from CompartmentalSystems.helpers_reservoir import igraph_plot
p=igraph_plot(mvs_mm.get_StateVariableTuple(),in_fluxes,internal_fluxes,out_fluxes)
p

type(p)
