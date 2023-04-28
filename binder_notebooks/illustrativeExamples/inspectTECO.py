# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

from IPython.display import HTML

display(HTML("<style>.container { width:100% !important; }</style>"))

# %load_ext autoreload
# %autoreload 2
import bgc_md2.helper as h

model_inspection = h.MvarSetInspectionBox()

from bgc_md2.resolve.mvars import CompartmentalMatrix, StateVariableTuple, VegetationCarbonInputPartitioningTuple,VegetationCarbonInputTuple
from ComputabilityGraphs import CMTVS

# +
from bgc_md2.models.TECOmm.source  import mvs as mvs_mm

h.compartmental_graph(mvs_mm)
# -

for key,fl in mvs_mm.get_InternalFluxesBySymbol().items():
    print(key);display(fl) 

# +
in_fluxes, internal_fluxes, out_fluxes = mvs_mm.get_InFluxesBySymbol(),mvs_mm.get_InternalFluxesBySymbol(),mvs_mm.get_OutFluxesBySymbol()

in_flux_targets, out_flux_sources = [[str(k) for k in d.keys()] for d in (in_fluxes, out_fluxes)] 

internal_connections = [(str(s),str(t)) for s,t in internal_fluxes.keys()]                                                                
# -


internal_connections

# +
#import networkx as nx
#import matplotlib.pyplot as plt
# -


import CompartmentalSystems.helpers_reservoir as hr
Gnx = hr.nxgraphs(mvs_mm.get_StateVariableTuple(),in_fluxes,internal_fluxes,out_fluxes)
#[Gnx.get_edge_data(s,t) for s,t in Gnx.edges]

# +
# hr.igraph_plot(mvs_mm.get_StateVariableTuple(),in_fluxes,internal_fluxes,out_fluxes)
# -

# # Let's inspect the vegetation part 

mvs_mm.get_VegetationCarbonStateVariableTuple()

mvs_mm.get_VegetationCarbonCompartmentalMatrix()

# what was formerly send to a soil pool is considered an output of the vegetation subsystem...
for k,v in mvs_mm.get_VegetationCarbonOutFluxesBySymbol().items():
    display(k,v)

for k,v in mvs_mm.get_VegetationCarbonInFluxesBySymbol().items():
    display(k,v)

for k,v in mvs_mm.get_VegetationCarbonInternalFluxesBySymbol().items():
    display(k,v)

combined = (
    set(mvs_mm.get_StateVariableTuple()),
    mvs_mm.get_InFluxesBySymbol(),
    mvs_mm.get_OutFluxesBySymbol(),
    mvs_mm.get_InternalFluxesBySymbol()
)
sv_set_veg = frozenset(mvs_mm.get_VegetationCarbonStateVariableTuple())

from sympy import Matrix,Symbol
state_vector_soil = Matrix([Symbol(s) for s in ["C_fastsom" ,"C_slowsom" ,"C_passsom"]])
# Probably the litter pools would be also  considered to be part of the soil subsystem.
# I just wanted to show that the division does not have tp be complete
# state_vector_soil = Matrix([C_metlit,C_stlit,C_fastsom,C_slowsom,C_passsom])
sv_set_soil = frozenset({sv for sv in state_vector_soil})


_,in_fluxes_veg,out_fluxes_veg,internal_fluxes_veg=hr.extract(combined,sv_set_veg) #obviously we do not need to return sv_set_veg, since it is an argument
_,in_fluxes_soil,out_fluxes_soil,internal_fluxes_soil=hr.extract(combined,sv_set_soil)

internal_fluxes_veg, in_fluxes_soil

part_dict =  {
    sv_set_veg:'green',
    sv_set_soil:'brown',
}


import igraph as ig
import matplotlib.pyplot as plt
#g = ig.Graph(n=10, edges=[[0, 1], [0, 5]])
Gnx = hr.nxgraphs(mvs_mm.get_StateVariableTuple(), in_fluxes, internal_fluxes, out_fluxes)
g=ig.Graph.from_networkx(Gnx)
type(Gnx)
l=g.layout_sugiyama()
#ig.plot(g,layout=g.layout_sugiyama())

hr.igraph_part_plot(
    mvs_mm.get_StateVariableTuple(),
    in_fluxes,
    internal_fluxes,
    out_fluxes,
    part_dict,
)

mvs_mm.provided_mvar_types


#Now we can compute the vegetation cycling matrix
hr.compartmental_matrix_2(
    out_fluxes_veg,
    internal_fluxes_veg,
    mvs_mm.get_VegetationCarbonStateVariableTuple()
)

#Now we can compute the soil cycling matrix
hr.compartmental_matrix_2(
    out_fluxes_soil,
    internal_fluxes_soil,
    state_vector_soil
)




