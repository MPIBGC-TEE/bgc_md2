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
from sympy import Matrix,Symbol,symbols
import networkx as nx

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

mvs_mm =  MVarSet.from_model_name('TECOmm')

combined = (
    set(mvs_mm.get_StateVariableTuple()),
    mvs_mm.get_InFluxesBySymbol(),
    mvs_mm.get_OutFluxesBySymbol(),
    mvs_mm.get_InternalFluxesBySymbol()
)

sv_set_veg = frozenset(mvs_mm.get_VegetationCarbonStateVariableTuple())

C_foliage,C_wood,C_roots,C_metlit,C_stlit,C_fastsom,C_slowsom,C_passsom=symbols("C_foliage,C_wood,C_roots,C_metlit,C_stlit,C_fastsom,C_slowsom,C_passsom")
# state_vector_soil = Matrix([C_fastsom,C_slowsom,C_passsom])
# Probably the litter pools would be also  considered to be part of the soil subsystem.
# I just wanted to show that the division does not have tp be complete
state_vector_soil = Matrix([C_metlit,C_stlit,C_fastsom,C_slowsom,C_passsom])
sv_set_soil = frozenset({sv for sv in state_vector_soil})


# %matplotlib notebook
import matplotlib.pyplot as plt


# +
# let us construct  seperate graphs for  vegetation and soil sub models but with several virtual in and out pools
#
# 'Ausschneiden mit Rand' 'cutout with boundary'
#
# The following function is not necessary. 
# The graphs for vegetation and soil models could be created completely independently.
# However in order to construct a round trip example it is useful to first extract the models from a bigger model
# and then reconstruct it from the parts.
# In a real application this function could be used to extract submodels from DIFFERENT bigger models.
    
def subgraph_with_virtual_nodes(mvs,sv_set):
    # Ausschneiden mit Rand...
    g=nx.DiGraph()
    for s in sv_set:
        g.add_node(s)

    for t,fl in mvs.get_InFluxesBySymbol().items() :
        if t in sv_set:
            s=Symbol("v_in_to_"+str(t))
            g.add_node(s,virtual=True,type='in')
            g.add_edge(s,t)

    for tup, fl in mvs.get_InternalFluxesBySymbol().items():
        s,t=tup
        if s in sv_set and t in sv_set:
            g.add_edge(s,t)
        elif s in sv_set:
            t=Symbol("v_out_from_"+str(s)+"_to_"+str(t))
            g.add_node(t,virtual=True,type='out')
            g.add_edge(s,t)
        elif t in sv_set:
            s=Symbol("v_in_from_"+str(s)+"_to_"+str(t))
            g.add_node(s,virtual=True,type='in')
            g.add_edge(s,t)
            
            

    for s,fl in mvs.get_OutFluxesBySymbol().items() :
        if s in sv_set:
            t=Symbol("v_out_from_"+str(s))
            g.add_node(t,virtual=True,type='out')
            g.add_edge(s,t)
    return g
# -

g_veg=subgraph_with_virtual_nodes(mvs_mm,sv_set_veg)
g_soil=subgraph_with_virtual_nodes(mvs_mm,sv_set_soil)

[str(n) for n in g_veg.nodes]

[str(n) for n in g_soil.nodes]

# +
# %matplotlib notebook
ax1=plt.subplot(2,1,1)
ax2=plt.subplot(2,1,2)

nx.draw(
    g_veg,
    ax=ax1,
    labels={n:str(n) for n in g_veg.nodes},
    nodelist=list(g_veg.nodes),
    node_size=[20 if 'virtual' in g_veg.nodes[n].keys() else 200 for n in g_veg.nodes],
    node_color='g',
    font_size=6,
    pos=nx.spring_layout(g_veg,1,seed=1)
)

#nx.draw(g_veg,ax=ax1,labels={n:str(n) for n in g_veg.nodes},node_size=200,node_color='g')
nx.draw(
    g_soil,
    ax=ax2,
    labels={n:str(n) for n in g_soil.nodes},
    nodelist=list(g_soil.nodes),
    node_size=[20 if 'virtual' in g_soil.nodes[n].keys() else 200 for n in g_soil.nodes],
    node_color='r',
    font_size=6,
    #pos=nx.kamada_kawai_layout(g_soil)
    pos=nx.spring_layout(g_soil,1,seed=1)

)
# -

#[n for n in g_soil.nodes if 'virtual' in g_soil.nodes[n].keys() and g_soil.nodes[n]['type']=='in']  
[n for n in g_soil.nodes if 'virtual' in g_soil.nodes[n].keys() and g_soil.nodes[n]['type']=='out']  


def plot(g,ax,pos,svs1,svs2):
    def node_color(n):
        if n in svs1:
            c='g'
        elif n in svs2:
            c='r'
        else:
            c='b'
        
        return c
    
    nx.draw(
    g,
    ax=ax,
    labels={n:str(n) for n in g.nodes},
    nodelist=list(g.nodes),
    node_size=[2 if 'virtual' in g.nodes[n].keys() else 400 for n in g.nodes],
    node_color=[node_color(n) for n in g.nodes],
    font_size=6,
    #pos=nx.kamada_kawai_layout(g_soil)
    pos=pos
)


# +
# In order to combine the two graphs 
# the first step is to build a composed graph (nodes= union of nodes , edges= union of edges) 
g0=nx.compose(g_veg,g_soil)

# %matplotlib notebook
ax1=plt.subplot(1,1,1)
pos=nx.spring_layout(g0,1,seed=1)
#pos=nx.kamada_kawai_layout(g0)
#pos=nx.spectral_layout(g0)
plot(g0,ax1,pos,sv_set_veg,sv_set_soil)


# -

# The second bunch of steps is to remove the virtual out and influxes and replace them by new internal fluxes
# this could be done graphically either 
# - by grabbing the arrow head of an outflux edge and move it to the new target node + removing the superflus input edge
#   or alternatively 
# - by grabbing a virtual output node and dragging it on the intended virtual input node, which would cause both virtual nodes to disappear and a direct 
#   egde to be created
#   here its done by a function that takes the edges that should be 'fused' together
def virtual_to_internal(
    g_old,
    e1, #outflux edge to be removed
    e2, #influx edge to be removed
):
    m,v_out=e1 
    v_in,n=e2
    g=g_old.copy()
    #expr=g.get_edge_data(C_roots,v_out)#['expr']# do to. this is the expression that we will later add to the new internal flux
    g.remove_edge(m,v_out) 
    g.remove_node(v_out)
    g.remove_edge(v_in,n) 
    g.remove_node(v_in)
    g.add_edge(m,n) # to do: add expression that was saved
    return g
    


g1=virtual_to_internal(g0,(C_roots,Symbol('v_out_from_C_roots_to_C_metlit')),(Symbol('v_in_from_C_roots_to_C_metlit'),C_metlit))
# %matplotlib notebook
ax1=plt.subplot(1,1,1)
plot(g1,ax1,pos,sv_set_veg,sv_set_soil)

g2=virtual_to_internal(g1,(C_roots,Symbol('v_out_from_C_roots_to_C_stlit')),(Symbol('v_in_from_C_roots_to_C_stlit'),C_stlit))
# %matplotlib notebook
ax1=plt.subplot(1,1,1)
plot(g2,ax1,pos,sv_set_veg,sv_set_soil)

g3=virtual_to_internal(g2,(C_wood,Symbol('v_out_from_C_wood_to_C_stlit')),(Symbol('v_in_from_C_wood_to_C_stlit'),C_stlit))
# %matplotlib notebook
ax1=plt.subplot(1,1,1)
plot(g3,ax1,pos,sv_set_veg,sv_set_soil)


g4=virtual_to_internal(g3,(C_foliage,Symbol('v_out_from_C_foliage_to_C_metlit')),(Symbol('v_in_from_C_foliage_to_C_metlit'),C_metlit))
# %matplotlib notebook
ax1=plt.subplot(1,1,1)
plot(g4,ax1,pos,sv_set_veg,sv_set_soil)


g5=virtual_to_internal(g4,(C_foliage,Symbol('v_out_from_C_foliage_to_C_stlit')),(Symbol('v_in_from_C_foliage_to_C_stlit'),C_stlit))
# %matplotlib notebook
ax1=plt.subplot(1,1,1)
plot(g5,ax1,pos,sv_set_veg,sv_set_soil)



