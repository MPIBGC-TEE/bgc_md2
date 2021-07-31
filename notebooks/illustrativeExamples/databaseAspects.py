# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# This illustrative notebook shows the interplay of all three packages on the collection of predefined models in bgc_md2 and ways to extend this collection by defining a new model that can be compared with respect to diagnostic variables that are computed using LAPM and CompartmentalSystems.
#

# +
# adjust the output
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# #%load_ext autoreload
# %autoreload 2
# -

# # Database aspects of bgc_md
# ## Collection
# First we demonstrate the collection aspect of bgc_md by listing all the models.<p>
#  **Fixme mm**: 
#   * *activate the thumbnail plots again*
#   
#   
#

# +
import bgc_md2.helper as h

model_inspection = h.MvarSetInspectionBox()
model_list = h.ModelListGridBox(
    inspection_box=model_inspection,
)
model_list
# -

# ## Inspection with respect to common diagnostics
# The variable in the next cell is connected to the list above and will show a summary of the model you click. This summary contains not only the data the author of the model put into the database but also all the results that can be computed from it. It also contains a button to create a model specific notebook.

model_inspection

# ## Queriyng
# One purpose behind the unified (and therefore restricted)  representation of data of a database is the ability to select subsets. In bgc_md a record typically describes what we know about a model or a particular simulations.
# We demonstrate the querying capability by selectiong models for which the variables
# `VegetationCarbonInputPartitioningTuple` and `NumericSolutionArray` are defined or computable which indicates models that include vegetation pools and for which we also have enough example data to run a simulation.<p>
# **Fixme mm** <p>
#
# * There are no examples for the computation of a transit time or age density computation. Not even a variable with that name exists (in bgc_md2) yet.

from bgc_md2.resolve.mvars import VegetationCarbonInputPartitioningTuple, NumericSolutionArray
from bgc_md2.resolve.MVarSet import MVarSet

li = h.list_target_models(
    target_classes=frozenset(
        {
            VegetationCarbonInputPartitioningTuple,
            NumericSolutionArray
        }
    ),
    # explicit_exclude_models=frozenset({'CARDAMOM'})
)
li    


# From these we chose two models to investigate more thoroughly.
# We load the first one with a helper function and the second by using standard python tools.
#
#

mvs_Luo =  MVarSet.from_model_name('Luo2012TE')
# or using a moduel import
from bgc_md2.models.TECO.source import mvs as mvs_TECO

# # ComputabilityGraphs
# Now that we actually have two records (MVarSets) we can exlore what we can comptute from them.
# Just add a dot "." behind the variable in the next cell and press the tab key!
#

mvs_TECO

# The options provided by the python interpreter are actually the result of a graph computation.
# To see all computable mvars of the TECO MVarSet execute the next cell!
#

mvs_TECO.computable_mvar_names


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

# Ok, not very informative. Lets investigate the additional information that the translator of the model provided about the meaning of these symbols

# this could be displayed much nicer...
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

hr.igraph_plot(mvs_mm.get_StateVariableTuple(),in_fluxes,internal_fluxes,out_fluxes)

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

state_vector_soil = Matrix([C_fastsom,C_slowsom,C_passsom])
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
hr.igraph_part_plot(
    mvs_mm.get_StateVariableTuple(),
    in_fluxes,
    internal_fluxes,
    out_fluxes,
    part_dict
)

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
