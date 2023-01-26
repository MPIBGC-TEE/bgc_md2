# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Luo2021 vs. Teco
# After finishing this notebook you will know:
# - how to compare models that look similar .. 
# - why you should avoid variable names like x1, x2.... 
#
#

from IPython.display import HTML

display(HTML("<style>.container { width:100% !important; }</style>"))

import bgc_md2.helper as h

model_inspection = h.MvarSetInspectionBox()

from bgc_md2.resolve.mvars import CompartmentalMatrix, StateVariableTuple, VegetationCarbonInputPartitioningTuple,VegetationCarbonInputTuple
from ComputabilityGraphs.CMTVS import CMTVS 
from ComputabilityGraphs.helpers import module_computers

# We look at two files that once were both part of the database but basically described the same model. 

from Luo2012TE.source import mvs as mvs_Luo #=  MVarSet.from_model_name('Luo2012TE')
from TECO.source import mvs as mvs_TECO  #MVarSet.from_model_name('TECO')

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

# So the two matrices have identical structure.
# However this could still be  misleading since  the state variables behind this description could be different.  
# Lets look at the state variables and their order as represented by the `StatevariableTuple`

mvs_Luo.get_StateVariableTuple()

mvs_TECO.get_StateVariableTuple()

# Ok, not very informative. Lets investigate the additional information that the translator of the model provided about the meaning of these symbols

# this could be displayed much nicer... room for a pullrequest 
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
mvs_TECO.provided_mvar_types


mvs_TECO.get_StateVariableTuple().subs(subs_dict)

import bgc_md2.resolve.computers
computers=module_computers(bgc_md2.resolve.computers)
type(computers)

mvs_subs = CMTVS(
    {
        var.subs(subs_dict) 
        for var in {
            mvs_TECO.get_CompartmentalMatrix(),
            mvs_TECO.get_StateVariableTuple(),
            mvs_TECO.get_InputTuple(),
            mvs_TECO.get_TimeSymbol()
        }
    },
    computers=module_computers(bgc_md2.resolve.computers)
)
# alternatively we could have created an update
#mvs_subs=mvs_TECO.update(
#    {
#        var.subs(subs_dict) 
#        for var in {mvs_TECO.get_CompartmentalMatrix(),mvs_TECO.get_StateVariableTuple(),mvs_TECO.get_InputTuple()}
#    }
#)
#mvs_subs

for v in mvs_subs.computable_mvar_types():
    display(mvs_subs._get_single_value_by_TypeTree(v))


# This description already shows that the CompartmentalMatrix and Statevector are probably not consistent.
# We can make this even more obvious by computing the outflux from the `C_root` pool.

 mvs_subs.get_OutFluxesBySymbol()[C_roots]

# and the internal flux from `C_roots` to `C_stlit`.

mvs_subs.get_InternalFluxesBySymbol()[(C_roots,C_stlit)]

# We have repaired this by exchanging the positions of `C_roots` and `C_woods` in the following version of the model which is part of the package (as you can see from the way it is imported).


from bgc_md2.models.TECOmm.source import mvs as mvs_mm

# Looking at the fluxes we can clearly see that they all co
for key,fl in mvs_mm.get_InternalFluxesBySymbol().items():
    display(key);display(fl) 

in_fluxes, internal_fluxes, out_fluxes = mvs_mm.get_InFluxesBySymbol(),mvs_mm.get_InternalFluxesBySymbol(),mvs_mm.get_OutFluxesBySymbol()
in_flux_targets, out_flux_sources = [[str(k) for k in d.keys()] for d in (in_fluxes, out_fluxes)] 
internal_connections = [(str(s),str(t)) for s,t in internal_fluxes.keys()]                                                                


internal_connections

# +
#import networkx as nx
#import matplotlib.pyplot as plt
# -


import CompartmentalSystems.helpers_reservoir as hr
Gnx = hr.nxgraphs(mvs_mm.get_StateVariableTuple(),in_fluxes,internal_fluxes,out_fluxes)
#[Gnx.get_edge_data(s,t) for s,t in Gnx.edges]

hr.igraph_plot(mvs_mm.get_StateVariableTuple(),in_fluxes,internal_fluxes,out_fluxes)
