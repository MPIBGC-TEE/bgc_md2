# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
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

li = h.list_target_models(
    target_classes=frozenset(
        {
            CompartmentalMatrix,
            StateVariableTuple,
            VegetationCarbonInputPartitioningTuple,
            VegetationCarbonInputTuple
            
        }
    ),
    # explicit_exclude_models=frozenset({'CARDAMOM'})
)
li    


# From these we chose two models to investigate more thoroughly.
#

mvs1 =  MVarSet.from_model_name('Luo2012TE')
mvs2 =  MVarSet.from_model_name('TECO')

mvs1.get_InputTuple()

mvs2.get_InputTuple()

B1 = mvs1.get_CompartmentalMatrix();B1

B2 = mvs2.get_CompartmentalMatrix();B2


# The matrices look structurally similar. Lets check if this is really the case.

# +
from sympy import Matrix
def MatrixStructure(M):
    return Matrix(M.rows,M.cols,lambda i, j : 1 if M[i,j]!=0 else 0)

S1 = MatrixStructure(B1)
S2 = MatrixStructure(B2)
for x in (S1,S2,S1 == S2):
    display(x)

# -

# So the two matrices look very similar.
# However this could still be very misleading since  the statevariables behind this discription could be very different.  
# Lets look at the statevariables and their order as represented by the StatevariableTuple

mvs1.get_StateVariableTuple()

mvs2.get_StateVariableTuple()

# Lets investigate the additional information that the translator of the model provided about the meaning of these symbols

bib1=mvs1.get_BibInfo();bib1.sym_dict


bib2=mvs2.get_BibInfo();bib2.sym_dict


# We can see that although the matrices are identical the ordering of the state variables differs between the two models!

[bib1.sym_dict[str(sym)] for sym in mvs1.get_StateVariableTuple()]


[bib2.sym_dict[str(sym)] for sym in mvs2.get_StateVariableTuple()]

# So the positions of roots and woods is exchanged

mvs1.


