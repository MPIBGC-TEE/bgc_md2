# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
#h.list_models()
# -

from ipywidgets import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h   
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    StateVariableTuple
)
dh.table(
    tups=[
        (name,h.CMTVS_from_model_name(name)) 
        for name in ["cable_general","TECOmm","TECO","Luo2012TE","TECO_general"]
        #for name in ["cable_yuanyuan","VISIT_Kostia"]
        #for name in h.list_models()
    ],
    types_compact = [],
    types_expanded = [InputTuple,CompartmentalMatrix]
    #types = [InputTuple,CompartmentalMatrix,StateVariableTuple]
    #types = [InputTuple,InternalFluxesBySymbol]
)
len(h.list_models())



