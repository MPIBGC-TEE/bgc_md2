# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    StateVariableTuple
)
from importlib import import_module
model_folders=['kv_visit2', 'jon_yib','Aneesh_SDGVM','cable-pop','cj_isam','yz_jules']#,'kv_ft_dlem']
#model_folders=['cj_isam']
name=model_folders[0]
import_module("{}.source".format(name))
def mvs(name):
    return import_module("{}.source".format(name)).mvs

tups=[(name,mvs(name))
      for name in model_folders
]


dh.table(
    tups=tups,
    types_compact = [],
    types_expanded= [InputTuple,CompartmentalMatrix]

)

# +
# #%ls
# -



