# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # VISIT_Kostia

import bgc_md2.helper as h
import importlib 

importlib.invalidate_caches()
mod = importlib.import_module('bgc_md2.models.VISIT_Kostia.source')
mvs = mod.mvs

h.compartmental_graph(mvs)

mvs.computable_mvar_types()

for t in mvs.computable_mvar_types():
    res = mvs._get_single_value(t)
    h.latex_render(t,res,capture=False) 

from IPython.display import Math


display(Math("\\text{InputTuple}"))


