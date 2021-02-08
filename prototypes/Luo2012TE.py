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

# # Luo2012TE

import bgc_md2.helper as h
import importlib 

importlib.invalidate_caches()
mod = importlib.import_module('bgc_md2.models.Luo2012TE.source')
mvs = mod.mvs

mvs.graph()

mvs.computable_mvar_names

for var in mvs.computable_mvar_types():
    mvs.render(var)
