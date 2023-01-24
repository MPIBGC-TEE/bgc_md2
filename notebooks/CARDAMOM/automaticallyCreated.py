# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # CARDAMOM

import bgc_md2.helper as h
import importlib 

importlib.invalidate_caches()
mod = importlib.import_module('bgc_md2.models.CARDAMOM.source')
mvs = mod.mvs

mvs.graph()

mvs.computable_mvar_names

for var in mvs.computable_mvar_types():
    mvs.render(var)
