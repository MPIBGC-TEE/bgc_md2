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

# # ElMasri2013AgricForMeteorol

import bgc_md2.helper as h
import importlib 

importlib.invalidate_caches()
mod = importlib.import_module('bgc_md2.models.ElMasri2013AgricForMeteorol.source')
mvs = mod.mvs

mvs.graph()

mvs.computable_mvar_names

mvs.get_InputTuple()
mvs.get_VegetationCarbonInputScalar()
mvs.get_VegetationCarbonInputPartitioningTuple()
mvs.get_VegetationCarbonInputTuple()
for var in mvs.computable_mvar_types():
    mvs.render(var)

mvs.computable_mvar_types()

mvs.get_OutFluxesBySymbol()

importlib.invalidate_caches()
mod2 = importlib.import_module('bgc_md2.models.Haverd2016Biogeosciences.source')
mvs2 = mod2.mvs

mvs2.computable_mvar_names

mvs2.get_VegetationCarbonInputPartitioningTuple()


