# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # ElMasri2013AgricForMeteorol

import bgc_md2.helper as h
import importlib 
from bgc_md2.resolve.mvars import VegetationCarbonInputScalar,VegetationCarbonInputPartitioningTuple

importlib.invalidate_caches()
mod = importlib.import_module('bgc_md2.models.ElMasri2013AgricForMeteorol.source')
mvs = mod.mvs

# +
#mvs.graph()
# -

mvs.get_InputTuple()

mvs.computable_mvar_types()

b=mvs.get_VegetationCarbonInputScalar()
u=mvs.get_VegetationCarbonInputPartitioningTuple()


from sympy import simplify
simplify(b)

it=mvs.get_VegetationCarbonInputTuple()
#mvs.get_InFluxesBySymbol()
#VegetationCarbonInputPartitioningTuple([tc/b for tc in it])

# +

VegetationCarbonInputScalar(sum(it))
# -

[t for t in mvs.provided_mvar_types]

# +
#for var in mvs.computable_mvar_types():
#    mvs.render(var)
