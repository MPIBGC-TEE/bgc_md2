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
from bgc_md2.resolve.mvars import VegetationCarbonInputScalar,VegetationCarbonInputPartitioningTuple

importlib.invalidate_caches()
mod = importlib.import_module('bgc_md2.models.ElMasri2013AgricForMeteorol.source')
mvs = mod.mvs

# +
#mvs.graph()
# -

mvs.get_InputTuple()

mvs.computable_mvar_types()

mvs.get_OutFluxesBySymbol()

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
# -

importlib.invalidate_caches()
mod2 = importlib.import_module('bgc_md2.models.Haverd2016Biogeosciences.source')
mvs2 = mod2.mvs

mvs2.computable_mvar_names

b = mvs2.get_VegetationCarbonInputPartitioningTuple()

from sympy import simplify
simplify(b)

b

mvs2.get_VegetationCarbonInputScalar()

mvs2.get_CompartmentalMatrix()


