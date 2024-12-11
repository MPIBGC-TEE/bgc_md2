# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
from time import time as now
from unittest import TestCase, skip
from ComputabilityGraphs.rec_graph_helpers import fast_graph
from ComputabilityGraphs.helpers import all_mvars
from ComputabilityGraphs.TypeSet import TypeSet

import bgc_md2.helper as h 
from bgc_md2.resolve.mvars import *


# -

mn='Aneesh_SDGVM'
mvs = h.CMTVS_from_model_name(mn)
mvars = mvs.computable_mvar_types()
mvars


mvs.get_NumericMeanAgeSolutionArray()





for var in mvars:
    print("########################################")
    print(str(var.__name__))
    #print(mvs._get_single_value(var))
    mvs._get_single_value(var)
    #print(mvs._get_single_value_by_depgraph(var))


