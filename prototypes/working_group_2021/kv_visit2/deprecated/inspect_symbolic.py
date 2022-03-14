# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # VISIT_Kostia

# +
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important;}</style>"))

import bgc_md2.helper as h
import importlib 
# -

importlib.invalidate_caches()
mod_org = importlib.import_module('bgc_md2.models.VISIT_Kostia.source_minus_1')
mod = importlib.import_module('bgc_md2.models.VISIT_Kostia.source')

mvs_org = mod_org.mvs
mvs = mod.mvs

h.compartmental_graph(mvs)

mvs_org.get_CompartmentalMatrix()

from IPython.display import Math
from sympy import latex
for k,v in mvs_org.get_InFluxesBySymbol().items():
    display( Math( "In_{"+str(k)+"} = " + latex(v) ))

# +
import sys
sys.path.insert(0,'..')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import json 

from model_specific_helpers import (
    get_example_site_vars,
    make_param_filter_func,
    make_weighted_cost_func,
    make_param2res_sym,
    UnEstimatedParameters,
    EstimatedParameters,
    Observables
)

from general_helpers import (
        make_uniform_proposer,
        make_multivariate_normal_proposer,
        #mcmc,
        make_feng_cost_func,
        mcmc,
        adaptive_mcmc,
        plot_solutions
)

with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

dataPath = Path(conf_dict['dataPath'])
# -

dataPath


# +
# We factor out xi automatically
A=mvs_org.get_CompartmentalMatrix()
from sympy import (factor,gcd,var,Piecewise,
gcd_list,gcd_terms,gcdex,factor_terms,prod,parallel_poly_from_expr,exp,Function)

#gcd(list(A))
var("x y z t k_x")
ex1 = (x**2+x)*k_x*Piecewise((t-1,t<0),(t**2,t>=0))
ex2 = (y**2+x**2+x)*3*Piecewise((t-1,t<0),(t**2,t>=0))
ex3 = (y**2+z**2+x)*3*Piecewise((t-1,t<0),(t**2,t>=0))
gcd_list([ex1,ex2,ex3],gens=(x,y,z))
# -

res1=ex1.cancel()
res2=ex2.cancel()

[ a for a in res1.args if len(a.free_symbols.intersection({x}))<1]


{x}.intersection


def xi(el):
    fs=[ f for f in el.factor().args if len(f.free_symbols.intersection(mvs.get_TimeSymbol()))>0]
    #return prod(fs) 
    return fs 


srm=mvs_org.get_SmoothReservoirModel()

_,T,N,x,u=srm.xi_T_N_u_representation()

N

var('xi')
simplify(N.diagonal()[3]/xi)

var ("E T0 KM t")
tsl=Function('tsl')
mrso=Function('mrso')
from sympy import simplify
xi= Piecewise(
    (exp(E * (1 / (10 - T0) - 1 / (tsl(t) - T0))) * mrso(t) / (KM + mrso(t)),tsl(t)>T0),
    #(0,tsl(t)<T0)
    (0,True) #catch all
    )

[simplify(el) for el in N.diagonal()]


