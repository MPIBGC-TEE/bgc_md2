# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# # yy_cable_translation_to_fluxrates

from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
import bgc_md2.helper as h
import importlib 

importlib.invalidate_caches()
mod = importlib.import_module('bgc_md2.models.cable_yuanyuan.source')
mvs = mod.mvs

h.compartmental_graph(mvs)

#mvs.computable_mvar_types()


mvs.get_CompartmentalMatrix()


# This matrix is equivalent to the following fluxes:

# +
def st(var):
    str= var.__str__()
    ps=str.split("_")
    p1=ps[0]
    rest=ps[1:]
    return p1+"_{" + "".join(rest) +"}"

from IPython.display import Math
for t,v in mvs.get_InternalFluxesBySymbol().items():
    donor, receiver = t
    display(Math(
        "F_{"+st(donor)+"2"+st(receiver)+"} = "+ st(v) 
    ))
    #display(Math( 
    #    "= k_{"+st(donor)+","+st(receiver)+"}" + st(donor) 
    #))
# -

# the outfluxes we choose the 0 to express the exterior 
for donor,v in mvs.get_OutFluxesBySymbol().items():
   display(Math(
       "F_{"+st(donor)+"}="+st(v)
   ))

# Now we want to create a matrix in terms of the flux rates of the single fluxes from pool to pool or to the exterior.
# We write every flux in the form $F_{x2y}=r_{x2y} x$ where x is the donor pool content $r_{x2y}$ is the flux rate. 

# +
from sympy import var, Symbol
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers    


from bgc_md2.resolve.mvars import (
        OutFluxesBySymbol,
        InternalFluxesBySymbol
)
mvs_F = CMTVS(
    {
        mvs.get_InFluxesBySymbol(),
        mvs.get_TimeSymbol(),
        mvs.get_StateVariableTuple(),
        OutFluxesBySymbol(
            {
                dp :
                Symbol("r_"+ str(dp))*dp
                for dp in mvs.get_OutFluxesBySymbol().keys()
            }
        ),        
        InternalFluxesBySymbol(
            {
                (dp,rp) : Symbol("r_"+ str(dp)+ "_2_" + str(rp))*dp
                for dp,rp in mvs.get_InternalFluxesBySymbol().keys()
            }
        )
    },
    bgc_md2_computers()
)

# -

# To get rid of the $f_{ij}$ and $k_i$ we have to express them in terms of the new variables $r_{i2j}$ 

# +
from sympy import simplify, Symbol
internal_fluxes_by_symbol = mvs_F.get_InternalFluxesBySymbol()
out_fluxes_by_symbol = mvs_F.get_OutFluxesBySymbol()
svt = mvs_F.get_StateVariableTuple()

def release_rate_by_symbol(donor):
    external = out_fluxes_by_symbol[donor] if donor in out_fluxes_by_symbol.keys() else 0
    internal = sum([flux for (d,r),flux in internal_fluxes_by_symbol.items() if d == donor])
    decomp_flux = external + internal                            
    return simplify(decomp_flux/donor)
    
release_rate_by_symbol(svt[0])


# +
def transfer_rate_by_symbol(tup):
    flux = internal_fluxes_by_symbol[tup] if tup in internal_fluxes_by_symbol.keys() else 0
    donor,receiver = tup
    return flux/donor/release_rate_by_symbol(donor)

transfer_rate_by_symbol((svt[0],svt[3]))


# -

# Now we can create a substitution dictionary to replace all the $f_{i,j}$ and $k_i$

# +
def cut(pool_name):
    return str(pool_name)[2:]
subsdict_f= {
    Symbol("f_"+cut(d)+"2"+cut(r)) : transfer_rate_by_symbol((d,r)) 
    for d,r in mvs.get_InternalFluxesBySymbol().keys()
}
subsdict_k= {
    Symbol("k_"+cut(d)) : release_rate_by_symbol(d) 
    for d in mvs.get_StateVariableTuple()
    
}
subsdict = {**subsdict_f,**subsdict_k}
# -

subsdict

mvs.get_CompartmentalMatrix().subs(subsdict) 

out_fs = {d: f.subs(subsdict) for d,f in mvs.get_OutFluxesBySymbol().items()}

internal_fs = {t: f.subs(subsdict) for t,f in mvs.get_InternalFluxesBySymbol().items()}

 #new_mvs = CMTVS()


