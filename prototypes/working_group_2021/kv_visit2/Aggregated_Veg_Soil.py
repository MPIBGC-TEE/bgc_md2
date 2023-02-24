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
# adjust the output to full width
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# make changes to imported files immidiately available 
# avoiding the need to reload (in most cases)
# %load_ext autoreload
# %autoreload 2

from pathlib import Path
import json 
from sympy import  Symbol, Function, simplify, Function, Symbol, diff, simplify, exp
from sympy.core.function import UndefinedFunction
import numpy as np
import matplotlib.pyplot as plt
from ComputabilityGraphs.CMTVS import CMTVS
import CompartmentalSystems.helpers_reservoir as hr
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c
import bgc_md2.helper as h
import bgc_md2.display_helpers as dh

# imports from new files 
import sys
sys.path.insert(0,'..')
from source import mvs 
import source_1 as s1
import model_specific_helpers_2 as msh
from general_helpers import day_2_month_index, make_B_u_funcs_2 
import general_helpers as gh

import bgc_md2.resolve.computers as bgc_c
from subs_1 import subs_dict
# +
#mvs.computable_mvar_types(),mvs.provided_mvar_types

# +
#mvs_abstract.computable_mvar_types(),mvs.provided_mvar_types
# -

mvs.get_AggregatedVegetation2SoilCarbonFlux()

mvs.get_AggregatedSoil2VegetationCarbonFlux()

simplify(mvs.get_AggregatedVegetationCarbonInFlux())

simplify(s1.mvs.get_AggregatedVegetationCarbonOutFlux())

C__Veg,C__Soil=map(Symbol,["C__Veg","C__Soil"])
mvs_vs=CMTVS(
    {
        mvs.get_TimeSymbol(),
        StateVariableTuple([C__Veg,C__Soil]),
        InFluxesBySymbol({C__Veg: mvs.get_AggregatedVegetationCarbonInFlux()}),
        OutFluxesBySymbol(
            {
                C__Veg: mvs.get_AggregatedVegetationCarbonOutFlux(),
                C__Soil: mvs.get_AggregatedSoilCarbonOutFlux()
            }
        ),
        InternalFluxesBySymbol({
            (C__Veg,C__Soil): mvs.get_AggregatedVegetation2SoilCarbonFlux(),
            (C__Soil, C__Veg): mvs.get_AggregatedSoil2VegetationCarbonFlux()
        })
    },
    computers=h.module_computers(bgc_c)
)

h.compartmental_graph(mvs_vs)

#luo_tau=mvs.get_LuoTau()
beta=s1.mvs.get_CarbonInputPartitioningTuple()
simplify(beta)

luo_rt=s1.mvs.get_LuoRT()
luo_rt.free_symbols
              

# +
#def deriv(cm,temp):
#    T=temp.args[0]
#    if isinstance(T,Symbol):
#        return diff(cm,T)
#    elif isinstance(T, UndefinedFunction):
#        dv=list(cm.atoms(T))[0]
#        return diff(cm,dv)
#    else:
#        raise(Exception("The Temperature has to be either of type Symbol or an UndefinedFunction"))
    
temp=s1.mvs.get_Temperature()
temp
T=temp.args[0]
d=simplify(diff(luo_rt,T))
# -
xi=subs_dict[s1.xi]
xi
dss=d.subs(subs_dict)

# +
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 
msh.download_my_TRENDY_output(conf_dict)

#     # Read NetCDF data  ******************************************************************************************************************************
svs,dvs=msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
svs_0=msh.Observables(*map(lambda v: v[0],svs))

cpa = msh.Constants(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 npp_0=dvs.npp[0],# * 86400,   # kg/m2/s kg/m2/day
 #xi_0=dvs.xi[0],
 rh_0=svs_0.rh,# * 86400,   # kg/m2/s kg/m2/day
 #ra_0=svs_0.ra,# * 86400,   # kg/m2/s kg/m2/day
 #r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
 #r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
 number_of_months=len(svs.rh)
 #number_of_months=120 # for testing and tuning mcmc
)

# provide an inital guess for the paramters to be estimated by the data assimilation
# this is the result of the elaborate procedures 
epa_0=msh.EstimatedParameters(
    beta_leaf=0.6,
    beta_wood=0.25,
    T_0=2,
    E=6.5,
    KM=10,
    #env_modifier=1,
    r_C_leaf_litter_rh=0.0004151100041511,
    r_C_wood_litter_rh=0.00012453300124533,
    r_C_root_litter_rh=0.000122042341220423,
    r_C_soil_fast_rh=0.00015220700152207,
    #r_C_soil_slow_rh=2.73972602739726e-05,
    r_C_soil_slow_rh=3.73972602739726e-05,
    #r_C_soil_passive_rh=7.82778864970646e-06,
    r_C_soil_passive_rh=8.82778864970646e-06,
    r_C_leaf_2_C_leaf_litter=0.00833333333333333,
    r_C_wood_2_C_wood_litter=9.1324200913242e-05,
    r_C_root_2_C_root_litter=0.00012453300124533,
    r_C_leaf_litter_2_C_soil_fast=0.000340390203403902,
    r_C_leaf_litter_2_C_soil_slow=5.8115400581154e-05,
    r_C_leaf_litter_2_C_soil_passive=1.6604400166044e-05,
    r_C_wood_litter_2_C_soil_fast=7.4719800747198e-05,
    r_C_wood_litter_2_C_soil_slow=2.98879202988792e-05,
    r_C_wood_litter_2_C_soil_passive=1.99252801992528e-05,
    r_C_root_litter_2_C_soil_fast=7.4719800747198e-05,
    r_C_root_litter_2_C_soil_slow=3.48692403486924e-05,
    r_C_root_litter_2_C_soil_passive=1.74346201743462e-05,
    C_leaf_0=0.051828761170322826,
    C_wood_0=1.970572690329994,
    C_leaf_litter_0=0.1202311902470766,
    C_wood_litter_0=0.2225433197876749,
    C_soil_fast_0=1.7309510511856925,
    C_soil_slow_0=2.4435101360092473
)



#it_sym_trace = msh.make_traceability_iterator(mvs,dvs,cpa,epa_0)#,epa_opt)
X0=msh.numeric_X_0(mvs,dvs,cpa,epa_0)

# +
epa=epa_0
delta_t_val=15
func_dict=msh.make_func_dict(dvs,cpa,epa_0)
traced_expressions={
    'AggregatedVegetation2SoilCarbonFlux': mvs.get_AggregatedVegetation2SoilCarbonFlux(),
    'dss': dss,
}
apa = {**cpa._asdict(), **epa._asdict()}
par_dict = gh.make_param_dict(mvs, cpa, epa)

B_func, I_func = make_B_u_funcs_2(mvs, par_dict, func_dict, delta_t_val)

# we produce functions of f(it,x_a,...,x_z) to compute the tracable expressions
#traced_funcs = {
#    key: gh.numfunc(
#        expr_cont=val,
#        mvs=mvs,
#        delta_t_val=delta_t_val,
#        par_dict=par_dict,
#        func_dict=func_dict
#    )
#    for key,val in traced_expressions.items()
#}

# -

itr = gh.traceability_iterator(
    X0,
    func_dict, 
    mvs,
    dvs,
    cpa,
    epa_0,
    delta_t_val=delta_t_val,
    traced_expressions=traced_expressions
)
vals=itr[0:150]

import matplotlib.pyplot as plt
n=len(mvs.get_StateVariableTuple())

vals._fields

vals.AggregatedVegetation2SoilCarbonFlux

vals.dss

from string import Template
svt=mvs.get_StateVariableTuple()
import matplotlib.pyplot as plt
fig=plt.figure(figsize=(16,32))
axs=fig.subplots(2)
for i in range(len(svt)):
    axs[0].plot(
        vals.dss[:,i],
        #label=Template('$$ \partial X_c,{$sv} / \partial TAS  $$').substitute
        label=Template('$$ \partial RT,{$sv} / \partial TAS  $$').substitute(sv=svt[i])
    )
axs[0].legend()    
axs[0].set_title("Temperature sensitivity of RT components")
axs[1].plot(vals.AggregatedVegetation2SoilCarbonFlux)
axs[1].set_title("Aggregated fluxes from Vegetation to Soil pools")
