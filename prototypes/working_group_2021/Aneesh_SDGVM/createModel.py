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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # SDGVM model

# +
import sys
sys.path.insert(0,'..') # necessary to import general_helpers
from sympy import Symbol, Function 
from sympy import var, exp
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

#var("""
#C_leaf C_wood C_root C_abvstrlit C_abvmetlit C_belowstrlit C_belowmetlit C_surface_microbe C_soil_microbe _slowsom C_passsom
#k_C_leaf k_C_wood k_C_root k_C_abvmetlit k_C_abvstrlit k_C_belowstrlit k_C_belowmetlit 
#k_C_surface_microbe k_C_soil_microbe k_C_slowsom k_C_passsom 
#beta1 beta2 beta3
#leached ls_aboveground 
#ls_belowground f41 f51 f52 f63 f73 f84 f104 f85 f96 f106 f97 f910 f911 f108 f119 f109 f1110 
#Ca clay silt_clay leachedwtr30""")
# we organize the new symbols a bit better by putting them in a dictionary
# together with some description that we can use later to display some metainformation
sym_dict = {
        'C_leaf': 'content of leaf pool',
        'C_wood': 'content of wood pool',
        'C_root': 'content of root pool',
        'C_abvstrlit': 'content of aboveground structural litter pool',
        'C_abvmetlit': 'content of aboveground metabolic litter pool',
        'C_belowstrlit': 'content of belowground structural litter pool',
        'C_belowmetlit': 'content of belowground metabolic litter pool',
        'C_surface_microbe': 'content of surface microbe pool',
        'C_soil_microbe': 'content of soil microbial pool',
        'C_slowsom': 'content of the slow soil pool',
        'C_passsom': 'content of the passive soil pool',
        'beta_leaf': 'NPP partitioning to leaf', 
        'beta_wood': 'NPP partitioning to wood',
        'beta_root': 'NPP partitioning to root',
        'k_C_leaf': 'Turn over rate of leaf pool',
        'k_C_wood':'Turn over rate of wood pool',
        'k_C_root': 'Turn over rate of root pool',
        'k_C_abvmetlit':'Turn over rate of aboveground metabolic litter pool',
        'k_C_abvstrlit' :'Turn over rate of aboveground structural litter pool',
        'k_C_belowstrlit' :'Turn over rate of belowground structural litter pool',
        'k_C_belowmetlit':'Turn over rate of belowground metabolic litter pool',
        'k_C_surface_microbe' :'Turn over rate of surface microbe pool',
        'k_C_soil_microbe' :'Turn over rate of soil microbial pool',
        'k_C_slowsom' :'Turn over rate of slow soil pool',
        'k_C_passsom':'Turn over rate of passive soil pool',
        'ls_aboveground': 'Fraction of structural C that is lignin in aboveground',
        'ls_belowground': 'Fraction of structural C that is lignin in belowground',
        'f_leaf2abvstrlit': 'Transfer coefficient of C from leaf to aboveground structural litter',
        'f_leaf2abvmetlit':'Transfer coefficient of C from leaf to aboveground metabolic litter',
        'f_wood2abvstrlit':'Transfer coefficient of C from wood to aboveground structural litter',
        'f_wood2abvmetlit':'Transfer coefficient of C from wood to aboveground metabolic litter',
        'f_root2belowstrlit':'Transfer coefficient of C from root to belowground structural litter',
        'f_root2belowmetlit':'Transfer coefficient of C from root to belowground metabolic litter',
        'f_abvstrlit2surface_microbe':'Transfer coefficient of C from aboveground structural litter to surface microbe',
        'f_abvstrlit2slowsom':'Transfer coefficient of C from aboveground structural litter to slow soil',
        'f_abvmetlit2surface_microbe':'Transfer coefficient of C from aboveground metabolic litter to surface microbe',
        'f_belowstrlit2soil_microbe':'Transfer coefficient of C from belowground structural litter to soil microbe',
        'f_belowstrlit2slowsom':'Transfer coefficient of C from belowground structural litter to slow soil',
        'f_belowmetlit2soil_microbe':'Transfer coefficient of C from belowground metabolic litter to soil microbe',
        'f_slowsom2soil_microbe':'Transfer coefficient of C from soil soil to soil microbe',
        'f_passsom2soil_microbe':'Transfer coefficient of C from passive soil to soil microbe',
        'f_surface_microbe2slowsom':'Transfer coefficient of C from surface microbe to slow soil',
        'f_soil_microbe2passsom':'Transfer coefficient of C from soil microbe to passive soil',
        'f_soil_microbe2slowsom':'Transfer coefficient of C from soil microbe to slow soil',
        'f_slowsom2passsom':'Transfer coefficient of C from slow soil to passive soil',
        'Ca': 'Need to find (unknown)',
        'clay': 'Clay content (fraction)',
        'silt_clay': 'Silt plus clay content (fraction)',
        'leachedwtr30': 'H20 leached below 30cm (cm/m)',
        'r_C_leaf2abvstrlit':'',
        'r_C_leaf2abvmetlit':'',
        'r_C_wood2abvstrlit':'', 
        'r_C_wood2abvmetlit' :'',
        'r_C_root2belowstrlit':'',  
        'r_C_root2belowmetlit':'',
        'r_C_abvstrlit2surface_microbe':'',
        'r_C_abvmetlit2surface_microbe':'', 
        'r_C_abvstrlit2slowsom':'',
        'r_C_belowstrlit2soil_microbe':'',  
        'r_C_belowmetlit2soil_microbe':'',
        'r_C_belowstrlit2slowsom':'' ,
        'r_C_surface_microbe2slowsom':'',
        'r_C_soil_microbe2slowsom' :'',''
        'r_C_slowsom2soil_microbe':'',
        'r_C_soil_microbe2passsom':'',
        'r_C_slowsom2passsom':'',
        'r_C_passsom2soil_microbe':'',
        'r_C_leached': ''
    
}

# for the moment only use the 
var(list(sym_dict.keys()))

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    'xi': 'a scalar function of temperature and moisture and thereby ultimately of time',
    'NPP': '',
}

for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

#ls_aboveground = 0.12
#ls_belowground = 0.35

#f_leaf2abvstrlit = 0.91986 + 0.00324* Ca
#f_leaf2abvmetlit = 1 - f_leaf2abvstrlit
#f_wood2abvstrlit = f_leaf2abvstrlit
#f_wood2abvmetlit = 1 - f_wood2abvstrlit
#f_root2belowstrlit = f_leaf2abvstrlit
#f_root2belowmetlit = 1 - f_root2belowstrlit
#f_abvstrlit2surface_microbe = (1 - exp(-3*ls_aboveground))*0.4
#f_abvstrlit2slowsom = exp(-3*ls_aboveground) * 0.7
#f_abvmetlit2surface_microbe = 0.4
#f_belowstrlit2soil_microbe = exp(1 - (-3*ls_belowground))*0.45
#f_belowstrlit2slowsom = exp(-3*ls_aboveground) * 0.7
#f_belowmetlit2soil_microbe = 0.45
#f_slowsom2soil_microbe = 0.447 + 0.009*clay
#f_passsom2soil_microbe = 0.45
#f_surface_microbe2slowsom = 0.4
#f_soil_microbe2passsom = 0.003 + 0.032 * clay
#f_soil_microbe2slowsom = 1 - f_soil_microbe2passsom - (leachedwtr30)/18 * (0.01 + 0.04* (1- silt_clay)) - 0.85 - 0.68*silt_clay
#f_slowsom2passsom = 0.003 + 0.009*clay
#leached = (leachedwtr30)/18 * (0.01 + 0.04* (1- silt_clay))

#r_C_leaf2abvstrlit=k_C_leaf *f_leaf2abvstrlit
#r_leaf2abvmetlit = k_C_leaf *f_leaf2abvmetlit
#r_wood2abvstrlit = k_C_wood *f_wood2abvstrlit
#r_wood2abvmetlit = k_C_wood *f_wood2abvmetlit
#r_root2belowstrlit = k_C_root * f_root2belowstrlit
#r_root2belowmetlit = k_C_root * f_root2belowmetlit
#r_abvstrlit2surface_microbe = k_C_abvstrlit *f_abvstrlit2surface_microbe
#r_abvmetlit2surface_microbe = k_C_abvmetlit *f_abvmetlit2surface_microbe
#r_abvstrlit2slowsom = k_C_abvstrlit*f_abvstrlit2slowsom
#r_belowstrlit2soil_microbe = k_C_belowstrlit * f_belowstrlit2soil_microbe
#r_belowmetlit2soil_microbe = k_C_belowmetlit  * f_belowmetlit2soil_microbe
#r_belowstrlit2slowsom = k_C_belowstrlit *f_belowstrlit2slowsom
#r_surface_microbe2slowsom = k_C_surface_microbe*f_surface_microbe2slowsom
#r_soil_microbe2slowsom = k_C_soil_microbe *f_soil_microbe2slowsom
#r_slowsom2soil_microbe = k_C_slowsom *f_slowsom2soil_microbe
#r_soil_microbe2passsom = k_C_soil_microbe*f_soil_microbe2passsom
#r_slowsom2passsom = k_C_slowsom*f_slowsom2passsom
#r_passsom2soil_microbe =k_C_passsom * f_passsom2soil_microbe
#r_leached = leached = (leachedwtr30)/18 * (0.01 + 0.04* (1- silt_clay))
t=TimeSymbol("t")
mvs = CMTVS(
    {
        StateVariableTuple((C_leaf, 
                            C_wood, 
                            C_root, 
                            C_abvstrlit,
                            C_abvmetlit, 
                            C_belowstrlit, 
                            C_belowmetlit, 
                            C_surface_microbe, 
                            C_soil_microbe, 
                            C_slowsom, 
                            C_passsom)),
        t,
        InFluxesBySymbol(
            {
                C_leaf: NPP(t)* beta_leaf, 
                C_wood: NPP(t)* beta_wood, 
                C_root:NPP(t)*  beta_root
            }),
        OutFluxesBySymbol(
            {
                C_soil_microbe: r_C_leached*C_soil_microbe
            }
        ),
        InternalFluxesBySymbol(
            {
                                (C_leaf, C_abvstrlit): r_C_leaf2abvstrlit* C_leaf, 
                                (C_leaf, C_abvmetlit): r_C_leaf2abvmetlit* C_leaf, 
                                (C_wood, C_abvstrlit): r_C_wood2abvstrlit* C_wood, 
                                (C_wood, C_abvmetlit): r_C_wood2abvmetlit* C_wood, 
                                (C_root, C_belowstrlit): r_C_root2belowstrlit* C_root, 
                                (C_root, C_belowmetlit): r_C_root2belowmetlit * C_root, 
                                (C_abvstrlit , C_surface_microbe ): r_C_abvstrlit2surface_microbe* C_abvstrlit, 
                                (C_abvstrlit , C_slowsom ): r_C_abvstrlit2slowsom*C_abvstrlit,
                                (C_abvmetlit, C_surface_microbe ): r_C_abvmetlit2surface_microbe* C_abvmetlit, 
                                (C_belowstrlit, C_soil_microbe): r_C_belowstrlit2soil_microbe*C_belowstrlit, 
                                (C_belowmetlit , C_soil_microbe): r_C_belowmetlit2soil_microbe*C_belowmetlit, 
                                (C_belowstrlit, C_slowsom): r_C_belowstrlit2slowsom* C_belowstrlit, 
                                (C_surface_microbe , C_slowsom): r_C_surface_microbe2slowsom*C_surface_microbe, 
                                (C_soil_microbe, C_slowsom): r_C_soil_microbe2slowsom*C_soil_microbe, 
                                (C_slowsom, C_soil_microbe): r_C_slowsom2soil_microbe*C_slowsom, 
                                (C_soil_microbe, C_passsom): r_C_soil_microbe2passsom*C_soil_microbe, 
                                (C_slowsom , C_passsom): r_C_slowsom2passsom*C_slowsom,
                               (C_passsom, C_soil_microbe): r_C_passsom2soil_microbe * C_passsom
            }
        )
    },

    computers=module_computers(bgc_c)
)
# -

h.compartmental_graph(mvs)

dh.mass_balance_equation(mvs)

# Nothing has changed in the model description but be have some more symbols to work with.
# We could type `C_leaf` somewhere without an error since it is now known as a variable.
# In the next step we replace all occurences of `vl` by `C_leaf` 

mvs.get_StateVariableTuple()

mvs.get_CompartmentalMatrix()


# +
srm = mvs.get_SmoothReservoirModel()

# The matrices T and N in the function refer to what in Yiqi Luo's group is usually called A and K
# and xi is a single scalar (not the diagonal matrix we are looking for here)
# The function has documentation which you can see by typing the following
# # ?srm.xi_T_N_u_representation
_,A,N,_,_=srm.xi_T_N_u_representation(factor_out_xi=False) 
# -


mvs.computable_mvar_types()

from sympy import diag
xi_d=diag([1,1,1]+[xi(t) for i in range(8)],unpack=True)
xi_d

# We can go on and decompose N =\xi K -> K=\xi^{-1}N
K=xi_d.inv()*N
K
# we now have the (symbolic) ingredients for the tracebility analysis.
#xi_d,K,A

# +
import sys
sys.path.insert(0,'..') # necessary to import general_helpers
from general_helpers import download_TRENDY_output
import json 
from pathlib import Path
from collections import namedtuple 

with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

# we will use the trendy output names directly in other parts of the output
Observables = namedtuple(
    'Observables',
    ["cLitter", "cSoil", "cVeg", "cRoot", "rh"]
)
Drivers=namedtuple(
    "Drivers",
    ["npp"]
)    
    
#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output():
    download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['SDGVM'],
        variables = Observables._fields + Drivers._fields
    )
#call it to test that the download works the data
#download_my_TRENDY_output()


# -

# Before we build a function to load the data lets look at it to get an idea.
#

# +
# this is just an analysis you dont need this code later
dataPath=Path(conf_dict['dataPath'])
import netCDF4 as nc
import numpy as np
from pathlib import Path
import json 

#  the reference point, values and number of recorded times in the data streams differ
# there are two gropus
ds_cLitter=nc.Dataset(dataPath.joinpath('SDGVM_S2_cLitter.nc'))
ds_cVeg=nc.Dataset(dataPath.joinpath('SDGVM_S2_cVeg.nc'))
print(ds_cLitter)
for ds in [ds_cLitter,ds_cVeg]:
    print(ds['time'].units)
    print(ds['time'].shape)

ds_rh=nc.Dataset(dataPath.joinpath('SDGVM_S2_rh.nc'))
ds_npp=nc.Dataset(dataPath.joinpath('SDGVM_S2_npp.nc'))
for ds in [ds_rh,ds_npp]:
    print(ds['time'].units)
    print(ds['time'].shape)

# refer them to a single reference point 01/01/1900 by substractin 200y in h's
h_from_1900=ds_cLitter.variables['time'][:]-200*12*30*24


# find the index after which we are after 01/01 1900
ind_start=0
while h_from_1900[ind_start]<0:
    ind_start+=1
print('ind_start={}'.format(ind_start))
#ind_start is 200 confirming that we have 1 value per year and have to remove the first 200 values
c_h_from_1900_after_1900=h_from_1900[ind_start:] 
c_h_from_1900_after_1900.shape

#ds=nc.Dataset(dataPath.joinpath('SDGVM_S2_cLitter.nc'))
#ds=nc.Dataset(dataPath.joinpath('SDGVM_S2_rh.nc'))
#ds.variables['rh']#[:,-56,26][20]
#print(ds)

# +
def get_example_site_vars(dataPath):
    (
        C_litter,
        C_soil,
        C_veg,
        C_root,
        npp,
        rh
    )= get_variables_from_files(dataPath)
    # pick up 1 site   wombat state forest
    s = slice(None, None, None)  # this is the same as :
    #t = s, 50, 33  # [t] = [:,49,325]
    t = s, -25, 16 
    npp = npp[t] * 86400   # kg/m2/s kg/m2/day;
    rh = rh[t]*86400  # per s to per day
    #tsl_mean = np.mean(tsl, axis=1)  # average soil temperature at different depth
    (
        C_litter,
        C_soil,
        C_veg,
        C_root
    ) = map(
        lambda var: var[t],
        (
            C_litter,
            C_soil,
            C_veg,
            C_root
        )
    )
    
    # the dataset is not uniform 
    # - npp and rh start at 16-01-1900 end at 03-28-2018 and are recorded every 30 days
    # - C_litter, C_soil, C_veg, C_root start at 06-30-1700 end at 01-12-2014 
    #
    # To make them uniform we will:
    # 1. Make C_litter, C_soil, C_veg, C_root start at time 1900 (cutting of the first 200y)
    # 2. Adapt the resolution of rh to yearly by averaging over the monthly values and
    #    also cutting them short to 2014
    # 3. We will cut short npp to 2014 bu will NOT change the resolution (since it is a driver)
    #    and does not affect the data assimilation functions 
    print(C_litter.shape)
    (
        C_litter,
        C_soil,
        C_veg,
        C_root
    ) = map(
        lambda var: var[204:],
        (
            C_litter,
            C_soil,
            C_veg,
            C_root
        )
    )
    
    return (Observables(C_litter, C_soil, C_veg, C_root, rh), Drivers(npp))

with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

dataPath = Path(conf_dict['dataPath'])

obs, dr =get_example_site_vars(dataPath)
obs.cLitter
# -

ds_cVeg.variables['cVeg'][ind_start:,-25,16], ds_cVeg.variables['cVeg'][ind_start:,-25,16].shape


# this is the important function that we will use later
def get_example_site_vars(dataPath):
    # According to the netcdf metadata the datasets are not uniform 
    # - npp and rh start at 360h (15 days) after 01-01-1900 and are recorded every 30 days
    #   these are interpreted as mid-monthly 
    # - C_litter, C_soil, C_veg, C_root start at 4320 h = 180 days = 6 months after 01-01-1700 
    #   These can be seen at midyearly values since there are 6 (midmonth) times of the npp and rh measurements after the last (midyear)
    #   measurement of C_litter, C_soil, C_veg, C_root
    
    # To combine these streams into a consistent array of observations we will:
    # 1. Make C_litter, C_soil, C_veg, C_root refer to hours after 01/01/1900 (as npp and rh) 
    # 
    # 2. cut away everything before 1900 from them (cutting of the first 200y)
    #    
    # Note:
    #    We will have to adapt the costfunction and param2res later to accommodate the different 
    #    resolution of the C_pool and rh observations.
    


    # 1.) 
    # pick one of the 1700 yearly example ds to get at the times 
    # convert time to refer to the same starting point (from h after 1700 to h after 1900) 
    hs_from_1900=nc.Dataset(dataPath.joinpath('SDGVM_S2_cLitter.nc')).variables['time'][:]-200*12*30*24
    
    #2.)
    # find the index after which we are after 01/01 1900
    ind_start = 200 
    c_h_from_1900_after_1900=hs_from_1900[ind_start:] 
    c_h_from_1900_after_1900.shape
    
    # pick up 1 site   wombat state forest for the spacial selection
    s_rh  = slice(None, None, None)  # this is the same as :
    s_c  = slice(ind_start, None, None)  # this is the same as ind_start:
    #t = s, 50, 33  # [t] = [:,49,325]
    loc=(-25,16)
    t_rh = s_rh,*loc
    t_c = s_c, *loc 
    print(t_c)
    
    # Read NetCDF data and slice out our site 
    arr_dict={
        **{
            vn:nc.Dataset(str(dataPath.joinpath(fn))).variables[vn][t_c] 
            for vn,fn in  {
                'cLitter': 'SDGVM_S2_cLitter.nc',
                'cSoil': 'SDGVM_S2_cSoil.nc',
                'cVeg': 'SDGVM_S2_cVeg.nc',
                'cRoot': 'SDGVM_S2_cRoot.nc',
            }.items()
        },
        **{
            vn:nc.Dataset(str(dataPath.joinpath(fn))).variables[vn][t_rh]*86400   # kg/m2/s kg/m2/day;
            for vn,fn in {
                'npp': 'SDGVM_S2_npp.nc',
                'rh': 'SDGVM_S2_rh.nc'
            }.items()
        }
    }
    
    return (
        Observables(*(arr_dict[k] for k in Observables._fields)),
        Drivers(*(arr_dict[k] for k in Drivers._fields))
    )
#    return (Observables(C_litter, C_soil, C_veg, C_root, rh), Drivers(npp))
#
#with Path('config.json').open(mode='r') as f:
#    conf_dict=json.load(f) 
#
#dataPath = Path(conf_dict['dataPath'])
#
#obs, dr =get_example_site_vars(dataPath)
#obs.cVeg.shape

# +
import sys 
sys.path.insert(0,'..')
from general_helpers import day_2_month_index
def NPP_fun(day ):
    return dr.npp[day_2_month_index(day)] 

func_dict={NPP: NPP_fun}
# -

svs,dvs=get_example_site_vars(dataPath=Path(conf_dict["dataPath"]))

svs.cLitter.shape

# +
from sympy import ImmutableMatrix
sv=mvs.get_StateVariableTuple()
n=len(sv)
# create new symbols for the f_{i,j}
for i in range(n):
    for j in range(n):
        if A[i,j]!=0 and i!=j:
            name="f_" + str(sv[j]) + "_2_" + str(sv[i])
            code="{0}=Symbol('{0}')".format(name)
            print(code)
            exec(code)
            
A_sym=ImmutableMatrix(
    n,n,
    lambda i,j:  -1 if i==j else (
        0 if A[i,j]==0 else Symbol("f_" + str(sv[j]) + "_2_" + str(sv[i]))
    )
)
A_sym

# +
# create new symbols for the k_{i}
for i in range(n):
    if K[i,i]!=0:
        name="k_{0}".format(sv[i])
        code="{0}=Symbol('{0}')".format(name)
        print(code)
        exec(code)
        
K_sym=ImmutableMatrix(
    n,n,
    lambda i,j: Symbol("k_" + str(sv[i])) if i==j else 0
)
K_sym
# -

M_sym=A_sym*K_sym
M_sym

import  CompartmentalSystems.helpers_reservoir as hr
hr.out_fluxes_by_symbol(sv,M_sym)

# +
# we create a dictionary for the outfluxrates (flux divided by dono pool content)
outflux_rates = {"r_"+str(key)+"_rh":value/key for key,value in hr.out_fluxes_by_symbol(sv,M_sym).items()}
internal_flux_rates = {"r_"+str(key[0])+"_2_"+str(key[1]):value/key[0] for key,value in hr.internal_fluxes_by_symbol(sv,M_sym).items()}
from copy import  deepcopy
all_rates=deepcopy(outflux_rates)
all_rates.update(internal_flux_rates)
all_rates

par_dict = {
     beta_leaf:0.15, 
     beta_wood:0.6, 
     beta_root:0.25,
     r_C_leaf2abvstrlit: 0.3,
     r_C_abvmetlit2surface_microbe:0.2,
     r_C_abvstrlit2slowsom:0.4,
     r_C_abvstrlit2surface_microbe:0.5,
     r_C_belowmetlit2soil_microbe:0.4,
     r_C_belowstrlit2slowsom:0.3,
     r_C_belowstrlit2soil_microbe:0.4,
     r_C_leached:0.1,
     r_C_leaf2abvmetlit:0.6,
     r_C_passsom2soil_microbe:0.2,
     r_C_root2belowmetlit:0.3,
     r_C_root2belowstrlit:0.4,
     r_C_slowsom2passsom:0.3,
     r_C_slowsom2soil_microbe:0.1,
     r_C_soil_microbe2passsom:0.05,
     r_C_soil_microbe2slowsom:0.1,
     r_C_surface_microbe2slowsom:0.02,
     r_C_wood2abvmetlit:0.4,
     r_C_wood2abvstrlit:0.6}

# +
from general_helpers import make_B_u_funcs_2, day_2_month_index
# check the numeric functions for B and u

def npp_func(day):
    month=day_2_month_index(day)
    return dvs.npp[month]

def xi_func(day):
    return 1.0 # preliminary fake for lack of better data... 

func_dict={
    'NPP':npp_func,
     'xi':xi_func
}
# for the next line to work the 
# two dictionaries par_dict and func_dict have to be complete.
# In sympy terms this means that the expressions for 
B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
# we create a numeric startvector for the compartmental system
# 
svs_0=Observables(*map(lambda v: v[0],svs))

X_0= np.array((
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cVeg/3,
    svs_0.cLitter/4,
    svs_0.cLitter/4,
    svs_0.cLitter/4,
    svs_0.cLitter/4,
    svs_0.cSoil/4,
    svs_0.cSoil/4,
    svs_0.cSoil/4,
    svs_0.cSoil/4
))#.reshape(11,)
# in general B and u are nonlinear and can depend on X, thats why we have to test them with index and X arguments
u_func(0,X_0),B_func(0,X_0)
# -

# to guard agains accidentally changed order we use a namedtuple again. Since B_func and u_func rely 
# on the correct ordering of the statevariables we build V dependent on this order 
StartVector=namedtuple(
    "StartVector",
    [str(v) for v in mvs.get_StateVariableTuple()]+
    ["rh"]
)
StartVector._fields

V_init= StartVector(
    C_leaf=svs_0.cVeg/3,
    C_wood=svs_0.cVeg/3,
    C_root=svs_0.cVeg/3,
    C_abvstrlit=svs_0.cLitter/4,
    C_abvmetlit=svs_0.cLitter/4,
    C_belowstrlit=svs_0.cLitter/4,
    C_belowmetlit=svs_0.cLitter/4,
    C_surface_microbe=svs_0.cSoil/4,
    C_soil_microbe=svs_0.cSoil/4,
    C_slowsom=svs_0.cSoil/4,
    C_passsom=svs_0.cSoil/4,
    rh=svs_0.rh        
)
V_init.__getattribute__("C_leaf")

# +
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
        TimeStepIterator2,
)
from copy import copy

def make_daily_iterator_sym(
        mvs,
        V_init: StartVector,
        par_dict,
        func_dict
    ):
    B_func, u_func = make_B_u_funcs_2(mvs,par_dict,func_dict)  
    sv=mvs.get_StateVariableTuple()
    n=len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct 
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr=np.array(
        [V_init.__getattribute__(str(v)) for v in sv]+
        #[V_init.ra,V_init.rh]
        [V_init.rh]
    ).reshape(n+1,1) #reshaping is neccessary for matmux
    
    def f(it,V):
        X = V[0:n]
        b = u_func(it,X)
        B = B_func(it,X)
        outfluxes = B @ X
        X_new = X + b + outfluxes
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
        #ra= -np.sum(outfluxes[0:3]) # check with  mvs.get_StateVariableTuple()[0:3]
        rh= -np.sum(outfluxes[3:n]) # check with  mvs.get_StateVariableTuple()[3:9]
        
        V_new = np.concatenate((X_new.reshape(n,1),np.array([rh]).reshape(1,1)), axis=0)
        return V_new
    
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


# -

V_init

np.array(V_init).shape

# +
# test the daily iterator
    
it_sym = make_daily_iterator_sym(
    mvs,
    V_init=V_init,
    par_dict=par_dict,
    func_dict=func_dict
)
# we will run the model for 15 steps
ns=15
res= np.zeros((ns,len(V_init)))
res_sym = copy(res)
for i in range(ns):
    res_sym[i,:]=it_sym.__next__().reshape(len(V_init),)
res_sym

# +
# As a safety measure we specify those parameters again as 'namedtuples', which are like a mixture of dictionaries and tuples
# They preserve order as numpy arrays which is great (and fast) for the numeric computations
# and they allow to access values by keys (like dictionaries) which makes it difficult to accidentally mix up values.

UnEstimatedParameters = namedtuple(
    "UnEstimatedParameters",
    [
        "cLitter_0",
        "cSoil_0",
        #"c_root_0",
        "cRoot_0",
        "cVeg_0",
        "npp_0",
        "rh_0",
        "number_of_months" # necessary to prepare the output in the correct lenght 
    ]
)
# Note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues and their length (number of months)
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 
# It is better to start with only a few

EstimatedParameters = namedtuple(
    "EstimatedParameters",[ 
         'beta_leaf', 
         'beta_wood', 
         #beta_root,
         'r_C_leaf2abvstrlit',
         'r_C_abvmetlit2surface_microbe',
         'r_C_abvstrlit2slowsom',
         'r_C_abvstrlit2surface_microbe',
         'r_C_belowmetlit2soil_microbe',
         'r_C_belowstrlit2slowsom',
         'r_C_belowstrlit2soil_microbe',
         'r_C_leached',
         'r_C_leaf2abvmetlit',
         'r_C_passsom2soil_microbe',
         'r_C_root2belowmetlit',
         'r_C_root2belowstrlit',
         'r_C_slowsom2passsom',
         'r_C_slowsom2soil_microbe',
         'r_C_soil_microbe2passsom',
         'r_C_soil_microbe2slowsom',
         'r_C_surface_microbe2slowsom',
         'r_C_wood2abvmetlit',
         'r_C_wood2abvstrlit',
        'C_leaf_0',
        #'C_root_0',
        'C_abvstrlit_0',
        'C_abvmetlit_0',
        'C_blwstrlit_0',
        'C_surfacemic_0',
        'C_soilmic_0',
        'C_slow_0'
    ]
)
# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues. 
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 
# -

EstimatedParameters._fields

UnEstimatedParameters._fields

cpa=UnEstimatedParameters(
 cVeg_0=svs_0.cVeg,
 cLitter_0=svs_0.cLitter,
 cSoil_0=svs_0.cSoil,
 cRoot_0 = svs_0.cRoot,
 npp_0=dvs.npp[0],
 rh_0=svs_0.rh,
 number_of_months=len(svs.rh)
)

cpa

epa_0=EstimatedParameters(
     beta_leaf=0.15, 
     beta_wood=0.6, 
     r_C_leaf2abvstrlit= 0.3,
     r_C_abvmetlit2surface_microbe=0.2,
     r_C_abvstrlit2slowsom=0.4,
     r_C_abvstrlit2surface_microbe=0.5,
     r_C_belowmetlit2soil_microbe=0.4,
     r_C_belowstrlit2slowsom=0.3,
     r_C_belowstrlit2soil_microbe=0.4,
     r_C_leached=0.1,
     r_C_leaf2abvmetlit=0.6,
     r_C_passsom2soil_microbe=0.2,
     r_C_root2belowmetlit=0.3,
     r_C_root2belowstrlit=0.4,
     r_C_slowsom2passsom=0.3,
     r_C_slowsom2soil_microbe=0.1,
     r_C_soil_microbe2passsom=0.05,
     r_C_soil_microbe2slowsom=0.1,
     r_C_surface_microbe2slowsom=0.02,
     r_C_wood2abvmetlit=0.4,
     r_C_wood2abvstrlit=0.6,
     C_leaf_0=svs_0.cVeg/3,
     #C_root_0=svs_0.cVeg/3,
     C_abvstrlit_0=svs_0.cLitter/4,
     C_abvmetlit_0=svs_0.cLitter/4,
     C_blwstrlit_0=svs_0.cLitter/4,
     C_surfacemic_0=svs_0.cSoil/4,
     C_soilmic_0=svs_0.cSoil/4,
     C_slow_0=svs_0.cSoil/4
)

# check initial values for the pool siz
V_init

# +
from typing import Callable
from general_helpers import month_2_day_index
from functools import reduce

def make_param2res_sym(
        cpa: UnEstimatedParameters
    ) -> Callable[[np.ndarray], np.ndarray]: 
    # To compute numeric solutions we will need to build and iterator 
    # as we did before. As it will need numeric values for all the parameters 
    # we will have to create a complete dictionary for all the symbols
    # exept those for the statevariables and time.
    # This set of symbols does not change during the mcmc procedure, since it only
    # depends on the symbolic model.
    # Therefore we create it outside the mcmc loop and bake the result into the 
    # param2res function.
    # The iterator does not care if they are estimated or not it only 
    # wants a dictionary containing key: value pairs for all
    # parameters that are not state variables or the symbol for time
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test
    def npp_func(day):
        month=day_2_month_index(day)
        return dvs.npp[month]   # kg/m2/s kg/m2/day;
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
        V_init = StartVector(
            C_leaf = epa.C_leaf_0,
            C_root= cpa.cRoot_0,
            C_wood = cpa.cVeg_0 - epa.C_leaf_0 - cpa.cRoot_0,
            C_abvstrlit = epa.C_abvstrlit_0,
            C_abvmetlit = epa.C_abvmetlit_0,
            C_belowstrlit = epa.C_blwstrlit_0,
            C_belowmetlit = cpa.cLitter_0 - epa.C_abvstrlit_0 - epa.C_abvmetlit_0 - epa.C_blwstrlit_0,
            C_surface_microbe = epa.C_surfacemic_0,
            C_soil_microbe = epa.C_soilmic_0,
            C_slowsom = epa.C_slow_0,
            C_passsom= cpa.cSoil_0 - epa.C_surfacemic_0 - epa.C_soilmic_0 - epa.C_slow_0,
            rh=cpa.rh_0
        )
        # next we create the parameter dict for the iterator
        # The iterator does not care if they are estimated or not so we look for them
        # in the combination
        apa = {**cpa._asdict(),**epa._asdict()}
        model_par_dict = {
            k:v for k,v in apa.items()
            if k in model_par_dict_keys
        }
        # Beside the par_dict the iterator also needs the python functions to replace the symbolic ones with
        # our fake xi_func could of course be defined outside of param2res but in general it
        # could be build from estimated parameters and would have to live here...
        def xi_func(day):
            return 1.0 # preliminary fake for lack of better data... 
    
        func_dict={
            'NPP':npp_func,
             'xi':xi_func
        }
    
        it_sym = make_daily_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=par_dict,
            func_dict=func_dict
        )
        
        # Now that we have the iterator we can start to compute.
        # the first thing to notice is that we don't need to store all daily values,
        # since the observations are recorded monthly while our iterator has a
        # daily timestep.
        # - for the pool contents we only need to record the last day of the month
        # - for the respiration(s) ra and rh we have to sum up the daily values 
        #   over a month
        # 
        # Note: check if TRENDY months are like this...
        days_per_month = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
        sols=[]
        for m in range(cpa.number_of_months):
            dpm = days_per_month[ m % 12]  
            #mra=0
            mrh=0
            for d in range(dpm):
                v = it_sym.__next__()
                #mra +=v[9,0]
                mrh +=v[11,0]
            V=StartVector(*v)
            o=Observables(
                cVeg=float(V.C_leaf+V.C_wood+V.C_root),
                cRoot = float(V.C_root),
                #cLitter=float(V.C_leaf_litter+V.C_wood_litter+V.C_root_litter),
                cLitter=float(V.C_abvstrlit + V.C_abvmetlit + V.C_belowstrlit + V.C_belowmetlit),
                #cSoil=float(V.C_soil_fast+V.C_soil_slow+V.C_soil_passive),
                cSoil=float(V.C_surface_microbe + V.C_soil_microbe + V.C_slowsom + V.C_passsom),
                #ra=mra,
                rh=mrh,
            )
            # equivalent
            #o=np.array([
            #    np.sum(v[0:3]),
            #    np.sum(v[3:6]),
            #    np.sum(v[6:9]),
            #    mra,
            #    mrh,
            #])
            sols.append(o)
            
        sol=np.stack(sols)       
        return sol
        
    return param2res
# -



# +
# now test it 
import matplotlib.pyplot as plt
from general_helpers import plot_solutions
const_params = cpa

param2res_sym = make_param2res_sym(const_params)
xs= param2res_sym(epa_0)

    
day_indices=month_2_day_index(range(cpa.number_of_months)),

fig = plt.figure()
plot_solutions(
        fig,
        #times=day_indices,
        times=range(cpa.number_of_months),
        var_names=Observables._fields,
        tup=(xs,)
)
fig.savefig('solutions_SDGVM.pdf')

# +
# now test it 
import matplotlib.pyplot as plt
from general_helpers import plot_solutions
const_params = cpa

param2res_sym = make_param2res_sym(const_params)
xs= param2res_sym(epa_0)

# convert model output to yearly
mod=np.zeros(obs.shape[0]*obs.shape[1]).reshape([obs.shape[0],obs.shape[1]])  
for i in range(obs.shape[1]):
    mod[:,i]=monthly_to_yearly(xs[:,i])

day_indices=month_2_day_index(range(cpa.number_of_months)),

fig = plt.figure()
plot_solutions(
        fig,
        #times=day_indices,
        times=range(int(cpa.number_of_months/12)),
        var_names=Observables._fields,
        tup=(mod,obs)
)
fig.savefig('solutions.pdf')

# -

ns=1400
days=list(range(ns))
npp=[ npp_func(d) for d in days]

import matplotlib.pyplot as plt
n=len(mvs.get_StateVariableTuple())
#fig=plt.figure(figsize=(10,(n+1)*10))
fig=plt.figure()
axs=fig.subplots(2,1)
axs[0].plot(days,npp)
days2=list(range(70))
axs[1].plot(days2,[day_2_month_index(d) for d in days2])



