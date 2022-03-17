import sys
import json 
from pathlib import Path
from collections import namedtuple 
import netCDF4 as nc
import numpy as np
from sympy import Symbol
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
        TimeStepIterator2,
)
from copy import copy
from typing import Callable
from general_helpers import month_2_day_index, monthly_to_yearly
from functools import reduce

sys.path.insert(0,'..') # necessary to import general_helpers
import general_helpers as gh


Constants = namedtuple(
    "Constants",
    [
        'C_leaf_0',
        'C_root_0',
        'C_wood_0',
        'clitter_0',
        'csoil_0',
        'rh_0',
        'clay',
        'silt',
        'lig_wood',
        'f_wood2CWD',
        'f_metlit2mic',
        'npp',
        'number_of_months'
    ]
)
#
## This set is used (as argument) by the functions that are called
## inside the mcmc
#EstimatedParameters = namedtuple(
#    "EstiamatedParameters",
#    [
#        "beta_leaf",    #  0 (indices uses in original code) 
#        "beta_root",    #  1
#        "lig_leaf",     #  2
#        "f_leaf2metlit",#  3
#        "f_root2metlit",#  4
#        "k_leaf",       #  5
#        "k_root",       #  6
#        "k_wood",       #  7
#        "k_metlit",	#  8
#        "k_mic",	#  9
#        "k_slowsom",	# 10
#        "k_passsom",	# 11
#        "C_metlit_0",	# 12
#        "C_CWD_0",	# 13
#        "C_mic_0",	# 14
#        "C_passom_0"    # 15
#    ]
#)
#
Observables = namedtuple(
    'Observables',
        ["cCwd","cLeaf", "cLitter", "cRoot", "cSoil", "cVeg", "cWood", "npp", "rh"]
)

Drivers=namedtuple(
    "Drivers",
    ["npp"] 
)    

def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['CABLE-POP'],
        variables = Observables._fields + Drivers._fields
    )

    
def get_example_site_vars(dataPath):
    # pick up 1 site
    s = slice(None, None, None)  # this is the same as :
    t = s, 50, 33  # [t] = [:,49,325]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them
        ds = nc.Dataset(str(path))
        res = ds.variables[vn][t]    
        # kg/m2/s kg/m2/day;
        return res*86400 if vn in ['npp','rh'] else res

    o_names=[(f,"CABLE-POP_S2_{}.nc".format(f)) for f in Observables._fields]
    d_names=[(f,"CABLE-POP_S2_{}.nc".format(f)) for f in Drivers._fields]

    dvs=Drivers(*map(f,d_names))

    obss=Observables(*map(f, o_names))


    return (obss, dvs)

    
def make_npp_func(dvs):
    def func(day):
        month=gh.day_2_month_index(day)
        return (dvs.npp[month]) 

    return func


def make_xi_func(dvs):
    def func(day):
        return 1.0 # preliminary fake for lack of better data... 
    return func


def make_func_dict(mvs,dvs):
    return {
        "NPP": make_npp_func(dvs),
        "xi": make_xi_func(dvs)
    }


