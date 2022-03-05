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
from general_helpers import download_TRENDY_output, day_2_month_index, make_B_u_funcs_2

observables_annual = namedtuple(
    'observables_annual',
    ["cSoil", "cVeg"]
)
observables_monthly = namedtuple(
    'observables_monthly',
    ["rh", "ra"]
)
observables = namedtuple(
    "observables",
    observables_annual._fields+observables_monthly._fields
)
Drivers=namedtuple(
    "Drivers",
    ["npp"]
)

#when downloading data make sure model names match TRENDY server names:
#"CABLE-POP","CLASSIC","CLM5","DLEM","IBIS","ISAM","ISBA_CTRIP",
#"JSBACH","JULES-ES-1.0","LPJ-GUESS","LPJwsl","LPX-Bern","OCN",
#"ORCHIDEE","ORCHIDEE-CNP","ORCHIDEEv3","ORCHIDEEv3_0.5deg"
#"SDGVM","VISIT","YIBs"

#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output(conf_dict):
    download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['YIBs'],
        variables =observables._fields + Drivers._fields
    )
#call it to test that the download works the data
#download_my_TRENDY_output()

def get_example_site_vars(dataPath):
    # pick up 1 site   
    s = slice(None, None, None)  # this is the same as :
    t = s, 74, 118  # [t] = [:,49,325]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them 
        ds = nc.Dataset(str(path))
        #check for npp/gpp/rh/ra to convert from kg/m2/s to kg/m2/day
        if vn in ["npp","gpp","rh","ra"]:
            return ds.variables[vn][t]*86400
        else:
            return ds.variables[vn][t]
    o_names=[(f,"YIBs_S2_Annual_{}.nc".format(f)) for f in observables_annual._fields]
    #YIBs has unique naming convention (Annual vs Monthly), create monthly names and extend 0_names vector
    monthly_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in observables_monthly._fields]
    o_names.extend(monthly_names)
    print(o_names)
    d_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Drivers._fields]
    return (observables(*map(f, o_names)),Drivers(*map(f,d_names)))


