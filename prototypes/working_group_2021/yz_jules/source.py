import sys
from collections import namedtuple
# Packages for symbolic code:
from sympy import Symbol, Function, diag, ImmutableMatrix 
from pathlib import Path
from copy import copy, deepcopy
from functools import reduce
from typing import Callable
from pprint import pprint
import json
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

from ComputabilityGraphs.CMTVS import CMTVS
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
import CompartmentalSystems.helpers_reservoir as hr

import bgc_md2.resolve.computers as bgc_c
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h
from bgc_md2.models.BibInfo import BibInfo
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

# Other packages

sys.path.insert(0, '..')  # necessary to import general_helpers
from general_helpers import (
    download_TRENDY_output,
    day_2_month_index,
    month_2_day_index,
    make_B_u_funcs_2,
    monthly_to_yearly,
    plot_solutions
)


# Make a small dictionary for the variables we will use
sym_dict={
    'mrso': 'Total Soil Moisture Content, in kg m-2',
    'tsl': 'Temperature of Soil - layer, four layers, in K',
    'Mw': 'soil moisture at wilting point as a fraction of saturation',
    'Ms': 'soil moisture content at saturation',
    'beta_leaf': 'NPP allocation fraction to leaf',
    'beta_wood': 'NPP allocation fraction to wood',
    #'beta_root': 'NPP allocation fraction to root',
    'c_leaf': '',  # Names: c_poolname
    'c_wood': '',
    'c_root': '',
    'c_DPM': '',  # full name: decomposable plant material
    'c_RPM': '',  # full name: resistant plant material
    'c_BIO': '',  # full name: microbial biomass
    'c_HUM': '',   # full name: long-lived humified
    'r_c_DPM_rh': '',  # Pools with loss from system will be listed here
    'r_c_RPM_rh': '',  # Names: r_c_poolname_rh
    'r_c_BIO_rh': '',
    'r_c_HUM_rh': '',
    'r_c_leaf_2_c_DPM': '',  # Pool transfer paths
    'r_c_leaf_2_c_RPM': '',  # Names: r_c_donorPool_2_recievingPool
    'r_c_wood_2_c_DPM': '',
    'r_c_wood_2_c_RPM': '',
    'r_c_root_2_c_DPM': '',
    'r_c_root_2_c_RPM': '',
    'r_c_DPM_2_c_BIO': '',
    'r_c_DPM_2_c_HUM': '',
    'r_c_RPM_2_c_BIO': '',
    'r_c_RPM_2_c_HUM': '',
    'r_c_BIO_2_c_HUM': '',
    'r_c_HUM_2_c_BIO': ''
}

for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# define beta wood from other allocation values
beta_root = 1.0 - (beta_leaf + beta_wood)

# create symbols for scaler and input functions
func_dict = {
    'xi': 'Environmental scaler as a function of time',
    'NPP': 'Inputs as a function of time',
}
for k in func_dict.keys():
    code = k + " = Function('{0}')".format(k)
    exec(code)

# define t as a symbol for time
t = TimeSymbol("t")
#
### Symbolic Model Description (Must Edit)
##Define your model using sympy:
#
## define model in sympy format
mvs = CMTVS(
    {
        t,
        StateVariableTuple(
            [
                c_leaf,  # Names: c_poolname
                c_wood,
                c_root,
                c_DPM,  # full name: decomposable plant material
                c_RPM,  # full name: resistant plant material
                c_BIO,  # full name: microbial biomass
                c_HUM   # full name: long-lived humified
            ]
        ),
        InFluxesBySymbol(
            {  # define input/allocation
                # RecievingPool: Input * Allocation
                c_leaf: NPP(t) * beta_leaf,
                c_root: NPP(t) * beta_root,
                c_wood: NPP(t) * beta_wood
            }
        ),
        OutFluxesBySymbol(
            {  # define fluxes leaving the system
                # Fluxes leaving the system: FluRate * DonorPool * EnvironmentalScaler
                c_DPM: r_c_DPM_rh * c_DPM * xi(t),
                c_RPM: r_c_RPM_rh * c_RPM * xi(t),
                c_BIO: r_c_BIO_rh * c_BIO * xi(t),
                c_HUM: r_c_HUM_rh * c_HUM * xi(t)
            }
        ),
        InternalFluxesBySymbol(
            {  # define fluxes between pools
                # (Donor pool, recieving pool): FluxRate * DonorPool
                (c_leaf, c_DPM): r_c_leaf_2_c_DPM * c_leaf,
                (c_leaf, c_RPM): r_c_leaf_2_c_RPM * c_leaf,
                (c_wood, c_DPM): r_c_wood_2_c_DPM * c_wood,
                (c_wood, c_RPM): r_c_wood_2_c_RPM * c_wood,
                (c_root, c_DPM): r_c_root_2_c_DPM * c_root,
                (c_root, c_RPM): r_c_root_2_c_RPM * c_root,
                (c_DPM,  c_BIO): r_c_DPM_2_c_BIO  * c_DPM  * xi(t),
                (c_DPM,  c_HUM): r_c_DPM_2_c_HUM  * c_DPM  * xi(t),
                (c_RPM,  c_BIO): r_c_RPM_2_c_BIO  * c_RPM  * xi(t),
                (c_RPM,  c_HUM): r_c_RPM_2_c_HUM  * c_RPM  * xi(t),
                (c_BIO,  c_HUM): r_c_BIO_2_c_HUM  * c_BIO  * xi(t),
                (c_HUM,  c_BIO): r_c_HUM_2_c_BIO  * c_HUM  * xi(t)
            }
        ),
        BibInfo(# Bibliographical Information
            name="Visit",
            longName="",
            version="1",
            entryAuthor="",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        ),
    },
    computers=module_computers(bgc_c)
)
