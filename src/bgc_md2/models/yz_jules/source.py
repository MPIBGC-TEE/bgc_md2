#import sys
from pathlib import Path
import numpy as np
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
from CompartmentalSystems import start_distributions as sd

#import bgc_md2.display_helpers as dh
from ... import helper as h
from ..BibInfo import BibInfo
from ...helper import module_computers
from ...resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonStateVariableTuple,
    SoilCarbonStateVariableTuple,
    NumericParameterization,
    NumericSimulationTimes,
    NumericStartValueArray,
    NumericStartMeanAgeVector,
    NumericParameterizedSmoothReservoirModel
)
from ...resolve import computers as bgc_c
from .CachedParameterization import CachedParameterization

# Other packages


# Make a small dictionary for the variables we will use
sym_dict={
    'mrsos': 'Moisture in top soil (10cm) layer, in kg m-2',
    'tsl': 'Temperature of Soil - layer, top layer extracted from a four-layer data, in K',
    'Mw': 'soil moisture at wilting point as a fraction of saturation',
    'Ms': 'soil moisture content at saturation',
    'Topt': 'optimal soil temperature for heterotrophic respiration (deg C)',
    'Tcons': 'constant value in temperature scaler',
    'beta_leaf': 'NPP allocation fraction to leaf',
    'beta_wood': 'NPP allocation fraction to wood',
    #'beta_root': 'NPP allocation fraction to root',
    'c_leaf': 'leaf carbon pool',  # Names: c_poolname
    'c_wood': 'wood carbon pool',
    'c_root': 'root carbon pool',

    'c_DPM': 'decomposable plant material carbon pool',  
    'c_RPM': 'resistant plant material carbon pool', 
    'c_BIO': 'microbial biomass carbon pool',  
    'c_HUM': 'long-lived humified carbon',  
    'r_c_DPM_rh': 'heterotrophic respiration from DPM',  # Pools with loss from system will be listed here
    'r_c_RPM_rh': 'heterotrophic respiration from RPM',  # Names: r_c_poolname_rh
    'r_c_BIO_rh': 'heterotrophic respiration from BIO',
    'r_c_HUM_rh': 'heterotrophic respiration from HUM',
    'r_c_leaf_2_c_DPM': 'internal C flux/transfer from leaf to DPM',  # Pool transfer paths
    'r_c_leaf_2_c_RPM': 'internal C flux/transfer from leaf to RPM',  # Names: r_c_donorPool_2_recievingPool
    'r_c_wood_2_c_DPM': 'internal C flux/transfer from wood to DPM',
    'r_c_wood_2_c_RPM': 'internal C flux/transfer from wood to RPM',
    'r_c_root_2_c_DPM': 'internal C flux/transfer from root to DPM',
    'r_c_root_2_c_RPM': 'internal C flux/transfer from root to RPM',
    'r_c_DPM_2_c_BIO':  'internal C flux/transfer from DPM to BIO',
    'r_c_DPM_2_c_HUM':  'internal C flux/transfer from DPM to HUM',
    'r_c_RPM_2_c_BIO':  'internal C flux/transfer from RPM to BIO',
    'r_c_RPM_2_c_HUM':  'internal C flux/transfer from RPM to HUM',
    'r_c_BIO_2_c_HUM':  'internal C flux/transfer from BIO to HUM',
    'r_c_HUM_2_c_BIO':  'internal C flux/transfer from HUM to BIO'
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
        VegetationCarbonStateVariableTuple((
                c_leaf,  # Names: c_poolname
                c_wood,
                c_root,
        )),
        SoilCarbonStateVariableTuple((
                c_DPM,  # full name: decomposable plant material
                c_RPM,  # full name: resistant plant material
                c_BIO,  # full name: microbial biomass
                c_HUM   # full name: long-lived humified
        )),
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
            name="JULES",
            longName="Joint UK Land Environment Simulator",
            version="1",
            entryAuthor="Yu Zhou",
            entryAuthorOrcid="0000-0002-5544-8342",
            entryCreationDate="2022-03-21",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        ),
    },
    computers=module_computers(bgc_c)
)
# we provide some example parameterization

dirPath = Path(__file__).parent

pp = Path(__file__).parent.joinpath("parameterization_from_test_args")
cp=CachedParameterization.from_path(pp)
par_dict = cp.parameter_dict
func_dict = cp.func_dict

t0 = 0  #3/2*np.pi
n_steps = 12  # 2881
number_of_months = cp.drivers.npp.shape[0] 
t_max = number_of_months*30 # time measured in days 
times = np.linspace(t0, t_max, n_steps)

# For this example we assume that the system was in steady state 
# at t_0 with X_fix given by X_fix = M^{-1} 
# since we know that the system is linear we can easily compute the steady state
# (in general (nonlinear fluxes) it is not clear that a steady state even exists
# let alone how to compute it
srm = mvs.get_SmoothReservoirModel()

start_mean_age_vec, X_fix = sd.start_mean_age_vector_from_steady_state_linear(
    srm,
    t0=t0,
    parameter_dict=par_dict,
    func_set=func_dict
)
mvs = mvs.update({
    NumericParameterization(
        par_dict=par_dict,
        func_dict=func_dict
    ),
    NumericStartValueArray(X_fix.reshape(-1)),
    NumericSimulationTimes(times),
    NumericStartMeanAgeVector(start_mean_age_vec)
})
