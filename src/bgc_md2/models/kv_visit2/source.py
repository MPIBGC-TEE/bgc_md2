from pathlib import Path
import numpy as np
from sympy import Symbol, Function, exp, Piecewise
import json
from CompartmentalSystems import start_distributions as sd
from ComputabilityGraphs.CMTVS import CMTVS
from ..BibInfo import BibInfo
from ... import helper as h
from ...resolve.mvars import (
    NumericParameterization,
    NumericSimulationTimes,
    NumericStartValueArray,
    NumericStartMeanAgeTuple,
    NumericParameterizedSmoothReservoirModel
)
#from bgc_md2 import general_helpers as gh
from . import source_1 as s1
from .CachedParameterization import CachedParameterization
from .subs_1 import subs_dict

def subs_xi(var):
    return var.subs(subs_dict)

mvs=CMTVS(
    {
        s1.mvs.get_TimeSymbol(),
        s1.mvs.get_Temperature(),
        s1.mvs.get_StateVariableTuple(),
        s1.mvs.get_CarbonStateVariableTuple(),
        s1.mvs.get_VegetationCarbonStateVariableTuple(),
        s1.mvs.get_SoilCarbonStateVariableTuple(),
        subs_xi(s1.mvs.get_InFluxesBySymbol()),
        subs_xi(s1.mvs.get_OutFluxesBySymbol()),
        subs_xi(s1.mvs.get_InternalFluxesBySymbol()),
        #NumericParameterization(
        #    par_dict=)
        s1.mvs.get_BibInfo(),
    },
    computers=s1.mvs.computers
)
# we provide some example parameterization
#
dirPath = Path(__file__).parent
#
pp = dirPath.joinpath("parameterization_from_test_args")
cp=CachedParameterization.from_path(pp)
par_dict = cp.parameter_dict
func_dict = cp.func_dict

days_per_month=h.date.days_per_month
start_shift_in_months = 4+27*12
#t0 = start_shift_in_months * days_per_month # lets start in spring
t0 = 120
#number_of_months = cp.drivers.npp.shape[0]-start_shift_in_months  #2881
number_of_months = 60 # this is only a test example so we save resources
t_max= t0 + number_of_months * days_per_month 
times = np.linspace(t0, t_max, number_of_months)

# For this example we assume that the system was in steady state 
# at t_0 with X_fix given by X_fix = -M(t_0)^{-1} I(t_0) 
# since we know that the system is linear we can easily compute the steady state
# (in general (nonlinear fluxes) it is not clear that a steady state even exists

# The following function computes consistent steady state start conditions for the whole system and all subsystems
def scm(
        npsrm #: NumericParameterizedSmoothReservoirModel or a subclass 
    ):
    srm=npsrm.srm # the Smoth reservoir model (or a subclass)
    nupa=npsrm.parameterization
    par_dict=nupa.par_dict
    func_dict=nupa.func_dict
    a_dens_function, X_fix = sd.start_age_distributions_from_steady_state(
        srm, 
        t0=t0, 
        parameter_dict= par_dict, 
        func_set=func_dict, 
        #x0=smr.start_values#x0=X_0
    )
    start_mean_age_vec = sd.start_age_moments_from_steady_state(
        srm,
        t0=t0,
        parameter_dict=par_dict,
        func_set=func_dict,
        max_order=1,
        x0=X_fix
    ).reshape(-1)
    return X_fix,start_mean_age_vec, a_dens_function

from bgc_md2.resolve.mvars import StartConditionMaker
mvs = mvs.update({
    NumericParameterization(
        par_dict=par_dict,
        func_dict=func_dict
    ),
    NumericSimulationTimes(times),
    StartConditionMaker(scm)
})
