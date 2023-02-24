from pathlib import Path
from sympy import Symbol, Function, exp, Piecewise
import numpy as np
import json
from CompartmentalSystems import start_distributions as sd
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.models.BibInfo import BibInfo
from bgc_md2.resolve.mvars import (
    NumericParameterization,
    NumericSimulationTimes,
    NumericStartValueArray,
    NumericStartMeanAgeVector,
    NumericParameterizedSmoothReservoirModel
)
from . import source_1 as s1
from . import model_specific_helpers_2 as msh
from .subs_1 import subs_dict
from bgc_md2 import general_helpers as gh

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
# we provide some example data

t0 = 0 #3/2*np.pi
n_steps = 12  # 2881
t_max = 144 
times = np.linspace(t0, t_max, n_steps)
delta_t_val = (t_max - t0)/n_steps
dirPath = Path(__file__).parent
dvs = msh.get_cached_global_mean_drivers(dirPath)
cpa = gh.load_named_tuple_from_json_path(
    msh.Constants,
    dirPath.joinpath("cpa.json")
)        
epa_opt = gh.load_named_tuple_from_json_path(
    msh.EstimatedParameters,
    dirPath.joinpath("epa_opt.json")
)        
func_dict = msh.make_func_dict(dvs,cpa,epa_opt)
par_dict = gh.make_param_dict(mvs, cpa, epa_opt)
srm = mvs.get_SmoothReservoirModel()
# For this example we assume that the system was in steady state 
# at t_0 with X_fix given by X_fix = M^{-1} 
# since we know that the system is linear we can easily compute the steady state
# (in general (nonlinear fluxes) it is not clear that a steady state even exists
# let alone how to compute it
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
