import numpy as np
from pathlib import Path
from sympy import Symbol, Function
from ComputabilityGraphs.CMTVS import CMTVS
from CompartmentalSystems import start_distributions as sd
from ...helper import module_computers
from ..BibInfo import BibInfo
from ...resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonStateVariableTuple,
    SoilCarbonStateVariableTuple,
    CarbonStateVariableTuple,
    NumericParameterization,
    NumericSimulationTimes,
    NumericStartValueArray,
    NumericStartMeanAgeTuple,
    NumericParameterizedSmoothReservoirModel,
    StartConditionMaker
)
from ...resolve import computers as bgc_c
from .CachedParameterization import CachedParameterization

# Make a dictionary for the variables we will use
# The variables should be described
sym_dict={
    'c_leaf_0': '',               #Names: c_poolname_0
    'c_root_0': '',               #Only initial pools that are estimated
    'c_lit_cwd_0': '',
    'c_lit_met_0': '',
    'c_lit_str_0': '',
    'c_lit_mic_0': '',
    'c_soil_met_0': '',
    'c_soil_str_0': '',
    'c_soil_mic_0': '',
    'c_soil_slow_0': '',
    'c_leaf': 'carbon content of leaf pool',
    'c_root': '',
    'c_wood': '',
    'c_lit_cwd': '',
    'c_lit_met': '',
    'c_lit_str': '',
    'c_lit_mic': '',
    'c_soil_met': '',
    'c_soil_str': '',
    'c_soil_mic': '',
    'c_soil_slow': '',
    'c_soil_passive': '',
    #'r_c_leaf_rh': '',
    #'r_c_root_rh': '',
    #'r_c_wood_rh': '',
    'r_c_lit_cwd_rh': '',
    'r_c_lit_met_rh': '',
    'r_c_lit_str_rh': '',
    'r_c_lit_mic_rh': '',
    'r_c_soil_met_rh': '',
    'r_c_soil_str_rh': '',
    'r_c_soil_mic_rh': '',
    'r_c_soil_slow_rh': '',
    'r_c_soil_passive_rh': '',
    'r_c_leaf_2_c_lit_met': '',
    'r_c_leaf_2_c_lit_str': '',
    'r_c_root_2_c_soil_met': '',    
    'r_c_root_2_c_soil_str': '',
    'r_c_wood_2_c_lit_cwd': '',
    'r_c_lit_cwd_2_c_lit_mic': '',
    'r_c_lit_cwd_2_c_soil_slow': '',
    'r_c_lit_met_2_c_lit_mic': '',
    'r_c_lit_str_2_c_lit_mic': '',
    'r_c_lit_str_2_c_soil_slow': '',
    'r_c_lit_mic_2_c_soil_slow': '',
    'r_c_soil_met_2_c_soil_mic': '',
    'r_c_soil_str_2_c_soil_mic': '',
    'r_c_soil_str_2_c_soil_slow': '',
    'r_c_soil_mic_2_c_soil_slow': '',
    'r_c_soil_mic_2_c_soil_passive': '',
    'r_c_soil_slow_2_c_soil_mic': '',
    'r_c_soil_slow_2_c_soil_passive': '',
    'r_c_soil_passive_2_c_soil_mic': '',
    'silt': '',
    'clay': '',
    'beta_leaf': '',
    'beta_root': '',
}
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

#create symbols for scaler and input functions
func_dict={
    'xi_leaf': 'Environmental scaler for leaf pool',
    'xi_soil': 'Environmental scaler for soil pool',
    'GPP': 'Aboveground gross input',
    'NPP': 'Aboveground net input'
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")

mvs = CMTVS(
    {
        t,
        StateVariableTuple(( 
            c_leaf,
            c_root,
            c_wood,
            c_lit_cwd,
            c_lit_met,
            c_lit_str,
            c_lit_mic,
            c_soil_met,
            c_soil_str,
            c_soil_mic,
            c_soil_slow,
            c_soil_passive,
        )),
        VegetationCarbonStateVariableTuple((
            c_leaf,
            c_root,
            c_wood,
        )),
        SoilCarbonStateVariableTuple((
            c_lit_cwd,
            c_lit_met,
            c_lit_str,
            c_lit_mic,
            c_soil_met,
            c_soil_str,
            c_soil_mic,
            c_soil_slow,
            c_soil_passive,
        )),
        InFluxesBySymbol(
            {   #define input/allocation
                #RecievingPool: Input * Allocation
                c_leaf: NPP(t) * beta_leaf,
                c_root: NPP(t) * beta_root,
                c_wood: NPP(t) * (1-beta_leaf-beta_root)
            }
        ),
        OutFluxesBySymbol(
            {   #define fluxes leaving the system
                #Fluxes leaving the system: FluRate * DonorPool * EnvironmentalScaler
                #c_leaf: r_c_leaf_rh * c_leaf * xi_leaf(t),
                #c_root: r_c_root_rh * c_root * xi_leaf(t),
                #c_wood: r_c_wood_rh * c_wood * xi_leaf(t),
                c_lit_cwd: r_c_lit_cwd_rh * c_lit_cwd * xi_soil(t),
                c_lit_met: r_c_lit_met_rh * c_lit_met * xi_soil(t),
                c_lit_str: r_c_lit_str_rh * c_lit_str * xi_soil(t),
                c_lit_mic: r_c_lit_mic_rh * c_lit_mic * xi_soil(t),
                c_soil_met: r_c_soil_met_rh * c_soil_met * xi_soil(t),
                c_soil_str: r_c_soil_str_rh * c_soil_str * xi_soil(t),
                c_soil_mic: r_c_soil_mic_rh * c_soil_mic * xi_soil(t),
                c_soil_slow: r_c_soil_slow_rh * c_soil_slow * xi_soil(t),
                c_soil_passive: r_c_soil_passive_rh * c_soil_passive * xi_soil(t),
            }
        ),
        InternalFluxesBySymbol(
            {   #define fluxes between pools
                #(Donor pool, recieving pool): FluxRate * DonorPool
                (c_leaf, c_lit_met): r_c_leaf_2_c_lit_met * c_leaf * xi_leaf(t),
                (c_leaf, c_lit_str): r_c_leaf_2_c_lit_str * c_leaf * xi_leaf(t),
                (c_root, c_soil_met): r_c_root_2_c_soil_met * c_root * xi_leaf(t),
                (c_root, c_soil_str): r_c_root_2_c_soil_str * c_root * xi_leaf(t),
                (c_wood, c_lit_cwd): r_c_wood_2_c_lit_cwd * c_wood * xi_leaf(t),
                (c_lit_cwd, c_lit_mic): r_c_lit_cwd_2_c_lit_mic * c_lit_cwd * xi_soil(t),
                (c_lit_cwd, c_soil_slow): r_c_lit_cwd_2_c_soil_slow * c_lit_cwd * xi_soil(t),
                (c_lit_met, c_lit_mic): r_c_lit_met_2_c_lit_mic * c_lit_met * xi_soil(t),
                (c_lit_str, c_lit_mic): r_c_lit_str_2_c_lit_mic * c_lit_str * xi_soil(t),
                (c_lit_str, c_soil_slow): r_c_lit_str_2_c_soil_slow * c_lit_str * xi_soil(t),
                (c_lit_mic, c_soil_slow): r_c_lit_mic_2_c_soil_slow * c_lit_mic * xi_soil(t),
                (c_soil_met, c_soil_mic): r_c_soil_met_2_c_soil_mic * c_soil_met * xi_soil(t),
                (c_soil_str, c_soil_mic): r_c_soil_str_2_c_soil_mic * c_soil_str * xi_soil(t),
                (c_soil_str, c_soil_slow): r_c_soil_str_2_c_soil_slow * c_soil_str * xi_soil(t),
                (c_soil_mic, c_soil_slow): r_c_soil_mic_2_c_soil_slow * c_soil_mic * xi_soil(t),
                (c_soil_mic, c_soil_passive): r_c_soil_mic_2_c_soil_passive * c_soil_mic * xi_soil(t),
                (c_soil_slow, c_soil_mic): r_c_soil_slow_2_c_soil_mic * c_soil_slow * xi_soil(t),
                (c_soil_slow, c_soil_passive): r_c_soil_slow_2_c_soil_passive * c_soil_slow * xi_soil(t),
                (c_soil_passive, c_soil_mic): r_c_soil_passive_2_c_soil_mic * c_soil_passive * xi_soil(t)
            }
        ),
        BibInfo(# Bibliographical Information
            name="YIBS",
            longName="",
            version="1",
            entryAuthor="Jon Wells",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        )
    }
    ,
    computers=module_computers(bgc_c)
)
# we provide some example parameterization

dirPath = Path(__file__).parent

pp = Path(__file__).parent.joinpath("parameterization_from_test_args")
cp=CachedParameterization.from_path(pp)
par_dict = cp.parameter_dict
func_dict = cp.func_dict

t0 = 0  #3/2*np.pi
number_of_months = cp.drivers.npp.shape[0] 
n_steps = number_of_months
t_max = number_of_months*30 # time measured in days 
times = np.linspace(t0, t_max, n_steps)

# For this example we assume that the system was in steady state 
# at t_0 with X_fix given by X_fix = M^{-1} 
# since we know that the system is linear we can easily compute the steady state
# (in general (nonlinear fluxes) it is not clear that a steady state even exists
# let alone how to compute it

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

mvs = mvs.update({
    NumericParameterization(
        par_dict=par_dict,
        func_dict=func_dict
    ),
    NumericSimulationTimes(times),
    StartConditionMaker(scm)
})
