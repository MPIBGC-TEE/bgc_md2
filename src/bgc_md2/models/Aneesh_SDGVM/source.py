# import sys
# sys.path.insert(0,'..') # necessary to import general_helpers
from pathlib import Path
from sympy import Symbol, Function
import numpy as np
from CompartmentalSystems import start_distributions as sd
from ComputabilityGraphs.CMTVS import CMTVS
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
    NumericStartMeanAgeTuple,
    NumericParameterizedSmoothReservoirModel,
    StartConditionMaker
)
from ...resolve import computers as bgc_c
from ..BibInfo import BibInfo
from .CachedParameterization import CachedParameterization
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
         #'k_C_leaf': 'Turn over rate of leaf pool',
         #'k_C_wood':'Turn over rate of wood pool',
         #'k_C_root': 'Turn over rate of root pool',
         #'k_C_abvmetlit':'Turn over rate of aboveground metabolic litter pool',
         #'k_C_abvstrlit' :'Turn over rate of aboveground structural litter pool',
         #'k_C_belowstrlit' :'Turn over rate of belowground structural litter pool',
         #'k_C_belowmetlit':'Turn over rate of belowground metabolic litter pool',
         #'k_C_surface_microbe' :'Turn over rate of surface microbe pool',
         #'k_C_soil_microbe' :'Turn over rate of soil microbial pool',
         #'k_C_slowsom' :'Turn over rate of slow soil pool',
         #'k_C_passsom':'Turn over rate of passive soil pool',
         #'ls_aboveground': 'Fraction of structural C that is lignin in aboveground',
         #'ls_belowground': 'Fraction of structural C that is lignin in belowground',
         #'f_leaf2abvstrlit': 'Transfer coefficient of C from leaf to aboveground structural litter',
         #'f_leaf2abvmetlit':'Transfer coefficient of C from leaf to aboveground metabolic litter',
         #'f_wood2abvstrlit':'Transfer coefficient of C from wood to aboveground structural litter',
         #'f_wood2abvmetlit':'Transfer coefficient of C from wood to aboveground metabolic litter',
         #'f_root2belowstrlit':'Transfer coefficient of C from root to belowground structural litter',
         #'f_root2belowmetlit':'Transfer coefficient of C from root to belowground metabolic litter',
         #'f_abvstrlit2surface_microbe':'Transfer coefficient of C from aboveground structural litter to surface microbe',
         #'f_abvstrlit2slowsom':'Transfer coefficient of C from aboveground structural litter to slow soil',
         #'f_abvmetlit2surface_microbe':'Transfer coefficient of C from aboveground metabolic litter to surface microbe',
         #'f_belowstrlit2soil_microbe':'Transfer coefficient of C from belowground structural litter to soil microbe',
         #'f_belowstrlit2slowsom':'Transfer coefficient of C from belowground structural litter to slow soil',
         #'f_belowmetlit2soil_microbe':'Transfer coefficient of C from belowground metabolic litter to soil microbe',
         #'f_slowsom2soil_microbe':'Transfer coefficient of C from soil soil to soil microbe',
         #'f_passsom2soil_microbe':'Transfer coefficient of C from passive soil to soil microbe',
         #'f_surface_microbe2slowsom':'Transfer coefficient of C from surface microbe to slow soil',
         #'f_soil_microbe2passsom':'Transfer coefficient of C from soil microbe to passive soil',
         #'f_soil_microbe2slowsom':'Transfer coefficient of C from soil microbe to slow soil',
         #'f_slowsom2passsom':'Transfer coefficient of C from slow soil to passive soil',
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
        #'r_C_leached': '',
        'r_C_abvstrlit_rh':'',
        'r_C_abvmetlit_rh':'',
        'r_C_belowstrlit_rh':'',
        'r_C_belowmetlit_rh':'',
        'r_C_surface_microbe_rh':'',
        'r_C_slowsom_rh':'',
        'r_C_passsom_rh':'',
        'r_C_soil_microbe_rh' :''
}

for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument)
func_dict={
    'xi': 'a scalar function of temperature and moisture and thereby ultimately of time',
    'NPP': 'NPP',
}

for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")
beta_root=1-(beta_leaf+beta_wood)
mvs = CMTVS(
    {
        StateVariableTuple((
            C_leaf, 
            C_wood, 
            C_root, 
            C_abvstrlit,
            C_abvmetlit, 
            C_belowstrlit, 
            C_belowmetlit, 
            C_surface_microbe, 
            C_soil_microbe, 
            C_slowsom, 
            C_passsom
        )),
        VegetationCarbonStateVariableTuple((
            C_leaf, 
            C_wood, 
            C_root, 
        )),
        SoilCarbonStateVariableTuple((
            C_abvstrlit,
            C_abvmetlit, 
            C_belowstrlit, 
            C_belowmetlit, 
            C_surface_microbe, 
            C_soil_microbe, 
            C_slowsom, 
            C_passsom
        )),
        t,
        InFluxesBySymbol(
            {
                C_leaf: NPP(t)* beta_leaf, 
                C_wood: NPP(t)* beta_wood, 
                C_root:NPP(t)*  beta_root
            }),
        OutFluxesBySymbol(
            {
                C_abvstrlit: r_C_abvstrlit_rh * C_abvstrlit,
                C_abvmetlit: r_C_abvmetlit_rh * C_abvmetlit, 
                C_belowstrlit: r_C_belowstrlit_rh * C_belowstrlit, 
                C_belowmetlit: r_C_belowmetlit_rh * C_belowmetlit, 
                C_surface_microbe: r_C_surface_microbe_rh * C_surface_microbe,  #r_C_leached* C_surface_microbe, 
                C_soil_microbe: r_C_soil_microbe_rh * C_soil_microbe, 
                C_slowsom: r_C_slowsom_rh * C_slowsom, 
                C_passsom: r_C_passsom_rh * C_passsom
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
        ),
        BibInfo(# Bibliographical Information
            name="SDGVM",
            longName="",
            version="1",
            entryAuthor="Aneesh",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        ),
    },

    computers=module_computers(bgc_c)
)

# provide example parameterization
pp = Path(__file__).parent.joinpath("parameterization_from_test_args")
cp = CachedParameterization.from_path(pp)
par_dict = cp.parameter_dict
func_dict = cp.func_dict
t0 = 0  #3/2*np.pi
#sneak at the underlying data to avoid extrapolation (
#from IPython import embed; embed()
days_per_month=365.25/12
number_of_months = cp.drivers.npp.shape[0] 
n_steps = number_of_months
t_max = number_of_months*days_per_month # time measured in days 
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
