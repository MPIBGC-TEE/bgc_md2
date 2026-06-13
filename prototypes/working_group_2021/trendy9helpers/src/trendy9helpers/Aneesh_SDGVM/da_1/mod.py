from importlib import import_module
from collections import namedtuple
from ComputabilityGraphs import CMTVS
from typing import Callable, Dict, Tuple
from pathlib import Path
import numpy as np

from CompartmentalSystems.ArrayDictResult import ArrayDictResult
import bgc_md2.helper as h

from  ... import general_helpers as gh
from .. import model_specific_helpers_2 as msh  

model_mod = 'bgc_md2.models.Aneesh_SDGVM'
cp_mod=import_module(f"{model_mod}.CachedParameterization")
mvs=import_module(f"{model_mod}.source").mvs
make_func_dict = cp_mod.make_func_dict
Drivers=cp_mod.Drivers

# this is what is read from cpa.json 
FreeConstants = namedtuple(
    "FreeConstants",
    []
    # funny you actually don't have any free constants
    # so your cpa.json is empty and everything is guessed, but you 
    # could transfer some of the FreeEstimatedParameters here if 
    # you have a good guess for their value
)    
DerivedConstants = namedtuple(
    "DerivedConstants",
    [
        "cVeg_0",
        "cRoot_0",
        "cLitter_0",
        "cSoil_0",
    ]
)
Constants = namedtuple(
    "Constants",
    [*FreeConstants._fields,*DerivedConstants._fields] 
)
def cpa(
        fcpa: FreeConstants,
        dvs: msh.Drivers,
        svs: msh.Observables
    ) -> Constants:
    dcpa = DerivedConstants(
        cVeg_0=svs.cVeg[0],
        cRoot_0=svs.cRoot[0],
        cLitter_0=svs.cLitter[0],
        cSoil_0=svs.cSoil[0],
    )
    print(f"############################## dcpa={dcpa}")
    return Constants(*fcpa, *dcpa)


FreeEstimatedParameters = namedtuple(
    "FreeEstimatedParameters",
    [
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
         #'r_C_leached',
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
         'r_C_abvstrlit_rh',
         'r_C_abvmetlit_rh',
         'r_C_belowstrlit_rh',
         'r_C_belowmetlit_rh',
         'r_C_surface_microbe_rh',
         'r_C_slowsom_rh',
         'r_C_passsom_rh',
         'r_C_soil_microbe_rh',
    ]
)
DerivedEstimatedParameters = namedtuple(
    "DerivedEstimatedParameters",
    [
        'C_leaf_0', 
        'C_abvstrlit_0',
        'C_abvmetlit_0',
        'C_blwstrlit_0',
        'C_surfacemic_0',
        'C_soilmic_0',
        'C_slow_0'
    ]
)
EstimatedParameters = namedtuple(
    "EstimatedParameters",
    [*FreeEstimatedParameters._fields, *DerivedEstimatedParameters._fields]
)
def epa_min(
    fepa: FreeEstimatedParameters,
    dvs: msh.Drivers,
    svs: msh.Observables
) -> EstimatedParameters:
    # this function can obviously be changed
    # it avoids contradictions between the startvalues and 
    # data 
    depa = DerivedEstimatedParameters(
         C_leaf_0=0,
         C_abvstrlit_0=0,
         C_abvmetlit_0=0,
         C_blwstrlit_0=0,
         C_surfacemic_0=0,
         C_soilmic_0=0,
         C_slow_0=0
    )
    return EstimatedParameters(*fepa, *depa)

def epa_max(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # this function can obviously be changed
    # it avoids contradictions between the startvalues and 
    # data 
    cv = svs.cVeg[0]
    cl = svs.cLitter[0]
    cs = svs.cSoil[0]
    depa = DerivedEstimatedParameters(
        C_leaf_0=cv, 
        C_abvstrlit_0=cl,
        C_abvmetlit_0=cl,
        C_blwstrlit_0=cl,
        C_surfacemic_0=cs,
        C_soilmic_0=cs,
        C_slow_0=cs
    )
    return EstimatedParameters(*fepa,*depa)

def epa_0(
    fepa: FreeEstimatedParameters,
    dvs: msh.Drivers,
    svs: msh.Observables
) -> EstimatedParameters:
    # this function can obviously be changed
    # it avoids contradictions between the startvalues and 
    # data 
    cv = svs.cVeg[0]/2
    cl = svs.cLitter[0]/4
    cs = svs.cSoil[0]/4
    depa = DerivedEstimatedParameters(
        C_leaf_0=cv, 
        C_abvstrlit_0=cl,
        C_abvmetlit_0=cl,
        C_blwstrlit_0=cl,
        C_surfacemic_0=cs,
        C_soilmic_0=cs,
        C_slow_0=cs
    )
    return EstimatedParameters(*fepa, *depa)


def make_proposer(
    c_max: EstimatedParameters,
    c_min: EstimatedParameters, 
    fcpa: FreeConstants,
    dvs: msh.Drivers,
    svs: msh.Observables
) -> Callable[[np.ndarray, float, bool], np.ndarray]:
    cpa_v = cpa(fcpa, dvs, svs)
    """Returns a function that will be used by the mcmc algorithm to propose
    a new parameter value tuple based on a given one.
    The two arrays c_max and c_min define the boundaries
    of the n-dimensional rectangular domain for the parameters and must be of
    the same shape.  After a possible parameter value has been sampled the
    filter_func will be applied to it to either accept or discard it.  So
    filter func must accept parameter array and return either True or False
    :param c_max: array of maximum parameter values
    :param c_min: array of minimum parameter values
    :param D: a damping parameter to regulate the proposer step 
     larger D means smaller step size. If the maximum distance for a new value is
     max_dist then the proposer will make a max_dist/ 
    :param filter_func: model-specific function to filter out impossible 
    parameter combinations
    """


    g = np.random.default_rng()
    # filter out the startvalues and use a dirichlet distribution 
    # for them  
    dirichlet_tups= [
        (["beta_leaf", "beta_wood"], 1),
        (["C_leaf_0" ], cpa_v.cVeg_0),
        (["C_abvstrlit_0", "C_abvmetlit_0", "C_blwstrlit_0"], cpa_v.cLitter_0),
        (["C_surfacemic_0", "C_soilmic_0", "C_slow_0"], cpa_v.cSoil_0),
    ]
    return gh.make_dirichlet_uniform_proposer(
        dirichlet_tups=dirichlet_tups,
        EstimatedParameters=EstimatedParameters,
        c_min=c_min,
        c_max=c_max
    )

def make_param_filter_func(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters, 
        fcpa: Constants,
        dvs: Drivers,
        svs: msh.Observables,
        ) -> Callable[[np.ndarray], bool]:

    # find position of beta_leaf and beta_wood
        
    
    def isQualified(
            c,
            print_conds=False
        )-> bool:
        def value(field_name):
            try:
                return c[EstimatedParameters._fields.index(field_name)]
            except Exception as e:
                print("###########################")
                print(e)
                print(field_name)
                raise e
        conds=[
            # not necessary any more since the
            # proposer does not choose values
            # outside the range except for a human error in epa_0...
            (c >= c_min).all(), 
            (c <= c_max).all(),

            sum(map(value, ["beta_leaf", "beta_wood"])) <= 0.99,
            sum(map(value, ["C_leaf_0"])) <= svs.cVeg[0],
            sum(
                map(
                    value,
                    ["C_abvstrlit_0", "C_abvmetlit_0", "C_blwstrlit_0"]
                )
            ) <= svs.cLitter[0],
            sum(
                map(
                    value,
                    ["C_surfacemic_0", "C_soilmic_0", "C_slow_0"]
                )
            ) <= svs.cSoil[0]  
        ]    
        res=all(conds)
        if print_conds:
            if not res:
                print(conds)
                #from IPython import embed; embed()
        return res
        
    return isQualified


def make_param2res_sym(
        mvs,
        fcpa: Constants,
        dvs: Drivers,
        svs: msh.Observables
) -> Callable[[np.ndarray], np.ndarray]: 
    
    cpa_v = cpa(fcpa, dvs, svs)
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        apa = {**cpa_v._asdict(), **epa._asdict()}
        X_0 = numeric_X_0(mvs, dvs, apa)
        par_dict = gh.make_param_dict2(mvs, apa)
        func_dict = make_func_dict(dvs , cpa=cpa_v, epa=epa)
        return msh.synthetic_observables(
            mvs,
            X_0,
            par_dict=par_dict,
            func_dict=func_dict,
            dvs=dvs
        )
    return param2res


def numeric_X_0(mvs,dvs,apa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    par_dict = gh.make_param_dict2(mvs, apa)
    X_0_dict = {
        "C_leaf": apa['C_leaf_0'],
        "C_root": apa['cRoot_0'],
        "C_wood": apa['cVeg_0'] - (apa['C_leaf_0'] ),
        "C_abvstrlit": apa['C_abvstrlit_0'],
        "C_abvmetlit": apa['C_abvmetlit_0'],
        "C_belowstrlit": apa["C_blwstrlit_0"],
        "C_belowmetlit": (
            apa["cLitter_0"]
            - apa["C_abvstrlit_0"]
            - apa["C_abvmetlit_0"]
            - apa["C_blwstrlit_0"]
        ),
        "C_surface_microbe": apa["C_surfacemic_0"],
        "C_soil_microbe": apa["C_soilmic_0"],
        "C_slowsom": apa["C_slow_0"],
        "C_passsom": (
            apa["cSoil_0"]
            - apa["C_surfacemic_0"]
            - apa["C_soilmic_0"]
            - apa["C_slow_0"]
        )
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict), 1)
    return X_0
