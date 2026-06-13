
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

model_mod = 'bgc_md2.models.yz_jules'
cp_mod=import_module(f"{model_mod}.CachedParameterization")
mvs=import_module(f"{model_mod}.source").mvs
make_func_dict = cp_mod.make_func_dict
Drivers=cp_mod.Drivers
FreeConstants = namedtuple(
    "FreeConstants",
    []
)    
DerivedConstants = namedtuple(
    "Constants",
    [
        'c_veg_0',
        'c_soil_0',
        'fVegSoil_0',  # Total carbon mass from vegetation directly into the soil
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
        # Total carbon mass from vegetation directly into the soil
        c_veg_0=svs.cVeg[0],
        c_soil_0=svs.cSoil[0],
        fVegSoil_0=svs.fVegSoil[0],  # Total carbon mass from vegetation directly into the soil
    )
    return Constants(*fcpa, *dcpa)


FreeEstimatedParameters = namedtuple(
    'FreeEstimatedParameters',
    [
        'c_leaf_0',
        'c_wood_0',
        'c_DPM_0',
        'c_RPM_0',
        'c_BIO_0',
        'beta_leaf',
        'beta_wood',
        'Mw',
        'Ms',
        'Topt',
        'Tcons',
        'r_c_DPM_rh',
        'r_c_RPM_rh',
        'r_c_BIO_rh',
        'r_c_HUM_rh',
        'r_c_leaf_2_c_DPM',
        'r_c_leaf_2_c_RPM',
        'r_c_wood_2_c_DPM',
        'r_c_wood_2_c_RPM',
        'r_c_root_2_c_DPM',
        'r_c_root_2_c_RPM',
        'r_c_DPM_2_c_BIO',
        'r_c_DPM_2_c_HUM',
        'r_c_RPM_2_c_BIO',
        'r_c_RPM_2_c_HUM',
        'r_c_BIO_2_c_HUM',
        'r_c_HUM_2_c_BIO',
    ]
)
DerivedEstimatedParameters = namedtuple(
    'DerivedEstimatedParameters',
    []
)
EstimatedParameters = namedtuple(
    "EstimatedParameters",
    [*FreeEstimatedParameters._fields, *DerivedEstimatedParameters._fields] 
)
def epa_0(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # in this scheme it does not do anything ...
    # in others some parameters are to difficult to set
    # in the epa_0.json file since they depend on other parameters...
    depa = DerivedEstimatedParameters()
    return EstimatedParameters(*fepa, *depa)


def epa_min(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # in this scheme it does not do anything ...
    # in others some parameters are to difficult to set
    # in the epa_min.json file since they depend on other parameters...
    depa = DerivedEstimatedParameters()
    return EstimatedParameters(*fepa, *depa)


def epa_max(
        fepa: FreeEstimatedParameters,
        dvs: msh.Drivers,
        svs: msh.Observables
    )->EstimatedParameters:
    # in this scheme it does not do anything ...
    # in others some parameters are to difficult to set
    # in the epa_max.json file since they depend on other parameters...
    depa = DerivedEstimatedParameters()
    return EstimatedParameters(*fepa, *depa)

def make_param2res_sym(
        mvs,
        fcpa: FreeConstants,
        dvs: Drivers,
        svs: msh.Observables
) -> Callable[[np.ndarray], np.ndarray]:
    def param2res(pa):
        epa = EstimatedParameters(*pa)

        apa = {
            **fcpa._asdict(),
            **epa._asdict(),
            'c_root_0': svs.cVeg[0] - (epa.c_leaf_0 + epa.c_wood_0),
            'c_veg_0': svs.cVeg[0],
            'c_soil_0': svs.cSoil[0],
            # Total carbon mass from vegetation directly into the soil
            'fVegSoil_0': svs.fVegSoil[0],
        }
        X_0 = numeric_X_0(mvs, dvs, apa)

        #from IPython import embed; embed()
        dpm=h.date.days_per_month
        steps_per_month = 2
        delta_t_val = dpm/steps_per_month 

        par_dict = gh.make_param_dict2(mvs, apa)
        func_dict = make_func_dict(
            dvs,
            epa.Mw,
            epa.Ms,
            epa.Topt,
            epa.Tcons,
        )
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
    par_dict=gh.make_param_dict2(mvs,apa)
    X_0_dict={
        "c_leaf": apa['c_leaf_0'],     
        "c_wood": apa['c_wood_0'],     
        "c_root": apa['c_veg_0'] - (apa['c_leaf_0'] +  apa['c_wood_0']),  
        "c_DPM": apa['c_DPM_0'],
        "c_RPM": apa['c_RPM_0'],
        "c_BIO": apa['c_BIO_0'],
        "c_HUM": apa['c_soil_0'] - (apa['c_DPM_0'] + apa['c_RPM_0'] + apa['c_BIO_0'])
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    return X_0

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
        (["beta_leaf", "beta_wood"],1),
        (["c_leaf_0", "c_wood_0"], cpa_v.c_veg_0),
        (['c_DPM_0', 'c_RPM_0', 'c_BIO_0'],cpa_v.c_soil_0),
    ]
    # Note mm 7-25-2023
    # we could also try to make r_c_root_2_c_RPM fit svs.fVegSoil[0]
    # c_root_0 = svs.cVeg[0] - (epa.c_leaf_0 + epa.c_wood_0)
    # r_c_root_2_c_RPM = (
    #    svs.fVegSoil[0]
    #    -
    #    (
    #        (epa.r_c_leaf_2_c_DPM + epa.r_c_leaf_2_c_RPM) * epa.c_leaf_0
    #        + 
    #        (epa.r_c_wood_2_c_DPM + epa.r_c_wood_2_c_RPM) * epa.c_wood_0
    #        + 
    #        epa.r_c_root_2_c_DPM * c_root_0
    #    )
    #)/c_root_0
    #
    # but this would also require to choose 
    # epa.r_c_leaf_2_c_DPM,
    # epa.r_c_leaf_2_c_RPM,
    # epa.r_c_wood_2_c_DPM,
    # epa.r_c_wood_2_c_RPM,
    # epa.r_c_root_2_c_DPM,
    # automatically 
    # to guarantee that: r_c_root_2_c_RPM > 0
    # which could be done via recursive dirichlet distributions 
    # If implemented we would move
    # r_c_leaf_2_c_DPM,
    # r_c_leaf_2_c_RPM,
    # r_c_wood_2_c_DPM,
    # r_c_wood_2_c_RPM,
    # r_c_root_2_c_DPM,
    # from the FreeEstimatedParameters to the DerivedEstimatedParameters
    # and remove 
    # r_c_root_2_c_RPM from EstimatedParameters alltogether 
    # (and compute it via the above expression 
    return gh.make_dirichlet_uniform_proposer(
        dirichlet_tups=dirichlet_tups,
        EstimatedParameters=EstimatedParameters,
        c_min=c_min,
        c_max=c_max
    )


def make_param_filter_func(
    c_max: EstimatedParameters,
    c_min: EstimatedParameters,
    fcpa: FreeConstants,
    dvs: Drivers,
    svs: msh.Observables
) -> Callable[[np.ndarray], bool]:

    def isQualified(
            c,
            print_conds=False
        ):
        epa=EstimatedParameters(*c)
        def value(field_name):
            try:
                return c[EstimatedParameters._fields.index(field_name)]
            except Exception as e:
                print("###########################")
                print(e)
                print(field_name)
                raise e

        conds=[
            (c >= c_min).all(),
            (c <= c_max).all(),
            sum(map(value, ["beta_leaf", "beta_wood"])) <= 0.99,
            sum(map(value, ["c_leaf_0", "c_wood_0"])) <= svs.cVeg[0],
            sum(map(value, ['c_DPM_0', 'c_RPM_0', 'c_BIO_0'])) <= svs.cSoil[0]
        ]

        # from IPython import embed; embed()
        print(conds)
        res=all(conds)
        if print_conds:
            if not res:
                print(conds)
        return res

    return isQualified
