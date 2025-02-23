import numpy as np
from scipy import interpolate, sparse
from tqdm import tqdm
from typing import Callable, Dict, Tuple, Iterable, List
from functools import reduce, lru_cache
from copy import copy, deepcopy
from itertools import islice
from time import time
from sympy.core.expr import Expr
from collections import namedtuple
from frozendict import frozendict
from importlib import import_module
from collections import OrderedDict
import matplotlib.pyplot as plt
import os
import datetime as dt
import pandas as pd
import inspect
from pathlib import Path
import json
import netCDF4 as nc
import pyresample  # new package: to update bgc_md2 run install_developer_conda script
from pyresample.bilinear import NumpyBilinearResampler
import paramiko
import tarfile
import gzip
import shutil
from importlib.resources import files as mod_files

from sympy import (
    Symbol,
    Function,
    sympify,
    simplify,
    lambdify,
    diff,
    exp,
    var,
    sin,
    Min,
    Max,
    pi,
    integrate,
    diag,
)
from ComputabilityGraphs.CMTVS import CMTVS
import CompartmentalSystems.helpers_reservoir as hr
from CompartmentalSystems.start_distributions import (
    # start_age_moments_from_empty_spinup,
    start_age_moments_from_steady_state,
    # start_age_moments_from_zero_initial_content,
    # compute_fixedpoint_numerically,
    start_age_distributions_from_steady_state,
    # start_age_distributions_from_empty_spinup,
    # start_age_distributions_from_zero_initial_content,
)

# from CompartmentalSystems.BlockArrayIterator import BlockArrayIterator
from CompartmentalSystems.BlockDictIterator import BlockDictIterator
from CompartmentalSystems.ArrayDict import ArrayDict
from CompartmentalSystems.ArrayDictResult import ArrayDictResult
from bgc_md2.resolve.mvars import (
    TimeSymbol,
    CompartmentalMatrix,
    InputTuple,
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    StateVariableTuple,
    NumericParameterization,
    NumericSimulationTimes,
    NumericStartValueArray,
)

# from bgc_md2.helper import bgc_md2_computers
import bgc_md2.display_helpers as dh
import bgc_md2.helper as h
from typing import TypeVar, Callable, Union, NewType

T = TypeVar("T")
# Array=Union(np.ma.core.MaskedArray,np.ndarray)
NewType("Array", np.ndarray)

boundaries = namedtuple("boundaries", ["min_lat", "max_lat", "min_lon", "max_lon"])
Transformers = namedtuple(
    "Transformers",
    [
        "i2lat",
        # "i2lat_min_max",
        # "lat2i",
        "i2lon",
        # "i2lon_min_max",
        # "lon2i",
    ],
)

CoordTransformers = namedtuple(
    "CoordTransformers", ["lat2LAT", "LAT2lat", "lon2LON", "LON2lon"]
)
# def __init__(self, year, month, day):
#    self.year=year
#    self.month=month
#    self.day=day

# @classmethod
# def from_trendy_days_since_AD(cls,tdsDA):
#    y = int(tdsDA/self.days_per_year)
#    m = int((tdsDA- self.days_per_year * y) / self.days_per_month)
#    d = tdsDA - y * self.days_per_year - m * self.days_per_month
#    #from IPython import embed; embed()
#    return cls(y,m+1,d+1)

# @property
# def trendy_days_since_AD(self) -> np.array:

#    return (
#        self.year * self.days_per_year
#        + self.days_per_month * (self.month - 1)
#        + self.day-1
#    )

# def __eq__(self,other):
#    return all(
#        [
#            self.year==other.year,
#            self.month==other.month,
#            self.day==other.day,
#        ]
#    )
# def __repr__(self):
#    return f"{self.month}-{self.day} {self.year})"


def compose_2(f: Callable, g: Callable) -> Callable:
    """Function composition"""
    return lambda arg: f(g(arg))


def globalMaskTransformers(mask: np.ndarray) -> "SymTransformers":
    n_lats, n_lons = mask.shape
    # this functions define how the indices of the combined (all models) mask
    # relates to our unified coord system LAT, LON
    # with -90< LAT < 90 , 0 at equator
    # and -180 < LON < 180, 0 at Greenwich
    step_lat = 180.0 / n_lats
    step_lon = 360.0 / n_lons

    lat_0 = -90 + step_lat / 2
    lon_0 = -180 + step_lon / 2
    itr = transform_maker(
        lat_0,
        lon_0,
        step_lat,
        step_lon,
    )
    # here we use the identical transformation
    ctr = identicalTransformers()
    return SymTransformers(itr=itr, ctr=ctr)


def globalMask(file_name="common_mask.nc") -> "CoordMask":
    # fixme mm 6-20-2022
    # this should actually be  package data
    # We create a target
    this_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ds = nc.Dataset(Path(this_dir).joinpath(file_name))
    mask = ds.variables["mask"][:, :].data
    return CoordMask(index_mask=mask, tr=globalMaskTransformers(mask))


def identicalTransformers():
    """This function can be used if the grid used by the model has
    - lat ==   0 at the equator with
    - lat == -90 at the south pole,
    - lat == +90 at the north pole,
    - lon ==   0 at Greenwich and
    - lon is counted positive eastwards from -180 to 180
    """
    return CoordTransformers(
        lat2LAT=lambda lat: lat,
        LAT2lat=lambda LAT: LAT,
        lon2LON=lambda lon: lon,
        LON2lon=lambda LON: LON,
    )


class SymTransformers:
    """as Transformers but with a fixed coord system
    with
    -90 <= lat <= 90
    -180 <= lon <= 180
    lat=0,Equator
    lon=0,Greenich
    """

    def __init__(self, itr: Transformers, ctr: CoordTransformers):
        self.i2lat = compose_2(ctr.lat2LAT, itr.i2lat)
        # self.i2lat_min_max=lambda i: map(ctr.lat2LAT,itr.i2lat_min_max(i))
        # self.lat2i=compose_2(itr.lat2i,ctr.LAT2lat)

        self.i2lon = compose_2(ctr.lon2LON, itr.i2lon)
        # self.i2lon_min_max=lambda i: map(ctr.lon2LON,itr.i2lon_min_max(i))
        # self.lon2i=compose_2(itr.lon2i,ctr.LON2lon)


# def pixel_boundaries(lat,lon,step_lat,step_lon):
#    return boundaries(
#        min_lat=lat-step_lat/2,
#        max_lat=lat+step_lat/2,
#        min_lon=lon-step_lon/2,
#        max_lon=lon+step_lon/2
#    )

# _TraceTuple = namedtuple(
#    "_TraceTuple",
#    [
#        "x",
#        "x_p",
#        "x_c",
#        "x_dot",
#        "rt",
#        #
#        "x",
#        "x_p",
#        "x_c",
#        "x_dot",
#        "rt",
#        "u",
#        "AggregatedVegetation2SoilCarbonFlux",
#    ],
# )


# some tiny helper functions for module loading
def mvs(mf):
    return import_module(f"{msh(mf).model_mod}.source").mvs


def msh(mf):
    return import_module(f"{__package__}.{mf}.model_specific_helpers_2")


def th(mf):
    return import_module(f"{__package__}.{mf}.test_helpers")


def confDict(mf):
    with mod_files(f"trendy9helpers.{mf}").joinpath("config.json").open("r") as f:
        confDict = frozendict(json.load(f))
    return confDict


def da_dir_path(p, mf, da_scheme, par_sub_dir):
    return p.joinpath(mf, da_scheme, par_sub_dir)


def data_path(mf):
    return Path(confDict(mf)["dataPath"])


def target_path(p, mf):
    return p.joinpath(mf, "global_means")


def da_mod_name(mf, da_scheme):
    return f"trendy9helpers.{mf}.{da_scheme}"


def da_mod(mf, da_scheme):
    return import_module(f"{da_mod_name(mf,da_scheme)}.mod")


def output_cache_path(p, mf, da_scheme, par_dir):
    return da_dir_path(p, mf, da_scheme, par_dir).joinpath("out")


def da_param_path(p, mf, da_scheme, par_dir):
    return da_dir_path(p, mf, da_scheme, par_dir).joinpath("in")


def gm_cp_from_folder_names(
    # just a convenience function implementing the folder hierarchy
    # just a convenience function implementing the folder hirarchy
    # p/
    #   model/
    #       da_scheme/
    #           param_dir/
    #               in/
    #                   cpa.json
    #                   epa_min.json
    #                   epa_max.json
    #                   epa_0.json
    #                   hyper.json
    #
    # p/
    #   model/
    #       da_scheme/
    #           param_dir/
    #               out/
    #                   epa_opt.json
    #
    # for the da modules and parameter dirs to be used in different
    # notebooks and scripts
    # It is not intended as a central funtion to build much on top of.
    # So please don't! ;-)
    p: Path,
    mf: str,
    da_scheme: str,
    param_dir: str,
):
    state_vector = mvs(mf).get_StateVariableTuple()
    # load global mean vars
    # data path not necessary unless global means are recomputed
    svs, dvs = msh(mf).get_global_mean_vars(
        data_path(mf), target_path(p, mf), flash_cache=False
    )
    #
    # The following lines perform data_assimilation

    cpa = da_mod(mf, da_scheme).Constants(
        **h.load_dict_from_json_path(
            da_param_path(p, mf, da_scheme, param_dir).joinpath("cpa.json")
        )
    )
    epa_opt = da_mod(mf, da_scheme).EstimatedParameters(
        **h.load_dict_from_json_path(
            output_cache_path(p, mf, da_scheme, param_dir).joinpath("epa_opt.json")
        )
    )
    apa = {**cpa._asdict(), **epa._asdict()}
    # create a parameterization for convinience and return it
    param_dict = make_param_dict(mvs(mf), cpa, epa_opt)
    X_0 = da_mod(mf, da_scheme).numeric_X_0(mvs(mf), dvs, cpa, epa_opt)
    X_0_dict = {
        str(sym): X_0[i, 0] for i, sym in enumerate(mvs(mf).get_StateVariableTuple())
    }
    # some models (e.g. yz_jules) need extra (not represented by symbols) parameters  to build
    # the func_dict for the parameterization
    apa = {**cpa._asdict(), **epa_opt._asdict()}
    func_dict_param_dict = {
        str(k): v
        for k, v in apa.items()
        if str(k) in msh(mf).CachedParameterization.func_dict_param_keys
    }
    cp = msh(mf).CachedParameterization(param_dict, dvs, X_0_dict, func_dict_param_dict)
    return cp


def gm_da_from_folder_names(
    # just a convenience function implementing the folder hierarchy
    # p/
    #   model/
    #       da_scheme/
    #           param_dir/
    #               in/
    #                   cpa.json
    #                   epa_min.json
    #                   epa_max.json
    #                   epa_0.json
    #                   hyper.json
    #
    # p/
    #   model/
    #       da_scheme/
    #           param_dir/
    #               out/
    #                   epa_opt.json
    #
    # where da_scheme is a submodule of trendy9helpers/model/ (e.g. da_1, da_2)
    # just for notebooks and scripts
    # It is not intended as a central funtion to build much on top of.
    # So please don't! ;-)
    p: Path,
    mf: str,
    da_scheme: str,
    param_dir: str,
):
    state_vector = mvs(mf).get_StateVariableTuple()
    # load global mean vars
    # data path not necessary unless global means are recomputed
    svs, dvs = msh(mf).get_global_mean_vars(
        data_path(mf), target_path(p, mf), flash_cache=False
    )
    #
    # The following lines perform data_assimilation
    # The following lines perform data_assimilation

    func = cached_da_res_1_maker(
        da_mod(mf, da_scheme).make_proposer,
        da_mod(mf, da_scheme).make_param_filter_func,
        da_mod(mf, da_scheme).make_param2res_sym,
        msh(mf).make_weighted_cost_func,
        da_mod(mf, da_scheme).numeric_X_0,
        da_mod(mf, da_scheme).EstimatedParameters,
    )
    # specify from where to read the start parameters ranges and constants
    sp = da_dir_path(p, mf, da_scheme, param_dir)
    if not sp.exists():
        # sp.mkdir(parents=True,exist_ok=True)
        # copy example files for the tests as part of the package
        # If the directory is empty we copy the one used for the da tests which
        # are part of the package.
        da_dir_ex_path = mod_files(da_mod_name(mf, da_scheme)).joinpath(param_dir)
        print(da_dir_ex_path)
        # from IPython import embed; embed()
        shutil.copytree(
            da_dir_ex_path, da_dir_path(p, mf, da_scheme, param_dir), dirs_exist_ok=True
        )

    cpa = da_mod(mf, da_scheme).FreeConstants(
        **h.load_dict_from_json_path(
            da_param_path(p, mf, da_scheme, param_dir).joinpath("cpa.json")
        )
    )
    fepa_min, fepa_max, fepa_0 = tuple(
        map(
            lambda p: da_mod(mf, da_scheme).FreeEstimatedParameters(
                **h.load_dict_from_json_path(p)
            ),
            [
                da_param_path(p, mf, da_scheme, param_dir).joinpath(f"{s}.json")
                for s in ["epa_min", "epa_max", "epa_0"]
            ],
        )
    )
    # check that we can read the hyperparameters
    hyper_dict = h.load_dict_from_json_path(
        da_param_path(p, mf, da_scheme, param_dir).joinpath("hyper.json")
    )
    # compute the chains and the optimal parameter tuple
    # (and automatically cache them)
    Cs, Js, epa_opt = func(
        output_cache_path(p, mf, da_scheme, param_dir),
        mvs(mf),
        svs,
        dvs,
        cpa,
        da_mod(mf, da_scheme).epa_min(fepa_min, dvs, svs),
        da_mod(mf, da_scheme).epa_max(fepa_max, dvs, svs),
        da_mod(mf, da_scheme).epa_0(fepa_0, dvs, svs),
        nsimu=hyper_dict["nsimu"],
        acceptance_rate=hyper_dict["acceptance_rate"],
        chunk_size=hyper_dict["chunk_size"],
        D_init=hyper_dict["D_init"],
        K=hyper_dict["K"],
    )
    # create a parameterization for convinience and return it
    param_dict = make_param_dict(mvs(mf), cpa, epa_opt)

    apa = {**cpa._asdict(), **epa_opt._asdict()}
    X_0 = da_mod(mf, da_scheme).numeric_X_0(mvs(mf), dvs, apa) #, cpa, epa_opt)

    X_0_dict = {
        str(sym): X_0[i, 0] for i, sym in enumerate(mvs(mf).get_StateVariableTuple())
    }
    # some models (e.g. yz_jules) need extra (not represented by symbols) parameters  to build
    # the func_dict for the parameterization
    apa = {**cpa._asdict(), **epa_opt._asdict()}
    func_dict_param_dict = {
        str(k): v
        for k, v in apa.items()
        if str(k) in msh(mf).CachedParameterization.func_dict_param_keys
    }
    cp = msh(mf).CachedParameterization(param_dict, dvs, X_0_dict, func_dict_param_dict)
    return Cs, Js, epa_opt, cp


# should be part  of CompartmentalSystems
def make_B_u_funcs(mvs, mpa, func_dict):
    model_params = {Symbol(k): v for k, v in mpa._asdict().items()}
    return make_B_u_funcs_2(mvs, model_params, func_dict)


def make_B_u_funcs_2(mvs, model_params, func_dict, delta_t_val=1):
    # symbol_names = mvs.get_BibInfo().sym_dict.keys()
    # for name in symbol_names:
    #    var(name)
    t = mvs.get_TimeSymbol()
    it = Symbol("it")
    delta_t = Symbol("delta_t")
    parameter_dict = {**model_params, delta_t: delta_t_val}
    state_vector = mvs.get_StateVariableTuple()

    sym_B = hr.discrete_time_sym(  # hr.euler_forward_B_sym(
        mvs.get_CompartmentalMatrix(), cont_time=t, delta_t=delta_t, iteration=it
    )
    sym_u = hr.discrete_time_sym(mvs.get_InputTuple(), t, delta_t, it)

    B_func = hr.numerical_array_func(
        state_vector=state_vector,
        time_symbol=it,
        expr=sym_B,
        parameter_dict=parameter_dict,
        func_dict=func_dict,
    )
    # u_func = hr.numerical_1d_vector_func(
    u_func = hr.numerical_array_func(
        state_vector=state_vector,
        time_symbol=it,
        expr=sym_u,
        parameter_dict=parameter_dict,
        func_dict=func_dict,
    )
    return (B_func, u_func)


# fixme mm 11-15-2022
# should be obsolete with the new iterator
def numfunc(expr_cont, mvs, delta_t_val, par_dict, func_dict):
    # This function
    # build the discrete expression (which depends on it,delta_t instead of
    # the continius one that depends on t (TimeSymbol))
    t = mvs.get_TimeSymbol()
    it = Symbol("it")  # arbitrary symbol for the step index )
    delta_t = Symbol("delta_t")
    expr_disc = expr_cont.subs({t: delta_t * it})
    expr_num = expr_disc.subs({delta_t: delta_t_val})
    # print(expr_cont,expr_disc,expr_num)
    return hr.numerical_function_from_expression(
        expr=expr_num,
        tup=(it, *mvs.get_StateVariableTuple()),
        parameter_dict=par_dict,
        func_set=func_dict,
    )


def make_param_dict(mvs, cpa, epa):
    # legacy wrapper
    # Parameter dictionary for the iterator
    apa = {**cpa._asdict(), **epa._asdict()}
    return make_param_dict2(mvs, apa)


def make_param_dict2(mvs: CMTVS, apa: Dict) -> Dict:
    srm = mvs.get_SmoothReservoirModel()
    model_par_dict_keys = srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))] + list(mvs.get_StateVariableTuple())
    )
    # Parameter dictionary for the iterator
    model_par_dict = {
        Symbol(k): v for k, v in apa.items() if Symbol(k) in model_par_dict_keys
    }
    return model_par_dict


def make_dirichlet_uniform_proposer(
    dirichlet_tups: Dict[List[str], float],
    EstimatedParameters: type,  # model and da_scheme specific
    c_min,
    c_max,
):
    # uses dirichlet distribution for the parameters that should add up
    # to a constant (like the betas that have to sum up to one
    # or start values of the veg pools that have to sum up to the first
    # observed cVeg (same for soil)
    dirichlet_name_lists = [names for names, _ in dirichlet_tups]  #
    dirichlet_p_sums = [p_sums for _, p_sums in dirichlet_tups]
    # compute the positions of the names in EstimatedParameters
    # in order to avoid this in every call of the proposer
    rest_keys = set(EstimatedParameters._fields).difference(
        sum(dirichlet_name_lists, [])
    )

    def ind(names):
        try:
            positions = [EstimatedParameters._fields.index(name) for name in names]
            return positions
        except ValueError:
            for name in names:
                try:
                    EstimatedParameters._fields.index(name)
                except ValueError:
                    raise ValueError(f"{name} is not in fields")

    pss = list(map(ind, [*dirichlet_name_lists, rest_keys]))
    dirichlet_pss = pss[:-1]
    rest_ps = pss[-1]

    paramNum = len(rest_keys)

    def GenerateParamValues(c_op, D, print_conds=False):
        p = 1.0 / D

        def dirichlet_update(positions, p_sum):
            v = np.array([c_op[p] for p in positions])
            # v_min=np.array(
            #    [
            #        c_min[p] for p in positions
            #    ]
            # )
            # v_max=np.array(
            #    [
            #        c_max[p] for p in positions
            #    ]
            # )
            x = np.random.dirichlet(alpha=(1,) * (len(v) + 1))[0:-1] * p_sum
            return (1 - p) * v + p * x

        valss = [
            dirichlet_update(ps, p_sum)
            for ps, p_sum in zip(dirichlet_pss, dirichlet_p_sums)
        ]
        dirichlet_val_pos_tups = zip(valss, dirichlet_pss)

        cu_op, cu_min, cu_max = map(lambda c: c[rest_ps], [c_op, c_min, c_max])
        cu_new = (1 - p) * cu_op + p * np.random.uniform(cu_min, cu_max, paramNum)

        ###########################
        cp_new = np.zeros_like(c_op)
        for vals, positions in [*dirichlet_val_pos_tups, (cu_new, rest_ps)]:
            cp_new[positions] = vals
        ###########################
        return cp_new

    return GenerateParamValues


def make_uniform_proposer_2(
    c_max: Iterable,
    c_min: Iterable,
    filter_func: Callable[[np.ndarray], bool],
) -> Callable[[Iterable], Iterable]:
    """Returns a function that will be used by the mcmc algorithm to propose
    a new parameter value tuple based on a given one.
    The two arrays c_max and c_min define the boundaries
    of the n-dimensional rectangular domain for the parameters and must be of
    the same shape.  After a possible parameter value has been sampled the
    filter_func will be applied to it to either accept or discard it.  So
    filter func must accept parameter array and return either True or False
    :param c_max: array of maximum parameter values
    :param c_min: array of minimum parameter values
    :param D: a parameter to regulate the proposer step .
    smaller D means smaller step size.
    :param filter_func: model-specific function to filter out impossible parameter combinations
    """

    g = np.random.default_rng()

    def GenerateParamValues(c_op, D, print_conds=False):
        paramNum = len(c_op)
        keep_searching = True
        counter = 0
        print_conds = False
        while keep_searching:
            counter += 1
            # c_new = c_op + (np.random.uniform(-0.5, 0.5, paramNum) * c_op * D) # does not stay between c_min,c_max
            c_new = c_op + (np.random.uniform(c_min, c_max, paramNum) - c_op) / D
            if counter % 10:
                print(
                    f"warning:{counter} unsuccsessful attempts to propose new parameters"
                )
                print_conds = True

            if filter_func(c_new, print_conds=print_conds):
                keep_searching = False
        return c_new

    return GenerateParamValues


def accept_costfunction(J_last: float, J_new: float, K=2):
    """Regulates how new cost functions are accepted or rejected. If the new cost function is lower than the old one,
    it is always accepted. If the the new cost function is higher than the old one, it has a random
    chance to be accepted based on percentage difference between the old and the new. The chance is defined
    by an exponential function and regulated by the K coefficient.
    :param J_last: old (last accepted) cost function
    :param J_new: new cost function
    :param K: regulates acceptance chance. Default 1 means that a 1% higher cost function has 37% chance to be accepted.
    Increase K to reduce the chance to accept higher cost functions
    """
    accept = False
    delta_J_percent = (
        (J_last - J_new) / J_last * 100
    )  # normalize delta_J as a percentage of current J
    randNum = np.random.uniform(0, 1)
    if (
        min(1.0, np.exp(delta_J_percent * K)) > randNum
    ):  # 1% higher cost function has 14% chance to be accepted
        accept = True
    return accept


# Autostep MCMC: modifying its step every 100 iterations
# depending on acceptance rate
# with uniform proposer
def autostep_mcmc_2(
    initial_parameters: Iterable,
    proposer: Callable,
    filter_func: Callable[[np.ndarray, bool], bool],
    param2res: Callable[[np.ndarray], np.ndarray],
    costfunction: Callable[[np.ndarray], np.float64],
    nsimu: int,
    c_max: np.ndarray,
    c_min: np.ndarray,
    acceptance_rate=0.23,
    chunk_size=100,
    D_init=5,
    K=1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    performs the Markov chain Monte Carlo simulation an returns a tuple of the array of sampled parameter(tuples)
    with shape (len(initial_parameters),nsimu) and the array of costfunction values with shape (q,nsimu)

    :param initial_parameters: The initial guess for the parameter (tuple) to be estimated
    :param filter_func: model-specific function to filter out impossible parameter combinations
    :param param2res: A function that given a parameter(tuple) returns
    the model output, which has to be an array of the same shape as the observations used to
    build the cost function.
    :param costfunction: A function that given a model output returns a real number. It is assumed to be created for a
    specific set of observations, which is why they do not appear as an argument.
    :param nsimu: The length of the chain
    :param c_max: Array of maximum values for each parameter
    :param c_min: Array of minimum values for each parameter
    :param acceptance_rate: Target acceptance rate in %
    :param chunk_size: number of iterations for which current acceptance ratio is assessed to modify the proposer step
    Set to 0 for constant step size. Default is 100.
    :param D_init: initial D value (decrease to get a bigger step size must be larger than 1 which is equivalent to using the full parameter range for the next guess)
    :param K: allowance coeffitient (increase K to reduce acceptance of higher cost functions)
    """
    # if D <1 then the proposer would produce values
    # outside the parameter range
    assert D_init >= 1

    np.random.seed(seed=10)

    # print(model_par_dict)
    C_op = initial_parameters

    # check if the initial parameters are valid (to avoid an infinite loop before we even start changing them)
    assert filter_func(
        C_op,
        print_conds=True,
    ), "initial conditions don't match the filter"

    tb = time()
    first_out = param2res(C_op)
    J_last = costfunction(first_out)
    J_min = J_last
    J_min_simu = 0
    print("First_iteration done after " + str(round(time() - tb, 2)) + "s")
    print("Status update every 10 iterations:")
    # J_last = 400 # original code

    # initialize the result arrays to the maximum length
    # Depending on many of the parameters will be accepted only
    # a part of them will be filled with real values
    paramNum = len(initial_parameters)
    C_upgraded = np.zeros((paramNum, nsimu))
    J_upgraded = np.zeros((2, nsimu))
    D = D_init
    st = time()
    if chunk_size == 0:
        chunk_size = (
            nsimu  # if chunk_size is set to 0 - proceed without updating step size.
        )

    upgraded = 0
    accepted_in_last_chunk = 0
    for simu in tqdm(range(nsimu)):
        if simu % 10 == 0 or simu == (nsimu - 1):
            print(
                f""" 
               #(upgraded): {upgraded}   
               D value: {D} 
               target acceptance rate: {acceptance_rate}%
               current chunk acceptance rate: {int(accepted_in_last_chunk/chunk_size  * 100)}%  
               overall acceptance rate: {int(upgraded / (simu + 1) * 100)}%  
               progress: {simu:05d}/{nsimu:05d} {int(simu / (nsimu - 1) * 100):02d}%
               time elapsed: {int((time() - st) / 60):02d}:{int((time() - st) % 60):02d}
               overall min cost: {round(J_min, 2)} achieved at {J_min_simu} iteration | last accepted cost: {round(J_last, 2)} 
               """
            )

        if (simu > 0) and (
            simu % chunk_size == 0
        ):  # every chunk size (e.g. 100 iterations) update the proposer step
            chunk_acceptance_rate = accepted_in_last_chunk / chunk_size
            D = D * 1.2 if chunk_acceptance_rate < acceptance_rate else max(1, D * 0.8)
            print(f"D={D}")
            accepted_in_last_chunk = 0
        if (
            simu % (chunk_size * 20) == 0
        ):  # every 20 chunks - return to the initial step size (to avoid local minimum)
            D = D_init

        keep_searching = True
        while keep_searching:
            c_new = proposer(C_op, D)
            if filter_func(c_new, print_conds=True):
                keep_searching = False
        out_simu = param2res(c_new)
        J_new = costfunction(out_simu)

        if accept_costfunction(J_last=J_last, J_new=J_new, K=K):
            C_op = c_new
            J_last = J_new
            if J_last < J_min:
                J_min = J_last
                J_min_simu = simu
            C_upgraded[:, upgraded] = C_op
            J_upgraded[1, upgraded] = J_last
            J_upgraded[0, upgraded] = simu
            upgraded += 1
            accepted_in_last_chunk += 1

        # print some metadata
        # (This could be added to the output file later)
    # remove the part of the arrays that is still filled with zeros
    useful_slice = slice(0, upgraded)
    return C_upgraded[:, useful_slice], J_upgraded[:, useful_slice]


# def adaptive_mcmc(
#    # Adaptive MCMC: with multivariate normal proposer based on adaptive covariance matrix
#    initial_parameters: Iterable,
#    covv: np.ndarray,
#    filter_func: Callable,
#    param2res: Callable[[np.ndarray], np.ndarray],
#    costfunction: Callable[[np.ndarray], np.float64],
#    nsimu: int,
#    sd_controlling_factor=10,
# ) -> Tuple[np.ndarray, np.ndarray]:
#    """
#    performs the Markov chain Monte Carlo simulation an returns a tuple of the array of sampled parameter(tuples) with
#    shape (len(initial_parameters),nsimu) and the array of cost function values with shape (q,nsimu)
#    :param initial_parameters: The initial guess for the parameter (tuple) to be estimated
#    :param covv: The covariance matrix (usually estimated from a previously run chain)
#    :param filter_func: function to remove impossible parameter combinations
#    :param param2res: A function that given a parameter(tuple) returns
#    the model output, which has to be an array of the same shape as the observations used to
#    build the cost function.
#    :param costfunction: A function that given a model output returns a real number. It is assumed to be created for
#    a specific set of observations, which is why they do not appear as an argument.
#    :param nsimu: The length of the chain
#    :param sd_controlling_factor: optional parameter to scale the covariance matrix. Increase to get a smaller step size
#    """
#
#    np.random.seed(seed=10)
#
#    paramNum = len(initial_parameters)
#
#    sd = 1 / sd_controlling_factor / paramNum
#    covv = covv * sd
#
#    proposer = make_multivariate_normal_proposer(covv, filter_func)
#
#    upgraded = 0
#    C_op = initial_parameters
#    tb = time()
#    first_out = param2res(C_op)
#
#    J_last = costfunction(first_out)
#    J_min = J_last
#    J_min_simu = 0
#    print("first_iteration done after " + str(time() - tb))
#    # J_last = 400 # original code
#
#    # initialize the result arrays to the maximum length
#    # Depending on many of the parameters will be accepted only
#    # a part of them will be filled with real values
#    C_upgraded = np.zeros((paramNum, nsimu))
#    J_upgraded = np.zeros((2, nsimu))
#
#    # for simu in tqdm(range(nsimu)):
#    st = time()
#    for simu in range(nsimu):
#        # if (upgraded%10 == 0) & (upgraded > nsimu/20):
#        if simu > nsimu / 10:
#            covv = sd * np.cov(C_accepted)
#            proposer = make_multivariate_normal_proposer(covv, filter_func)
#        c_new = proposer(C_op)
#        out_simu = param2res(c_new)
#        J_new = costfunction(out_simu)
#
#        if accept_costfunction(J_last=J_last, J_new=J_new):
#            C_op = c_new
#            J_last = J_new
#            if J_last < J_min:
#                J_min = J_last
#                J_min_simu = simu
#            C_upgraded[:, upgraded] = C_op
#            C_accepted = C_upgraded[:, 0:upgraded]
#            J_upgraded[1, upgraded] = J_last
#            J_upgraded[0, upgraded] = simu
#            upgraded = upgraded + 1
#        # print some metadata
#        # (This could be added to the output file later)
#
#        if simu % 10 == 0 or simu == (nsimu - 1):
#            print(
#                """
##(upgraded): {n}
# overall acceptance ratio till now: {r}%
# progress: {simu:05d}/{nsimu:05d} {pbs} {p:02d}%
# time elapsed: {minutes:02d}:{sec:02d}
# overall minimum cost: {cost} achieved at {s} iteration | last accepted cost: {cost2}
# """.format(
#                    n=upgraded,
#                    r=int(upgraded / (simu + 1) * 100),
#                    simu=simu,
#                    nsimu=nsimu,
#                    pbs="|"
#                    + int(50 * simu / (nsimu - 1)) * "#"
#                    + int((1 - simu / (nsimu - 1)) * 50) * " "
#                    + "|",
#                    p=int(simu / (nsimu - 1) * 100),
#                    minutes=int((time() - st) / 60),
#                    sec=int((time() - st) % 60),
#                    cost=round(J_min, 2),
#                    cost2=round(J_last, 2),
#                    s=J_min_simu,
#                ),
#                end="\033[5A",  # print always on the same spot of the screen...
#            )
#
#    # remove the part of the arrays that is still filled with zeros
#    useful_slice = slice(0, upgraded)
#    return C_upgraded[:, useful_slice], J_upgraded[:, useful_slice]


def mcmc(
    initial_parameters: Iterable,
    proposer: Callable[[Iterable], Iterable],
    # param2res: Callable[[np.ndarray], np.ndarray],
    # costfunction: Callable[[np.ndarray], np.float64],
    param2res: Callable,
    costfunction: Callable,
    nsimu: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    performs the Markov chain Monte Carlo simulation and returns a tuple of the array of sampled parameter(tuples) with
    shape (len(initial_parameters),nsimu) and the array of costfunction values with shape (q,nsimu)

    :param initial_parameters: The initial guess for the parameter (tuple) to be estimated
    :param proposer: A function that proposes a new parameter(tuple) from a given parameter (tuple).
    :param param2res: A function that given a parameter(tuple) returns
    the model output, which has to be an array of the same shape as the observations used to
    build the costfunction.
    :param costfunction: A function that given a model output returns a real number. It is assumed to be created for
    a specific set of observations, which is why they do not appear as an argument.
    :param nsimu: The length of the chain
    """
    np.random.seed(seed=10)

    paramNum = len(initial_parameters)

    upgraded = 0
    C_op = initial_parameters
    tb = time()
    first_out = param2res(C_op)
    J_last = costfunction(first_out)
    J_min = J_last
    J_min_simu = 0
    print("first_iteration done after " + str(time() - tb))
    # J_last = 400 # original code

    # initialize the result arrays to the maximum length
    # Depending on many of the parameters will be accepted only
    # a part of them will be filled with real values
    C_upgraded = np.zeros((paramNum, nsimu))
    J_upgraded = np.zeros((2, nsimu))

    # for simu in tqdm(range(nsimu)):
    st = time()
    for simu in range(nsimu):
        c_new = proposer(C_op)
        out_simu = param2res(c_new)
        J_new = costfunction(out_simu)

        if accept_costfunction(J_last=J_last, J_new=J_new):
            C_op = c_new
            J_last = J_new
            if J_last < J_min:
                J_min = J_last
                J_min_simu = simu
            C_upgraded[:, upgraded] = C_op
            J_upgraded[1, upgraded] = J_last
            J_upgraded[0, upgraded] = simu
            upgraded = upgraded + 1
        # print some metadata
        # (This could be added to the output file later)

        if simu % 10 == 0 or simu == (nsimu - 1):
            print(
                """ 
#(upgraded): {n}
overall acceptance ratio till now: {r}% 
progress: {simu:05d}/{nsimu:05d} {pbs} {p:02d}%
time elapsed: {minutes:02d}:{sec:02d}
overall minimum cost: {cost} achieved at {s} iteration | last accepted cost: {cost2} 
""".format(
                    n=upgraded,
                    r=int(upgraded / (simu + 1) * 100),
                    simu=simu,
                    nsimu=nsimu,
                    pbs="|"
                    + int(50 * simu / (nsimu - 1)) * "#"
                    + int((1 - simu / (nsimu - 1)) * 50) * " "
                    + "|",
                    p=int(simu / (nsimu - 1) * 100),
                    minutes=int((time() - st) / 60),
                    sec=int((time() - st) % 60),
                    cost=round(J_min, 2),
                    cost2=round(J_last, 2),
                    s=J_min_simu,
                ),
                end="\033[5A",  # print always on the same spot of the screen...
            )

    # remove the part of the arryas that is still filled with zeros
    useful_slice = slice(0, upgraded)
    return C_upgraded[:, useful_slice], J_upgraded[:, useful_slice]


# def make_feng_cost_func(
#    obs: np.ndarray,
# ) -> Callable[[np.ndarray], np.float64]:
#    # Note:
#    # in our code the dimension 0 is the time
#    # and dimension 1 the pool index
#    means = obs.mean(axis=0)
#    mean_centered_obs = obs - means
#    # now we compute a scaling factor per observable stream
#    # fixme mm 10-28-2021
#    #   The denominators in this case are actually the TEMPORAL variances of the data streams
#    denominators = np.sum(mean_centered_obs**2, axis=0)
#
#    #   The desired effect of automatically adjusting weight could be achieved
#    #   by the mean itself.
#    # dominators = means
#    def costfunction(mod: np.ndarray) -> np.float64:
#        cost = np.mean(np.sum((obs - mod) ** 2, axis=0) / denominators * 100)
#        return cost
#
#    return costfunction


# def make_jon_cost_func(
#    obs: np.ndarray,
# ) -> Callable[[np.ndarray], np.float64]:
#    # Note:
#    # in our code the dimension 0 is the time
#    # and dimension 1 the pool index
#    n = obs.shape[0]
#    means = obs.mean(axis=0)
#    denominators = means**2
#
#    def costfunction(mod: np.ndarray) -> np.float64:
#        cost = (100 / n) * np.sum(100 * np.sum((obs - mod) ** 2, axis=0) / denominators)
#        return cost
#
#    return costfunction
#
#
# def respiration_from_compartmental_matrix(B, X):
#    """This function computes the combined respiration from all pools"""
#    return -np.sum(B @ X)
#
#
# def plot_solutions(fig, times, var_names, tup, names=None):
#    if names is None:
#        names = tuple(str(i) for i in range(len(tup)))
#
#    assert all([tup[0].shape == el.shape for el in tup])
#
#    if tup[0].ndim == 1:
#        n_times = tup[0].shape[0]
#        ax = fig.subplots(1, 1)
#        for i, sol in enumerate(tup):
#            ax.plot(
#                np.array(times).reshape(
#                    n_times,
#                ),
#                sol,
#                marker="o",
#                label=names[i],
#            )
#            ax.set_title(var_names[0])
#            ax.legend()
#    else:
#        n_times, n_vars = tup[0].shape
#
#        fig.set_figheight(n_vars * fig.get_figwidth())
#        axs = fig.subplots(n_vars, 1)
#        colors = ("red", "blue", "green", "orange")
#        for j in range(n_vars):
#            for i, sol in enumerate(tup):
#                axs[j].plot(
#                    np.array(times).reshape(
#                        n_times,
#                    ),
#                    sol[:, j],
#                    marker="+",
#                    label=names[i],
#                )
#                axs[j].set_title(var_names[j])
#                axs[j].legend()
#
#
# def plot_observations_vs_simulations(fig, svs_cut, obs_simu):
#    # svn and obs_simu are both
#    n_plots = len(svs_cut)
#    fig.set_figheight(n_plots * fig.get_figwidth())
#    axs = fig.subplots(n_plots)
#    for i, name in enumerate(svs_cut.__class__._fields):
#        var = svs_cut[i]
#        var_simu = obs_simu[i]
#        axs[i].plot(range(len(var_simu)), var_simu, label="simulation", zorder=2)
#        axs[i].plot(range(len(var)), var, label="observation", zorder=1)
#        axs[i].legend()
#        axs[i].set_title(name)


def is_infinite_rec(chunk):
    return reduce_or_rec(np.logical_not(np.isfinite(chunk)))


def reduce_or_rec(bool_arr: np.ma.core.MaskedArray) -> np.ma.core.MaskedArray:
    # works recursively for any dimension (some variables are 4, most are 3 dimensional)
    # 2d array (lat,lon are always the last dimensions in our datasets)
    if bool_arr.ndim == 2:
        return bool_arr
    else:
        return reduce_or_rec(
            reduce(
                np.logical_or,
                map(lambda i: bool_arr[i], range(bool_arr.shape[0])),
            )
        )


def get_nan_pixel_mask(
    var: Union[nc._netCDF4.Variable, np.ma.core.MaskedArray]
) -> np.ndarray:
    try:
        # only works if var is a netCDF4.Variable
        print(var.name)
    except:
        AttributeError()

    # partition along the the chunks of the netcdf file (measured about a 100
    # times faster) and much less memory consuming than reading the complete
    # array which requiers more than 16GB RAM

    try:
        chunksizes = var.chunking()
        chunksize = 1 if chunksizes == "contiguous" else chunksizes[0]
    except:
        # chunksizes could still be None
        TypeError
        chunksize = 1
    # Note: no list comprehension but a lazy map object to save memory
    subarrs = map(
        lambda t: var[t[0] : t[1]],
        partitions(start=0, stop=var.shape[0], nr_acc=chunksize),
    )

    res = reduce(
        lambda acc, el: np.logical_or(acc, el), tqdm(map(is_infinite_rec, subarrs))
    )
    # res is a masked array where both res.data and res.mask
    # are boolean arrays. res.mask is the logical_or for any masked
    # array operation and in our case logical_or is also applied to the
    # data (in the reduce call in this function and is_finite_rec)

    # the ultimate result is a mask that masks every lat,lon pixel
    # that is already masked for any time or any layer
    # or has a nan entry at any time or any layer
    mask = np.logical_or(res.data, res.mask)

    # we check if the original mask was appropriate
    org_mask_ind = (
        tuple(0 for i in range(var.ndim - 2)) + (slice(None),) + (slice(None),)
    )
    org_mask = var[org_mask_ind].mask

    unmasked_nan_pixels = np.argwhere(np.logical_and(mask, np.logical_not(org_mask)))
    if len(unmasked_nan_pixels) > 0:
        print(
            "found {} previously unmasked NaN pixels {}".format(
                len(unmasked_nan_pixels), unmasked_nan_pixels
            )
        )
    return mask


def get_nan_pixels(var: nc._netCDF4.Variable):
    """We use a netCDF4.Variable to avoid having to load the whole array into memory"""
    Nn_t, N_lat, N_lon = var.shape
    cs = 30

    def f(I_lat, I_lon):
        n_lat = min(cs, N_lat - I_lat)
        n_lon = min(cs, N_lon - I_lon)
        chunk = var[:, I_lat : I_lat + n_lat, I_lon : I_lon + n_lon]
        return tuple(
            (
                (I_lat + i_lat, I_lon + i_lon)
                for i_lat in range(n_lat)
                for i_lon in range(n_lon)
                if np.isnan(chunk[:, i_lat, i_lon]).any()
            )
        )

    l = (
        f(I_lat, I_lon)
        for I_lat in range(0, N_lat, cs)
        for I_lon in range(0, N_lon, cs)
    )
    return reduce(lambda x, y: x + y, l)


def reduce_is_infinite_over_first_dim(arr):
    def f(acc, el):
        new_acc = np.logical_or(acc, el)
        new_ind = np.argwhere(np.logical_not(new_acc == acc))
        # if len(new_ind)>0:
        #    print(new_ind)
        return new_acc

    def make_sub_arr_mask(i):
        # print(i)
        return np.logical_not(np.isfinite(arr[i]))

    return reduce(f, map(make_sub_arr_mask, range(arr.shape[0])))


def get_weight_mat(lats: np.ma.core.MaskedArray, lons: np.ma.core.MaskedArray):
    # assuming an equidistant grid.
    delta_lat = (lats.max() - lats.min()) / (len(lats) - 1)
    delta_lon = (lons.max() - lons.min()) / (len(lons) - 1)
    pixel_area = make_pixel_area_on_unit_spehre(delta_lat, delta_lon)

    return np.array(
        [
            [pixel_area(lats[lat_ind]) for lon_ind in range(len(lons))]
            for lat_ind in range(len(lats))
        ]
    )


# note mm 4-13 2023:
# Only use the function if you can not
# use global_mean_var because it requires np.array arguments
# that have to be in memory while the global_mean_var allows
# netcdf variables and thus can work on objects bigger than
# memory.
def global_mean(
    lats: np.ma.core.MaskedArray,
    lons: np.ma.core.MaskedArray,
    arr: np.ma.core.MaskedArray,
) -> np.array:
    """As the signature shows this function expects masked arrays.
    These occure naturaly if netCDF4.Variables are sliced.
    e.g.
    ds = nc.Dataset("example.nc")
    var=ds.variables['npp'] #->type(var)=netCDF4._netCDF4.Variable
    arr=var[:,:,;] # type(arr)=np.ma.core.MaskedArray
    # or if we don't know the shape
    arr=var.__array__()
    The important thing is not to call this functions with the netcdfc.Variables but the arrays sliced from them.
    If the mean of a whole variable is requiered it is
    better to use the sister function global_mean_vars which works also for variables that don't fit into memory.
    """

    # copy the mask from the array (first time step)
    weight_mask = arr.mask[0, :, :] if arr.mask.any() else False

    weight_mat = np.ma.array(get_weight_mat(lats, lons), mask=weight_mask)
    wms = weight_mat.sum()
    # to compute the sum of weights we add only those weights that
    # do not correspond to an unmasked grid cell
    return ((weight_mat * arr).sum(axis=(1, 2)) / wms).data


def global_mean_var_with_resampled_mask(
    template: np.ndarray,
    ctr: CoordTransformers,
    itr: Transformers,
    lats: np.ma.core.MaskedArray,
    lons: np.ma.core.MaskedArray,
    var: nc._netCDF4.Variable,
    time_slice: slice = slice(None, None, None),  # equals [:]
):
    gcm = global_coord_mask_resampled(
        template=template,
        ctr=ctr,
        itr=itr,
    )
    mask = gcm.index_mask.__array__()
    return global_mean_var(lats, lons, mask, var, time_slice)


def global_mean_var(
    lats: np.ma.core.MaskedArray,
    lons: np.ma.core.MaskedArray,
    mask: np.array,
    var: nc._netCDF4.Variable,
    time_slice: slice = slice(None, None, None),  # equals [:]
) -> np.array:
    """As the signature shows this function expects a netCDF4.Variable
    This is basically metadata which allows us to compute the maean even
    if the whole array would not fit into memory.

    ds = nc.Dataset("example.nc")
    var=ds.variables['npp'] #->type(var)=netCDF4._netCDF4.Variable

    the mask array is used to block out extra pixels that are not
    masked in var
    """
    weight_mat = np.ma.array(get_weight_mat(lats, lons), mask=mask)

    # to compute the sum of weights we add only those weights that
    # do not correspond to an unmasked grid cell
    wms = weight_mat.sum()
    n_t = var.shape[0]  # we assume that time is the first dimension of every var
    ra = range(
        0 if time_slice.start is None else max(0, time_slice.start),
        n_t if time_slice.stop is None else min(n_t, time_slice.stop),
        1 if time_slice.step is None else time_slice.step,
    )

    # len(ra) takes care of the corner case:
    # ra.max-ra.min not a multiple of ra.step
    res = np.zeros(len(ra))
    print(ra)
    res_ind = 0
    for it in tqdm(ra):
        el = (weight_mat * var[it, :, :]).sum() / wms
        res[res_ind] = el
        res_ind += 1
    return res


def grad2rad(alpha_in_grad):
    return np.pi / 180 * alpha_in_grad


def make_pixel_area_on_unit_spehre(delta_lat, delta_lon, sym=False):
    # we compute the are of a delta_phi * delta_theta patch
    # on the unit ball centered around phi,theta
    # (which depends on theta but not
    # on phi)
    # the infinitesimal area element dA = sin(theta)*d_phi * d_theta
    # we have to integrate it from phi_min to phi_max
    # and from theta_min to theta_max
    if sym:
        # we can do this with sympy (for testing)
        for v in ("theta", "phi", "theta_min", "theta_max", "phi_min", "phi_max"):
            var(v)

        # We can do this symbolicaly with sympy just for testing...
        A_sym = integrate(
            integrate(sin(theta), (theta, theta_min, theta_max)),
            (phi, phi_min, phi_max),
        )
        # translate this to a numeric function
        A_num = lambdify(
            (theta_min, theta_max, phi_min, phi_max), A_sym, modules=["numpy"]
        )
    else:
        # or manually solve the integral since it is very simple
        def A_num(theta_min, theta_max, phi_min, phi_max):
            return (phi_max - phi_min) * (-np.cos(theta_max) + np.cos(theta_min))

    delta_theta, delta_phi = map(grad2rad, (delta_lat, delta_lon))
    dth = delta_theta / 2.0
    dph = delta_phi / 2.0

    def A_patch(theta):
        # computes the area of a pixel on the unitsphere
        if np.abs(theta < dth / 100):  # (==0)
            # pixel centered at north pole only extends northwards
            # print("##################### north pole ##########")
            theta_min_v = 0.0
            theta_max_v = dth
        elif np.abs(theta > np.pi - dth / 100):  # ==pi)
            # pixel centered at south pole only extends northwards
            # print("##################### south pole ##########")
            theta_min_v = np.pi - dth
            theta_max_v = np.pi
        else:
            # normal pixel extends south and north-wards
            theta_min_v = theta - dth
            theta_max_v = theta + dth

        phi_min_v = -dph
        phi_max_v = +dph
        res = A_num(theta_min_v, theta_max_v, phi_min_v, phi_max_v)
        # print(res)
        return res

    def pixel_area_on_unit_sphere(lat):
        # computes the fraction of the area of the sphere covered by this pixel
        theta_grad = lat + 90
        theta = grad2rad(theta_grad)
        # the area of the unitsphere is 4 * pi
        return A_patch(theta)

    return pixel_area_on_unit_sphere


def download_TRENDY_output(
    username: str,
    password: str,
    dataPath: Path,
    models: List[str],
    variables: List[str],
    experiments=["S2"],  # We are using s2 data, can add more
):
    def unzip_shutil(source_filepath, dest_filepath, model):
        if model == "YIBs":
            f = tarfile.open(source_filepath, "r:gz")
            f.extractall(path=dataPath)
            f.close()
        else:
            with gzip.open(source_filepath, "rb") as s_file, open(
                dest_filepath, "wb"
            ) as d_file:
                shutil.copyfileobj(s_file, d_file)

    # open a transport
    host = "trendy.ex.ac.uk"
    port = 22
    transport = paramiko.Transport(host)

    # authentication
    transport.connect(None, username=username, password=password)
    sftp = paramiko.SFTPClient.from_transport(transport)

    # We are using s2 data, can add more
    experiments = experiments

    # Loop through models, experiments, and variables to download
    for model in models:
        print("downloading data for", model, "model")
        for experiment in experiments:
            for variable in variables:
                modelname = model
                modelname_file = model
                ext = "nc"
                extra = ""

                if model == "CLM5":
                    modelname_file = "CLM5.0"
                elif model == "ORCHIDEEv3" or model == "ORCHIDEEv3_0.5deg":
                    modelname_file = "ORCHIDEEv3"
                elif model == "ISBA_CTRIP":
                    modelname_file = "ISBA-CTRIP"
                elif model == "JULES-ES" or model == "JULES-ES-1.0":
                    modelname = "JULES-ES-1.0"
                    modelname_file = "JULES-ES-1p0"
                elif model == "SDGVM" or model == "VISIT":
                    ext = "nc.gz"
                elif model == "YIBs":
                    ext = "nc.tar.gz"
                    if (
                        variable == "cSoil"
                        or variable == "cVeg"
                        or variable == "landCoverFrac"
                    ):
                        extra = "Annual_"
                    else:
                        extra = "Monthly_"
                elif model == "LPJwsl":
                    modelname_file = "LPJ"
                    ext = "nc.gz"

                filename = (
                    modelname_file
                    + "_"
                    + experiment
                    + "_"
                    + extra
                    + variable
                    + "."
                    + ext
                )

                try:
                    dataPath.mkdir(exist_ok=True)
                    complete_path = (
                        "output/" + modelname + "/" + experiment + "/" + filename
                    )
                    zipped_path = dataPath.joinpath(filename)
                    unzipped_filename = (
                        modelname_file
                        + "_"
                        + experiment
                        + "_"
                        + extra
                        + variable
                        + ".nc"
                    )
                    unzipped_path = dataPath.joinpath(unzipped_filename)
                    try:
                        unzipped_path.resolve(strict=True)
                    except FileNotFoundError:
                        try:
                            zipped_path.resolve(strict=True)
                        except FileNotFoundError:
                            print("downloading missing data:", variable)
                            sftp.get(remotepath=complete_path, localpath=zipped_path)
                            if zipped_path != unzipped_path:
                                print("unzipping", zipped_path)
                                unzip_shutil(zipped_path, unzipped_path, model)
                        else:
                            print("unzipping", zipped_path)
                            unzip_shutil(zipped_path, unzipped_path, model)
                    else:
                        print(unzipped_path, "exists, skipping")
                except FileNotFoundError as e:
                    print(e)
                    print(complete_path)
                    print(zipped_path)
    print("finished!")


# def make_feng_cost_func_2(svs):  #: Observables
#    # now we compute a scaling factor per observable stream
#    # fixme mm 10-28-2021
#    # The denominators in this case are actually the TEMPORAL variances of the data streams
#    # why would one want to punish something for
#    # changing????
#    obs_arr = np.stack([arr for arr in svs], axis=1)
#    means = obs_arr.mean(axis=0)
#    mean_centered_obs = obs_arr - means
#    denominators = np.sum(mean_centered_obs**2, axis=0)
#
#    def feng_cost_func_2(simu):  #: Observables
#        def f(i):
#            arr = simu[i]
#            obs = obs_arr[:, i]
#            diff = ((arr - obs) ** 2).sum() / denominators[i] * 100
#            return diff
#
#        return np.array([f(i) for i in range(len(simu))]).mean()
#
#    return feng_cost_func_2


# def make_param_filter_func(
#    c_max, c_min, betas: List[str] = []
# ) -> Callable[[np.ndarray], bool]:
#    positions = [c_max.__class__._fields.index(beta) for beta in betas]
#
#    def isQualified(c):
#        cond1 = (c >= c_min).all()
#        cond2 = (c <= c_max).all()
#
#        cond3 = np.sum([c[p] for p in positions]) <= 1
#        return cond1 and cond2 and cond3
#
#    return isQualified


# def make_StartVectorTrace(mvs):
#    # deprecated
#    svt = mvs.get_StateVariableTuple()
#    return namedtuple(
#        "StartVectorTrace",
#        [str(v) for v in svt]
#        + [str(v) + "_p" for v in svt]
#        + [str(v) + "_c" for v in svt]
#        + [str(v) + "_RT" for v in svt],
#    )
#
#
# def make_InitialStartVectorTrace(X_0, mvs, par_dict, func_dict):
#    # test the new iterator
#
#    # we make the X_p and X_c  parts compatible with the ones computed by the iterator
#    # for following timesteps (to keep it)
#    # As you can see in the definition of the iterator these values have no impact on further results
#    B_func, u_func = make_B_u_funcs_2(mvs, par_dict, func_dict)
#    I = u_func(0, X_0)
#    u = I.sum()
#    b = I / u
#    B = B_func(0, X_0)
#    B_inf = np.linalg.inv(B)
#    X_p_0 = B_inf @ I
#    X_c_0 = X_0 + X_p_0
#    RT_0 = B_inf @ b
#    # combine the three
#    # here we rely on order to be consistent
#    # (although we could use also the names of the namedtuple)
#    V_arr = np.concatenate((X_0, X_p_0, X_c_0, RT_0), axis=0)
#    StartVectorTrace = make_StartVectorTrace(mvs)
#    V_init = StartVectorTrace(*V_arr)
#    return V_init


def minimal_iterator(
    X_0,
    func_dict,
    mvs,  #: CMTVS,
    dvs,  #: Drivers,
    cpa,  #: Constants,
    epa,  #: EstimatedParameters
    t_0=0,
    delta_t_val: int = 1,  # defaults to 1day timestep
    # traced_expressions: Dict[str, Expr] = dict(),
    # extra_functions: Dict[str, Callable] = dict(),
):
    # apa = {**cpa._asdict(), **epa._asdict()}
    par_dict = make_param_dict(mvs, cpa, epa)
    return minimal_iterator_internal(
        mvs,
        X_0,
        par_dict,
        func_dict=func_dict,
        t_0=0,
        delta_t_val=delta_t_val,
    )


# fixme mm: 3-10-2023
# could have even more succint arguments like mvs
def minimal_iterator_internal(
    mvs,
    X_0,
    par_dict,  #: EstimatedParameters
    func_dict,
    t_0=0,
    delta_t_val=1,  # defaults to 1day timestep
    # traced_expressions: Dict[str, Expr] = dict(),
    # extra_functions: Dict[str, Callable] = dict(),
):
    t = mvs.get_TimeSymbol()
    state_vector = mvs.get_StateVariableTuple()
    delta_t = Symbol("delta_t")
    disc_par_dict = {**par_dict, delta_t: delta_t_val}
    B_func = hr.numerical_array_func(
        state_vector=state_vector,
        time_symbol=t,
        expr=mvs.get_CompartmentalMatrix(),
        parameter_dict=par_dict,
        func_dict=func_dict,
    )
    I_func = hr.numerical_1d_vector_func(
        state_vector=state_vector,
        time_symbol=t,
        expr=mvs.get_InputTuple(),
        parameter_dict=par_dict,
        func_dict=func_dict,
    )
    # make sure that the startvector X_0 is a ONE dimensional  vector
    # since it is important that  X_0 and the result of I(t,X) have
    # the same dimension
    X_0 = X_0.reshape(-1)
    assert I_func(0, X_0).ndim == 1

    bit = BlockDictIterator(
        iteration_str="it",  # access the inner counter
        start_seed_dict=ArrayDict({"X": X_0, "t": t_0}),
        present_step_funcs=OrderedDict(
            {
                # these are functions that are  applied in order
                # on the start_seed_dict
                # they might compute variables that are purely
                # diagnostic or those that are necessary for the
                # next step
                #
                # The first 2 are essential for any compartmental system
                "B": lambda t, X: -B_func(t, X),
                "I": lambda t, X: I_func(t, X),
                #
                # **extra_functions
            }
        ),
        next_step_funcs=OrderedDict(
            {
                # these functions have to compute the seed for the next timestep
                "X": lambda X, B, I: X + (I - B @ X) * delta_t_val,
                "t": lambda t: t + delta_t_val,
            }
        ),
    )
    return bit


# fixme mm 12-2 2022
# should be better called tracebility Iterator result
def traceability_iterator(
    X_0,
    func_dict,
    mvs,  #: CMTVS,
    dvs,  #: Drivers,
    cpa,  #: Constants,
    epa,  #: EstimatedParameters
    t_0,
    delta_t_val,
):
    print(
        """
    #######################################################
    Deprecation Warning:
    rather use traceability_iterator_internal directly
    for which is called by this  wrapper.
    """
    )
    par_dict = make_param_dict(mvs, cpa, epa)
    mit = traceability_iterator_internal(
        mvs,
        X_0,
        par_dict,
        func_dict,
        delta_t_val,
        t_0=0,
    )
    return ArrayDictResult(mit)


def traceability_iterator_internal(
    mvs: CMTVS,
    X_0,
    par_dict,
    func_dict,
    delta_t_val: int = 1,  # defaults to 1day timestep
    t_0=0,
):
    mit = minimal_iterator_internal(
        mvs=mvs,  #: CMTVS,
        X_0=X_0,
        par_dict=par_dict,
        func_dict=func_dict,
        t_0=t_0,
        delta_t_val=delta_t_val,
    )

    (
        veg_x_func,
        soil_x_func,
        out_2_veg_func,
        veg_2_soil_func,
        veg_2_out_func,
        soil_2_out_func,
    ) = map(
        lambda expr: hr.numerical_func_of_t_and_Xvec(
            state_vector=mvs.get_StateVariableTuple(),
            time_symbol=mvs.get_TimeSymbol(),
            expr=expr,
            parameter_dict=par_dict,
            func_dict=func_dict,
        ),
        [
            mvs.get_AggregatedVegetationCarbon(),
            mvs.get_AggregatedSoilCarbon(),
            mvs.get_AggregatedVegetationCarbonInFlux(),
            mvs.get_AggregatedVegetation2SoilCarbonFlux(),
            mvs.get_AggregatedVegetationCarbonOutFlux(),
            mvs.get_AggregatedSoilCarbonOutFlux(),
        ],
    )
    present_step_funcs = OrderedDict(
        {
            "u": lambda I: I.sum(),
            "b": lambda I, u: I / u,
            "B_inv": lambda B: np.linalg.inv(B),
            "X_c": lambda B_inv, I: B_inv @ I,
            "X_p": lambda X_c, X: X_c - X,
            "X_dot": lambda I, B, X: I - B @ X,
            "system_RT": lambda X_c, u: X_c / u,  # =B_inv@b but cheeper to compute
            "system_RT_sum": lambda system_RT: system_RT.sum(),
            "x": lambda X: X.sum(),
            "x_dot": lambda X_dot: X_dot.sum(),
            "system_m": lambda B, X, x: (B @ X).sum()
            / x,  # rate of the surrogate system
            "system_tot": lambda system_m: 1 / system_m,
            "x_c": lambda system_m, u: 1 / system_m * u,  # x_c of the surrogate system
            "x_p": lambda x_c, x: x_c - x,
            "rt": lambda x_c, u: x_c / u,
            "veg_x": lambda t, X: veg_x_func(t, X),
            "soil_x": lambda t, X: soil_x_func(t, X),
            "out_2_veg": lambda t, X: out_2_veg_func(t, X),
            "veg_2_soil": lambda t, X: veg_2_soil_func(t, X),
            "veg_2_out": lambda t, X: veg_2_out_func(t, X),
            "soil_2_out": lambda t, X: soil_2_out_func(t, X),
            "veg_m": lambda veg_2_out, veg_2_soil, veg_x: (veg_2_out + veg_2_soil)
            / veg_x,  # veg rate
            "veg_tot": lambda veg_m: 1 / veg_m,
            "soil_m": lambda soil_2_out, soil_x: soil_2_out / soil_x,  # soil rate
            "soil_tot": lambda soil_m: 1 / soil_m,
            # **extra_functions
        }
    )
    mit.add_present_step_funcs(present_step_funcs)
    return mit


def all_timelines_starting_at_steady_state(
    # X_0, is computed at time tmin
    mvs,  #: CMTVS,
    func_dict,
    par_dict,
    t_min,  # time of equilibrium computation in days after the start of the simulation
    index_slice: slice,  # making up the iterator slice to be output start_ind =0 refers to the
    delta_t_val: int = 1,  # iterator timestep defaults to 1day timestep
):
    # We want to compute a steady state from driver values in the
    # spring and start our computation from there

    a_dens_function, X_fix = start_age_distributions_from_steady_state(
        t0=t_min,  # determines B_0=B(t_0) and u_0=B(t_0) from which the equilibrium is
        # computed  (via the functions in func_dict which are called at the appropriate
        # time)
        srm=mvs.get_SmoothReservoirModel(),
        parameter_dict=par_dict,
        func_set=func_dict,
    )
    bit = ArrayDictResult(
        traceability_iterator_internal(
            mvs, X_fix, par_dict, func_dict, delta_t_val=delta_t_val, t_0=t_min
        )
    )
    vals = bit[index_slice]

    times = vals.t
    mvs = mvs.update(
        {
            NumericParameterization(par_dict=par_dict, func_dict=func_dict),
            NumericStartValueArray(X_fix),
            NumericSimulationTimes(times),
        }
    )
    sv = mvs.get_StateVariableTuple()
    n_pools = len(sv)
    smr = mvs.get_SmoothModelRun()
    # now create two submodels
    start_mean_age_vec = start_age_moments_from_steady_state(
        smr.model,
        t0=t_min,
        parameter_dict=par_dict,
        func_set=func_dict,
        max_order=1,
    ).reshape(-1)
    # compute solutions for the mean age system starting at X_fix and star
    order = 1
    s_arr, s_func = smr._solve_age_moment_system(
        order,
        start_mean_age_vec.reshape(1, n_pools),
    )
    # the first n colums are the solution
    # solutions = smr.solve()
    solutions = s_arr
    t = mvs.get_TimeSymbol()
    mean_btts = smr.backward_transit_time_moment(
        order=1, start_age_moments=start_mean_age_vec.reshape(1, n_pools)
    )

    def sub_mr_smav(mvs, svt):
        sv = mvs.get_StateVariableTuple()
        svl = list(sv)
        svs = set(sv)
        combined = (
            set(sv),
            mvs.get_InFluxesBySymbol(),
            mvs.get_OutFluxesBySymbol(),
            mvs.get_InternalFluxesBySymbol(),
        )

        (
            sub_sv_set,
            sub_in_fluxes,
            sub_out_fluxes,
            sub_internal_fluxes,
        ) = hr.extract(combined, set(svt))

        # we have to provide solutions (functions of time)  for the
        # state variables of the outer system ( In our case soil
        # influxes depend on the poolvalues of the veg system. In
        # general every subsystem could depend on pools of every other
        # subsystem
        outer_pools = svs.difference(sub_sv_set)
        # first replace the symbols by functions
        subs_dict = {sym: Function(str(sym))(t) for sym in outer_pools}
        sub_in_fluxes_f, sub_internal_fluxes_f, sub_out_fluxes_f = map(
            lambda d: {k: sympify(flux_ex).subs(subs_dict) for k, flux_ex in d.items()},
            (sub_in_fluxes, sub_internal_fluxes, sub_out_fluxes),
        )

        # now prepare the extension of the func_dict by the solutions
        # for the variables
        def s_func_maker(sym):
            return lambda t: s_func(t)[svl.index(sym)]

        outer_pool_funcs = {
            # note that we create the function in another
            # function since the normal dictionary comprehension
            # fails due to pythons lazy evaluation not
            # protecting the function specific variable
            # (in this case sym)
            Function(str(sym)): s_func_maker(sym)
            for sym in outer_pools
        }
        sub_func_dict = {**func_dict, **outer_pool_funcs}
        sub_X_fix = X_fix[[svl.index(w) for w in svt]]
        sub_mvs = CMTVS(
            {
                t,
                StateVariableTuple(svt),
                InFluxesBySymbol(sub_in_fluxes_f),
                InternalFluxesBySymbol(sub_internal_fluxes_f),
                OutFluxesBySymbol(sub_out_fluxes_f),
                NumericParameterization(
                    par_dict=par_dict,
                    func_dict=sub_func_dict,
                ),
                NumericStartValueArray(sub_X_fix),
                NumericSimulationTimes(times),
            },
            mvs.computers,
        )
        sub_smr = sub_mvs.get_SmoothModelRun()
        # for the start mean ages we can not just take the respective parts of the system mean ages
        # since they refer to the times stuff has spent in the SYSTEM and not in the SUB system.
        # Although the pools in the subsystem are pools in the
        sub_start_mean_age_vec = start_age_moments_from_steady_state(
            sub_smr.model,
            t0=t_min,
            parameter_dict=par_dict,
            func_set=sub_func_dict,
            max_order=1,
        )
        return sub_mvs, sub_smr, sub_start_mean_age_vec

    veg_sv = mvs.get_VegetationCarbonStateVariableTuple()
    veg_mvs, veg_smr, veg_smav = sub_mr_smav(mvs, veg_sv)
    veg_btt = veg_smr.backward_transit_time_moment(order=1, start_age_moments=veg_smav)
    veg_arr, veg_func = veg_smr._solve_age_moment_system(
        order, start_age_moments=veg_smav
    )
    soil_sv = mvs.get_SoilCarbonStateVariableTuple()
    soil_mvs, soil_smr, soil_smav = sub_mr_smav(mvs, soil_sv)
    # soil_smr.initialize_state_transition_operator_cache(lru_maxsize=None)
    soil_arr, soil_func = soil_smr._solve_age_moment_system(
        order, start_age_moments=soil_smav
    )
    soil_btt = soil_smr.backward_transit_time_moment(
        order=1, start_age_moments=soil_smav
    )
    vnp = veg_smr.nr_pools
    snp = soil_smr.nr_pools
    vals.update(
        {
            "system_continuous_solution": s_arr[:, 0:n_pools],
            "system_continuous_mean_age": s_arr[:, n_pools : 2 * n_pools],
            "system_continuous_mean_btt": mean_btts,
            "veg_continuous_solution": veg_arr[:, 0:vnp],
            "veg_continuous_mean_age": veg_arr[:, vnp : 2 * vnp],
            "veg_continuous_mean_btt": veg_smr.backward_transit_time_moment(
                order=1, start_age_moments=veg_smav
            ),
            "soil_continuous_solution": soil_arr[:, 0:snp],
            "soil_continuous_mean_age": soil_arr[:, snp : 2 * snp],
            "soil_continuous_mean_btt": soil_smr.backward_transit_time_moment(
                order=1, start_age_moments=soil_smav
            ),
            #'continuous_times': np.array(times)
        }
    )
    return vals


# deprecated
# def write_global_mean_cache(gm_path, gm: np.array, var_name: str):
#    # fixme deprecated
#    # just a safety wrapper until every reference to this function
#    # is replaced
#    write_timeline_to_nc_file(gm_path, gm, var_name)
#
## deprecated
def write_timeline_to_nc_file(gm_path, gm: np.array, var_name: str):
    # var=ds.variables[var_name]
    if gm_path.exists():
        print("removing old cache file{}")
        os.remove(gm_path)

    n_t = gm.shape[0]
    time_dim_name = "time"
    ds_gm = nc.Dataset(str(gm_path), "w", persist=True)
    time = ds_gm.createDimension(time_dim_name, size=n_t)
    var_gm = ds_gm.createVariable(var_name, np.float64, [time_dim_name])
    gm_ma = np.ma.array(gm, mask=np.zeros(gm.shape, dtype=np.bool_))
    var_gm[:] = gm_ma
    ds_gm.close()


# deprecated
# just for savety
# def get_cached_global_mean(gm_path, vn):
#    return get_nc_array(gm_path, vn)


def get_nc_array(gm_path, vn):
    return nc.Dataset(str(gm_path)).variables[vn].__array__()


def read_var_dict(
    targetPath: Path,
    varnames: List[str],
    varname2filename: Callable,
):
    def get_cached_global_mean(vn):
        res = get_nc_array(targetPath.joinpath(varname2filename(vn)), vn)
        print(
            """ Found cached global mean files. If you want to recompute the global means
            remove the following files: """
        )
        print(targetPath.joinpath(varname2filename(vn)))
        return res

    return {vn: get_cached_global_mean(vn) for vn in varnames}


def write_var_dict(
    arr_dict: Dict,
    targetPath: Path,
    varname2filename: Callable,
):
    if not targetPath.exists():
        targetPath.mkdir(parents=True)
    for vn, val in arr_dict.items():
        write_timeline_to_nc_file(targetPath.joinpath(varname2filename(vn)), val, vn)


# fixme: possibly obsolete - see combined_masks_2 using nearest neighbor resampling
def combine_masks(coord_masks: List["CoordMask"]):
    def k(cm):
        arr = cm.index_mask
        return arr.shape[0] + arr.shape[1]

    return reduce(
        project_2,  # (acc,el)
        sorted(coord_masks, key=k),  # we want the finest grid as the last
    )


# def combine_masks(masks,i2cs):
#    m1,m2 = masks
#    i2c_1,i2c_2 = i2cs
#
#    def g(m,i2c):
#        s=m.shape
#        print(s[0])
#        bounds=[]
#        lat_0,lon_0=i2c(0,0)
#        lat_1,lon_1=i2c(1,1)
#        step_lat=lat_1-lat_0
#        step_lon=lon_1-lon_0
#        for i in range(s[0]):
#            for j in range(s[1]):
#                print(i,j)
#                print(m[i,j])
#                if  m[i,j]:
#                    lat,lon=i2c(i,j)
#                    res=boundaries(
#                        min_lat=lat-step_lat/2,
#                        max_lat=lat+step_lat/2,
#                        min_lon=lon-step_lon/2,
#                        max_lon=lon+step_lon/2
#                    )
#                    bounds.append(res)
#        return bounds
#
#    #return g(m1,i2c_1)
#    return reduce(
#        lambda acc,el:acc+el,
#        [g(m,i2c) for m,i2c in zip(masks,i2cs)]
#    )
#
#
# def open_interval_intersect(i1,i2):
#    min1,max1=i1
#    min2,max2=i2
#    mid2=(min2+max2)/2
#    return (
#        min1 < min2 and min2< max1
#        or
#        min1 < max2 and max2< max1
#        or
#        min1 < mid2 and mid2< max1
#    )
#
# def pixel_intersect(b1,b2):
#    return (
#        open_interval_intersect(
#            (b1.min_lat,b1.max_lat),
#            (b2.min_lat,b2.max_lat)
#        )
#        and
#        open_interval_intersect(
#            (b1.min_lon,b1.max_lon),
#            (b2.min_lon,b2.max_lon)
#        )
#    )
#
# def project(s,i2c,common_mask):
#    mask=np.zeros(s)
#
#    lat_0,lon_0=i2c(0,0)
#    lat_1,lon_1=i2c(1,1)
#    step_lat=lat_1-lat_0
#    step_lon=lon_1-lon_0
#
#    def in_cm(i,j):
#        lat,lon=i2c(i,j)
#        p_b=boundaries(
#            min_lat=lat-step_lat/2,
#            max_lat=lat+step_lat/2,
#            min_lon=lon-step_lon/2,
#            max_lon=lon+step_lon/2
#        )
#        return reduce(
#            lambda acc,mpb: acc or pixel_intersect(p_b,mpb),
#            common_mask,
#            False
#        )
#
#    for i in range(s[0]):
#        for j in range(s[1]):
#            mask[i,j] = True if in_cm(i,j) else False
#
#    return mask


def permutation(unordered_vec):
    n = len(unordered_vec)
    tups = [(i, unordered_vec[i]) for i in range(n)]
    ordered_tups = sorted(tups, key=lambda t: t[1])

    data = np.ones(n)
    col = np.array([tup[0] for tup in ordered_tups])
    row = np.arange(n)
    p = sparse.bsr_matrix((data, (row, col)), dtype=np.int64)
    ordered_vec = np.array([tup[1] for tup in ordered_tups])
    t = p @ row
    assert (t == col).all()
    return ordered_vec, p, p.transpose()


class CoordMask:
    """A CoordMask
    with -90 <= lat <=90
    with -180 <= lon <=  180
    lat=0,Equator
    lon=0,Greenich
    """

    def __init__(self, index_mask: np.ndarray, tr: SymTransformers):
        s = index_mask.shape

        self.index_mask = index_mask.astype(np.float64, casting="safe")
        self.tr = tr

    # def plot_cell(
    #        ax,
    #        b: boundaries,
    #        color: str
    #    ):
    #    xs=[
    #        b.min_lon,
    #        b.min_lon,
    #        b.max_lon,
    #        b.max_lon,
    #        b.min_lon
    #    ]
    #    ys=[
    #        b.min_lat,
    #        b.max_lat,
    #        b.max_lat,
    #        b.min_lat,
    #        b.min_lat
    #    ]
    #    ax.plot(xs,ys,color=color)
    #    ax.fill(xs,ys,color=color,alpha=0.4)

    # def plot(self,ax,color="black"):
    #    mask = self.index_mask
    #
    #    s = mask.shape
    #    tr = self.tr

    #    for i in range(s[0]):
    #        for j in range(s[1]):
    #            min_lat, max_lat = tr.i2lat_min_max(i)
    #            min_lon, max_lon = tr.i2lon_min_max(j)
    #            self.__class__.plot_cell(
    #                ax,
    #                boundaries(
    #                    min_lat=min_lat,
    #                    max_lat=max_lat,
    #                    min_lon=min_lon,
    #                    max_lon=max_lon,
    #                ),
    #                color="red" if mask[i,j] == 1 else color
    #            )
    @property
    def n_lats(self):
        return self.index_mask.shape[0]

    @property
    def n_lons(self):
        return self.index_mask.shape[1]

    @property
    def lats(self):
        return np.array([self.tr.i2lat(i) for i in range(self.n_lats)])

    @property
    def lons(self):
        return np.array([self.tr.i2lon(i) for i in range(self.n_lons)])

    def ordered_lats(self):
        return permutation(self.lats)

    def ordered_lons(self):
        return permutation(self.lons)

    def ordered_mask(self):
        ordered_lats, p_lats, p_lats_inv = self.ordered_lats()
        ordered_lons, p_lons, p_lons_inv = self.ordered_lons()
        mask = self.index_mask
        new_mask = p_lats @ mask @ p_lons
        return ordered_lats, ordered_lons, new_mask

    def plot_dots(self, ax, color="black"):
        ordered_lats, ordered_lons, new_mask = self.ordered_mask()
        new_mask
        X, Y = np.meshgrid(ordered_lons, ordered_lats)
        ax.pcolormesh(X, Y, new_mask)

    def write_netCDF4(self, path):
        # for visualization we write out the lat and lon arrays too
        # although we don't use them at the moment to reconstruct
        # the transformers.
        index_mask = self.index_mask
        s = index_mask.shape
        n_lats, n_lons = s
        ds = nc.Dataset(str(path), "w", persist=True)
        # mask = ds.createDimension('mask',size=s)
        lat = ds.createDimension("lat", size=n_lats)
        lon = ds.createDimension("lon", size=n_lons)

        test = ds.createVariable("mask", "i1", ["lat", "lon"])
        test[:, :] = index_mask

        lats = ds.createVariable("lat", "float64", ["lat"])
        lats[:] = list(map(self.tr.i2lat, range(n_lats)))
        lons = ds.createVariable("lon", "float64", ["lon"])
        lons[:] = list(map(self.tr.i2lon, range(n_lons)))
        ds.close()


def project_2(source: CoordMask, target: CoordMask):
    s_mask = source.index_mask  # .astype(np.float64,casting='safe')
    # print(s_mask.mean())
    s_n_lat, s_n_lon = s_mask.shape

    t_mask = target.index_mask
    t_n_lat, t_n_lon = t_mask.shape

    t_tr = target.tr
    s_tr = source.tr
    # we first create the (unified) lat lon coordinates of the source array
    # and interpolate the source mask on them
    lats = np.array(list(map(s_tr.i2lat, range(s_n_lat))))
    lons = np.array(list(map(s_tr.i2lon, range(s_n_lon))))
    kind = "cubic"
    # f_old = interpolate.interp2d(x=lons, y=lats, z=s_mask,kind=kind) # deprecated by scipy
    # now we apply this interpolating function to the target grid
    # points
    target_lats = np.array(list(map(t_tr.i2lat, range(t_n_lat))))
    target_lons = np.array(list(map(t_tr.i2lon, range(t_n_lon))))
    # order lats and lons since f returns an array that assumes this anyway
    otlats, p_lat, p_lat_inv = target.ordered_lats()
    otlons, p_lon, p_lon_inv = target.ordered_lons()
    # val_old = f_old(otlons, otlats)
    # float_grid_old = p_lat_inv @ f_old(otlons, otlats) @ p_lon_inv
    from scipy.interpolate import RegularGridInterpolator as RGI

    f = RGI((lons, lats), s_mask.T, method=kind, bounds_error=False)
    ootlons, ootlats = np.meshgrid(otlons, otlats, indexing="ij", sparse=True)
    val_new = f((ootlons, ootlats)).T
    float_grid = p_lat_inv @ val_new @ p_lon_inv
    #    try:
    #        assert(np.allclose(float_grid_old, float_grid))
    #    except AssertionError:
    #        print(float_grid_old, float_grid)
    #
    # print(float_grid.mean())
    projected_mask = float_grid > 0.5
    return CoordMask(index_mask=np.logical_or(projected_mask, t_mask), tr=t_tr)


def resample_grid(
    source_coord_mask,
    target_coord_mask,
    var,
    method="nearest",
    radius_of_influence=500000,
    neighbours=10,
):
    lon2d, lat2d = np.meshgrid(source_coord_mask.lons, source_coord_mask.lats)
    lon2d_t, lat2d_t = np.meshgrid(target_coord_mask.lons, target_coord_mask.lats)

    lats_source = source_coord_mask.lats
    lons_source = source_coord_mask.lons

    lats = target_coord_mask.lats
    lons = target_coord_mask.lons
    orig_def = pyresample.geometry.SwathDefinition(lons=lon2d, lats=lat2d)

    targ_def = pyresample.geometry.SwathDefinition(lons=lon2d_t, lats=lat2d_t)

    if method == "nearest":
        target_var = pyresample.kd_tree.resample_nearest(
            orig_def,
            var,
            targ_def,
            radius_of_influence=radius_of_influence,
            fill_value=None,
        )

    elif method == "idw":
        wf = lambda r: 1 / r**2
        target_var = pyresample.kd_tree.resample_custom(
            orig_def,
            var,
            targ_def,
            radius_of_influence=radius_of_influence,
            neighbours=neighbours,
            weight_funcs=wf,
            fill_value=None,
        )
    elif method == "gauss":
        target_var = pyresample.kd_tree.resample_gauss(
            orig_def,
            var,
            targ_def,
            radius_of_influence=radius_of_influence,
            neighbours=neighbours,
            sigmas=250000,
            fill_value=None,
        )
    else:
        raise Exception(
            "Invalid resample method. Valid options are: 'nearest', 'idw' and 'gauss'"
        )

    return CoordMask(index_mask=target_var, tr=target_coord_mask.tr)


def global_coord_mask_resampled(template: np.ndarray, ctr, itr):
    # common wrapper to be used by all the models
    # to guarantee the same projection procedure
    gm = common_global_mask_expanded
    gcm = resample_grid(
        source_coord_mask=gm,
        target_coord_mask=CoordMask(
            index_mask=np.zeros_like(template),
            tr=SymTransformers(
                ctr=ctr,
                itr=itr,
            ),
        ),
        var=gm.index_mask,
        method="nearest",
    )
    return gcm


def resample_nc(
    model_names,  # dictionary e.g. "ab_classic":"CLASSIC"
    experiment_names,  # e.g. ['S2', 'S3']
    target_mask,
    method="nearest",
    radius_of_influence=500000,
):
    for experiment in experiment_names:
        print(
            "\033[1m" + ". . . Resampling data for " + experiment + " experiment . . ."
        )
        model_folders = [(m) for m in model_names]
        m_names = list(model_names.values())
        g_mask = target_mask.index_mask
        k = 0  # model counter
        for mf in model_folders:
            print("\033[1m" + m_names[k])
            print("\033[0m")
            experiment_name = m_names[k] + "_" + experiment + "_"
            conf_dict = confDict(mf)
            dataPath = Path(conf_dict["dataPath"])
            model_mask = msh(mf).spatial_mask(dataPath=Path(conf_dict["dataPath"]))
            for vn in msh(mf).data_str._fields:
                print("Resampling " + vn)
                file_path = dataPath.joinpath(
                    msh(mf).nc_file_name(vn, experiment_name=experiment_name)
                )
                ds = nc.Dataset(str(file_path))
                var = ds.variables[vn][:, :, :].data
                # preparing a narray to store results
                zero_array = np.zeros((var.shape[0], g_mask.shape[0], g_mask.shape[1]))
                gm = zero_array.copy()
                for i in range(gm.shape[0]):
                    gm[i, :, :] = g_mask
                var_array = zero_array
                # procesing all time steps
                for i in range(var.shape[0]):
                    var_current = var[i, :, :]
                    # initial masking
                    var_masked = np.ma.array(var_current, mask=model_mask.index_mask)
                    # resampling
                    var_resampled = resample_grid(
                        source_coord_mask=model_mask,
                        target_coord_mask=target_mask,
                        var=var_masked,
                        method=method,
                        radius_of_influence=radius_of_influence,
                    )
                    var_array[i, :, :] = var_resampled.index_mask
                    if i // 100 == i / 100:
                        print(
                            str(i + 1)
                            + " out of "
                            + str(var.shape[0])
                            + " time steps completed"
                        )
                # final masking
                var_final = np.ma.array(var_array, mask=gm)
                # creating and writing a new NetCDF file
                s = g_mask.shape
                n_lats, n_lons = s
                new_path = dataPath.joinpath(dataPath, experiment_name + vn + "_res.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)
                source_times = ds.variables["time"][:].data
                time = ds_new.createDimension("time", size=len(source_times))
                # creating variables
                nc_var = ds_new.createVariable(vn, "float32", ["time", "lat", "lon"])
                nc_var[:, :, :] = var_final
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(target_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(target_mask.tr.i2lon, range(n_lons)))
                times = ds_new.createVariable("time", "float32", ["time"])
                times[:] = source_times
                # closing NetCDF files
                ds.close()
                ds_new.close()
            k += 1  # model counter
        print("Done!")


def combine_masks_2(coord_masks: List["CoordMask"]):
    def combine(source, target):
        resampled_mask = resample_grid(
            source_coord_mask=source,
            target_coord_mask=target,
            var=source.index_mask,
            method="nearest",
            radius_of_influence=500000,
            neighbours=10,
        )
        combined_mask = CoordMask(
            index_mask=np.logical_or(resampled_mask.index_mask, target.index_mask),
            tr=target.tr,
        )
        return combined_mask

    target_mask = coord_masks[-1]  # we use last mask in the list as template grid
    for i in range(len(coord_masks) - 1):
        target_mask = combine(coord_masks[i], target_mask)
    return target_mask


# outputs a table with flow diagrams, compartmental matrices and allocation vectors
def model_table(
    model_names,  # dictionary (folder name : model name)
):
    model_folders = [(k) for k in model_names]
    mf = model_folders[0]
    import_module("{}.source".format(mf))

    tups = [(model_names[mf], mvs(mf)) for mf in model_folders]

    return dh.table(
        tups=tups, types_compact=[], types_expanded=[InputTuple, CompartmentalMatrix]
    )


def transform_maker(lat_0, lon_0, step_lat, step_lon):
    # fixme mm: 8-27-2022
    # since we do not need the inverse transformation anymore
    # we could use the lat lon variables of the users data sets
    # if we cache them it might be fast enough.
    """
    The result of this function is a tuple of functions
    To translate from indices (i,j)  in the array  to (lat,lon) with the following
    properties
    values of the center of the indexed pixel
    and back
    These functions do not necessarily work for every model
    You have to check them on the latitude longitude variables of your netcdf files.
    """
    n_lat = 180.0 / abs(step_lat)
    if int(n_lat) != n_lat:
        raise Exception("step_lat has to be a divisor of 180")
    n_lon = 360.0 / abs(step_lon)
    if int(n_lon) != n_lon:
        raise Exception("step_lon has to be a divisor of 360")

    def i2lat(i_lat):
        if i_lat > (n_lat - 1):
            raise IndexError(
                "i_lat > n_lat-1; with i_lat={}, n_lat={}".format(i_lat, n_lat)
            )
        return lat_0 + (step_lat * i_lat)

    #
    #    def i2lat_min_max(i):
    #        #compute the lat boundaries of pixel i
    #        center=i2lat(i)
    #        lat_min = center - step_lat/2
    #        lat_max=  center + step_lat/2
    #        return lat_min,lat_max
    #
    #    # the inverse finds the indices of the pixel containing
    #    # the point with the given coordinates
    #    def lat2i(lat):
    #        ir=(lat-lat_0)/step_lat
    #        ii=int(ir)
    #        d=ir-ii
    #        if ii == (n_lat -1):
    #            print("ii",ii)
    #
    #        return ii if (d<0.5 or ii == (n_lat - 1))  else ii+1
    #
    def i2lon(i_lon):
        if i_lon > (n_lon - 1):
            raise IndexError(
                "i_lon > n_lon; with i_lon={0}, n_lon={1}".format(i_lon, n_lon)
            )
        return lon_0 + (step_lon * i_lon)

    #
    #    def i2lon_min_max(i):
    #        #compute the lon boundaries of pixel i
    #        center=i2lon(i)
    #        lon_min = center - step_lon/2
    #        lon_max=  center + step_lon/2
    #        return lon_min,lon_max
    #
    #    def lon2i(lon):
    #        # we cant use round since we want ir=3.5 to be already in pixel 4
    #        ir=(lon-lon_0)/step_lon
    #        ii=int(ir)
    #        d=ir-ii
    #        return ii if d<0.5 else ii+1

    return Transformers(
        i2lat=i2lat,
        # i2lat_min_max=i2lat_min_max,
        # lat2i=lat2i,
        i2lon=i2lon,
        # i2lon_min_max=i2lon_min_max,
        # lon2i=lon2i,
    )


def make_cached_func(create_and_write: Callable, read: Callable):
    """
    basic cache functionality:
    param: create_and_write: The funciton that creates the  object to be cached.
    param: The first argument must be a path
    read: The function that extracts the cached object. It has only one argument <path>.
    """

    def read_or_create(*args, **kwargs):
        path = args[0]
        if path.exists():
            print(
                "Found cache path {}. If you want to recompute the result remove it.".format(
                    str(path)
                )
            )
            return read(path)
        else:
            pp = path.parent
            if not (pp.exists()):
                pp.mkdir(parents=True)
            return create_and_write(*args, **kwargs)

    return read_or_create


# +
# For harmonizing timelines and plotting
# -


def package_name():
    return __package__


def t_min_max_overlap_gm(
    model_folders: List[str],
    delta_t_val: float,  # iterator time step
    start_shift=0,  # time b etween the start.date and the time of the iterator's first step in
):
    tl = [
        h.date.timestep_dates_in_days_since_AD(
            msh(mf).start_dt(), msh(mf).n_months(), delta_t_val, start_shift
        )
        for mf in model_folders
    ]
    t_min = max([t.min() for t in tl])
    t_max = min([t.max() for t in tl])
    # from IPython import embed;embed()
    return (t_min, t_max)


def min_max_index_2(
    mf: str,  #
    delta_t_val: float,  # iterator time step
    t_min,  # output of t_min_tmax_overlap or t_min_tmax_full
    t_max,  # output of t_min_tmax_overlap or t_min_tmax_full
    start_shift=0,
):
    ts = h.date.timestep_dates_in_days_since_AD(
        msh(mf).start_dt(), msh(mf).n_months(), delta_t_val, start_shift
    )
    # ts = times_in_days_aD(test_args, delta_t_val, start_shift)

    def count(acc, i):
        min_i, max_i = acc
        t = ts[i]
        min_i = min_i + 1 if t < t_min else min_i
        max_i = max_i + 1 if t < t_max else max_i
        return (min_i, max_i)

    return reduce(count, range(len(ts)), (0, 0))


# fixme: do we need this partition?
# partitions is only used in averaged_1d_array and in Trace_Tuple class
# we can build this partition by a little function
def partitions(start, stop, nr_acc=1):
    diff = stop - start
    step = nr_acc
    number_of_steps = int(diff / step)
    last_start = start + number_of_steps * step
    last_tup_l = [(last_start, stop)] if last_start < stop else []
    return [
        (start + step * i, start + step * (i + 1)) for i in range(number_of_steps)
    ] + last_tup_l


# fixme: merge with avg_timeline?
# avg_timeline has a better interface - 2nd argument is a number rather than an array of parts
# so avg_timeline does not need a partitions function
# which is very useful in graphs when ve average times as well (partitions make it very inconvenient)
# averaged_1d_array and partiotions are used in Trace_Tuple class. Can they be changed to avg_timeline?
def averaged_1d_array(arr, partitions):
    """this function also works for multidimensional arrays
    It assumes that the first dimension is time/iteration
    """
    return np.array([arr[p[0] : p[1]].sum(axis=0) / (p[1] - p[0]) for p in partitions])


# The iterator allows us to compute and easily access the desired timelines with python index notation [2:5]
# it will return the values from position 2 to 5 of the solution (and desired variables).

# for attribution using variance decomposition (Sha Zhou et al., 2018)


# fixme: possible obsolete = gh.averaged_1d_array
def avg_timeline(timeline, averaging):  # array  # number of steps over which to average
    if averaging < 1:
        raise Exception("Invalid averaging in avg_timeline: should be >=1")
    output = timeline
    if averaging > 1:
        n = len(timeline) // averaging
        if len(timeline) % averaging > 0:
            n += 1
        output = np.zeros(n)
        counter = 0
        i = 0
        while i < (len(timeline)):
            x = 0
            sum = 0
            while x < (averaging):
                if i + x > (len(timeline) - 1):
                    break
                sum += timeline[i + x]
                x += 1
            output[counter] = sum / (x)
            counter += 1
            i += x
    return output


def cov_mat(a, b):
    return np.cov(a, b, bias=True)


def sum_attribution(y, x1, x2):
    cm1 = cov_mat(x1, y)
    var_y = cm1[1, 1]
    cov_x1_y = cm1[0, 1]

    cm2 = cov_mat(x2, y)
    cov_x2_y = cm2[0, 1]

    return tuple((c / var_y for c in (cov_x1_y, cov_x2_y)))


def product_attribution(v, z1, z2):
    y, x1, x2 = map(np.log, (v, z1, z2))
    return sum_attribution(y, x1, x2)


# fixme mm 4-7-2023
# deprecated
def write_data_streams_cache(gm_path, gm):
    # var=ds.variables[var_name]
    if gm_path.exists():
        print("removing old cache file{}")
        os.remove(gm_path)
    names = gm._fields

    n_t = gm.cVeg.shape[0]
    time_dim_name = "time"
    ds_gm = nc.Dataset(str(gm_path), "w", persist=True)
    time = ds_gm.createDimension(time_dim_name, size=n_t)
    var_gm0 = ds_gm.createVariable(names[0], np.float64, [time_dim_name])
    var_gm1 = ds_gm.createVariable(names[1], np.float64, [time_dim_name])
    var_gm2 = ds_gm.createVariable(names[2], np.float64, [time_dim_name])
    var_gm3 = ds_gm.createVariable(names[3], np.float64, [time_dim_name])
    var_gm4 = ds_gm.createVariable(names[4], np.float64, [time_dim_name])
    var_gm5 = ds_gm.createVariable(names[5], np.float64, [time_dim_name])
    gm_ma0 = np.ma.array(gm[0], mask=np.zeros(gm[0].shape, dtype=np.bool_))
    gm_ma1 = np.ma.array(gm[1], mask=np.zeros(gm[1].shape, dtype=np.bool_))
    gm_ma2 = np.ma.array(gm[2], mask=np.zeros(gm[2].shape, dtype=np.bool_))
    gm_ma3 = np.ma.array(gm[3], mask=np.zeros(gm[3].shape, dtype=np.bool_))
    gm_ma4 = np.ma.array(gm[4], mask=np.zeros(gm[4].shape, dtype=np.bool_))
    gm_ma5 = np.ma.array(gm[5], mask=np.zeros(gm[5].shape, dtype=np.bool_))
    var_gm0[:] = gm_ma0
    var_gm1[:] = gm_ma1
    var_gm2[:] = gm_ma2
    var_gm3[:] = gm_ma3
    var_gm4[:] = gm_ma4
    var_gm5[:] = gm_ma5
    ds_gm.close()


# fixme mm 4-7-2023:
# A lot of duplication with get_global_mean_vars
def get_global_mean_vars_all(
    model_folder,  # string e.g. "ab_classic"
    experiment_name,  # string, e.g. "CLASSIC_S2_"
    lat_var,
    lon_var,
):

    def nc_global_mean_file_name(experiment_name):
        return experiment_name + "gm_all_vars.nc"

    data_str = msh(model_folder).data_str
    names = data_str._fields
    conf_dict = confDict(model_folder)
    dataPath = Path(conf_dict["dataPath"])

    if dataPath.joinpath(
        nc_global_mean_file_name(experiment_name=experiment_name)
    ).exists():
        print(
            """ Found cached global mean files. If you want to recompute the global means
            remove the following files: """
        )
        print(
            dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))
        )

        def get_cached_global_mean(vn):
            gm_path = dataPath.joinpath(
                nc_global_mean_file_name(experiment_name=experiment_name)
            )
            return nc.Dataset(str(gm_path)).variables[vn].__array__()

        output = data_streams(*map(get_cached_global_mean, data_streams._fields))
        return output

    else:
        #gm = globalMask(file_name="common_mask_expanded.nc")
        gm = common_global_mask_expanded 
        msh_mod = msh(model_folder)
        tvm = msh_mod.template_var_name
        template = (
            nc.Dataset(
                dataPath.joinpath(
                    msh(model_folder).nc_file_name(
                        tvm, experiment_name=experiment_name
                    )
                )
            )
            .variables[tvm][0, :, :]
            .mask
        )
        # gcm=project_2(
        # source=gm,
        # target=CoordMask(
        # index_mask=np.zeros_like(template),
        # tr=SymTransformers(
        # ctr=msh(model_folder).make_model_coord_transforms(),
        # itr=msh(model_folder).make_model_index_transforms()
        # )
        # )
        # )
        gcm = resample_grid(
            source_coord_mask=gm,
            target_coord_mask=CoordMask(
                index_mask=np.zeros_like(template),
                tr=SymTransformers(
                    ctr=msh(model_folder).make_model_coord_transforms(),
                    itr=msh(model_folder).make_model_index_transforms(),
                ),
            ),
            var=gm.index_mask,
            method="nearest",
        )

        def compute_and_cache_global_mean(vn):
            path = dataPath.joinpath(
                msh(model_folder).nc_file_name(vn, experiment_name=experiment_name)
            )
            print(path)
            ds = nc.Dataset(str(path))
            vs = ds.variables
            lats = vs[lat_var].__array__()
            lons = vs[lon_var].__array__()
            print(vn)
            var = ds.variables[vn]
            # check if we have a cached version (which is much faster)
            gm_path = dataPath.joinpath(
                nc_global_mean_file_name(experiment_name=experiment_name)
            )

            # model_mask = gh.msh(model_folder).spatial_mask(dataPath=Path(conf_dict["dataPath"]))
            # combined_mask = combine_masks ([model_mask,gcm])
            gm = global_mean_var(
                lats,
                lons,
                # combined_mask.index_mask,
                gcm.index_mask,
                var,
            )
            # project to yearly averages for some fields 
            avg_candidates=["cVeg","cLitter", "cSoil", "gpp", "npp", "npp_nlim", "ra"]
            gm_y = gm if (var in avg_candidates) & gm.shape[0] < 500 else avg_timeline(gm ,12)
            return gm_y * 86400 if vn in ["gpp", "npp", "npp_nlim", "rh", "ra"] else gm_y

        # map variables to data
        print(data_str._fields)
        print("computing means, this may take some minutes...")
        output = data_str(*map(compute_and_cache_global_mean, data_str._fields))
        
        # fixme mm: 
        # hidden desing flaw (this is a model specific thing)
        # for models like SDGVM where pool data starts earlier than gpp data
        # this should be handled in model_specific helpers
        time_cut_candidates=frozenset(["cVeg","cLitter", "cSoil"]).intersection(
                data_str._fields
        )
        for var in time_cut_candidates:
             gm=output.__getattribute__(var)
             gms=gm.shape[0]
             gpps= output.gpp.shape[0]
             if gms > gpps:
                gm = gm[(gms - gpps): ]
        #if cVeg.shape[0] > gpp.shape[0]:
        #    cVeg = cVeg[cVeg.shape[0] - gpp.shape[0] :]
        #if "cLitter" in names and cLitter.shape[0] > gpp.shape[0]:
        #    cLitter = cLitter[cLitter.shape[0] - gpp.shape[0] :]
        #if cSoil.shape[0] > gpp.shape[0]:
        #    cSoil = cSoil[cSoil.shape[0] - gpp.shape[0] :]

        output_final = data_streams(
            cVeg=output.cVeg,
            cSoil=output.cLitter + output.cSoil if "cLitter" in names else output.cSoil,
            gpp=output.gpp,
            npp=output.npp if ("npp" in names) or ("npp_nlim" in names) else output.gpp - output.ra,
            ra=output.ra if "ra" in names else output.gpp - output.npp,
            rh=output.rh,
        )
        gm_path = dataPath.joinpath(
            nc_global_mean_file_name(experiment_name=experiment_name)
        )
        write_data_streams_cache(gm_path, output_final)
        return output_final


def nc_classes_2_masks(FilePath, var_name, classes, global_mask):
    g_mask = global_mask.index_mask
    ds = nc.Dataset(FilePath)
    var = ds.variables[var_name][:]
    ds.close()
    for nclass in classes:
        var_0 = var.copy()
        var_0[var_0 != nclass] = -9999
        var_0[var_0 == nclass] = 0
        var_0[var_0 == -9999] = 1
        var_0_corrected = var_0.copy()
        for i in range(var_0.shape[0]):
            var_0_corrected[i, :] = var_0[var_0.shape[0] - 1 - i, :]
        var_0_masked = np.logical_or(var_0_corrected, g_mask)
        cm_0 = CoordMask(var_0_masked, global_mask.tr)
        FileName = FilePath + "_mask_" + str(nclass) + ".nc"
        cm_0.write_netCDF4(FileName)
        print("File written as " + FileName)


data_streams = namedtuple("data_streams", ["cVeg", "cSoil", "gpp", "npp", "ra", "rh"])


def cached_var_dict(
    dataPath: Path,
    targetPath: Path,
    nc_global_mean_file_name: Callable,
    compute_arr_var_dict: Callable,
    names,
    flash_cache=False,
):
    if targetPath is None:
        targetPath = dataPath

    if flash_cache:
        for n in names:
            p = targetPath.joinpath(nc_global_mean_file_name(n))
            try:
                os.remove(p)
                print(f"removing {p}")

            except FileNotFoundError as e:
                print("could not find {p} to remove")

    try:
        return read_var_dict(
            targetPath,
            names,
            nc_global_mean_file_name,
        )
    except:
        arr_dict = compute_arr_var_dict(
            dataPath,
        )
        write_var_dict(arr_dict, targetPath, nc_global_mean_file_name)
        return arr_dict


def da_res_1_maker(
    # msh,#: module,
    make_proposer,
    make_param_filter_func,
    make_param2res_sym,
    make_weighted_cost_func,
    numeric_X_0,
    EstimatedParameters: type,
):
    def compute_func(
        mvs: CMTVS,
        svs,
        dvs,
        fcpa,
        epa_min,
        epa_max,
        epa_0,
        nsimu,
        acceptance_rate,
        chunk_size,
        D_init,
        K,
    ) -> Tuple[Dict, Dict, np.array]:
        """one of possibly many functions to reproduce the optimal parameters
        and startvalues to run the model forword. It replaces the model
        specific data assimilation procedures (different models used a
        different MCMC implementations and none worked for all combinations )
        this one does.
        """
        print(f"D_init={D_init}")
        c_max = np.array(epa_max)
        c_min = np.array(epa_min)
        c_0 = np.array(epa_0)
        epa_0 = EstimatedParameters(*c_0)
        C_autostep, J_autostep = autostep_mcmc_2(
            initial_parameters=np.array(epa_0),
            proposer=make_proposer(
                c_max=c_max, c_min=c_min, fcpa=fcpa, dvs=dvs, svs=svs
            ),
            filter_func=make_param_filter_func(
                c_max=c_max, c_min=c_min, fcpa=fcpa, dvs=dvs, svs=svs
            ),
            param2res=make_param2res_sym(mvs, fcpa, dvs, svs),
            costfunction=make_weighted_cost_func(svs),
            nsimu=nsimu,
            c_max=c_max,
            c_min=c_min,
            acceptance_rate=acceptance_rate,
            chunk_size=chunk_size,
            D_init=D_init,
            K=K,
        )
        # optimized parameter set (lowest cost function)
        J_vals = J_autostep[1]
        par_opt = C_autostep[:, np.argmin(J_vals)]
        epa_opt = EstimatedParameters(*par_opt)
        print("Data assimilation finished!")

        return (C_autostep, J_autostep, epa_opt)

    return compute_func


def cached_da_res_1_maker(
    make_proposer,
    make_param_filter_func,
    make_param2res_sym,
    make_weighted_cost_func,
    numeric_X_0,
    EstimatedParameters: type,
):
    compute_func = da_res_1_maker(
        make_proposer,
        make_param_filter_func,
        make_param2res_sym,
        make_weighted_cost_func,
        numeric_X_0,
        EstimatedParameters,
    )

    cached_func = mcdaf(compute_func, EstimatedParameters)
    return cached_func


def mcdaf(
    func: Callable,
    EstimatedParameters: type,
) -> Callable:
    # create a function with two more parameters: output_cache_path
    # than the original function and pipe through all the other parameters.
    # sig=inspect.signature(func)

    def cached_func(
        output_cache_path, *args, **kwargs  # add an argument
    ) -> Tuple[Dict, Dict, np.array]:
        # check for cached files
        sep = ","
        Cs_fn = "da_aa.csv"
        Js_fn = "da_j_aa.csv"
        epa_opt_fn = "epa_opt.json"

        def read_arr(p):
            return (
                pd.read_csv(
                    p,
                    sep=sep,
                    dtype=float,
                )
                .to_numpy()
                .transpose()
            )

        try:
            Cs = read_arr(output_cache_path.joinpath(Cs_fn))
            Js = read_arr(output_cache_path.joinpath(Js_fn))
            epa_opt = h.load_named_tuple_from_json_path(
                EstimatedParameters, output_cache_path.joinpath(epa_opt_fn)
            )
            fs = "\n".join([Cs_fn, Js_fn, epa_opt_fn])
            print(
                f"""
            Found cache files in {output_cache_path}:
            {fs}
            remove at least one of them to recompute
            """
            )

        # except Exception as e:
        except FileNotFoundError:
            Cs, Js, epa_opt = func(*args, **kwargs)
            if not output_cache_path.exists():
                output_cache_path.mkdir(parents=True)
            h.dump_named_tuple_to_json_path(
                epa_opt,
                output_cache_path.joinpath(epa_opt_fn),
            )
            pd.DataFrame(Cs.transpose(), columns=EstimatedParameters._fields).to_csv(
                output_cache_path.joinpath(Cs_fn), sep=sep, index=None
            )
            pd.DataFrame(Js.transpose()).to_csv(
                output_cache_path.joinpath(Js_fn), sep=sep, index=None
            )
        return Cs, Js, epa_opt

    return cached_func


def rh_iterator(
    mvs,
    X_0,  #: StartVector,
    par_dict,
    func_dict,
    delta_t_val=1,  # defaults to 1day timestep
):
    mit = minimal_iterator_internal(mvs, X_0, par_dict, func_dict, delta_t_val)

    def numfunc(expr):
        return hr.numerical_func_of_t_and_Xvec(
            state_vector=mvs.get_StateVariableTuple(),
            time_symbol=mvs.get_TimeSymbol(),
            expr=expr,
            parameter_dict=par_dict,
            func_dict=func_dict,
        )

    fd = {
        "rh": numfunc(mvs.get_AggregatedSoilCarbonOutFlux()),
        "cVeg": numfunc(mvs.get_AggregatedVegetationCarbon()),
        "cSoil": numfunc(mvs.get_AggregatedSoilCarbon()),
    }

    def make_func(key):
        return lambda t, X: fd[key](t, X)

    present_step_funcs = OrderedDict({key: make_func(key) for key in fd.keys()})
    mit.add_present_step_funcs(present_step_funcs)
    return mit


def single_stream_cost(obs, out_simu, key):
    # calculate costs for each data stream
    # take into account that model output might be shorter
    # than the observation (mainly for testing)
    simu = out_simu.__getattribute__(key)
    n_simu = simu.shape[0]
    _obs = obs.__getattribute__(key)[0:n_simu]  # cut obs to length
    cost = (n_simu) * np.sum((simu - _obs) ** 2, axis=0) / (_obs.mean(axis=0) ** 2)
    print(f"J_{key}  = {cost}")
    return cost


# constants
common_global_mask_expanded = globalMask("common_mask_expanded.nc")
