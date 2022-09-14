import numpy as np
from numpy.lib.function_base import meshgrid
from scipy import interpolate, sparse
from tqdm import tqdm
from typing import Callable, Dict, Tuple, Iterable, List
from functools import reduce, lru_cache
from copy import copy, deepcopy
from itertools import islice
from time import time
from sympy import var, Symbol, sin, Min, Max, pi, integrate, lambdify
from sympy.core.expr import Expr
from collections import namedtuple
from frozendict import frozendict
from importlib import import_module
from collections import OrderedDict
import matplotlib.pyplot as plt
import os
import inspect
import datetime as dt

from pathlib import Path
import json
import netCDF4 as nc

from ComputabilityGraphs.CMTVS import CMTVS
import CompartmentalSystems.helpers_reservoir as hr
from bgc_md2.resolve.mvars import (
    #    OutFluxesBySymbol,
    #    InternalFluxesBySymbol
    CompartmentalMatrix,
    InputTuple,
    StateVariableTuple,
)

# from bgc_md2.helper import bgc_md2_computers
import bgc_md2.display_helpers as dh
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
date = namedtuple("date", ["year", "month", "day"])


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


def globalMask() -> "CoordMask":
    # fixme mm 6-20-2022
    # this should actually be  package data
    # We create a target
    this_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ds = nc.Dataset(Path(this_dir).joinpath("common_mask.nc"))
    mask = ds.variables["mask"][:, :].data
    return CoordMask(index_mask=mask, tr=globalMaskTransformers(mask))


def identicalTransformers():
    """This function can be used if the grid used by the model has
    - lat ==   0 at the equator with
    - lat == -90 at the south pole,
    - lat == +90 at the north pole,
    - lon ==   0 at Greenich and
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

#_TraceTuple = namedtuple(
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
#)

class TraceTuple():
    def __init__(
        self,
        fd: dict
        ):    
        for k,v in fd.items():
            self.fd=fd
            self.__setattr__(k,v)

    @property     
    def _fields(self):
        return self.fd.keys()

    def averages(self, partitions):
        return self.__class__(
            {
                name:
                averaged_1d_array(self.__getattribute__(name), partitions)
                for name in self._fields
            }
        )

    def __add__(self, other):
        """overload + which is useful for averaging"""
        return self.__class__(
            {   
                name:
                self.__getattribute__(name) + other.__getattribute__(name)
                for name in self._fields
            }
        )

    def __truediv__(self, number):
        """overload / for scalars  which is useful for averaging"""

        return self.__class__(
            { 
                name: 
                self.__getattribute__(name) / number 
                for name in self._fields
            }
        )

    def __eq__(self, other):
        """overload == which is useful for tests"""
        return np.all(
            self.__getattribute__(name) == other.__getattribute__(name)
            for name in self._fields
        )


# some tiny helper functions for module loading
def mvs(mf):
    return import_module("{}.source".format(mf)).mvs


def msh(mf):
    return import_module("{}.model_specific_helpers_2".format(mf))


def th(mf):
    return import_module("{}.test_helpers".format(mf))


def confDict(mf):
    with Path(mf).joinpath("config.json").open(mode="r") as f:
        confDict = frozendict(json.load(f))
    return confDict


def test_args(mf):
    print("Loading data and parameters for " + mf + " model...")
    return th(mf).make_test_args(conf_dict=confDict(mf), msh=msh(mf), mvs=mvs(mf))


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

    sym_B = hr.euler_forward_B_sym(
        mvs.get_CompartmentalMatrix(), cont_time=t, delta_t=delta_t, iteration=it
    )
    # from IPython import embed;embed()
    sym_u = hr.discrete_time_sym(mvs.get_InputTuple(), t, delta_t, it)

    B_func = hr.numerical_array_func(
        state_vector=state_vector,
        time_symbol=it,
        expr=sym_B,
        parameter_dict=parameter_dict,
        func_dict=func_dict,
    )
    u_func = hr.numerical_array_func(
        state_vector=state_vector,
        time_symbol=it,
        expr=sym_u,
        parameter_dict=parameter_dict,
        func_dict=func_dict,
    )
    return (B_func, u_func)


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
    srm = mvs.get_SmoothReservoirModel()
    model_par_dict_keys = srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))] + list(mvs.get_StateVariableTuple())
    )
    # Parameter dictionary for the iterator
    apa = {**cpa._asdict(), **epa._asdict()}
    model_par_dict = {
        Symbol(k): v for k, v in apa.items() if Symbol(k) in model_par_dict_keys
    }
    return model_par_dict


def make_uniform_proposer(
    c_max: Iterable,
    c_min: Iterable,
    D: float,
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
    :param D: a parameter to regulate the proposer step. Higher D means smaller step size
    :param filter_func: model-specific function to filter out impossible parameter combinations
    """

    g = np.random.default_rng()

    def GenerateParamValues(c_op):
        paramNum = len(c_op)
        keep_searching = True
        while keep_searching:
            c_new = c_op + np.random.uniform(-0.5, 0.5, paramNum) * (
                (c_max - c_min) / D
            )
            if filter_func(c_new):
                keep_searching = False
        return c_new

    return GenerateParamValues


def make_uniform_proposer_2(
    c_max: Iterable,
    c_min: Iterable,
    D: float,
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
    :param D: a parameter to regulate the proposer step. Higher D means smaller step size
    :param filter_func: model-specific function to filter out impossible parameter combinations
    """

    g = np.random.default_rng()

    def GenerateParamValues(c_op, D):
        paramNum = len(c_op)
        keep_searching = True
        while keep_searching:
            c_new = c_op + (np.random.uniform(-0.5, 0.5, paramNum) * c_op * D)
            if filter_func(c_new):
                keep_searching = False
        return c_new

    return GenerateParamValues


def make_multivariate_normal_proposer(
    covv: np.ndarray,
    filter_func: Callable[[Iterable], bool],
) -> Callable[[Iterable], Iterable]:
    """Returns a function that will be used by mcmc algorithm to propose
    a new parameter(tuple) based on a given one.
    :param covv: The covariance matrix (usually estimated from a previously run chain)
    :param filter_func: model-specific function to filter out impossible parameter combinations
    """

    def GenerateParamValues(c_op):
        flag = True
        while flag:
            c_new = c_op + np.random.multivariate_normal(np.zeros(len(c_op)), covv)
            if filter_func(c_new):
                flag = False
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


# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
def autostep_mcmc(
    initial_parameters: Iterable,
    filter_func: Callable,
    param2res: Callable[[np.ndarray], np.ndarray],
    costfunction: Callable[[np.ndarray], np.float64],
    nsimu: int,
    c_max: np.ndarray,
    c_min: np.ndarray,
    acceptance_rate=15,
    chunk_size=100,
    D_init=1,
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
    :param acceptance_rate: Target acceptance rate in %, default is 20%
    :param chunk_size: number of iterations for which current acceptance ratio is assessed to modify the proposer step
    Set to 0 for constant step size. Default is 100.
    :param D_init: initial D value (increase to get a smaller step size), default = 1
    :param K: allowance coeffitient (increase K to reduce acceptance of higher cost functions), default = 1
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
    print("First_iteration done after " + str(round(time() - tb, 2)) + "s")
    print("Status update every 10 iterations:")
    # J_last = 400 # original code

    # initialize the result arrays to the maximum length
    # Depending on many of the parameters will be accepted only
    # a part of them will be filled with real values
    C_upgraded = np.zeros((paramNum, nsimu))
    J_upgraded = np.zeros((2, nsimu))
    D = D_init
    proposer = make_uniform_proposer(
        c_max=c_max, c_min=c_min, D=D * paramNum, filter_func=filter_func
    )
    # for simu in tqdm(range(nsimu)):
    st = time()
    accepted_current = 0
    if chunk_size == 0:
        chunk_size = (
            nsimu  # if chunk_size is set to 0 - proceed without updating step size.
        )
    for simu in range(nsimu):
        if (simu > 0) and (
            simu % chunk_size == 0
        ):  # every chunk size (e.g. 100 iterations) update the proposer step
            if accepted_current == 0:
                accepted_current = 1  # to avoid division by 0
            D = D * np.sqrt(
                acceptance_rate / (accepted_current / chunk_size * 100)
            )  # compare acceptance and update step
            if D < (
                1 / paramNum
            ):  # to avoid being stuck in too large steps that will always fail the filter.
                D = 1 / paramNum
            accepted_current = 0
            proposer = make_uniform_proposer(
                c_max=c_max, c_min=c_min, D=D * paramNum, filter_func=filter_func
            )
        if (
            simu % (chunk_size * 20) == 0
        ):  # every 20 chunks - return to the initial step size (to avoid local minimum)
            D = D_init

        c_new = proposer(C_op)
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
            upgraded = upgraded + 1
            accepted_current = accepted_current + 1

        # print some metadata
        # (This could be added to the output file later)
        if simu % 10 == 0 or simu == (nsimu - 1):
            print(
                """ 
               #(upgraded): {n}  | D value: {d} | overall acceptance rate: {r}%  
               progress: {simu:05d}/{nsimu:05d} {pbs} {p:02d}%
               time elapsed: {minutes:02d}:{sec:02d}
               overall min cost: {cost} achieved at {s} iteration | last accepted cost: {cost2} 
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
                    ac=accepted_current / chunk_size * 100,
                    # rr=int(accepted_current / chunk_size * 100),
                    ch=chunk_size,
                    d=round(D, 3),
                    s=J_min_simu,
                ),
                end="\033[5A",  # print always on the same spot of the screen...
            )

    # remove the part of the arrays that is still filled with zeros
    useful_slice = slice(0, upgraded)
    return C_upgraded[:, useful_slice], J_upgraded[:, useful_slice]


# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
def autostep_mcmc_2(
    initial_parameters: Iterable,
    filter_func: Callable,
    param2res: Callable[[np.ndarray], np.ndarray],
    costfunction: Callable[[np.ndarray], np.float64],
    nsimu: int,
    c_max: np.ndarray,
    c_min: np.ndarray,
    acceptance_rate=0.23,
    chunk_size=100,
    D_init=0.10,
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
    :param acceptance_rate: Target acceptance rate in %, default is 20%
    :param chunk_size: number of iterations for which current acceptance ratio is assessed to modify the proposer step
    Set to 0 for constant step size. Default is 100.
    :param D_init: initial D value (increase to get a smaller step size), default = 1
    :param K: allowance coeffitient (increase K to reduce acceptance of higher cost functions), default = 1
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
    print("First_iteration done after " + str(round(time() - tb, 2)) + "s")
    print("Status update every 10 iterations:")
    # J_last = 400 # original code

    # initialize the result arrays to the maximum length
    # Depending on many of the parameters will be accepted only
    # a part of them will be filled with real values
    C_upgraded = np.zeros((paramNum, nsimu))
    J_upgraded = np.zeros((2, nsimu))
    D = D_init
    proposer = make_uniform_proposer_2(
        c_max=c_max, c_min=c_min, D=D, filter_func=filter_func
    )
    # for simu in tqdm(range(nsimu)):
    st = time()
    accepted_current = 0
    if chunk_size == 0:
        chunk_size = (
            nsimu  # if chunk_size is set to 0 - proceed without updating step size.
        )
    for simu in range(nsimu):
        if (simu > 0) and (
            simu % chunk_size == 0
        ):  # every chunk size (e.g. 100 iterations) update the proposer step
            if accepted_current == 0:
                accepted_current = 1  # to avoid division by 0
            if accepted_current / chunk_size > acceptance_rate:
                D = D * (1 + 0.2)
            else:
                D = D * (1 - 0.2)
            # D = D * np.sqrt(
            #    acceptance_rate / (accepted_current / chunk_size * 100))  # compare acceptance and update step
            # if D < (1 / paramNum):  # to avoid being stuck in too large steps that will always fail the filter.
            #    D = (1 / paramNum)
            accepted_current = 0
            # proposer = make_uniform_proposer(c_max=c_max, c_min=c_min, D=D, filter_func=filter_func)
        if (
            simu % (chunk_size * 20) == 0
        ):  # every 20 chunks - return to the initial step size (to avoid local minimum)
            D = D_init

        c_new = proposer(C_op, D)
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
            upgraded = upgraded + 1
            accepted_current = accepted_current + 1

        # print some metadata
        # (This could be added to the output file later)
        if simu % 10 == 0 or simu == (nsimu - 1):
            print(
                """ 
               #(upgraded): {n}  | D value: {d} | overall acceptance rate: {r}%  
               progress: {simu:05d}/{nsimu:05d} {pbs} {p:02d}%
               time elapsed: {minutes:02d}:{sec:02d}
               overall min cost: {cost} achieved at {s} iteration | last accepted cost: {cost2} 
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
                    ac=accepted_current / chunk_size * 100,
                    # rr=int(accepted_current / chunk_size * 100),
                    ch=chunk_size,
                    d=round(D, 3),
                    s=J_min_simu,
                ),
                end="\033[5A",  # print always on the same spot of the screen...
            )

    # remove the part of the arrays that is still filled with zeros
    useful_slice = slice(0, upgraded)
    return C_upgraded[:, useful_slice], J_upgraded[:, useful_slice]


# Adaptive MCMC: with multivariate normal proposer based on adaptive covariance matrix
def adaptive_mcmc(
    initial_parameters: Iterable,
    covv: np.ndarray,
    filter_func: Callable,
    param2res: Callable[[np.ndarray], np.ndarray],
    costfunction: Callable[[np.ndarray], np.float64],
    nsimu: int,
    sd_controlling_factor=10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    performs the Markov chain Monte Carlo simulation an returns a tuple of the array of sampled parameter(tuples) with
    shape (len(initial_parameters),nsimu) and the array of cost function values with shape (q,nsimu)
    :param initial_parameters: The initial guess for the parameter (tuple) to be estimated
    :param covv: The covariance matrix (usually estimated from a previously run chain)
    :param filter_func: function to remove impossible parameter combinations
    :param param2res: A function that given a parameter(tuple) returns
    the model output, which has to be an array of the same shape as the observations used to
    build the cost function.
    :param costfunction: A function that given a model output returns a real number. It is assumed to be created for
    a specific set of observations, which is why they do not appear as an argument.
    :param nsimu: The length of the chain
    :param sd_controlling_factor: optional parameter to scale the covariance matrix. Increase to get a smaller step size
    """

    np.random.seed(seed=10)

    paramNum = len(initial_parameters)

    sd = 1 / sd_controlling_factor / paramNum
    covv = covv * sd

    proposer = make_multivariate_normal_proposer(covv, filter_func)

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
    #  from IPython import embed;embed()
    for simu in range(nsimu):
        # if (upgraded%10 == 0) & (upgraded > nsimu/20):
        if simu > nsimu / 10:
            covv = sd * np.cov(C_accepted)
            proposer = make_multivariate_normal_proposer(covv, filter_func)
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
            C_accepted = C_upgraded[:, 0:upgraded]
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

    # remove the part of the arrays that is still filled with zeros
    useful_slice = slice(0, upgraded)
    return C_upgraded[:, useful_slice], J_upgraded[:, useful_slice]


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
    performs the Markov chain Monte Carlo simulation an returns a tuple of the array of sampled parameter(tuples) with
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


def make_feng_cost_func(
    obs: np.ndarray,
) -> Callable[[np.ndarray], np.float64]:
    # Note:
    # in our code the dimension 0 is the time
    # and dimension 1 the pool index
    means = obs.mean(axis=0)
    mean_centered_obs = obs - means
    # now we compute a scaling factor per observable stream
    # fixme mm 10-28-2021
    #   The denominators in this case are actually the TEMPORAL variances of the data streams
    denominators = np.sum(mean_centered_obs ** 2, axis=0)

    #   The desired effect of automatically adjusting weight could be achieved
    #   by the mean itself.
    # dominators = means
    def costfunction(mod: np.ndarray) -> np.float64:
        cost = np.mean(np.sum((obs - mod) ** 2, axis=0) / denominators * 100)
        return cost

    return costfunction


def make_jon_cost_func(
    obs: np.ndarray,
) -> Callable[[np.ndarray], np.float64]:
    # Note:
    # in our code the dimension 0 is the time
    # and dimension 1 the pool index
    n = obs.shape[0]
    means = obs.mean(axis=0)
    denominators = means ** 2

    def costfunction(mod: np.ndarray) -> np.float64:
        cost = (100 / n) * np.sum(100 * np.sum((obs - mod) ** 2, axis=0) / denominators)
        return cost

    return costfunction


def days_per_month():
    # dpm= [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    dpm = [30 for i in range(12)]
    return dpm


def days_per_year():
    return sum(days_per_month())


def day_2_month_index(day):
    months_by_day_arr = np.concatenate(
        tuple(
            map(lambda m: m * np.ones(days_per_month()[m], dtype=np.int64), range(12))
        )
    )
    #  for variable months
    dpy = days_per_year()
    return months_by_day_arr[(day % dpy)] + int(day / dpy) * 12


def days_since_AD(iteration, delta_t_val, start_date):
    # To get the startdate we used the model specific function start_date()
    # where it can be derivde using the original calendar (365 day or 360 day)
    # from then on we might use a different counting (see days_per_month())
    start_year, start_month, start_day = start_date
    td_AD = (
        start_year * days_per_year()
        + sum(days_per_month()[0 : (start_month - 1)])
        + (start_day - 1)
    )
    return td_AD + iteration * delta_t_val


def month_2_day_index(ns, start_date):
    start_month = start_date.month
    """ computes the index of the day at the end of the month n in ns
    this works on vectors and is faster than a recursive version working
    on a single index (since the smaller indices are handled anyway)
    """

    # We first compute the sequence of day indices up to the highest month in ns
    # and then select from this sequence the day indices for the months in ns
    d = days_per_month()
    dpm = (d[i % len(d)] for i in range(max(ns)))

    # compute indices for which we want to store the results which is the
    # list of partial sums of the above list  (repeated)

    def f(acc, el):
        if len(acc) < 1:
            res = (el,)
        else:
            last = acc[-1]
            res = acc + (el + last,)
        return res

    day_indices_for_continuous_moths = reduce(f, dpm, (0,))
    day_indices = reduce(
        lambda acc, n: acc + [day_indices_for_continuous_moths[n]],  # for n=0 we want 0
        ns,
        [],
    )
    return day_indices


# +
# def year_2_day_index(ns):
#     """ computes the index of the day at the end of the year n in ns
#     this works on vectors
#     """
#     return np.array(list(map(lambda n:days_per_year*n,ns)))
# -


class TimeStepIterator2:
    """iterator for looping forward over the results of a difference equation
    X_{i+1}=f(X_{i},i)"""

    def __init__(
        self,
        initial_values,  # a tuple of values that will be
        f,  # the function to compute the next ts
        max_it=False,
    ):
        self.initial_values = initial_values
        self.f = f
        self.reset()
        self.max_it = max_it

    def reset(self):
        self.i = 0
        self.ts = self.initial_values

    def __iter__(self):
        self.reset()
        return self

    def __next__(self):
        if self.max_it:
            if self.i == self.max_it:
                raise StopIteration

        ts = copy(self.ts)
        ts_new = self.f(self.i, ts)
        self.ts = ts_new
        self.i += 1
        return ts

    def values(self, day_indices):
        # we traverse the iterator to the highest index and
        # collect the results we want to keep in a list (acc)
        tsi = copy(self)
        tsi.reset()

        def g(acc, i):
            v = tsi.__next__()
            if i in day_indices:
                acc += [v]
            return acc

        xs = reduce(g, range(max(day_indices) + 1), ([]))
        return xs


def respiration_from_compartmental_matrix(B, X):
    """This function computes the combined respiration from all pools"""
    return -np.sum(B @ X)


def plot_solutions(fig, times, var_names, tup, names=None):
    if names is None:
        names = tuple(str(i) for i in range(len(tup)))

    # from IPython import embed; embed()
    assert all([tup[0].shape == el.shape for el in tup])

    if tup[0].ndim == 1:
        n_times = tup[0].shape[0]
        ax = fig.subplots(1, 1)
        for i, sol in enumerate(tup):
            ax.plot(
                np.array(times).reshape(
                    n_times,
                ),
                sol,
                marker="o",
                label=names[i],
            )
            ax.set_title(var_names[0])
            ax.legend()
    else:
        n_times, n_vars = tup[0].shape

        fig.set_figheight(n_vars * fig.get_figwidth())
        axs = fig.subplots(n_vars, 1)
        colors = ("red", "blue", "green", "orange")
        for j in range(n_vars):
            for i, sol in enumerate(tup):
                axs[j].plot(
                    np.array(times).reshape(
                        n_times,
                    ),
                    sol[:, j],
                    marker="+",
                    label=names[i],
                )
                axs[j].set_title(var_names[j])
                axs[j].legend()


def plot_observations_vs_simulations(fig, svs_cut, obs_simu):
    # svn and obs_simu are both
    n_plots = len(svs_cut)
    fig.set_figheight(n_plots * fig.get_figwidth())
    axs = fig.subplots(n_plots)
    for i, name in enumerate(svs_cut.__class__._fields):
        var = svs_cut[i]
        var_simu = obs_simu[i]
        axs[i].plot(range(len(var_simu)), var_simu, label="simulation", zorder=2)
        axs[i].plot(range(len(var)), var, label="observation", zorder=1)
        axs[i].legend()
        axs[i].set_title(name)


def is_infinite_rec(chunk):
    return reduce_or_rec(
        np.logical_not(
            np.isfinite(chunk)
        )
    )


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
        # from IPython import embed;embed()
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
    The important thing is not to call this functions with the variables but the arrays.
    """

    # copy the mask from the array (first time step)
    weight_mask = arr.mask[0, :, :] if arr.mask.any() else False

    weight_mat = np.ma.array(get_weight_mat(lats, lons), mask=weight_mask)
    wms = weight_mat.sum()
    # to compute the sum of weights we add only those weights that
    # do not correspond to an unmasked grid cell
    return ((weight_mat * arr).sum(axis=(1, 2)) / wms).data


def global_mean_var(
    lats: np.ma.core.MaskedArray,
    lons: np.ma.core.MaskedArray,
    mask: np.array,
    var: nc._netCDF4.Variable,
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

    n_t = var.shape[0]
    res = np.zeros(n_t)
    for it in tqdm(range(n_t)):
        el = (weight_mat * var[it, :, :]).sum() / wms
        res[it] = el
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
    experiments = ["S2"] # We are using s2 data, can add more
):
    import paramiko
    import tarfile
    import gzip
    import shutil

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
                elif model == "JULES-ES"or model == "JULES-ES-1.0":
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


def monthly_to_yearly(monthly):
    # TRENDY specific - months weighted like all months are 30 days
    if len(monthly.shape) > 1:
        sub_arrays = [
            monthly[i * 12 : (i + 1) * 12, :, :]
            for i in range(int(monthly.shape[0] / 12))
        ]
    else:
        sub_arrays = [
            monthly[
                i * 12 : (i + 1) * 12,
            ]
            for i in range(int(monthly.shape[0] / 12))
        ]
    return np.stack(list(map(lambda sa: sa.mean(axis=0), sub_arrays)), axis=0)


def pseudo_daily_to_yearly(daily):
    # compute a yearly average from pseudo daily data
    # for one data point
    pseudo_days_per_year = pseudo_days_per_month * 12
    sub_arrays = [
        daily[i * pseudo_days_per_year : (i + 1) * pseudo_days_per_year, :]
        for i in range(int(daily.shape[0] / pseudo_days_per_year))
    ]
    return np.stack(list(map(lambda sa: sa.mean(axis=0), sub_arrays)), axis=0)


def make_feng_cost_func_2(svs):  #: Observables
    # now we compute a scaling factor per observable stream
    # fixme mm 10-28-2021
    # The denominators in this case are actually the TEMPORAL variances of the data streams
    obs_arr = np.stack([arr for arr in svs], axis=1)
    means = obs_arr.mean(axis=0)
    mean_centered_obs = obs_arr - means
    denominators = np.sum(mean_centered_obs ** 2, axis=0)

    def feng_cost_func_2(simu):  #: Observables
        def f(i):
            arr = simu[i]
            obs = obs_arr[:, i]
            diff = ((arr - obs) ** 2).sum() / denominators[i] * 100
            return diff

        return np.array([f(i) for i in range(len(simu))]).mean()

    return feng_cost_func_2


def make_param_filter_func(
    c_max, c_min, betas: List[str] = []
) -> Callable[[np.ndarray], bool]:

    positions = [c_max.__class__._fields.index(beta) for beta in betas]

    def isQualified(c):
        cond1 = (c >= c_min).all()
        cond2 = (c <= c_max).all()

        cond3 = np.sum([c[p] for p in positions]) <= 1
        return cond1 and cond2 and cond3

    return isQualified


def make_StartVectorTrace(mvs):
    # deprecated
    svt = mvs.get_StateVariableTuple()
    return namedtuple(
        "StartVectorTrace",
        [str(v) for v in svt]
        + [str(v) + "_p" for v in svt]
        + [str(v) + "_c" for v in svt]
        + [str(v) + "_RT" for v in svt],
    )


def make_InitialStartVectorTrace(X_0, mvs, par_dict, func_dict):
    # test the new iterator

    # we make the X_p and X_c  parts compatible with the ones computed by the iterator
    # for following timesteps (to keep it)
    # As you can see in the definition of the iterator these values have no impact on further results
    B_func, u_func = make_B_u_funcs_2(mvs, par_dict, func_dict)
    I = u_func(0, X_0)
    u = I.sum()
    b = I / u
    B = B_func(0, X_0)
    B_inf = np.linalg.inv(B)
    X_p_0 = B_inf @ I
    X_c_0 = X_0 + X_p_0
    RT_0 = B_inf @ b
    # combine the three
    # here we rely on order to be consistent
    # (although we could use also the names of the namedtuple)
    V_arr = np.concatenate((X_0, X_p_0, X_c_0, RT_0), axis=0)
    StartVectorTrace = make_StartVectorTrace(mvs)
    V_init = StartVectorTrace(*V_arr)
    return V_init


# fixme: mm 04-22-2022
# this function is deprecated rather use traceability_iterator
# def make_daily_iterator_sym_trace(
#        mvs,
#        V_init, #: StartVectorTrace,
#        par_dict,
#        func_dict
#    ):
#    B_func, I_func = make_B_u_funcs_2(mvs,par_dict,func_dict)
#    V_arr=np.array(V_init).reshape(-1,1) #reshaping for matmul which expects a one column vector (nr,1)
#
#    n=len(mvs.get_StateVariableTuple())
#    def f(it,V):
#        #the pools are the first n values
#        X = V[0:n]
#        I = I_func(it,X)
#        # we decompose I
#        u=I.sum()
#        b=I/u
#        B = B_func(it,X)
#        B_inf = np.linalg.inv(B)
#        X_new = X + I + B @ X
#        X_p = B_inf @ I
#        X_c = X_new+X_p
#        RT = B_inf @ b
#        V_new = np.concatenate(
#            (
#                X_new.reshape(n,1),
#                X_p.reshape(n,1),
#                X_c.reshape(n,1),
#                RT.reshape(n,1),
#            ),
#            axis=0
#        )
#        return V_new
#
#    return TimeStepIterator2(
#        initial_values=V_arr,
#        f=f,
#    )
#


class InfiniteIterator:
    def __init__(self, x0, func):  # ,n):
        self.x0 = x0
        self.func = func

        self.cur = x0
        self.pos = 0

    def __iter__(self):
        # return a fresh instance that starts from the first step)
        # This is important for the __getitem__ to work
        # as expected and not have side effects
        # for
        # res1=itr[0:10]
        # res2=itr[0:10]

        c = self.__class__(self.x0, self.func)
        return c
        # return self

    def __next__(self):
        # print(self.pos, self.cur)
        val = self.func(self.pos, self.cur)
        self.cur = val
        self.pos += 1
        return val
        # raise StopIteration()

    # @lru_cache
    def value_at(self, it_max):
        I = self.__iter__()

        def f_i(acc, i):
            return I.__next__()

        return reduce(f_i, range(it_max), I.x0)

    def __getitem__(self, arg):
        # this functions implements the python index notation itr[start:stop:step]
        # fixme mm 4-26-2022
        # we could use the cache for value_at if we dont use
        if isinstance(arg, slice):
            start = arg.start
            stop = arg.stop
            step = arg.step
            return tuple(islice(self, start, stop, step))

        elif isinstance(arg, int):
            return (self.value_at(it_max=arg),)
        else:
            raise IndexError(
                """arguments to __getitem__ have to be either
                indeces or slices."""
            )


def values_2_TraceTuple(tups):
    # instead of the tuple of TraceTuples that the InfiniteIterator returns
    # we want a TraceTuple of arrays whith time (iterations)  added as the first dimension
    return TraceTuple(
        {
            name:
            np.stack(tuple((tup.__getattribute__(name) for tup in tups)))
            for name in tups[0]._fields
        }
    )


class TraceTupleIterator(InfiniteIterator):
    # overload methods specific to the TraceTupleIterator
    def __getitem__(self, arg):
        # we call the [] method of the superclass
        # which returns a tuple of TraceTuples
        tups = super().__getitem__(arg)

        # But we want a TraceTuple of arrays
        return values_2_TraceTuple(tups)

    def averaged_values(self, partitions):
        start = partitions[0][0]

        def tt_avg(tups):
            # note that this involves overloaded + and /
            # for TraceTuple
            l = len(tups)
            return reduce(lambda acc, el: acc + el, tups) / l

        # make a copy to avoid side effects on self
        itr = self.__iter__()
        # move to the start
        for i in range(start):
            itr.__next__()
        # from IPython import embed;embed()

        tts = [
            tt_avg([itr.__next__() for i in range(stop_p - start_p)])
            for (start_p, stop_p) in partitions
        ]
        return values_2_TraceTuple(tts)


def make_trace_tuple_func(
        traced_functions: Dict[str, Callable]
    ):

    # create a function that produces all the values 
    def f(X, B, I, it):
        # These are values that are computable from the momentary values of X and B
        u = I.sum()
        b = I / u
        B_inv = np.linalg.inv(B)
        X_c = B_inv @ I
        X_p = X_c - X
        X_dot = I - B @ X
        RT = X_c / u  # =B_inv@b but cheeper to compute
        # we now compute the system X_c and X_p
        # This is in general not equal to the sum of the component, but rather a property
        # of a surrogate system.
        x = X.sum()
        x_dot = X_dot.sum()
        m_s = (B @ X).sum() / x
        x_c = 1 / m_s * u
        x_p = x_c - x
        rt = x_c / u
            
        static={ 
            "X": X,
            "X_p": X_p,
            "X_c": X_c,
            "X_dot": X_dot,
            "RT": RT,
            "x": x,
            "x_p": x_p,
            "x_c": x_c,
            "x_dot": x_dot,
            "rt": rt,
            "u": u,
        }
        dynamic={
            k: f(it,*X)
            for k,f in traced_functions.items()
        }
        def dict_merge(d1,d2):
            d=deepcopy(d1)
            d.update(d2)
            return d

        return TraceTuple(
            dict_merge(static,dynamic)
        )

    return f 

# def trace_tuple_instance(X, B, I):
#     # These are values that are computable from the momentary values of X and B
#     u = I.sum()
#     b = I / u
#     B_inv = np.linalg.inv(B)
#     X_c = B_inv @ I
#     X_p = X_c - X
#     X_dot = I - B @ X
#     RT = X_c / u  # =B_inv@b but cheeper to compute
#     # we now compute the system X_c and X_p
#     # This is in general not equal to the sum of the component, but rather a property
#     # of a surrogate system.
#     x = X.sum()
#     x_dot = X_dot.sum()
#     m_s = (B @ X).sum() / x
#     x_c = 1 / m_s * u
#     x_p = x_c - x
#     rt = x_c / u
#
#
#     return TraceTuple(
#         X=X,
#         X_p=X_p,
#         X_c=X_c,
#         X_dot=X_dot,
#         RT=RT,
#         x=x,
#         x_p=x_p,
#         x_c=x_c,
#         x_dot=x_dot,
#         rt=rt,
#         u=u,
#         AggregatedVegetation2SoilCarbonFlux=0
#     )
#

def traceability_iterator(
        X_0,
        func_dict,
        mvs,  #: CMTVS,
        dvs,  #: Drivers,
        cpa,  #: Constants,
        epa,  #: EstimatedParameters
        delta_t_val: int = 1,  # defaults to 1day timestep
        traced_expressions: Dict[str, Expr] = dict()
    ):
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict = make_param_dict(mvs, cpa, epa)

    B_func, I_func = make_B_u_funcs_2(mvs, par_dict, func_dict, delta_t_val)

    # we produce functions of f(it,x_a,...,x_z) to compute the tracable expressions
    traced_funcs = {
        key: numfunc(
            expr_cont=val,
            mvs=mvs,
            delta_t_val=delta_t_val,
            par_dict=par_dict,
            func_dict=func_dict
        )
        for key,val in traced_expressions.items()
    }
    #from IPython import embed; embed()
    trace_tuple_instance = make_trace_tuple_func(traced_funcs)
    # define the start tuple for the iterator
    V_init = trace_tuple_instance(
        X_0,
        # in Yiqi's nomenclature: dx/dt=I-Bx
        # instead of           : dx/dt=I+Bx
        # as assumed by B_u_func
        -B_func(0, X_0),
        I_func(0, X_0),
        it=0
    )

    # define the function with V_{i+1}=f(i,V_i)
    def f(it: int, V: TraceTuple) -> TraceTuple:
        X = V.X
        I = I_func(it, X)
        B = B_func(it, X)
        X_new = X + I + B @ X
        return trace_tuple_instance(X_new, -B, I, it)

    return TraceTupleIterator(x0=V_init, func=f)


def write_global_mean_cache(gm_path, gm: np.array, var_name: str):
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


def get_cached_global_mean(gm_path, vn):
    return nc.Dataset(str(gm_path)).variables[vn].__array__()


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
    # from IPython import embed; embed()
    p = sparse.bsr_array((data, (row, col)), dtype=np.int64)
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


def project_2(source: CoordMask, target: CoordMask):

    s_mask = source.index_mask  # .astype(np.float64,casting='safe')
    print(s_mask.mean())
    s_n_lat, s_n_lon = s_mask.shape

    t_mask = target.index_mask
    t_n_lat, t_n_lon = t_mask.shape

    t_tr = target.tr
    s_tr = source.tr
    # we first create the (unified) lat lon coordinates of the source array
    # and interpolate the source mask on them
    lats = np.array(list(map(s_tr.i2lat, range(s_n_lat))))
    lons = np.array(list(map(s_tr.i2lon, range(s_n_lon))))
    f = interpolate.interp2d(x=lons, y=lats, z=s_mask)
    # now we apply this interpolating function to the target grid
    # points
    target_lats = np.array(list(map(t_tr.i2lat, range(t_n_lat))))
    target_lons = np.array(list(map(t_tr.i2lon, range(t_n_lon))))
    # order lats and lons since f returns an array that assumes this anyway
    otlats, p_lat, p_lat_inv = target.ordered_lats()
    otlons, p_lon, p_lon_inv = target.ordered_lons()
    float_grid = p_lat_inv @ f(otlons, otlats) @ p_lon_inv
    print(float_grid.mean())
    projected_mask = float_grid > 0.5
    # from IPython import embed; embed()
    return CoordMask(index_mask=np.logical_or(projected_mask, t_mask), tr=t_tr)


# outputs a table with flow diagrams, compartmental matrices and allocation vectors
def model_table(
    model_names,  # dictionary (folder name : model name)
):
    model_folders = [(k) for k in model_names]
    mf = model_folders[0]
    import_module("{}.source".format(mf))

    def mvs(mf):
        return import_module("{}.source".format(mf)).mvs

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
            # from IPython import embed; embed()
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


def read_or_create(
    path: Path, create_and_write: Callable[[Path], T], read: Callable[[Path], T]
) -> T:
    """
    A the basic cache functionality:
    path: the path of the cache file
    create_and_write: The funciton that creates the  object to be cached and writes it to <path>.
    read: The function that extracts the cached object from <path>.
    """
    print("type", type(path))
    if path.exists():
        print(
            "Found cache file {}. If you want to recompute the result remove it.".format(
                str(path)
            )
        )
        return read(path)
    else:
        pp = path.parent
        if not (pp.exists()):
            pp.mkdir(parents=True)
        return create_and_write(path)


# +
# For harmonizing timelines and plotting
# -

# a function to get a list of test_args from all models involved in the comparison
def get_test_arg_list(model_folders):
    test_arg_list = []
    for mf in model_folders:
        current_ta = test_args(mf)
        test_arg_list.append(current_ta)
    print("Done!")
    return test_arg_list
# fixme mm 8-12: 
# it would make sense to create a dictionary indexed by the model name 
# so it can be used for model comparisons by name like this one

def get_test_arg_dict(model_folders):
    return { mf:  test_args(mf)
        for mf in model_folders
    }    


def times_in_days_aD(
    test_args,  # a tuple defined in model-specific test_helpers.py
    delta_t_val,  # iterator time step
):
    n_months = len(test_args.dvs[0])
    n_days = month_2_day_index([n_months], test_args.start_date)[0]
    n_iter = int(n_days / delta_t_val)
    return np.array(
        tuple(
            (
                days_since_AD(i, delta_t_val, test_args.start_date)
                for i in np.arange(n_iter)  # days_after_sim_start
            )
        )
    )


# function to determine overlapping time frames for models simulations
def t_min_tmax_overlap(
    test_arg_list,  # a list of test_args from all models involved
    delta_t_val,  # iterator time step
):
    td = {
        i: times_in_days_aD(test_arg_list[i], delta_t_val)
        for i in range(len(test_arg_list))
    }
    t_min = max([t.min() for t in td.values()])
    t_max = min([t.max() for t in td.values()])
    return (t_min, t_max)


def t_min_tmax_full(
    test_arg_list,  # a list of test_args from all models involved
    delta_t_val,  # iterator time step
):
    td = {
        i: times_in_days_aD(test_arg_list[i], delta_t_val)
        for i in range(len(test_arg_list))
    }
    t_min = min([t.min() for t in td.values()])
    t_max = max([t.max() for t in td.values()])
    return (t_min, t_max)


# function find the timesteps corresponding to shared times
def min_max_index(
    test_args,  # a tuple defined in model-specific test_helpers.py
    delta_t_val,  # iterator time step
    t_min,  # output of t_min_tmax_overlap or t_min_tmax_full
    t_max,  # output of t_min_tmax_overlap or t_min_tmax_full
):
    ts = times_in_days_aD(test_args, delta_t_val)

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


def traceability_iterator_instance(
    mf,  # name of the model folder
    test_args,  # a tuple defined in model-specific test_helpers.py
    delta_t_val,  # model time step
):
    """Wrapper that just needs the folder name and the step
    size.
    """
    mvs_t = mvs(mf)
    dvs_t = test_args.dvs
    cpa_t = test_args.cpa
    epa_t = test_args.epa_opt
    X_0 = msh(mf).numeric_X_0(mvs_t, dvs_t, cpa_t, epa_t)
    func_dict = msh(mf).make_func_dict(mvs_t, dvs_t, cpa_t, epa_t)

    return traceability_iterator(
        X_0,
        func_dict,
        mvs=mvs_t,
        dvs=dvs_t,
        cpa=cpa_t,
        epa=epa_t,
        delta_t_val=delta_t_val,
    )


# fixme: possible obsolete = averaged_1d_array
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

# this is to export yearly traceable components for use outside of framework (experimenting with Enqing's code)
import pandas as pd
def write_yearly_components(
    model_names,  # dictionary (folder name : model name)
    test_arg_list,  # a list of test_args from all models involved
    var_names,  # dictionary (trace_tuple name : descriptive name)
    delta_t_val,  # model time step
    model_cols,  # dictionary (folder name :color)
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging=12,  # number of iterator steps over which to average results. 1 for no averaging
):
    if (part < 0) | (part > 1):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        )
    model_folders = [(k) for k in model_names]
    variables = [(k) for k in var_names]
    n = len(variables)

    Xs=np.array((0,))
    X_cs=np.array((0,))
    X_ps=np.array((0,))
    us=np.array((0,))
    RTs=np.array((0,))
 
    for i, name in enumerate(variables):
        k = 0
        for mf in model_folders:
            itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
            start_min, stop_max = min_max_index(
                test_arg_list[k],
                delta_t_val,
                *t_min_tmax_overlap(test_arg_list, delta_t_val)
            )
            # if we do not want the whole interval but look at a smaller part to observe the dynamics
            # start,stop = start_min, int(start_min+(stop_max-start_min)*part)
            start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
            times = (
                times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
                / days_per_year()
            )
            print("Gathering " + str(name) + " data for " + str(mf) + " model")
            vals = itr[start:stop]
            vals_output = avg_timeline(vals.__getattribute__(name), averaging)
            
            if name == "x":  
                Xs=np.concatenate((Xs,vals_output))
            if name == "x_c":  
                X_cs=np.concatenate((X_cs,vals_output))  
            if name == "x_p":  
                X_ps=np.concatenate((X_ps,vals_output)) 
            if name == "u":  
                us=np.concatenate((us,vals_output))                  
            if name == "rt":  
                RTs=np.concatenate((RTs,vals_output)) 
            k += 1       
        
    # Saving yearly components to a csv for further analysis
    test_keys = ["Model", "Year", "X", "X_c", "X_p", "u", "RT"]
        
    model = np.repeat(model_folders[0], 160).reshape(160,1) 
    models=np.array((model)).reshape(160,1)
    for i in range (len(model_folders)-1):
        model=np.repeat(model_folders[i+1], 160).reshape(160,1) 
        models=np.concatenate((models,model))
        
    print("Model")
    print(models.shape)
    #print(models)  
        
    year = np.array(range(1860,2020)).reshape(160,1)
    years=np.array((year)).reshape(160,1)
    for i in range (len(model_folders)-1):
        years=np.concatenate((years,year))
        
    print("Year")
    print(years.shape)
    
    test_values = [models.flatten(), years.flatten(), Xs.flatten()[1:], X_cs.flatten()[1:], X_ps.flatten()[1:], us.flatten()[1:], RTs.flatten()[1:]]
    
    template = {test_keys[i]: test_values[i] for i in range(len(test_keys))}
    output = pd.DataFrame(template)    
    
    pd.DataFrame(output).to_csv('traceable_components_output.csv', sep=',')
    
    return output    

def get_traceable_components(
    model_names,  # dictionary (folder name : model name)
    test_arg_list,  # a list of test_args from all models involved
    delta_t_val,  # model time step
    model_cols,  # dictionary (folder name :color)
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging,  # number of iterator steps over which to average results. 1 for no averaging
    overlap=True,  # compute overlapping timeframe or plot whole duration for all models
):
    if (part < -1) | (part > 1) | (part == 0):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between -1 and 1 excluding 0"
        )
    model_folders = [(k) for k in model_names]
    k = 0
    sum_x = np.array(0)
    sum_x_c = np.array(0)
    sum_x_p = np.array(0)
    sum_u = np.array(0)
    sum_rt = np.array(0)
    all_components = list ()
    for mf in model_folders:
        print ("Computing traceable components for "+mf+"...")
        itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        if overlap == True:
            start_min, stop_max = min_max_index(
                test_arg_list[k],
                delta_t_val,
                *t_min_tmax_overlap(test_arg_list, delta_t_val)
            )
        else:
            start_min, stop_max = min_max_index(
                test_arg_list[k],
                delta_t_val,
                *t_min_tmax_full(test_arg_list, delta_t_val)
            )
        # if we do not want the whole interval but look at a smaller part to observe the dynamics
        if part < 0:
            start, stop = int(stop_max - (stop_max - start_min) * abs(part)), stop_max
        else:
            start, stop = start_min, int(start_min + (stop_max - start_min) * part)
        times = (
            times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            / days_per_year()
        )
        vals = itr[start:stop]
        k += 1       
        
        times_avg = avg_timeline(times, averaging)     
        
        comp_dict = {
            "x": avg_timeline(vals.x, averaging),
            "x_c": avg_timeline(vals.x_c, averaging),
            "x_p": avg_timeline(vals.x_p, averaging),
            "u": avg_timeline(vals.u, averaging),
            "rt": avg_timeline(vals.rt, averaging),
            }
        all_components.append(comp_dict) 
        if sum_x.all()==0:
            sum_x = np.append(sum_x,vals.x)[1:]
            sum_x_c = np.append(sum_x_c,vals.x_c)[1:]
            sum_x_p = np.append(sum_x_p,vals.x_p)[1:]
            sum_u = np.append(sum_u,vals.u)[1:]
            sum_rt = np.append(sum_rt,vals.rt)[1:]
        else:
            sum_x = sum_x + vals.x
            sum_x_c = sum_x_c + vals.x_c
            sum_x_p = sum_x_p + vals.x_p
            sum_u = sum_u + vals.u
            sum_rt = sum_rt + vals.rt   
                
    
    ave_dict = {
            "x": avg_timeline(sum_x / len(model_folders), averaging),
            "x_c": avg_timeline(sum_x_c / len(model_folders), averaging),
            "x_p": avg_timeline(sum_x_p / len(model_folders), averaging),
            "u": avg_timeline(sum_u / len(model_folders), averaging),
            "rt": avg_timeline(sum_rt / len(model_folders), averaging),
            }       
    
    all_components.append(ave_dict)
    all_components.append(times_avg)       
    
    mods=list(model_names.values())
    mods.append("Mean")
    mods.append("Times")
    all_comp_dict = {mods[i]: all_components[i] for i in range(len(mods))}      

    return all_comp_dict 
    
def plot_traceable_component(
    all_comp_dict,
    comp_name,
    model_cols,
    delta=False,
):
    models=list(all_comp_dict.keys())[:-2]
    fig = plt.figure(figsize=(17, 8))
    ax = fig.subplots(1, 1)
    vals_mean=0
    st_dev=0
    if comp_name == "x_x_c":
        ax.plot(
                all_comp_dict["Times"],
                all_comp_dict["Mean"]["x"]
                * 148940000
                * 1000000
                * 0.000000000001,  # convert to global C in Gt
                label="Multi-model mean - X",
                color="black",
                linewidth=3,
        )
        ax.plot(
                all_comp_dict["Times"],
                all_comp_dict["Mean"]["x_c"]
                * 148940000
                * 1000000
                * 0.000000000001,  # convert to global C in Gt
                label="Multi-model mean - X_c",
                color="black",
                linestyle="dashed",
                linewidth=2,
        )
        for m in models:      
            
            ax.plot(
                    all_comp_dict["Times"],
                    all_comp_dict[m]["x"]
                    * 148940000
                    * 1000000
                    * 0.000000000001,  # convert to global C in Gt
                    label=m + " - X",
                    color=model_cols[m],
                    linewidth=0.8,
                )
            ax.plot(
                    all_comp_dict["Times"],
                    all_comp_dict[m]["x_c"]
                    * 148940000
                    * 1000000
                    * 0.000000000001,  # convert to global C in Gt
                    label=m + " - X_c",
                    color=model_cols[m],
                    linestyle="dashed",
                    linewidth=0.8,
            )
    else:
        if delta:
            vals_array_mean = all_comp_dict["Mean"][comp_name]
            vals_mean = (vals_array_mean-vals_array_mean[0]) #/ vals_array_mean[0] * 100
        else:
            vals_mean = all_comp_dict["Mean"][comp_name]
                
        ax.plot(
                all_comp_dict["Times"],
                vals_mean,
                label="Multi-model mean",
                color="black",
                linewidth=3, 
                )                
        
        diff_sqr = np.zeros_like(all_comp_dict[models[0]][comp_name])             
        
        for m in models:
            if delta:
                vals_array = all_comp_dict[m][comp_name]
                vals=(vals_array-vals_array[0])#/vals_array[0] * 100
            else:
                vals = all_comp_dict[m][comp_name]  
                      
            diff_sqr=diff_sqr + (vals-vals_mean)**2
            
            ax.plot(
                    all_comp_dict["Times"],
                    vals,
                    label=m,
                    color=model_cols[m],
                    linewidth=0.8, 
                )
        variance = diff_sqr / (len(models)-1)
        st_dev=np.sqrt(variance)
      
        ax.fill_between(
                all_comp_dict["Times"],
                vals_mean-st_dev*2, 
                vals_mean+st_dev*2,
                label="\u00B12$\sigma$ confidence interval",
                color="grey",
                alpha=0.2                
                )                          
    ax.grid()
    if delta:
        ax.set_title("$\Delta$ "+str(comp_name))
    else:
        ax.set_title(comp_name)
    #ax.set_ylabel("Gt C")
    ax.legend(bbox_to_anchor =(1.3, 1))   
    return(vals_mean, st_dev)
         
def plot_attribution (
    times,
    delta_x_pos,
    delta_x_neg,
    rt_contrib_pos,
    rt_contrib_neg,
    u_contrib_pos,
    u_contrib_neg,
    rt_u_inter_pos,
    rt_u_inter_neg,
    x_p_contrib_pos,
    x_p_contrib_neg,
    percent=False,
):
        x_p_contrib=x_p_contrib_pos-x_p_contrib_neg
        rt_contrib=rt_contrib_pos-rt_contrib_neg
        u_contrib=u_contrib_pos-u_contrib_neg
        rt_u_inter=rt_u_inter_pos-rt_u_inter_neg
        
        fig4=plt.figure(figsize=(15,5))
        axs=fig4.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]})
        
        # contributions timeline
        ax=axs[0]

        # quadratic trends
            
        if abs(np.mean(x_p_contrib_pos))>0.01:
            z = np.polyfit(times,  x_p_contrib_pos, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ X_p positive",
                color="blue",
            ) 
            ax.fill_between(
                    times,
                    x_p_contrib_pos, 
                    p(times),
                    color="blue",
                    alpha=0.2                
                    )         
        if abs(np.mean(x_p_contrib_neg))>0.01:
            z = np.polyfit(times,  x_p_contrib_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ X_p negative",
                color="blue",
                linestyle="dashed",            
            ) 
            ax.fill_between(
                    times,
                    x_p_contrib_neg, 
                    p(times),
                    color="blue",
                    #linewidth=0.1,
                    alpha=0.2                
                    )          
        if abs(np.mean(u_contrib_pos))>0.01:        
            z = np.polyfit(times,  u_contrib_pos, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ NPP positive",
                color="green",
            ) 
            ax.fill_between(
                    times,
                    u_contrib_pos, 
                    p(times),
                    color="green",
                    alpha=0.2                
                    )
        if abs(np.mean(u_contrib_neg))>0.01:        
            z = np.polyfit(times,  u_contrib_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ NPP negative",
                color="green",
                linestyle="dashed",            
            ) 
            ax.fill_between(
                    times,
                    u_contrib_neg, 
                    p(times),
                    color="green",
                    alpha=0.2                
                    )           
        if abs(np.mean(rt_contrib_pos))>0.01:                
            z = np.polyfit(times,  rt_contrib_pos, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ RT positive",
                color="darkorange",
            )              
            ax.fill_between(
                    times,
                    rt_contrib_pos, 
                    p(times),
                    color="darkorange",
                    alpha=0.2                
                    )   
        if abs(np.mean(rt_contrib_neg))>0.01:        
            z = np.polyfit(times,  rt_contrib_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ RT negative",
                color="darkorange",
                linestyle="dashed",
            )              
            ax.fill_between(
                    times,
                    rt_contrib_neg, 
                    p(times),
                    color="darkorange",
                    alpha=0.2                
                    )      
        if abs(np.mean(rt_u_inter_pos))>0.01:          
            z = np.polyfit(times,  rt_u_inter_pos, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$NPP*$\Delta$RT positive",
                color="grey",
            )         
            ax.fill_between(
                    times,
                    rt_u_inter_pos, 
                    p(times),
                    color="grey",
                    alpha=0.2                
                    )  
        if abs(np.mean(rt_u_inter_neg))>0.01:         
            z = np.polyfit(times,  rt_u_inter_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$NPP*$\Delta$RT negative",
                color="grey",
                linestyle="dashed",            
            )         
            ax.fill_between(
                    times,
                    rt_u_inter_neg, 
                    p(times),
                    color="grey",
                    alpha=0.2                
                    )                         
        if abs(np.mean(delta_x_pos+delta_x_neg))>0.01:        
            z = np.polyfit(times,  delta_x_pos+delta_x_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="Net $\Delta$ X",
                color="black",
            )         
            ax.fill_between(
                    times,
                    delta_x_pos+delta_x_neg, 
                    p(times),
                    color="black",
                    alpha=0.2                
                    )                                   
        ax.set_title('Contributions over time')
        ax.grid()                 
                
        # bar charts       
        
        positive_delta_X=sum(delta_x_pos)/len(delta_x_pos)
        negative_delta_X=sum(delta_x_neg)/len(delta_x_neg)
        positive_x_p_contrib=sum(x_p_contrib_pos)/len(x_p_contrib_pos)
        negative_x_p_contrib=sum(x_p_contrib_neg)/len(x_p_contrib_neg)        
        positive_rt_contrib=sum(rt_contrib_pos)/len(rt_contrib_pos)
        negative_rt_contrib=sum(rt_contrib_neg)/len(rt_contrib_neg)
        positive_u_contrib=sum(u_contrib_pos)/len(u_contrib_pos)
        negative_u_contrib=sum(u_contrib_neg)/len(u_contrib_neg)      
        positive_rt_u_inter=sum(rt_u_inter_pos)/len(rt_u_inter_pos)
        negative_rt_u_inter=sum(rt_u_inter_neg)/len(rt_u_inter_neg)        
        
        ax0=axs[1]  
        ax0.set_title('Temporal average of contributions')       
        ax0.axhline(0, color='black', ls='dashed')
        
        if abs(np.mean(positive_delta_X+negative_delta_X))>0.01:
            ax0.bar ('Net $\Delta$ X', positive_delta_X+negative_delta_X, width=0.4, color="black", label='Net $\Delta$ X')        
        ax0.bar ('$\Delta$ RT', positive_rt_contrib, color="darkorange", label='$\Delta$ RT')
        ax0.bar ('$\Delta$ RT', negative_rt_contrib, color="darkorange")        
        ax0.bar ('$\Delta$ NPP', positive_u_contrib, color="green", label='$\Delta$ NPP')
        ax0.bar ('$\Delta$ NPP', negative_u_contrib, color="green")
        ax0.bar ('$\Delta$npp*$\Delta$rt', positive_rt_u_inter, color="lightgrey", label='$\Delta$npp*$\Delta$rt')
        ax0.bar ('$\Delta$npp*$\Delta$rt', negative_rt_u_inter, color="lightgrey")        
        ax0.bar ('$\Delta$ X_p', positive_x_p_contrib, color="blue", label='$\Delta$ X_p')
        ax0.bar ('$\Delta$ X_p', negative_x_p_contrib, color="blue")   
        ax0.legend(bbox_to_anchor =(1,1))  
        ax0.grid() 

        abs_total=x_p_contrib+rt_contrib+u_contrib+rt_u_inter
        percent_x_p=x_p_contrib/abs_total*100
        percent_rt=rt_contrib/abs_total*100
        percent_u=u_contrib/abs_total*100
        percent_inter=rt_u_inter/abs_total*100
        
        ######### Percents
        if percent==True:   
            fig3=plt.figure(figsize=(15,5))
            axs=fig3.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]})         
            # if delta==False: 
            ax=axs[0]      
            
            # % timeline
       
            # quadratic trends
            z = np.polyfit(times,  percent_x_p, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ X_p",
                color="blue",
            ) 
            ax.fill_between(
                    times,
                    percent_x_p, 
                    p(times),
                    color="blue",
                    alpha=0.2                
                    )  
                    
            z = np.polyfit(times,  percent_u, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ NPP",
                color="green",
            ) 
            ax.fill_between(
                    times,
                    percent_u, 
                    p(times),
                    color="green",
                    alpha=0.2                
                    )  

            z = np.polyfit(times,  percent_rt, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ RT",
                color="darkorange",
            )         
            ax.fill_between(
                    times,
                    percent_rt, 
                    p(times),
                    color="darkorange",
                    alpha=0.2                
                    )  

            z = np.polyfit(times,  percent_inter, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$NPP*$\Delta$RT",
                color="grey",
            ) 
            ax.fill_between(
                    times,
                    percent_inter, 
                    p(times),
                    color="grey",
                    alpha=0.2                
                    )  
            
            ax.set_title('% Contributions over time')
            ax.set_ylabel('%')
            ax.grid()
            
            # pie charts
            ax1=axs[1]
            ax1.set_title('Temporal average % contributions')
              
            if np.mean(percent_inter) > 0.001:
                labels = '$\Delta$ RT', '$\Delta$ NPP', '$\Delta$NPP*$\Delta$RT', '$\Delta$ X_p'
                sizes = [np.mean(percent_rt), np.mean(percent_u), 
                    np.mean(percent_inter), np.mean(percent_x_p)]
                ax1.pie(sizes, autopct='%1.1f%%', 
                    startangle=90, counterclock=False, colors=("darkorange", "green", "lightgrey", "blue"))
            else:
                labels = '$\Delta$ RT','$\Delta$ NPP', '$\Delta$ X_p'
                sizes = [np.mean(percent_rt), np.mean(percent_u), np.mean(percent_x_p)]
                ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                    startangle=90, counterclock=False, colors=("darkorange", "green", "blue"))        
            ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
            ax1.legend(labels, bbox_to_anchor =(1,1))
        
        plt.show()


def plot_attribution_per_model (
    all_comp_dict,
    percent=False,
):
    models=list(all_comp_dict.keys())[:-2]        
    
    for m in models: 
    
        print ('\033[1m'+m) 
        print ('\033[1m'+'Deviation from the mean')
           
        x=all_comp_dict[m]["x"]
        x_c=all_comp_dict[m]["x_c"]
        x_p=all_comp_dict[m]["x_p"]
        u=all_comp_dict[m]["u"]
        rt=all_comp_dict[m]["rt"]   
        
        delta_x=x-all_comp_dict["Mean"]["x"]
        delta_x_c=x_c-all_comp_dict["Mean"]["x_c"]
        delta_x_p=x_p-all_comp_dict["Mean"]["x_p"]
        delta_u=u-all_comp_dict["Mean"]["u"]
        delta_rt=rt-all_comp_dict["Mean"]["rt"] 
        
        # attribution of delta X to delta X_c and delta X_p
        x_c_contrib=delta_x_c
        x_p_contrib=-delta_x_p
         
        # attribution of delta X_c to delta u and delta RT
        rt_contrib=delta_rt*(u-delta_u/2)
        u_contrib=delta_u*(rt-delta_rt/2) 
        rt_u_inter=delta_x_c-rt_contrib-u_contrib           

        # separation of positive and negative contributions
          
        pos_delta_x=delta_x.copy(); pos_delta_x[pos_delta_x<0]=0
        neg_delta_x=delta_x.copy(); neg_delta_x[neg_delta_x>0]=0
 
        pos_cont_rt=rt_contrib.copy(); pos_cont_rt[pos_cont_rt<0]=0
        neg_cont_rt=rt_contrib.copy(); neg_cont_rt[neg_cont_rt>0]=0

        pos_cont_u=u_contrib.copy(); pos_cont_u[pos_cont_u<0]=0
        neg_cont_u=u_contrib.copy(); neg_cont_u[neg_cont_u>0]=0

        pos_cont_rt_u_inter=rt_u_inter.copy(); pos_cont_rt_u_inter[pos_cont_rt_u_inter<0]=0
        neg_cont_rt_u_inter=rt_u_inter.copy(); neg_cont_rt_u_inter[neg_cont_rt_u_inter>0]=0

        pos_cont_x_p=x_p_contrib.copy(); pos_cont_x_p[pos_cont_x_p<0]=0
        neg_cont_x_p=x_p_contrib.copy(); neg_cont_x_p[neg_cont_x_p>0]=0
        
        plot_attribution (
            times=all_comp_dict["Times"],
            delta_x_pos=pos_delta_x,
            delta_x_neg=neg_delta_x,
            rt_contrib_pos=pos_cont_rt,
            rt_contrib_neg=neg_cont_rt,
            u_contrib_pos=pos_cont_u,
            u_contrib_neg=neg_cont_u,
            rt_u_inter_pos=pos_cont_rt_u_inter,
            rt_u_inter_neg=neg_cont_rt_u_inter,
            x_p_contrib_pos=pos_cont_x_p,
            x_p_contrib_neg=neg_cont_x_p,
        )
        
        print ('\033[1m'+'Deviation from other models')
        print("\033[0;0m")
        for m2 in models:       
            if m2!=m:
                print ('Deviation from '+m2)
                delta_x=x-all_comp_dict[m2]["x"]
                delta_x_c=x_c-all_comp_dict[m2]["x_c"]
                delta_x_p=x_p-all_comp_dict[m2]["x_p"]
                delta_u=u-all_comp_dict[m2]["u"]
                delta_rt=rt-all_comp_dict[m2]["rt"]
                
                # attribution of delta X to delta X_c and delta X_p
                x_c_contrib=delta_x_c
                x_p_contrib=-delta_x_p
                 
                # attribution of delta X_c to delta u and delta RT
                rt_contrib=delta_rt*(u-delta_u/2)
                u_contrib=delta_u*(rt-delta_rt/2) 
                rt_u_inter=delta_x_c-rt_contrib-u_contrib           

                # separation of positive and negative contributions                  
                pos_delta_x=delta_x.copy(); pos_delta_x[pos_delta_x<0]=0
                neg_delta_x=delta_x.copy(); neg_delta_x[neg_delta_x>0]=0
         
                pos_cont_rt=rt_contrib.copy(); pos_cont_rt[pos_cont_rt<0]=0
                neg_cont_rt=rt_contrib.copy(); neg_cont_rt[neg_cont_rt>0]=0

                pos_cont_u=u_contrib.copy(); pos_cont_u[pos_cont_u<0]=0
                neg_cont_u=u_contrib.copy(); neg_cont_u[neg_cont_u>0]=0

                pos_cont_rt_u_inter=rt_u_inter.copy(); pos_cont_rt_u_inter[pos_cont_rt_u_inter<0]=0
                neg_cont_rt_u_inter=rt_u_inter.copy(); neg_cont_rt_u_inter[neg_cont_rt_u_inter>0]=0

                pos_cont_x_p=x_p_contrib.copy(); pos_cont_x_p[pos_cont_x_p<0]=0
                neg_cont_x_p=x_p_contrib.copy(); neg_cont_x_p[neg_cont_x_p>0]=0
                
                plot_attribution (
                    times=all_comp_dict["Times"],
                    delta_x_pos=pos_delta_x,
                    delta_x_neg=neg_delta_x,
                    rt_contrib_pos=pos_cont_rt,
                    rt_contrib_neg=neg_cont_rt,
                    u_contrib_pos=pos_cont_u,
                    u_contrib_neg=neg_cont_u,
                    rt_u_inter_pos=pos_cont_rt_u_inter,
                    rt_u_inter_neg=neg_cont_rt_u_inter,
                    x_p_contrib_pos=pos_cont_x_p,
                    x_p_contrib_neg=neg_cont_x_p,
                    percent=percent,
                )

def plot_attribution_sum (
    all_comp_dict,
    percent=False,
    part=1,
):
    models=list(all_comp_dict.keys())[:-2]  

    # if we do not want the whole interval but look at a smaller part to observe the dynamics
    start_min=0
    stop_max=len(all_comp_dict["Times"])
    
    # if we do not want the whole interval but look at a smaller part to observe the dynamics
    if part < 0:
        start, stop = int(stop_max - (stop_max - start_min) * abs(part)), stop_max
    else:
        start, stop = start_min, int(start_min + (stop_max - start_min) * part)
    
    times=all_comp_dict["Times"][start:stop]   
    
    sum_pos_diff_x=0
    sum_pos_cont_rt=0
    sum_pos_cont_u=0
    sum_pos_cont_rt_u_inter=0
    sum_pos_cont_x_p=0
    
    sum_neg_diff_x=0
    sum_neg_cont_rt=0
    sum_neg_cont_u=0
    sum_neg_cont_rt_u_inter=0
    sum_neg_cont_x_p=0  
    
    print ('\033[1m'+'Attribution of summed deviations from the mean for all models ' +
        'to the differences in traceable components')     
    for m in models: 
           
        x=all_comp_dict[m]["x"][start:stop]
        x_c=all_comp_dict[m]["x_c"][start:stop]
        x_p=all_comp_dict[m]["x_p"][start:stop]
        u=all_comp_dict[m]["u"][start:stop]
        rt=all_comp_dict[m]["rt"][start:stop]
        x_mean=all_comp_dict["Mean"]["x"][start:stop]
        x_c_mean=all_comp_dict["Mean"]["x_c"][start:stop]
        x_p_mean=all_comp_dict["Mean"]["x_p"][start:stop]
        u_mean=all_comp_dict["Mean"]["u"][start:stop]
        rt_mean=all_comp_dict["Mean"]["rt"][start:stop]           
            
        delta_x=x-x_mean
        delta_x_c=x_c-x_c_mean
        delta_x_p=x_p-x_p_mean
        delta_u=u-u_mean
        delta_rt=rt-rt_mean
        
        # attribution of delta X to delta X_c and delta X_p
        x_c_contrib=delta_x_c
        x_p_contrib=-delta_x_p
         
        # attribution of delta X_c to delta u and delta RT
        rt_contrib=delta_rt*(u-delta_u/2)
        u_contrib=delta_u*(rt-delta_rt/2) 
        rt_u_inter=delta_x_c-rt_contrib-u_contrib         

        # summation of positive and negative contributions separately         
        
        pos_delta_x=delta_x.copy(); pos_delta_x[pos_delta_x<0]=0
        neg_delta_x=delta_x.copy(); neg_delta_x[neg_delta_x>0]=0
        sum_pos_diff_x+=pos_delta_x
        sum_neg_diff_x+=neg_delta_x       

        pos_cont_rt=rt_contrib.copy(); pos_cont_rt[pos_cont_rt<0]=0
        neg_cont_rt=rt_contrib.copy(); neg_cont_rt[neg_cont_rt>0]=0
        sum_pos_cont_rt+=pos_cont_rt
        sum_neg_cont_rt+=neg_cont_rt

        pos_cont_u=u_contrib.copy(); pos_cont_u[pos_cont_u<0]=0
        neg_cont_u=u_contrib.copy(); neg_cont_u[neg_cont_u>0]=0
        sum_pos_cont_u+=pos_cont_u
        sum_neg_cont_u+=neg_cont_u

        pos_cont_rt_u_inter=rt_u_inter.copy(); pos_cont_rt_u_inter[pos_cont_rt_u_inter<0]=0
        neg_cont_rt_u_inter=rt_u_inter.copy(); neg_cont_rt_u_inter[neg_cont_rt_u_inter>0]=0
        sum_pos_cont_rt_u_inter+=pos_cont_rt_u_inter
        sum_neg_cont_rt_u_inter+=neg_cont_rt_u_inter

        pos_cont_x_p=x_p_contrib.copy(); pos_cont_x_p[pos_cont_x_p<0]=0
        neg_cont_x_p=x_p_contrib.copy(); neg_cont_x_p[neg_cont_x_p>0]=0
        sum_pos_cont_x_p+=pos_cont_x_p
        sum_neg_cont_x_p+=neg_cont_x_p 
        
    plot_attribution (
        times=times,
        delta_x_pos=sum_pos_diff_x,
        delta_x_neg=sum_neg_diff_x,
        rt_contrib_pos=sum_pos_cont_rt,
        rt_contrib_neg=sum_neg_cont_rt,
        u_contrib_pos=sum_pos_cont_u,
        u_contrib_neg=sum_neg_cont_u,
        rt_u_inter_pos=sum_pos_cont_rt_u_inter,
        rt_u_inter_neg=sum_neg_cont_rt_u_inter,
        x_p_contrib_pos=sum_pos_cont_x_p,
        x_p_contrib_neg=sum_neg_cont_x_p,
        percent=percent,
    )    

def plot_single_trend(var, times, polynom_order,title):
    fig=plt.figure(figsize=(15,5))
    ax=fig.subplots()
    z = np.polyfit(times, var, polynom_order)
    p = np.poly1d(z)  
    ax.fill_between(
            times,
            var, 
            p(times),
            color="black",
            label="original value",
            alpha=0.2                
            ) 
    ax.plot(        
        times, 
        p(times),
        label="polynomial trend",
        color="black",
        ) 
    ax.grid()
    ax.set_title(title)
    ax.legend()
    plt.show()
    
def get_components_from_output(
    model_names,  # dictionary (folder name : model name)
    test_arg_list,  # a list of test_args from all models involved
    delta_t_val,  # model time step
    model_cols,  # dictionary (folder name :color)
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging,  # number of iterator steps over which to average results. 1 for no averaging
    overlap=True,  # compute overlapping timeframe or plot whole duration for all models
):
    if (part < -1) | (part > 1) | (part == 0):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between -1 and 1 excluding 0"
        )
    model_folders = [(k) for k in model_names]
    k = 0
    sum_cVeg = np.array(0)
    sum_cSoil = np.array(0)
    sum_npp = np.array(0)
    sum_rh = np.array(0)
    sum_x = np.array(0)
    sum_rt = np.array(0)
    sum_x_c = np.array(0)
    sum_x_p = np.array(0)
    
    all_components = list ()
    for mf in model_folders:
        print ("Getting traceable components for "+mf+"...")
        if overlap == True:
            start_min, stop_max = min_max_index(
                test_arg_list[k],
                delta_t_val,
                *t_min_tmax_overlap(test_arg_list, delta_t_val)
            )
        else:
            start_min, stop_max = min_max_index(
                test_arg_list[k],
                delta_t_val,
                *t_min_tmax_full(test_arg_list, delta_t_val)
            )
        # if we do not want the whole interval but look at a smaller part to observe the dynamics
        if part < 0:
            start, stop = int(stop_max - (stop_max - start_min) * abs(part)), stop_max
        else:
            start, stop = start_min, int(start_min + (stop_max - start_min) * part)
        times = (
            times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            / days_per_year()
        )
        # if len(test_arg_list[k].svs.cVeg)==320:
            # start_pool=160
            # stop_pool=320
        # elif len(test_arg_list[k].svs.cVeg)==3840:
            # start_pool=1920
            # stop_pool=3839
        # else:    
            # start_pool=0
            # stop_pool=1919        
        # if len(test_arg_list[k].svs.rh)==320:
            # start_flux=160
            # stop_flux=320
        # elif len(test_arg_list[k].svs.rh)==3840:
            # start_flux=1920
            # stop_flux=3839            
        # else:
            # start_flux=0
            # stop_flux=1919  
        # times = np.array(range(1860,2020))
        # times +=1
        # harmonising model outputs
        start_flux=start
        stop_flux=stop
        print(start)
        print(stop)
        print(len(test_arg_list[k].svs.cVeg)<1000)
        if len(test_arg_list[k].svs.cVeg)<1000:
            start_pool=start//12
            stop_pool=stop//12+1
        else:
            start_pool=start
            stop_pool=stop        
        cVeg=test_arg_list[k].svs.cVeg[start_pool:stop_pool]
        if cVeg.shape[0]>500: cVeg=avg_timeline(cVeg, averaging)
        if "cLitter" in test_arg_list[k].svs._fields: 
            cSoil=(test_arg_list[k].svs.cLitter[start_pool:stop_pool]+
                test_arg_list[k].svs.cSoil[start_pool:stop_pool])
        else: cSoil=test_arg_list[k].svs.cSoil[start_pool:stop_pool]
        if cSoil.shape[0]>500: cSoil=avg_timeline(cSoil, averaging)
        npp=test_arg_list[k].dvs.npp[start_flux:stop_flux]
        if npp.shape[0]>500: npp=avg_timeline(npp, averaging)
        rh=test_arg_list[k].svs.rh[start_flux:stop_flux]
        if rh.shape[0]>500: rh=avg_timeline(rh, averaging)
        # print(times)
        # print(cVeg.shape)
        # print(cSoil.shape)
        # print(npp.shape)
        # print(rh.shape)
        # print(start_pool, stop_pool)
        # print(start_flux, stop_flux)        
        # traditional traceable components
        x=cVeg+cSoil
        rt=x/rh
        x_c=rt*npp
        x_p=x_c-x
        # pool-wise
        #rt_soil=rh/cSoil
        comp_dict = {
            "cVeg": cVeg,
            "cSoil": cSoil,
            "u":npp,
            "rh":rh,
            "x":x,
            "rt":rt,
            "x_c":x_c,
            "x_p":x_p
            #"rt_soil":rt_soil,           
            }                
        all_components.append(comp_dict)
        if sum_cVeg.all()==0:
            sum_cVeg = np.append(sum_cVeg,cVeg)[1:]
            sum_cSoil = np.append(sum_cSoil,cSoil)[1:]
            sum_npp = np.append(sum_npp,npp)[1:]
            sum_rh = np.append(sum_rh,rh)[1:]
            sum_x = np.append(sum_x,x)[1:]
            sum_rt = np.append(sum_rt,rt)[1:]
            sum_x_c = np.append(sum_x_c,x_c)[1:]
            sum_x_p = np.append(sum_x_p,x_p)[1:]
        else:
            sum_cVeg = sum_cVeg + cVeg
            sum_cSoil = sum_cSoil + cSoil
            sum_npp = sum_npp + npp
            sum_rh = sum_rh + rh
            sum_x = sum_x + x
            sum_rt = sum_rt + rt
            sum_x_c = sum_x_c + x_c
            sum_x_p = sum_x_p + x_p
        k += 1
    ave_dict = {
            "cVeg": sum_cVeg / len(model_folders),
            "cSoil": sum_cSoil / len(model_folders), 
            "u": sum_npp / len(model_folders), 
            "rh": sum_rh / len(model_folders), 
            "x": sum_x / len(model_folders), 
            "rt": sum_rt / len(model_folders),
            "x_c": sum_x_c / len(model_folders),
            "x_p": sum_x_p / len(model_folders),           
            }    
    all_components.append(ave_dict) 
    
    #times_avg = times
    times_avg = avg_timeline(times, averaging)
    all_components.append(times_avg)  
    mods=list(model_names.values())
    mods.append("Mean")
    mods.append("Times")
    all_comp_dict = {mods[i]: all_components[i] for i in range(len(mods))}             
        
    return all_comp_dict
        