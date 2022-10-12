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

#fixme: possibly obsolete - see combined_masks_2 using nearest neighbor resampling
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
        ds.close()

def project_2(source: CoordMask, target: CoordMask):

    s_mask = source.index_mask  # .astype(np.float64,casting='safe')
    #print(s_mask.mean())
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
    #print(float_grid.mean())
    projected_mask = float_grid > 0.5
    # from IPython import embed; embed()
    return CoordMask(index_mask=np.logical_or(projected_mask, t_mask), tr=t_tr)

def resample_grid (source_coord_mask, target_coord_mask, var, 
            method="nearest", radius_of_influence=500000, neighbours=10):
    # NEED TO INSTALL PACKAGE: conda install -c conda-forge pyresample
    import pyresample
    from pyresample.bilinear import NumpyBilinearResampler    
    lon2d, lat2d = np.meshgrid(source_coord_mask.lons, source_coord_mask.lats)  
    lon2d_t, lat2d_t = np.meshgrid(target_coord_mask.lons, target_coord_mask.lats)     
    
    lats_source=source_coord_mask.lats
    lons_source=source_coord_mask.lons
    
    lats=target_coord_mask.lats
    lons=target_coord_mask.lons
    orig_def = pyresample.geometry.SwathDefinition(lons=lon2d, lats=lat2d)
    
    # orig_def = pyresample.AreaDefinition(area_id="world", description="CLASSIC mask", 
                          # proj_id="lat_lon", 
                          # projection='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                          # width=len(lons_source), height=len(lats_source), 
                          # area_extent=(min(lons_source),min(lats_source),max(lons_source),max(lats_source)),
                          # )      
    #orig_def = pyresample.geometry.GridDefinition(lons=lons_source, lats=lats_source)

    targ_def = pyresample.geometry.SwathDefinition(lons=lon2d_t, lats=lat2d_t)
    # targ_def = pyresample.AreaDefinition(area_id="world", description="global mask", 
                          # proj_id="lat_lon", 
                          # projection='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                          # width=len(lons), height=len(lats), 
                          # area_extent=(min(lons),min(lats),max(lons),max(lats)),
                          # )   
    # orig_con = pyresample.image.ImageContainerNearest(var, orig_def, radius_of_influence=5000)
    # target_con = orig_con.resample(targ_def)
    # result = target_con.image_data
                  
    if method=="nearest":
      target_var = pyresample.kd_tree.resample_nearest(orig_def, var, 
          targ_def, radius_of_influence=radius_of_influence, fill_value=None) 
          
    elif method=="idw":
        wf = lambda r: 1/r**2
        target_var = pyresample.kd_tree.resample_custom(orig_def, var, 
                              targ_def, radius_of_influence=radius_of_influence, 
                              neighbours=neighbours,
                              weight_funcs=wf, fill_value=None)
    elif method=="gauss":
        target_var = pyresample.kd_tree.resample_gauss(orig_def, var, 
                               targ_def, radius_of_influence=radius_of_influence, 
                               neighbours=neighbours,
                               sigmas=250000, fill_value=None)                               
    else:
              raise Exception(
           "Invalid resample method. Valid options are: 'nearest', 'idw' and 'gauss'"
       )

    return(
        CoordMask (
          index_mask=target_var,
          tr=target_coord_mask.tr
        )
    ) 
    
def resample_nc(
            model_names, # dictionary e.g. "ab_classic":"CLASSIC"
            experiment_names, # e.g. ['S2', 'S3']
            target_mask,
            method="nearest",
            radius_of_influence=500000,
            ):
    for experiment in experiment_names:
        print('\033[1m'+'. . . Resampling data for '+experiment+' experiment . . .')
        model_folders=[(m) for m in model_names] 
        m_names=list(model_names.values())  
        g_mask = target_mask.index_mask
        k=0 # model counter
        for mf in model_folders:
            print('\033[1m'+m_names[k])
            print('\033[0m')
            experiment_name=m_names[k]+"_"+experiment+"_"
            conf_dict = confDict(mf)
            dataPath=Path(conf_dict["dataPath"])
            model_mask=msh(mf).spatial_mask(dataPath=Path(conf_dict["dataPath"]))    
            for vn in msh(mf).data_str._fields:
                print("Resampling "+vn)        
                file_path = dataPath.joinpath(msh(mf).nc_file_name(vn, experiment_name=experiment_name))
                ds = nc.Dataset(str(file_path))
                var=ds.variables[vn][:, :, :].data             
                # preparing a narray to store results
                zero_array=np.zeros((var.shape[0],g_mask.shape[0],g_mask.shape[1]))
                gm=zero_array.copy()
                for i in range(gm.shape[0]):
                    gm[i,:,:]=g_mask
                var_array=zero_array
                # procesing all time steps
                for i in range(var.shape[0]):
                    var_current=var[i,:,:]
                    # initial masking
                    var_masked = np.ma.array(var_current, mask = model_mask.index_mask)    
                    # resampling
                    var_resampled=resample_grid (
                        source_coord_mask=model_mask, 
                        target_coord_mask=target_mask, 
                        var=var_masked, 
                        method=method,
                        radius_of_influence=radius_of_influence,               
                    )
                    var_array[i,:,:]=var_resampled.index_mask
                    if i//100 == i/100:
                        print(str(i+1)+" out of "+str(var.shape[0])+" time steps completed")
                # final masking
                var_final = np.ma.array(var_array,mask = gm)
                # creating and writing a new NetCDF file 
                s = g_mask.shape
                n_lats, n_lons = s
                new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_res.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)
                source_times=ds.variables["time"][:].data        
                time = ds_new.createDimension("time", size=len(source_times))               
                # creating variables                         
                nc_var = ds_new.createVariable(vn, "float32", ["time", "lat", "lon"])
                nc_var[:, :, :] = var_final
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(target_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(target_mask.tr.i2lon, range(n_lons)))               
                times = ds_new.createVariable ("time", "float32", ["time"])
                times[:] = source_times
                # closing NetCDF files
                ds.close()        
                ds_new.close() 
            k+=1 # model counter
        print("Done!")
        
def average_and_resample_nc(
            model_names, # dictionary e.g. "ab_classic":"CLASSIC"
            experiment_names, # e.g. ['S2', 'S3']
            target_mask,
            method="nearest",
            radius_of_influence=500000,
            ):
    for experiment in experiment_names:
        print('\033[1m'+'Resampling data for '+experiment+' experiment...')
        k=0 # model counter
        model_folders=[(m) for m in model_names] 
        m_names=list(model_names.values())  
        for mf in model_folders:
            print('\033[1m'+m_names[k])
            print('\033[0m')
            experiment_name=m_names[k]+"_"+experiment+"_"
            conf_dict = confDict(mf)
            dataPath=Path(conf_dict["dataPath"])
            model_mask=msh(mf).spatial_mask(dataPath=Path(conf_dict["dataPath"]))       
            for vn in msh(mf).data_str._fields:      
                if vn=="npp_nlim": file_path = dataPath.joinpath(msh(mf).nc_file_name("npp", experiment_name=experiment_name))
                else: file_path = dataPath.joinpath(msh(mf).nc_file_name(vn, experiment_name=experiment_name))
                print(file_path)
                ds = nc.Dataset(str(file_path))
                var=ds.variables[vn][:, :, :].data
                
                # finding start value (differs for different models)
                if var.shape[0]>1000: # monthly data
                    start=var.shape[0]-1440 # smallest time range - SDGVM fluxes 
                else: # yearly data
                    start=var.shape[0]-120 # smallest time range - SDGVM fluxes
                # temporal average + difference
                if vn in ['gpp', 'npp', 'ra', 'rh', 'npp_nlim']: # fluxes per sec to total per 120 years
                    var_avg=np.ma.mean(var[start:,:,:], axis=0)*86400*365*120                  
                else:  # for pools we compute mean              
                    var_avg=np.ma.mean(var[start:,:,:], axis=0)    
                #var_avg=np.ma.mean(var[start:,:,:], axis=0)                
                var_diff=var[-1,:,:]-var[start,:,:]
                # masking
                var_masked_avg = np.ma.array(var_avg, mask = model_mask.index_mask)
                var_masked_diff = np.ma.array(var_diff, mask = model_mask.index_mask) 
                # resampling to target grid
                var_resampled_avg=resample_grid (
                    source_coord_mask=model_mask, 
                    target_coord_mask=target_mask, 
                    var=var_masked_avg, 
                    method=method,
                    radius_of_influence=radius_of_influence,
                )
                var_resampled_diff=resample_grid (
                    source_coord_mask=model_mask, 
                    target_coord_mask=target_mask, 
                    var=var_masked_diff, 
                    method=method,
                    radius_of_influence=radius_of_influence,
                )
                # final masking
                var_final_avg = np.ma.array(var_resampled_avg.index_mask, 
                                            mask = target_mask.index_mask)
                var_final_diff = np.ma.array(var_resampled_diff.index_mask, 
                                            mask = target_mask.index_mask)            
                # creating and writing a new NetCDF file               
                s = target_mask.index_mask.shape
                n_lats, n_lons = s
                
                # if vn in ['gpp', 'npp', 'ra', 'rh', 'npp_nlim']: # for fluxes               
                    # new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_sum_res.nc")
                # else: # for pools 
                    # new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_ave_res.nc")
                new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_ave_res.nc")                    
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)
                # creating variables            
                avg = ds_new.createVariable(vn, "float32", ["lat", "lon"])     
                avg[:, :] = var_final_avg                
                diff = ds_new.createVariable(str(vn)+"_diff", "float32", ["lat", "lon"])                 
                diff[:, :] = var_final_diff
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(target_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(target_mask.tr.i2lon, range(n_lons)))        
                # closing NetCDF files
                ds.close()        
                ds_new.close()
            k+=1 # model counter
        print("Done!")
    
def combine_masks_2(coord_masks: List["CoordMask"]):
   
    def combine (source, target):
        resampled_mask=resample_grid(
            source_coord_mask=source,
            target_coord_mask=target,
            var=source.index_mask,
            method="nearest",
            radius_of_influence=500000, 
            neighbours=10
            )
        combined_mask=CoordMask(
            index_mask=np.logical_or(resampled_mask.index_mask, target.index_mask), 
            tr=target.tr)
        return combined_mask
    target_mask = coord_masks[-1] # we use last mask in the list as template grid
    for i in range(len(coord_masks)-1):
        target_mask = combine (coord_masks[i], target_mask)
    return target_mask    

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
    fig = plt.figure(figsize=(20, 10))
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
                label="Ensemble mean - X",
                color="black",
                linewidth=3,
        )
        ax.plot(
                all_comp_dict["Times"],
                all_comp_dict["Mean"]["x_c"]
                * 148940000
                * 1000000
                * 0.000000000001,  # convert to global C in Gt
                label="Ensemble mean - X_c",
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
                    linewidth=1,
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
                    linewidth=1,
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
                label="Ensemble mean",
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
                    linewidth=1.2, 
                )
        variance = diff_sqr / (len(models)-1)
        st_dev=np.sqrt(variance)
      
        ax.fill_between(
                all_comp_dict["Times"],
                vals_mean-st_dev*2, 
                vals_mean+st_dev*2,
                label="\u00B12$\sigma$ interval",
                color="grey",
                alpha=0.2                
                )                          
    ax.grid()
    if delta:
        ax.set_title("$\Delta$ "+str(comp_name))
    else:
        ax.set_title(comp_name)
    #ax.set_ylabel("Gt C")
    ax.legend(bbox_to_anchor =(1.2, 1))   
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



    
def get_vars_all_list(model_folders, experiment_names):
    vars_all_list = []
    i=0
    for mf in model_folders:
        print('\033[1m'+mf)
        print ('\033[0m')
        current_var_list = msh(mf).get_global_mean_vars_all(experiment_name=experiment_names[i])
        vars_all_list.append(current_var_list)
        i+=1
    print("Done!")
    return vars_all_list        
    
def get_components_from_output(
    model_names,  # dictionary (folder name : model name)
    vars_all_list,  # a list of test_args from all models involved
    delta_t_val,  # model time step
    #model_cols,  # dictionary (folder name :color)
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging=12,  # number of iterator steps over which to average results. 12 - avarage monthly to yearly 
    overlap=True,  # compute overlapping timeframe or plot whole duration for all models
):

    def times_in_days_aD(
        vars_all,  # a tuple defined in model-specific test_helpers.py
        delta_t_val,  # iterator time step
        start_date,
    ):
        n_months = len(vars_all.gpp)*12
        #n_days = month_2_day_index([n_months], start_date)[0]
        n_iter = n_months #int(n_days / delta_t_val)
        #print("n_months: "+ str(n_months))
        #print("n_days: "+ str(n_days))
        #print("n_iter: "+ str(n_iter))
        return np.array(
            tuple(
                (
                    days_since_AD(i, delta_t_val, start_date)
                    for i in np.arange(n_iter)  # days_after_sim_start
                )
            )
        )
    # function to determine overlapping time frames for models simulations
    def t_min_tmax_overlap(
        vars_all_list,  # a list of test_args from all models involved
        delta_t_val,  # iterator time step
        start_date,
    ):
        delta_t_val=1
        td = {
            i: times_in_days_aD(vars_all_list[i], delta_t_val, start_date)
            for i in range(len(vars_all_list))
        }
        t_min = max([t.min() for t in td.values()])
        t_max = min([t.max() for t in td.values()])
        return (t_min, t_max)
    
    def min_max_index(
    vars_all,  # a tuple defined in model-specific test_helpers.py
    delta_t_val,  # iterator time step
    t_min,  # output of t_min_tmax_overlap or t_min_tmax_full
    t_max,  # output of t_min_tmax_overlap or t_min_tmax_full
    start_date,
    ):
        ts = times_in_days_aD(vars_all, delta_t_val, start_date)

        def count(acc, i):
            min_i, max_i = acc
            t = ts[i]
            min_i = min_i + 1 if t < t_min else min_i
            max_i = max_i + 1 if t < t_max else max_i
            return (min_i, max_i)

        return reduce(count, range(len(ts)), (0, 0))

    if (part < -1) | (part > 1) | (part == 0):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between -1 and 1 excluding 0"
        )
    model_folders = [(k) for k in model_names]
    k = 0
    sum_cVeg = np.array(0)
    sum_cSoil = np.array(0)
    sum_gpp = np.array(0)    
    sum_npp = np.array(0)
    sum_ra = np.array(0)
    sum_rh = np.array(0)
    sum_x = np.array(0)
    sum_rt = np.array(0)
    sum_x_c = np.array(0)
    sum_x_p = np.array(0)
    sum_nep = np.array(0)    
    sum_rt_soil = np.array(0)
    sum_f_v2s = np.array(0)
    sum_f_v2s_2 = np.array(0)    
    sum_rt_soil = np.array(0)
    sum_rt_soil_2 = np.array(0)    
    sum_rt_veg = np.array(0)
    sum_rt_veg_2 = np.array(0)    
    sum_x_c_soil = np.array(0)
    sum_x_c_veg = np.array(0)
    sum_x_p_veg = np.array(0)
    sum_x_p_soil = np.array(0)
    sum_ra_rate = np.array(0)
    sum_rh_rate = np.array(0)
    sum_delta_cVeg = np.array(0)
    sum_delta_cSoil = np.array(0)
    sum_delta_dist = np.array(0)
    sum_dist = np.array(0)
    all_components = list ()
    
    def annual_delta (data_stream):
        delta=np.zeros_like(data_stream)
        for i in range(len(data_stream)-1):
            delta[i]=data_stream[i+1]-data_stream[i]  
        return delta

    min_len=1000
    for i in range (len(model_folders)):
        if len(vars_all_list[i].cVeg)<min_len: min_len=len(vars_all_list[i].cVeg)
    print(min_len)
    
    for mf in model_folders:
        print ("Getting traceable components for "+mf+"...")
        start_date=msh(mf).start_date()
        # start_min, stop_max = min_max_index(
            # vars_all_list[k],
            # delta_t_val,
            # *t_min_tmax_overlap(vars_all_list, delta_t_val, start_date),
            # start_date
        # )
        # # if we do not want the whole interval but look at a smaller part to observe the dynamics
        # if part < 0:
            # start, stop = int(stop_max - (stop_max - start_min) * abs(part)), stop_max
        # else:
            # start, stop = start_min, int(start_min + (stop_max - start_min) * part)
        # # times = (
            # # times_in_days_aD(vars_all_list[k], delta_t_val, start_date)[start:stop]
            # # / days_per_year()
        # # )
        
        
        
        stop=len(vars_all_list[k].cVeg)#-20
        start=len(vars_all_list[k].cVeg)-min_len
        
        times=np.array(range(start,stop))+1700
        
        #print(times)
        # harmonising model outputs
        start_flux=start
        stop_flux=stop
        #print(len(vars_all_list[k].cVeg)>1000)
        # if len(vars_all_list[k].cVeg)>1000:
            # start_pool=start//12
            # stop_pool=stop//12+1
        # else:
        start_pool=start
        stop_pool=stop       
           
        cVeg=vars_all_list[k].cVeg[start_pool:stop_pool]
        if cVeg.shape[0]>500: cVeg=avg_timeline(cVeg, averaging)
        cSoil=vars_all_list[k].cSoil[start_pool:stop_pool]
        if cSoil.shape[0]>500: cSoil=avg_timeline(cSoil, averaging)
        gpp=vars_all_list[k].gpp[start_flux:stop_flux]
        if gpp.shape[0]>500: gpp=avg_timeline(gpp, averaging)        
        npp=vars_all_list[k].npp[start_flux:stop_flux]
        if npp.shape[0]>500: npp=avg_timeline(npp, averaging)
        rh=vars_all_list[k].rh[start_flux:stop_flux]
        if rh.shape[0]>500: rh=avg_timeline(rh, averaging)
        ra=vars_all_list[k].ra[start_flux:stop_flux]
        if ra.shape[0]>500: ra=avg_timeline(ra, averaging)        
        # print(times.shape)
        # print(cVeg.shape)
        # print(cSoil.shape)
        # print(npp.shape)
        # print(rh.shape)
        # print(start_pool, stop_pool)
        # print(start_flux, stop_flux)        
        
        # correction for ISAM data issue
        if mf=="cj_isam": 
            # print (cVeg[-9])
            # print(cVeg[-8])
            # print(cVeg[-10])
            cVeg[-9]=np.mean((cVeg[-8], cVeg[-10]))
        
        # traditional traceable components 
        n_years=len(cSoil)
        #print(n_years)
        x=cVeg+cSoil
        nep=annual_delta(x)
        delta_cVeg=annual_delta(cVeg)
        delta_cSoil=annual_delta(cSoil)
        rt=x/(rh+ra)
        x_c=rt*gpp
        x_p=x_c-x
        # flux-based                
        dist=gpp-nep-ra-rh
        f_v2s=delta_cSoil+rh
        f_v2s_2=gpp-delta_cVeg-ra-dist
        rt_soil=cSoil/rh
        rt_soil_2=cSoil/(f_v2s-delta_cSoil)
        rt_veg=cVeg/(ra+dist+f_v2s)
        rt_veg_2=cVeg/(gpp-delta_cVeg)
        x_c_soil=f_v2s*rt_soil
        x_c_veg=gpp*rt_veg
        x_p_soil=x_c_soil-cSoil
        x_p_veg=x_c_veg-cVeg
        ra_rate=ra/cVeg
        rh_rate=rh/cSoil
        
        comp_dict = {
            "cVeg": cVeg,
            "cSoil": cSoil,
            "delta_cVeg": delta_cVeg,
            "delta_cSoil": delta_cSoil,
            "gpp":gpp,
            "npp":npp,
            "nep":nep,            
            "rh":rh,
            "ra":ra,
            "f_v2s":f_v2s,
            "f_v2s_2":f_v2s_2,
            "dist":dist,
            "x":x,
            "rt":rt,
            "x_c":x_c,
            "x_p":x_p,
            "rt_veg":rt_veg,
            "rt_veg_2":rt_veg_2,            
            "rt_soil":rt_soil,
            "rt_soil_2":rt_soil_2,
            "x_c_soil":x_c_soil,
            "x_c_veg":x_c_veg,
            #"x_c_veg2":x_c_veg2,
            "x_p_veg":x_p_veg,
            "x_p_soil":x_p_soil,
            #"x_p_veg2":x_p_veg2,
            "ra_rate":ra_rate,
            "rh_rate":rh_rate,
            }  
                                           
        all_components.append(comp_dict)
        if sum_cVeg.all()==0:
            sum_cVeg = np.append(sum_cVeg,cVeg)[1:]
            sum_cSoil = np.append(sum_cSoil,cSoil)[1:]
            sum_gpp = np.append(sum_gpp,gpp)[1:]
            sum_npp = np.append(sum_npp,npp)[1:]
            sum_rh = np.append(sum_rh,rh)[1:]
            sum_ra = np.append(sum_ra,ra)[1:]            
            sum_x = np.append(sum_x,x)[1:]
            sum_rt = np.append(sum_rt,rt)[1:]
            sum_x_c = np.append(sum_x_c,x_c)[1:]
            sum_x_p = np.append(sum_x_p,x_p)[1:]
            #sum_nep = np.append(sum_nep,nep)[1:]  
            sum_f_v2s = np.append(sum_f_v2s,f_v2s)[1:]
            sum_f_v2s_2 = np.append(sum_f_v2s_2,f_v2s_2)[1:]            
            sum_rt_soil = np.append(sum_rt_soil,rt_soil)[1:]
            sum_rt_soil_2 = np.append(sum_rt_soil_2,rt_soil_2)[1:]            
            sum_rt_veg = np.append(sum_rt_veg,rt_veg)[1:]
            sum_rt_veg_2 = np.append(sum_rt_veg_2,rt_veg_2)[1:]            
            sum_x_c_soil = np.append(sum_x_c_soil,x_c_soil)[1:]
            sum_x_c_veg = np.append(sum_x_c_veg,x_c_veg)[1:]
            sum_x_p_veg = np.append(sum_x_p_veg,x_p_veg)[1:]
            sum_x_p_soil = np.append(sum_x_p_soil,x_p_soil)[1:]            
            sum_ra_rate = np.append(sum_ra_rate,ra_rate)[1:]
            sum_rh_rate = np.append(sum_rh_rate,rh_rate)[1:] 
            sum_dist = np.append(sum_dist,dist)[1:] 
            sum_delta_cVeg = np.append(sum_delta_cVeg,delta_cVeg)[1:] 
            sum_delta_cSoil = np.append(sum_delta_cSoil,delta_cSoil)[1:]             
        else:
            sum_cVeg = sum_cVeg + cVeg
            sum_cSoil = sum_cSoil + cSoil
            sum_gpp = sum_gpp + gpp
            sum_npp = sum_npp + npp
            sum_rh = sum_rh + rh
            sum_ra = sum_ra + ra            
            sum_x = sum_x + x
            sum_rt = sum_rt + rt
            sum_x_c = sum_x_c + x_c
            sum_x_p = sum_x_p + x_p
            sum_nep = sum_nep + nep  
            sum_f_v2s = sum_f_v2s + f_v2s
            sum_f_v2s_2 = sum_f_v2s_2 + f_v2s_2            
            sum_rt_soil = sum_rt_soil + rt_soil
            sum_rt_soil_2 = sum_rt_soil_2 + rt_soil_2            
            sum_rt_veg = sum_rt_veg + rt_veg
            sum_rt_veg_2 = sum_rt_veg_2 + rt_veg_2            
            sum_x_c_soil = sum_x_c_soil + x_c_soil
            sum_x_c_veg = sum_x_c_veg + x_c_veg
            sum_x_p_veg = sum_x_p_veg + x_p_veg
            sum_x_p_soil = sum_x_p_soil + x_p_soil
            sum_ra_rate = sum_ra_rate + ra_rate
            sum_rh_rate = sum_rh_rate + rh_rate 
            sum_dist = sum_dist + dist
            sum_delta_cVeg = sum_delta_cVeg + delta_cVeg
            sum_delta_cSoil = sum_delta_cSoil + delta_cSoil
        k += 1
    ave_dict = {
            "cVeg": sum_cVeg / len(model_folders),
            "cSoil": sum_cSoil / len(model_folders), 
            "delta_cVeg": sum_delta_cVeg / len(model_folders),
            "delta_cSoil": sum_delta_cSoil / len(model_folders),
            "nep": sum_nep / len(model_folders),
            "f_v2s": sum_f_v2s / len(model_folders),  
            "f_v2s_2": sum_f_v2s_2 / len(model_folders),            
            "gpp": sum_gpp / len(model_folders), 
            "rh": sum_rh / len(model_folders),
            "ra": sum_ra / len(model_folders),
            "x": sum_x / len(model_folders), 
            "rt": sum_rt / len(model_folders),
            "x_c": sum_x_c / len(model_folders),
            "x_p": sum_x_p / len(model_folders), 
            #"nep": sum_nep,# / len(model_folders), 
            "dist": sum_dist / len(model_folders),
            "rt_veg": sum_rt_veg / len(model_folders),
            "rt_veg_2": sum_rt_veg_2 / len(model_folders),              
            "rt_soil": sum_rt_soil / len(model_folders), 
            "rt_soil_2": sum_rt_soil_2 / len(model_folders),
            "x_c_veg": sum_x_c_veg / len(model_folders),            
            "x_c_soil": sum_x_c_soil / len(model_folders),
            "x_p_veg": sum_x_p_veg / len(model_folders),
            "x_p_soil": sum_x_p_soil / len(model_folders),
            "ra_rate": sum_ra_rate / len(model_folders),
            "rh_rate": sum_rh_rate / len(model_folders),
            }    
    all_components.append(ave_dict) 
    
    times_avg = times
    #times_avg = avg_timeline(times, averaging)
    all_components.append(times_avg)  
    mods=list(model_names.values())
    mods.append("Mean")
    mods.append("Times")
    all_comp_dict = {mods[i]: all_components[i] for i in range(len(mods))}             
    return all_comp_dict

data_streams = namedtuple(
    'data_streams',
    ["cVeg", "cSoil", "gpp", "npp", "ra", "rh"]
    )        

def write_data_streams_cache(gm_path, gm):
    # var=ds.variables[var_name]
    if gm_path.exists():
        print("removing old cache file{}")
        os.remove(gm_path)
    names=gm._fields    

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
    
def get_global_mean_vars_all(model_folder,   # string e.g. "ab_classic"
                            experiment_name, # string, e.g. "CLASSIC_S2_"
                            lat_var,
                            lon_var,
                            ):
    # def nc_file_name(nc_var_name, experiment_name):
        # return experiment_name+"{}.nc".format(nc_var_name) if nc_var_name!="npp_nlim" else experiment_name+"npp.nc"

    def nc_global_mean_file_name(experiment_name):
        return experiment_name+"gm_all_vars.nc"

    data_str = msh(model_folder).data_str   
    names = data_str._fields
    conf_dict = confDict(model_folder)
    dataPath=Path(conf_dict["dataPath"])     
    
    if dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name)).exists():
        print(""" Found cached global mean files. If you want to recompute the global means
            remove the following files: """
        )
        print( dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name)))

        def get_cached_global_mean(vn):
            gm_path=dataPath.joinpath(
                nc_global_mean_file_name(experiment_name=experiment_name))               
            return nc.Dataset(str(gm_path)).variables[vn].__array__() 

        output=data_streams(*map(get_cached_global_mean, data_streams._fields))      
        return (
            output
        )

    else:
        gm=globalMask(file_name="common_mask_all_models.nc")
        # load an example file with mask
        # special case for YIBS that doesn't have a mask in all files ecept tas
        if model_folder=="jon_yib":
            template = nc.Dataset(dataPath.joinpath(msh(model_folder).nc_file_name("tas", 
                experiment_name=experiment_name))).variables['tas'][0,:,:].mask
        else:
            template = nc.Dataset(dataPath.joinpath(msh(model_folder).nc_file_name("cSoil", 
                experiment_name=experiment_name))).variables['cSoil'][0,:,:].mask
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
        gcm=resample_grid(
            source_coord_mask=gm, 
            target_coord_mask=CoordMask(
                    index_mask=np.zeros_like(template),
                    tr=SymTransformers(
                        ctr=msh(model_folder).make_model_coord_transforms(),
                        itr=msh(model_folder).make_model_index_transforms(), 
                        ),
                    ),    
            var=gm.index_mask, 
            method="nearest"
            )

        print("computing means, this may take some minutes...")

        def compute_and_cache_global_mean(vn):
            path = dataPath.joinpath(msh(model_folder).nc_file_name(vn, experiment_name=experiment_name))
            if vn=="npp_nlim": path=dataPath.joinpath(msh(model_folder).nc_file_name("npp", experiment_name=experiment_name))
            print(path)
            ds = nc.Dataset(str(path))
            vs=ds.variables
            lats= vs[lat_var].__array__()
            lons= vs[lon_var].__array__()
            var=ds.variables[vn]
            # check if we have a cached version (which is much faster)
            gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))
            
            #model_mask = msh(model_folder).spatial_mask(dataPath=Path(conf_dict["dataPath"]))
            #combined_mask = combine_masks ([model_mask,gcm])
            gm=global_mean_var(
                    lats,
                    lons,
                    #combined_mask.index_mask,
                    gcm.index_mask,
                    var      
            )     
            return gm * 86400 if vn in ["gpp", "npp", "npp_nlim", "rh", "ra"] else gm
        
        #map variables to data
        print(data_str._fields)
        output=data_str(*map(compute_and_cache_global_mean, data_str._fields)) 
        cVeg=output.cVeg if output.cVeg.shape[0]<500 else avg_timeline(output.cVeg, 12)
        if "cLitter" in names:
            cLitter=output.cLitter if output.cLitter.shape[0]<500 else avg_timeline(output.cLitter, 12)
        cSoil=output.cSoil if output.cSoil.shape[0]<500 else avg_timeline(output.cSoil, 12)        
        gpp=output.gpp if output.gpp.shape[0]<500 else avg_timeline(output.gpp, 12)
        if "npp" in names:
            npp=output.npp if output.npp.shape[0]<500 else avg_timeline(output.npp, 12)
        if "npp_nlim" in names:
            npp=output.npp_nlim if output.npp_nlim.shape[0]<500 else avg_timeline(output.npp_nlim, 12)            
        if "ra" in names:
            ra=output.ra if output.ra.shape[0]<500 else avg_timeline(output.ra, 12)    
        rh=output.rh if output.rh.shape[0]<500 else avg_timeline(output.rh, 12)          
        # for models like SDGVM where pool data starts earlier than gpp data
        if cVeg.shape[0]>gpp.shape[0]: cVeg=cVeg[cVeg.shape[0]-gpp.shape[0]:]        
        if "cLitter" in names and cLitter.shape[0]>gpp.shape[0]: 
            cLitter=cLitter[cLitter.shape[0]-gpp.shape[0]:]
        if cSoil.shape[0]>gpp.shape[0]: cSoil=cSoil[cSoil.shape[0]-gpp.shape[0]:]
        
        output_final=data_streams(
            cVeg=cVeg,
            cSoil=cLitter+cSoil if "cLitter" in names else cSoil,
            gpp=gpp, 
            npp=npp if ("npp" in names) or ("npp_nlim" in names) else gpp-ra,
            ra=ra if "ra" in names else gpp-npp,
            rh=rh,
            )      
        gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))
        write_data_streams_cache(gm_path, output_final)        
        return (
            output_final
        )    

def cov_mat(a, b):
    return np.cov(a, b, bias=True)


def sum_attribution(y, x1, x2):
    cm1 = cov_mat(x1, y)
    var_y = cm1[1, 1]
    cov_x1_y = cm1[0, 1]

    cm2 = cov_mat(x2, y)
    cov_x2_y = cm2[0, 1]

    return tuple((c/var_y for c in (cov_x1_y, cov_x2_y)))

def product_attribution(v, z1, z2):
    y, x1, x2 = map(np.log, (v, z1, z2))
    return sum_attribution(y, x1, x2)

def add_gridded_vars (
        model_names,
        experiment_names,
        global_mask,       
        ):
    model_folders=[(m) for m in model_names] 
    m_names=list(model_names.values())      
    g_mask=global_mask.index_mask    
    for experiment in experiment_names:
        print('\033[1m'+'. . . Adding variables for '+experiment+' experiment . . .')   
        print('\033[0m')                         
        k=0 # model counter         
        for mf in model_folders:
            print('\033[1m'+model_names[mf]+'\033[0m')
            #print('cSoil_total')
            experiment_name=m_names[k]+"_"+experiment+"_"
            conf_dict = confDict(mf)
            dataPath=Path(conf_dict["dataPath"])              
            file_path = dataPath.joinpath(experiment_name+"cSoil_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_soil_avg=ds.variables["cSoil"][:, :].data
            var_soil_diff=ds.variables["cSoil_diff"][:, :].data
            ds.close()
            if "cLitter" in msh(mf).data_str._fields:           
                file_path = dataPath.joinpath(experiment_name+"cLitter_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_litter_avg=ds.variables["cLitter"][:, :].data
                var_litter_diff=ds.variables["cLitter_diff"][:, :].data                
                ds.close()
                var_soil_total_avg=np.ma.array(var_soil_avg+var_litter_avg, mask=g_mask)
                var_soil_total_diff=np.ma.array(var_soil_diff+var_litter_diff, mask=g_mask)                
            else:
                var_soil_total_avg=np.ma.array(var_soil_avg, mask=g_mask)
                var_soil_total_diff=np.ma.array(var_soil_diff, mask=g_mask)     
            var_soil_total_avg[var_soil_total_avg<0.00001]=0.00001
            #print(np.ma.mean(var_soil_total_avg))                
            #print(np.ma.mean(var_soil_total_diff))                
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"cSoil_total_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var_avg = ds_new.createVariable('cSoil_total', "float32", ["lat", "lon"])
            var_avg[:, :] = var_soil_total_avg
            var_diff = ds_new.createVariable('cSoil_total_diff', "float32", ["lat", "lon"])
            var_diff[:, :] = var_soil_total_diff            
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()                      
            if "npp_nlim" in msh(mf).data_str._fields:  
                #print('npp_nlim')            
                file_path = dataPath.joinpath(experiment_name+"npp_nlim_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_npp_nlim_avg=ds.variables["npp_nlim"][:, :].data
                var_npp_nlim_diff=ds.variables["npp_nlim_diff"][:, :].data                
                ds.close()
                var_npp_avg=np.ma.array(var_npp_nlim_avg, mask=g_mask)
                var_npp_diff=np.ma.array(var_npp_nlim_diff, mask=g_mask)
                #print(np.ma.mean(var_npp_avg))                
                #print(np.ma.mean(var_npp_diff))
                # creating and writing a new NetCDF file 
                s = g_mask.shape
                n_lats, n_lons = s
                new_path=dataPath.joinpath(experiment_name+"npp_ave_res.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)         
                # creating variables                         
                var_avg = ds_new.createVariable('npp', "float32", ["lat", "lon"])
                var_avg[:, :] = var_npp_avg
                var_diff = ds_new.createVariable('npp_diff', "float32", ["lat", "lon"])
                var_diff[:, :] = var_npp_diff            
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
                # closing NetCDF file      
                ds_new.close()                       
            elif not ("npp" in msh(mf).data_str._fields):  
                #print('npp')            
                file_path = dataPath.joinpath(experiment_name+"gpp_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_gpp_avg=ds.variables["gpp"][:, :].data
                var_gpp_diff=ds.variables["gpp_diff"][:, :].data                
                ds.close()
                file_path = dataPath.joinpath(experiment_name+"ra_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_ra_avg=ds.variables["ra"][:, :].data
                var_ra_diff=ds.variables["ra_diff"][:, :].data                
                ds.close()                                              
                var_npp_avg=np.ma.array(var_gpp_avg-var_ra_avg, mask=g_mask)
                var_npp_diff=np.ma.array(var_gpp_diff-var_ra_diff, mask=g_mask)            
                # creating and writing a new NetCDF file 
                s = g_mask.shape
                n_lats, n_lons = s
                new_path=dataPath.joinpath(experiment_name+"npp_ave_res.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)         
                # creating variables                         
                var_avg = ds_new.createVariable('npp', "float32", ["lat", "lon"])
                var_avg[:, :] = var_npp_avg
                var_diff = ds_new.createVariable('npp_diff', "float32", ["lat", "lon"])
                var_diff[:, :] = var_npp_diff            
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
                # closing NetCDF file      
                ds_new.close() 
            if not ("ra" in msh(mf).data_str._fields):  
                #print('ra')             
                file_path = dataPath.joinpath(experiment_name+"gpp_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_gpp_avg=ds.variables["gpp"][:, :].data
                var_gpp_diff=ds.variables["gpp_diff"][:, :].data                
                ds.close()
                file_path = dataPath.joinpath(experiment_name+"npp_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_npp_avg=ds.variables["npp"][:, :].data
                var_npp_diff=ds.variables["npp_diff"][:, :].data                
                ds.close()                                              
                var_ra_avg=np.ma.array(var_gpp_avg-var_npp_avg, mask=g_mask)
                var_ra_diff=np.ma.array(var_gpp_diff-var_npp_diff, mask=g_mask)    
                #print(np.ma.mean(var_ra_avg))                
                #print(np.ma.mean(var_ra_diff))                  
                # creating and writing a new NetCDF file 
                s = g_mask.shape
                n_lats, n_lons = s
                new_path=dataPath.joinpath(experiment_name+"ra_ave_res.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)         
                # creating variables                         
                var_avg = ds_new.createVariable('ra', "float32", ["lat", "lon"])
                var_avg[:, :] = var_ra_avg
                var_diff = ds_new.createVariable('ra_diff', "float32", ["lat", "lon"])
                var_diff[:, :] = var_ra_diff            
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
                # closing NetCDF file      
                ds_new.close()
                
            # adding variables for uncertianty propagation
            
            # adding NEP = cVeg_diff + cSoil_total_diff  
            #print('nep')           
            file_path = dataPath.joinpath(experiment_name+"cVeg_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_cVeg_diff_array=ds.variables["cVeg_diff"][:, :].data          
            var_cVeg_diff=np.ma.array(var_cVeg_diff_array, mask=g_mask)
            var_cVeg_array=ds.variables["cVeg"][:, :].data 
            var_cVeg_array[var_cVeg_array<=0.00001]=0.00001           
            var_cVeg=np.ma.array(var_cVeg_array,mask=g_mask)
            #print('cVEG!!! '+star(np.ma.mean(var_cVeg)))
            ds.close()                                        
            var_nep=np.ma.array(var_cVeg_diff+var_soil_total_diff, mask=g_mask)  
            #print(np.ma.mean(var_nep))                                                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"nep_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('nep', "float32", ["lat", "lon"])
            var[:, :] = var_nep   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()
            
            # adding C loss from disturbances: dist = gpp - nep - ra - rh 
            #print('disturbance')           
            file_path = dataPath.joinpath(experiment_name+"gpp_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_gpp_array=ds.variables["gpp"][:, :].data
            var_gpp=np.ma.array(var_gpp_array, mask = g_mask)
            ds.close()                                        
            file_path = dataPath.joinpath(experiment_name+"ra_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_ra_array=ds.variables["ra"][:, :].data 
            var_ra = np.ma.array(var_ra_array, mask=g_mask) 
            var_ra [var_ra<=0.00001]=0.00001            
            ds.close()  
            file_path = dataPath.joinpath(experiment_name+"rh_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_rh_array=ds.variables["rh"][:, :].data 
            var_rh = np.ma.array(var_rh_array, mask=g_mask) 
            var_rh [var_rh<=0.00001]=0.00001  
            ds.close()  
            var_dist = np.ma.array(var_gpp - var_nep - var_ra - var_rh, mask = g_mask)
            var_dist[var_dist<0]=0  
            #print(np.ma.mean(var_dist))            
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"dist_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('dist', "float32", ["lat", "lon"])
            var[:, :] = var_dist   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding flux from vegetation to soil (cSoil_diff + rh) or (gpp - cVeg_diff - ra - dist)
            #print ('adding f_v2s')  
            #veg_dist_frac=1            
            var_f_v2s=np.ma.array(var_soil_diff+var_rh, mask = g_mask)
            var_f_v2s[var_f_v2s<=0]=0.00001  
            var_f_v2s_2=np.ma.array(var_gpp-var_cVeg_diff-var_ra-var_dist, mask = g_mask)
            var_f_v2s_2[var_f_v2s_2<=0]=0.00001 
            #print(np.ma.mean(var_f_v2s_1))
            #print(np.ma.mean(var_f_v2s_2))
            #print(np.ma.mean(var_f_v2s_1-var_f_v2s_2))
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"f_v2s_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('f_v2s', "float32", ["lat", "lon"])
            var[:, :] = var_f_v2s  
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding vegetation residence time (RT_veg = cVeg / (f_veg2soil + ra + dist) )
            #print ('adding RT_veg')                      
            var_RT_veg=np.ma.array(var_cVeg / (var_ra + var_dist + var_f_v2s) * 120, mask = g_mask)
            var_RT_veg2=np.ma.array(var_cVeg / (var_gpp - var_cVeg_diff) * 120, mask = g_mask)              
            #print(np.ma.mean(var_RT_veg))            
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"RT_veg_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('RT_veg', "float32", ["lat", "lon"])
            var[:, :] = var_RT_veg   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding soil residence time (RT_soil = cSoil / rh )
            #print ('adding RT_soil')           
            var_RT_soil=np.ma.array(var_soil_total_avg / var_rh * 120, mask = g_mask)            
            #var_RT_soil2=np.ma.array(var_soil_total_avg / (var_f_v2s - var_soil_total_diff) * 120, mask = g_mask)
            var_RT_soil[var_RT_soil>10000]=10000
            #print(np.ma.mean(var_RT_soil))            
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"RT_soil_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('RT_soil', "float32", ["lat", "lon"])
            var[:, :] = var_RT_soil   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()
            
            # adding total ecosystem carbon ( X = cVeg + cSoil_total )
            var_X=np.ma.array(var_cVeg + var_soil_total_avg, mask = g_mask)                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X', "float32", ["lat", "lon"])
            var[:, :] = var_X   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding ecosystem residence time ( RT = X / (ra + rh + dist) )
            var_RT=np.ma.array(var_X / (var_ra + var_rh + var_dist)* 120, mask = g_mask)                     
            var_RT[var_RT>10000]=10000
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"RT_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('RT', "float32", ["lat", "lon"])
            var[:, :] = var_RT   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()
            
            # adding equilibrium carbon storage capacity ( X_c = gpp * RT )
            var_X_c=np.ma.array(var_gpp * var_RT / 120, mask = g_mask)                     
            # print('X_c=gpp*RT/120: '+str(np.ma.mean(var_X_c))+' = '+str(np.ma.mean(var_gpp))+
                # ' * '+str(np.ma.mean(var_RT))+ ' / 120')
            var_X_c[var_X_c<0]=0.000001  
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_c_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_c', "float32", ["lat", "lon"])
            var[:, :] = var_X_c   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding carbon storage potential ( X_p = X_c - X )
            var_X_p=np.ma.array(var_X_c - var_X, mask = g_mask)                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_p_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_p', "float32", ["lat", "lon"])
            var[:, :] = var_X_p   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()
            
            # adding equilibrium carbon storage capacity for Vegetation ( X_c_veg = gpp * RT_veg )
            var_X_c_veg=np.ma.array(var_gpp * var_RT_veg / 120, mask = g_mask)                     
            # print('X_c=gpp*RT/120: '+str(np.ma.mean(var_X_c))+' = '+str(np.ma.mean(var_gpp))+
                # ' * '+str(np.ma.mean(var_RT))+ ' / 120')
            var_X_c_veg[var_X_c_veg<0]=0.000001  
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_c_veg_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_c_veg', "float32", ["lat", "lon"])
            var[:, :] = var_X_c_veg   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding carbon storage potential for Vegetation ( X_p_veg = X_c_veg - C_Veg )
            var_X_p_veg=np.ma.array(var_X_c_veg - var_cVeg, mask = g_mask)                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_p_veg_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_p_veg', "float32", ["lat", "lon"])
            var[:, :] = var_X_p_veg   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()            
            
            # adding equilibrium carbon storage capacity for soil ( X_c_soil = f_v2s * RT_soil )
            var_X_c_soil=np.ma.array(var_f_v2s * var_RT_soil / 120, mask = g_mask)                     
            # print('X_c=gpp*RT/120: '+str(np.ma.mean(var_X_c))+' = '+str(np.ma.mean(var_gpp))+
                # ' * '+str(np.ma.mean(var_RT))+ ' / 120')
            var_X_c_soil[var_X_c_soil<0]=0.000001  
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_c_soil_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_c_soil', "float32", ["lat", "lon"])
            var[:, :] = var_X_c_soil   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding carbon storage potential for soiletation ( X_p_soil = X_c_soil - C_soil )
            var_X_p_soil=np.ma.array(var_X_c_soil - var_soil_total_avg, mask = g_mask)                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_p_soil_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_p_soil', "float32", ["lat", "lon"])
            var[:, :] = var_X_p_soil   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()             
            
            k+=1 # model counter

            print('cVeg: '+ str(np.ma.mean(var_cVeg)))
            print('cVeg_max: '+ str(np.ma.max(var_cVeg)))            
            print('cVeg_min: '+ str(np.ma.min(var_cVeg)))
            print('cVeg_diff: '+ str(np.ma.mean(var_cVeg_diff)))
            print('cSoil: '+ str(np.ma.mean(var_soil_total_avg)))
            print('cSoil_min: '+ str(np.ma.min(var_soil_total_avg)))            
            #print('cSoil_max: '+ str(np.ma.max(var_soil_total_avg))) 
            #print('cSoil_min: '+ str(np.ma.min(var_soil_total_avg)))             
            print('cSoil_diff: '+ str(np.ma.mean(var_soil_total_diff)))
            print('X: '+str(np.ma.mean(var_X)))    
            print('X_max: '+str(np.ma.max(var_X))) 
            print('X_min: '+str(np.ma.min(var_X)))             
            print('nep: '+str(np.ma.mean(var_nep)))
            print('gpp: '+str(np.ma.mean(var_gpp)))
            #print('gpp_max: '+str(np.ma.max(var_gpp)))
            #print('gpp_min: '+str(np.ma.min(var_gpp)))
            #print('npp: '+str(np.mean(var_npp)))            
            print('ra: '+str(np.ma.mean(var_ra)))
            print('ra_max: '+str(np.ma.max(var_ra)))
            print('ra_min: '+str(np.ma.min(var_ra)))
            print('rh: '+str(np.ma.mean(var_rh)))
            print('rh_max: '+str(np.ma.max(var_rh)))
            print('rh_min: '+str(np.ma.min(var_rh)))            
            print('dist: '+str(np.ma.mean(var_dist)))
            print('dist_max: '+str(np.ma.max(var_dist)))
            print('dist_min: '+str(np.ma.min(var_dist)))            
            print('f_v2s: '+str(np.ma.mean(var_f_v2s)))
            #print('f_v2s_max: '+str(np.ma.max(var_f_v2s)))
            #print('f_v2s_min: '+str(np.ma.min(var_f_v2s)))             
            print('f_v2s2: '+str(np.ma.mean(var_f_v2s_2)))
            
            print('RT_veg: '+str(np.ma.mean(var_RT_veg)))
            print('RT_veg_max: '+str(np.ma.max(var_RT_veg)))
            print('RT_veg_min: '+str(np.ma.min(var_RT_veg)))            
            print('RT_veg2: '+str(np.ma.mean(var_RT_veg2)))
            print('RT_soil: '+str(np.ma.mean(var_RT_soil))) 
            print('RT_soil_max: '+str(np.ma.max(var_RT_soil)))   
            print('RT_soil_min: '+str(np.ma.min(var_RT_soil)))               
            #print('RT_soil2: '+str(np.ma.mean(var_RT_soil2)))               
            #print('dist_veg_error: '+str(np.ma.mean( var_gpp - var_cVeg_diff - var_ra - var_dist - var_f_v2s_1 )))    
            print('RT: '+str(np.ma.mean(var_RT))) 
            print('RT_max: '+str(np.ma.max(var_RT)))  
            print('RT_min: '+str(np.ma.min(var_RT)))  
            print('X_c: '+str(np.ma.mean(var_X_c))) 
            print('X_c_max: '+str(np.ma.max(var_X_c)))  
            print('X_c_min: '+str(np.ma.min(var_X_c))) 
            print('X_p: '+str(np.ma.mean(var_X_p))) 
            print('X_p_max: '+str(np.ma.max(var_X_p)))  
            print('X_p_min: '+str(np.ma.min(var_X_p)))  

            print('X_c_veg '+str(np.ma.max(var_X_c_veg)))  
            print('X_p_veg '+str(np.ma.max(var_X_p_veg)))  
            print('X_c_soil '+str(np.ma.max(var_X_c_soil)))  
            print('X_p_soil '+str(np.ma.max(var_X_p_soil)))              
            
    print('\033[1m'+'Done!')
    
def uncertainty_grids (
        model_names,
        experiment_names,
        global_mask,
        output_path,
        ):
    model_folders=[(m) for m in model_names]
    m_names=list(model_names.values())
    g_mask=global_mask.index_mask
    for experiment in experiment_names:
        print('\033[1m'+'. . . Computing uncertainty for '+experiment+' experiment . . .')
        print('\033[0m')
                  
        print('computing uncertainties...')
        for vn in ['cVeg', 'cSoil_total', 'gpp','ra','rh','f_v2s','nep',
            'RT_veg', 'RT_soil', 'dist', 'X', 'RT', 'X_c', 'X_p', 'X_c_veg', 'X_p_veg', 'X_c_soil', 'X_p_soil']:
            #print(vn)
            var_sum_zero=np.zeros_like(g_mask)
            var_sum=np.ma.array(var_sum_zero,mask = g_mask)
            var_diff_sqr_zero=np.zeros_like(g_mask)
            var_diff_sqr=np.ma.array(var_diff_sqr_zero,mask = g_mask)
            var_abs_diff_zero=np.zeros_like(g_mask)
            var_abs_diff=np.ma.array(var_abs_diff_zero,mask = g_mask)            
            delta_sum_zero=np.zeros_like(g_mask)
            delta_sum=np.ma.array(delta_sum_zero,mask = g_mask)
            delta_diff_sqr_zero=np.zeros_like(g_mask)
            delta_diff_sqr=np.ma.array(delta_diff_sqr_zero,mask = g_mask)   
            delta_abs_diff_zero=np.zeros_like(g_mask)
            delta_abs_diff=np.ma.array(delta_abs_diff_zero,mask = g_mask)            
            # computing ensemble mean       
            k=0 # model counter        
            for mf in model_folders:
                experiment_name=m_names[k]+"_"+experiment+"_"
                conf_dict = confDict(mf)
                dataPath=Path(conf_dict["dataPath"])       
                # if vn in ['gpp','ra','rh','cVeg','cSoil_total']: # for fluxes we use sum
                    # file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc")                
                # else: # for pools we use mean
                    # file_path = dataPath.joinpath(experiment_name+vn+"_res.nc")
                file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc") 
                ds = nc.Dataset(str(file_path))
                var=ds.variables[vn][:, :].data
                var_sum=var_sum+var
                if vn=="cVeg" or vn=="cSoil_total":
                    delta=ds.variables[vn+"_diff"][:, :].data
                    delta_sum=delta_sum+delta
                k+=1 # model counter
                ds.close() 
            mean=var_sum/len(model_folders) 
            delta_mean=delta_sum/len(model_folders)
            # computing uncertainty measures: standard deviation and average deviation                   
            k=0 # model counter        
            for mf in model_folders:
                experiment_name=m_names[k]+"_"+experiment+"_"
                conf_dict = confDict(mf)
                dataPath=Path(conf_dict["dataPath"])             
                # if vn in ['gpp','ra','rh','cVeg','cSoil_total']: # for fluxes we use sum    
                    # file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc")                
                # else:
                    # file_path = dataPath.joinpath(experiment_name+vn+"_res.nc")    
                file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc")    
                ds = nc.Dataset(str(file_path))
                var=ds.variables[vn][:, :].data           
                var_diff_sqr = var_diff_sqr + (var-mean)**2
                var_abs_diff = var_abs_diff + np.abs(var-mean)
                if vn=="cVeg" or vn=="cSoil_total":
                    delta=ds.variables[vn+"_diff"][:, :].data
                    delta_diff_sqr=delta_diff_sqr+(delta-delta_mean)**2 
                    delta_abs_diff=delta_abs_diff+np.abs(delta-delta_mean)
                k+=1 # model counter 
                ds.close()             
            variance = var_diff_sqr / (len(model_folders)-1)
            st_dev=np.sqrt(variance) 
            ave_dev=var_abs_diff / len(model_folders)
            # final masking
            var_mean_final = np.ma.array(mean,mask = g_mask)
            var_sd_final = np.ma.array(st_dev,mask = g_mask)
            var_avd_final = np.ma.array(ave_dev,mask = g_mask)
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=Path(output_path).joinpath(experiment+"_"+vn+"_uncertainty.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var_mean = ds_new.createVariable(vn+'_mean', "float32", ["lat", "lon"])
            var_mean[:, :] = var_mean_final
            var_sd = ds_new.createVariable(vn+'_sd', "float32", ["lat", "lon"])
            var_sd[:, :] = var_sd_final            
            var_avd = ds_new.createVariable(vn+'_avd', "float32", ["lat", "lon"])
            var_avd[:, :] = var_avd_final
            var_avd_relative = ds_new.createVariable(vn+'_avd_relative', "float32", ["lat", "lon"])
            relative_uncertainty = var_avd_final / abs(var_mean_final) *100
            relative_uncertainty [relative_uncertainty>300]=300
            var_avd_relative[:, :] = relative_uncertainty              
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close() 
            print(vn+': '+str(np.ma.mean(var_mean_final)))
            if vn=="cVeg" or vn=="cSoil_total":
                delta_variance = delta_diff_sqr / (len(model_folders)-1)
                delta_st_dev=np.sqrt(delta_variance)
                delta_ave_dev=delta_abs_diff / len(model_folders)
                # final masking
                var_mean_final = np.ma.array(delta_mean,mask = g_mask)
                var_sd_final = np.ma.array(delta_st_dev,mask = g_mask)
                var_avd_final = np.ma.array(delta_st_dev,mask = g_mask)                
                # creating and writing a new NetCDF file 
                s = g_mask.shape
                n_lats, n_lons = s
                new_path=Path(output_path).joinpath(experiment+"_"+vn+"_diff_uncertainty.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)         
                # creating variables                         
                var_mean = ds_new.createVariable(vn+'_diff_mean', "float32", ["lat", "lon"])
                var_mean[:, :] = var_mean_final
                var_sd = ds_new.createVariable(vn+'_diff_sd', "float32", ["lat", "lon"])
                var_sd[:, :] = var_sd_final            
                var_avd = ds_new.createVariable(vn+'_diff_avd', "float32", ["lat", "lon"])
                var_avd[:, :] = var_avd_final
                var_avd_relative = ds_new.createVariable(vn+'_diff_avd_relative', "float32", ["lat", "lon"])
                var_avd_relative[:, :] = var_avd_final / var_mean_final               
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
                # closing NetCDF file      
                ds_new.close()                 
            k=0 # model counter 
        
    print('\033[1m'+'Done!')    

def grid_attribution (
        experiment_names,
        global_mask,    
        data_path,    
        ):
    g_mask=global_mask.index_mask
    for experiment in experiment_names:
        print('\033[1m'+experiment+'\033[0m')  
        # attribution of uncertainty in X to GPP, RT and X_p 
        file_path = Path(data_path).joinpath(experiment+"_X_c_uncertainty.nc")
        ds = nc.Dataset(str(file_path))        
        X_c_mean_array=ds.variables["X_c_mean"][:, :].data
        X_c_mean=np.ma.array(X_c_mean_array, mask=g_mask)                
        X_c_avd_array=ds.variables["X_c_avd"][:, :].data
        X_c_avd=np.ma.array(X_c_avd_array, mask=g_mask)
        ds.close()        
        file_path = Path(data_path).joinpath(experiment+"_gpp_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        gpp_avd_array=ds.variables["gpp_avd"][:, :].data/120
        gpp_avd=np.ma.array(gpp_avd_array, mask=g_mask)
        gpp_mean_array=ds.variables["gpp_mean"][:, :].data/120
        gpp_mean=np.ma.array(gpp_mean_array, mask=g_mask)
        ds.close() 
        file_path = Path(data_path).joinpath(experiment+"_RT_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        rt_avd_array=ds.variables["RT_avd"][:, :].data
        rt_avd=np.ma.array(rt_avd_array, mask=g_mask)
        rt_mean_array=ds.variables["RT_mean"][:, :].data  
        rt_mean=np.ma.array(rt_mean_array, mask=g_mask)
        ds.close() 
        file_path = Path(data_path).joinpath(experiment+"_X_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_avd_array=ds.variables["X_avd"][:, :].data
        X_avd=np.ma.array(X_avd_array, mask=g_mask)
        X_mean_array=ds.variables["X_mean"][:, :].data  
        X_mean=np.ma.array(X_mean_array, mask=g_mask)
        ds.close()
        file_path = Path(data_path).joinpath(experiment+"_X_p_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_p_avd_array=ds.variables["X_p_avd"][:, :].data
        X_p_avd=np.ma.array(X_p_avd_array, mask=g_mask)
        X_p_mean_array=ds.variables["X_p_mean"][:, :].data  
        X_p_mean=np.ma.array(X_p_mean_array, mask=g_mask)
        ds.close()
        
        gpp_contrib_ecosystem = np.ma.array(gpp_avd * rt_mean, mask=g_mask)
        rt_contrib = np.ma.array(rt_avd * gpp_mean, mask=g_mask)
        #gpp_rt_contrib = np.ma.array(X_c_avd - (gpp_contrib + rt_contrib), mask=g_mask)
        X_c=np.ma.array(X_c_avd, mask=g_mask)
        
        # print('RT_contrib: '+str(np.ma.mean(rt_contrib))) 
        # print('gpp_contrib: '+str(np.ma.mean(gpp_contrib_ecosystem)))         
        # print('X_c_avd: '+str(np.ma.mean(X_c_avd))) 
        # print('X_p_contrib: '+str(np.ma.mean(X_p_avd)))
       
        # print('gpp: '+str(np.ma.mean(gpp_mean)))   
        # print('RT: '+str(np.ma.mean(rt_mean))) 
        # print('X_c: '+str(np.ma.mean(X_c_mean))) 

        # print('X_c reconstr: '+str(np.ma.mean(rt_mean*gpp_mean)))
        # print('X_p: '+str(np.ma.mean(X_p_mean)))        
        # print('X: '+str(np.ma.mean(X_mean)))         
        # print('X_p reconstr: '+str(np.ma.mean(X_c_mean - X_mean)))   

        file_path = Path(data_path).joinpath(experiment+"_RT_veg_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        RT_veg_avd_array=ds.variables["RT_veg_avd"][:, :].data
        RT_veg_mean_array=ds.variables["RT_veg_mean"][:, :].data
        RT_veg_avd=np.ma.array(RT_veg_avd_array, mask=g_mask)
        RT_veg_mean=np.ma.array(RT_veg_mean_array, mask=g_mask)        
        ds.close()          
        file_path = Path(data_path).joinpath(experiment+"_RT_soil_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        RT_soil_avd_array=ds.variables["RT_soil_avd"][:, :].data
        RT_soil_mean_array=ds.variables["RT_soil_mean"][:, :].data
        RT_soil_avd=np.ma.array(RT_soil_avd_array, mask=g_mask)
        RT_soil_mean=np.ma.array(RT_soil_mean_array, mask=g_mask)         
        ds.close()          
        file_path = Path(data_path).joinpath(experiment+"_f_v2s_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        f_v2s_avd_array=ds.variables["f_v2s_avd"][:, :].data/120
        f_v2s_mean_array=ds.variables["f_v2s_mean"][:, :].data/120
        f_v2s_avd=np.ma.array(f_v2s_avd_array, mask=g_mask)
        f_v2s_mean=np.ma.array(f_v2s_mean_array, mask=g_mask)         
        ds.close()          
        
        gpp_contrib = np.ma.array(gpp_avd * RT_veg_mean, mask=g_mask)
        RT_veg_contrib = np.ma.array(RT_veg_avd * gpp_mean, mask=g_mask)
        X_c_veg_mean = np.ma.array(gpp_mean*RT_veg_mean, mask=g_mask)    
     
        
        f_v2s_contrib = np.ma.array(f_v2s_avd * RT_soil_mean, mask=g_mask)
        RT_soil_contrib = np.ma.array(RT_soil_avd * f_v2s_mean, mask=g_mask)
        X_c_soil_mean = np.ma.array(f_v2s_mean*RT_soil_mean, mask=g_mask)           


        file_path = Path(data_path).joinpath(experiment+"_cVeg_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        cVeg_avd_array=ds.variables["cVeg_avd"][:, :].data
        cVeg_mean_array=ds.variables["cVeg_mean"][:, :].data
        cVeg_avd=np.ma.array(cVeg_avd_array, mask=g_mask)
        cVeg_mean=np.ma.array(cVeg_mean_array, mask=g_mask)        
        ds.close()          
        file_path = Path(data_path).joinpath(experiment+"_cSoil_total_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        cSoil_avd_array=ds.variables["cSoil_total_avd"][:, :].data
        cSoil_mean_array=ds.variables["cSoil_total_mean"][:, :].data
        cSoil_avd=np.ma.array(cSoil_avd_array, mask=g_mask)
        cSoil_mean=np.ma.array(cSoil_mean_array, mask=g_mask)         
        ds.close() 
        file_path = Path(data_path).joinpath(experiment+"_X_c_veg_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_c_veg_avd_array=ds.variables["X_c_veg_avd"][:, :].data
        X_c_veg_mean_array=ds.variables["X_c_veg_mean"][:, :].data
        X_c_veg_avd=np.ma.array(X_c_veg_avd_array, mask=g_mask)
        X_c_veg_mean=np.ma.array(X_c_veg_mean_array, mask=g_mask)        
        ds.close()
        file_path = Path(data_path).joinpath(experiment+"_X_p_veg_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_p_veg_avd_array=ds.variables["X_p_veg_avd"][:, :].data
        X_p_veg_mean_array=ds.variables["X_p_veg_mean"][:, :].data
        X_p_veg_avd=np.ma.array(X_p_veg_avd_array, mask=g_mask)
        X_p_veg_mean=np.ma.array(X_p_veg_mean_array, mask=g_mask)        
        ds.close() 
        file_path = Path(data_path).joinpath(experiment+"_X_c_soil_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_c_soil_avd_array=ds.variables["X_c_soil_avd"][:, :].data
        X_c_soil_mean_array=ds.variables["X_c_soil_mean"][:, :].data
        X_c_soil_avd=np.ma.array(X_c_soil_avd_array, mask=g_mask)
        X_c_soil_mean=np.ma.array(X_c_soil_mean_array, mask=g_mask)        
        ds.close()
        file_path = Path(data_path).joinpath(experiment+"_X_p_soil_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_p_soil_avd_array=ds.variables["X_p_soil_avd"][:, :].data
        X_p_soil_mean_array=ds.variables["X_p_soil_mean"][:, :].data
        X_p_soil_avd=np.ma.array(X_p_soil_avd_array, mask=g_mask)
        X_p_soil_mean=np.ma.array(X_p_soil_mean_array, mask=g_mask)        
        ds.close()          
        
        #X_p_veg_contrib = np.ma.array(RT_veg_contrib+gpp_contrib-cVeg_avd)
        #X_p_veg_mean = np.ma.array(RT_veg_mean+gpp_mean-cVeg_mean)   
        
        #X_p_soil_contrib = np.ma.array(RT_soil_contrib+f_v2s_contrib-cSoil_avd)
        #X_p_soil_mean = np.ma.array(RT_soil_mean+f_v2s_mean-cSoil_mean)
                
        # print('veg:')
        # print('gpp_contrib: '+str(np.ma.mean(gpp_contrib)))         
        # print('RT_contrib: '+str(np.ma.mean(RT_veg_contrib))) 
        # print('cVeg_avd: '+str(np.ma.mean(cVeg_avd))) 
        # print('X_p_veg_contrib: '+str(np.ma.mean(X_p_veg_avd))) 
        # print('gpp_mean: '+str(np.ma.mean(gpp_mean)))         
        # print('RT_mean: '+str(np.ma.mean(RT_veg_mean)))         
        # print('X_c_mean: '+str(np.ma.mean(X_c_veg_mean))) 
        # print('cVeg_mean:'+str(np.ma.mean(cVeg_mean)))  
        # print('X_p_veg_mean: '+str(np.ma.mean(X_p_veg_mean)))         
            
        # print('soil:')
        # print('f_v2s_contrib: '+str(np.ma.mean(f_v2s_contrib)))   
        # print('RT_contrib: '+str(np.ma.mean(RT_soil_contrib))) 
        # print('cSoil_avd: '+str(np.ma.mean(cSoil_avd)))         
        # print('X_p_soil_contrib: '+str(np.ma.mean(X_p_soil_avd)))         
        # print('f_v2s_mean: '+str(np.ma.mean(f_v2s_mean)))   
        # print('RT_mean: '+str(np.ma.mean(RT_soil_mean)))         
        # print('X_c_mean: '+str(np.ma.mean(X_c_soil_mean)))  
        # print('cSoil_mean:'+str(np.ma.mean(cSoil_mean))) 
        # print('X_p_soil_mean: '+str(np.ma.mean(X_p_soil_mean)))  
        
        
        gcm=globalMask(file_name="common_mask_all_models.nc")

        print("computing means, this may take some minutes...")

        # data_str = namedtuple(
        # 'data_str',
        # ["cVeg", "cSoil_total","cVeg_diff", "cSoil_total_diff","nep", "gpp", "ra", "rh",  "dist", "f_v2s", 
        # "X", "X_c", "X_p","RT", "RT_veg", "RT_soil"]
        # )  
                
        def global_mean_var(
            lats: np.ma.core.MaskedArray,
            lons: np.ma.core.MaskedArray,
            mask: np.array,
            var: nc._netCDF4.Variable,
        ) -> np.array:

            weight_mat = np.ma.array(get_weight_mat(lats, lons), mask=mask)
            wms = weight_mat.sum()
            res = (weight_mat * var[:, :]).sum() / wms #np.zeros(n_t)
            return res
        
        
        
        rt_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=rt_contrib)
        print("rt_contrib_global: "+str(rt_contrib_global))        
        gpp_contrib_eco_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=gpp_contrib_ecosystem)
        print("gpp_contrib_eco_global: "+str(gpp_contrib_eco_global)) 
        X_p_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=X_p_avd)
        print("X_p_contrib_global: "+str(X_p_contrib_global))         


        cVeg_avd_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=cVeg_avd)
        print("cVeg_avd_global: "+str(cVeg_avd_global)) 

        cSoil_avd_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=cSoil_avd)
        print("cSoil_total_avd_global: "+str(cSoil_avd_global)) 

        rt_veg_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=RT_veg_contrib)
        print("rt_veg_contrib_global: "+str(rt_veg_contrib_global))        
        gpp_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=gpp_contrib)
        print("gpp_contrib_global: "+str(gpp_contrib_global)) 
        X_p_veg_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=X_p_veg_avd)
        print("X_p_veg_contrib_global: "+str(X_p_veg_contrib_global))

        rt_soil_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=RT_soil_contrib)
        print("rt_soil_contrib_global: "+str(rt_soil_contrib_global))        
        f_v2s_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=f_v2s_contrib)
        print("f_v2s_contrib_global: "+str(f_v2s_contrib_global)) 
        X_p_soil_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                mask=gcm.index_mask,
                var=X_p_soil_avd)
        print("rt_soil_contrib_global: "+str(X_p_soil_contrib_global))        
        
        total_eco=rt_contrib_global+gpp_contrib_eco_global+X_p_contrib_global
        percent_rt=rt_contrib_global/total_eco*100
        percent_gpp_eco=gpp_contrib_eco_global/total_eco*100
        percent_X_p=X_p_contrib_global/total_eco*100        
        
        
        total_C=cVeg_avd_global+cSoil_avd_global
        cVeg_avd_part=cVeg_avd_global/total_C
        cSoil_avd_part=cSoil_avd_global/total_C
        print('total_C: '+str(total_C))
        print('cVeg_avd_part: '+str(cVeg_avd_part))
        print('cSoil_avd_part: '+str(cSoil_avd_part))
        
        # total_veg_soil=(rt_veg_contrib_global+gpp_contrib_global+X_p_veg_contrib_global+
                            # rt_soil_contrib_global+f_v2s_contrib_global+X_p_soil_contrib_global)        
        total_veg=rt_veg_contrib_global+gpp_contrib_global+X_p_veg_contrib_global
        total_soil=rt_soil_contrib_global+f_v2s_contrib_global+X_p_soil_contrib_global
        
        percent_rt_veg=rt_veg_contrib_global/total_veg*100*cVeg_avd_part
        percent_gpp=gpp_contrib_global/total_veg*100*cVeg_avd_part
        percent_X_p_veg=X_p_veg_contrib_global/total_veg*100*cVeg_avd_part      
        percent_rt_soil=rt_soil_contrib_global/total_soil*100*cSoil_avd_part
        percent_f_v2s=f_v2s_contrib_global/total_soil*100*cSoil_avd_part
        percent_X_p_soil=X_p_soil_contrib_global/total_soil*100*cSoil_avd_part       
        
        # pie charts
        
        fig=plt.figure(figsize=(7,7))
        ax1=fig.subplots()#axs[1]
        ax1.set_title('Global average % attribution of uncertainty in C storage')
                        
        #labels = '$\Delta$ GPP', '$\Delta$ RA', '$\Delta$NPP*$\Delta$ RH'#, '$\Delta$ Dist'
        labels = 'GPP', 'RT', 'X_p'#, '$\Delta$ Dist'
        sizes = [percent_gpp_eco, percent_rt, percent_X_p]
        ax1.pie(sizes, autopct='%1.1f%%', 
            startangle=90, counterclock=False, colors=("green", "darkorange", "blue"))#, "purple"))
   
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
        ax1.legend(labels, bbox_to_anchor =(1,1))
    
        plt.show() 

        # pie charts
        
        fig=plt.figure(figsize=(7,7))
        ax1=fig.subplots()#axs[1]
        ax1.set_title('Global average % attribution of uncertainty in C storage: veg + soil')
                        
        #labels = '$\Delta$ GPP', '$\Delta$ RA', '$\Delta$NPP*$\Delta$ RH'#, '$\Delta$ Dist'
        labels = 'GPP', 'RT_veg', 'X_p_veg', 'f_v2s', 'RT_soil', 'X_p_soil'#, '$\Delta$ Dist'
        sizes = [percent_gpp, percent_rt_veg, percent_X_p_veg, percent_f_v2s, percent_rt_soil, percent_X_p_soil ]
        ax1.pie(sizes, autopct='%1.1f%%', 
            startangle=90, counterclock=False, colors=("green", "darkorange", "blue", "cyan", "brown", "purple"))#, "purple"))
   
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
        ax1.legend(labels, bbox_to_anchor =(1,1))
    
        plt.show()
        
        #rt_contrib_global
        #gpp_contrib_eco_global
        #X_p_contrib_global
        #RGB = np.array(list(zip(rt_contrib_global, gpp_contrib_eco_global, X_p_contrib_global)))
        
        # height = len(gcm.lats) # can be any value
        # width = len(gcm.lons) # can be any value
        # channel = 3 # should be 3 for RGB images

        # rt_2d = np.reshape(rt_contrib, (-1, width))
        # gpp_2d=np.reshape(gpp_contrib_ecosystem, (-1, width))
        # X_p_2d=np.reshape(X_p_avd, (-1, width))

        # array = np.zeros((height, width, channel))
        # array[:,:,0] = rt_2d
        # array[:,:,1] = gpp_2d
        # array[:,:,2] = X_p_2d
        
        # image_data = array.astype(np.uint8)
        # image_filename = "3Dmap_image.jpeg"
        # image_data.save(image_filename)
        ############ flux approach ####################
        
        
    
        # if avd: suffix = '_avd' 
        # else: suffix = '_mean'
        
        # def compute_and_cache_global_mean(vn):
            # path = Path(dataPath).joinpath(experiment_name+'_'+vn+'_uncertainty.nc')
            # #if vn=="npp_nlim": path=dataPath.joinpath(msh(model_folder).nc_file_name("npp", experiment_name=experiment_name))
            # print(path)
            # ds = nc.Dataset(str(path))
            # vs=ds.variables
            # lats= vs['lat'].__array__()
            # lons= vs['lon'].__array__()
            # var=ds.variables[vn+suffix]
            # # check if we have a cached version (which is much faster)
            # #gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))
            
            # #model_mask = msh(model_folder).spatial_mask(dataPath=Path(conf_dict["dataPath"]))
            # #combined_mask = combine_masks ([model_mask,gcm])
            # gm=global_mean_var(
                    # lats,
                    # lons,
                    # #combined_mask.index_mask,
                    # gcm.index_mask,
                    # var      
            # )     
            # return gm #* 86400 if vn in ["gpp", "npp", "npp_nlim", "rh", "ra"] else gm
            
        # #map variables to data
        # #print(data_str._fields)
        # output=data_str(*map(compute_and_cache_global_mean, data_str._fields))         
        
       
        
        
        
        file_path = Path(data_path).joinpath(experiment+"_ra_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        ra_avd=ds.variables["ra_avd"][:, :].data
        ds.close()         
        file_path = Path(data_path).joinpath(experiment+"_rh_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        rh_avd=ds.variables["rh_avd"][:, :].data
        ds.close() 
        
def plot_attribution_output (
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
        u=all_comp_dict[m]["gpp"][start:stop]
        rt=all_comp_dict[m]["rt"][start:stop]
        x_mean=all_comp_dict["Mean"]["x"][start:stop]
        x_c_mean=all_comp_dict["Mean"]["x_c"][start:stop]
        x_p_mean=all_comp_dict["Mean"]["x_p"][start:stop]
        u_mean=all_comp_dict["Mean"]["gpp"][start:stop]
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
    
    print ('\033[1m'+'Attribution of summed deviations from the mean for all models ' +
        'to the differences in fluxes')     

    sum_pos_diff_nep=0
    sum_neg_diff_nep=0
    sum_gpp_pos=0
    sum_gpp_neg=0        
    sum_rh_pos=0
    sum_rh_neg=0
    sum_ra_pos=0
    sum_ra_neg=0
    sum_dist_pos=0
    sum_dist_neg=0     
    sum_cVeg_contrib_pos=0
    sum_cVeg_contrib_neg=0
    sum_ra_rate_contrib_pos=0
    sum_ra_rate_contrib_neg=0
    sum_ra_cVeg_inter_pos=0
    sum_ra_cVeg_inter_neg=0      
    sum_cSoil_contrib_pos=0
    sum_cSoil_contrib_neg=0   
    sum_rh_rate_contrib_pos=0
    sum_rh_rate_contrib_neg=0   
    sum_rh_cSoil_inter_pos=0
    sum_rh_cSoil_inter_neg=0        
    sum_f_v2s_pos=0
    sum_f_v2s_neg=0    

    for m in models: 
           
        nep=all_comp_dict[m]["nep"][start:stop]
        gpp=all_comp_dict[m]["gpp"][start:stop]        
        ra=all_comp_dict[m]["ra"][start:stop]
        rh=all_comp_dict[m]["rh"][start:stop]
        dist=all_comp_dict[m]["dist"][start:stop]
        cVeg=all_comp_dict[m]["cVeg"][start:stop]
        cSoil=all_comp_dict[m]["cSoil"][start:stop]
        ra_rate=all_comp_dict[m]["ra_rate"][start:stop]
        rh_rate=all_comp_dict[m]["rh_rate"][start:stop] 
        f_v2s=all_comp_dict[m]["f_v2s"][start:stop] 
         
        nep_mean=all_comp_dict["Mean"]["nep"][start:stop]
        gpp_mean=all_comp_dict["Mean"]["gpp"][start:stop]        
        ra_mean=all_comp_dict["Mean"]["ra"][start:stop]
        rh_mean=all_comp_dict["Mean"]["rh"][start:stop]
        dist_mean=all_comp_dict["Mean"]["dist"][start:stop]        
        cVeg_mean=all_comp_dict["Mean"]["cVeg"][start:stop]
        cSoil_mean=all_comp_dict["Mean"]["cSoil"][start:stop]
        ra_rate_mean=all_comp_dict["Mean"]["ra_rate"][start:stop]
        rh_rate_mean=all_comp_dict["Mean"]["rh_rate"][start:stop] 
        f_v2s_mean=all_comp_dict["Mean"]["f_v2s"][start:stop]         
        
        delta_nep=nep-nep_mean
        delta_gpp=gpp-gpp_mean
        delta_ra_contrib=ra-ra_mean
        delta_rh_contrib=rh-rh_mean
        delta_dist_contrib=dist-dist_mean        
        delta_cVeg=cVeg-cVeg_mean
        delta_cSoil=cSoil-cSoil_mean
        delta_ra_rate=ra_rate-ra_rate_mean
        delta_rh_rate=rh_rate-rh_rate_mean
        delta_f_v2s_contrib=f_v2s-f_v2s_mean
        
        # attribution of delta nep to delta gpp, ra, rh, disturbance
        #x_c_contrib=delta_x_c
        #x_p_contrib=-delta_x_p
         
        # attribution of product of pools and rates
        delta_cVeg_contrib=delta_cVeg*(ra_rate-delta_ra_rate/2)
        delta_ra_rate_contrib=delta_ra_rate*(cVeg-delta_cVeg/2) 
        delta_ra_cVeg_inter=delta_ra_contrib-delta_ra_rate_contrib-delta_cVeg_contrib         

        delta_cSoil_contrib=delta_cSoil*(rh_rate-delta_rh_rate/2)
        delta_rh_rate_contrib=delta_rh_rate*(cSoil-delta_cSoil/2) 
        delta_rh_cSoil_inter=delta_rh_contrib-delta_rh_rate_contrib-delta_cSoil_contrib  

        # summation of positive and negative contributions separately         
        
        pos_delta_nep=delta_nep.copy(); pos_delta_nep[pos_delta_nep<0]=0
        neg_delta_nep=delta_nep.copy(); neg_delta_nep[neg_delta_nep>0]=0
        sum_pos_diff_nep+=pos_delta_nep
        sum_neg_diff_nep+=neg_delta_nep     

        pos_delta_gpp=delta_gpp.copy(); pos_delta_gpp[pos_delta_gpp<0]=0
        neg_delta_gpp=delta_gpp.copy(); neg_delta_gpp[neg_delta_gpp>0]=0
        sum_gpp_pos+=pos_delta_gpp
        sum_gpp_neg+=neg_delta_gpp  
        
        pos_delta_ra=delta_ra_contrib.copy(); pos_delta_ra[pos_delta_ra<0]=0
        neg_delta_ra=delta_ra_contrib.copy(); neg_delta_ra[neg_delta_ra>0]=0
        sum_ra_pos+=pos_delta_ra
        sum_ra_neg+=neg_delta_ra    

        pos_delta_rh=delta_rh_contrib.copy(); pos_delta_rh[pos_delta_rh<0]=0
        neg_delta_rh=delta_rh_contrib.copy(); neg_delta_rh[neg_delta_rh>0]=0
        sum_rh_pos+=pos_delta_rh
        sum_rh_neg+=neg_delta_rh    

        pos_delta_dist=delta_dist_contrib.copy(); pos_delta_dist[pos_delta_dist<0]=0
        neg_delta_dist=delta_dist_contrib.copy(); neg_delta_dist[neg_delta_dist>0]=0
        sum_dist_pos+=pos_delta_dist
        sum_dist_neg+=neg_delta_dist    

        pos_delta_cVeg_contrib=delta_cVeg_contrib.copy(); pos_delta_cVeg_contrib[pos_delta_cVeg_contrib<0]=0
        neg_delta_cVeg_contrib=delta_cVeg_contrib.copy(); neg_delta_cVeg_contrib[neg_delta_cVeg_contrib>0]=0
        sum_cVeg_contrib_pos+=pos_delta_cVeg_contrib
        sum_cVeg_contrib_neg+=neg_delta_cVeg_contrib 

        pos_delta_cSoil_contrib=delta_cSoil_contrib.copy(); pos_delta_cSoil_contrib[pos_delta_cSoil_contrib<0]=0
        neg_delta_cSoil_contrib=delta_cSoil_contrib.copy(); neg_delta_cSoil_contrib[neg_delta_cSoil_contrib>0]=0
        sum_cSoil_contrib_pos+=pos_delta_cSoil_contrib
        sum_cSoil_contrib_neg+=neg_delta_cSoil_contrib 
        
        pos_delta_ra_rate_contrib=delta_ra_rate_contrib.copy(); pos_delta_ra_rate_contrib[pos_delta_ra_rate_contrib<0]=0
        neg_delta_ra_rate_contrib=delta_ra_rate_contrib.copy(); neg_delta_ra_rate_contrib[neg_delta_ra_rate_contrib>0]=0
        sum_ra_rate_contrib_pos+=pos_delta_ra_rate_contrib
        sum_ra_rate_contrib_neg+=neg_delta_ra_rate_contrib  
        
        pos_delta_ra_cVeg_inter=delta_ra_cVeg_inter.copy(); pos_delta_ra_cVeg_inter[pos_delta_ra_cVeg_inter<0]=0
        neg_delta_ra_cVeg_inter=delta_ra_cVeg_inter.copy(); neg_delta_ra_cVeg_inter[neg_delta_ra_cVeg_inter>0]=0
        sum_ra_cVeg_inter_pos+=pos_delta_ra_cVeg_inter
        sum_ra_cVeg_inter_neg+=neg_delta_ra_cVeg_inter          

        pos_delta_rh_rate_contrib=delta_rh_rate_contrib.copy(); pos_delta_rh_rate_contrib[pos_delta_rh_rate_contrib<0]=0
        neg_delta_rh_rate_contrib=delta_rh_rate_contrib.copy(); neg_delta_rh_rate_contrib[neg_delta_rh_rate_contrib>0]=0
        sum_rh_rate_contrib_pos+=pos_delta_rh_rate_contrib
        sum_rh_rate_contrib_neg+=neg_delta_rh_rate_contrib    

        pos_delta_rh_cSoil_inter=delta_rh_cSoil_inter.copy(); pos_delta_rh_cSoil_inter[pos_delta_rh_cSoil_inter<0]=0
        neg_delta_rh_cSoil_inter=delta_rh_cSoil_inter.copy(); neg_delta_rh_cSoil_inter[neg_delta_rh_cSoil_inter>0]=0
        sum_rh_cSoil_inter_pos+=pos_delta_rh_cSoil_inter
        sum_rh_cSoil_inter_neg+=neg_delta_rh_cSoil_inter         

        pos_delta_f_v2s=delta_f_v2s_contrib.copy(); pos_delta_f_v2s[pos_delta_f_v2s<0]=0
        neg_delta_f_v2s=delta_f_v2s_contrib.copy(); neg_delta_f_v2s[neg_delta_f_v2s>0]=0
        sum_f_v2s_pos+=pos_delta_f_v2s
        sum_f_v2s_neg+=neg_delta_f_v2s         
        
        
    plot_attribution_2 (
        times=times,
        delta_nep_pos=sum_pos_diff_nep,
        delta_nep_neg=sum_neg_diff_nep,
        gpp_pos=sum_gpp_pos,
        gpp_neg=sum_gpp_neg,        
        rh_pos=sum_rh_pos,
        rh_neg=sum_rh_neg,
        ra_pos=sum_ra_pos,
        ra_neg=sum_ra_neg,
        dist_pos=sum_dist_pos,
        dist_neg=sum_dist_neg,     
        cVeg_contrib_pos=sum_cVeg_contrib_pos,
        cVeg_contrib_neg=sum_cVeg_contrib_neg,
        ra_rate_contrib_pos=sum_ra_rate_contrib_pos,
        ra_rate_contrib_neg=sum_ra_rate_contrib_neg,
        ra_cVeg_inter_pos=sum_ra_cVeg_inter_pos,
        ra_cVeg_inter_neg=sum_ra_cVeg_inter_neg,      
        cSoil_contrib_pos=sum_cSoil_contrib_pos,
        cSoil_contrib_neg=sum_cSoil_contrib_neg,   
        rh_rate_contrib_pos=sum_rh_rate_contrib_pos,
        rh_rate_contrib_neg=sum_rh_rate_contrib_neg,   
        rh_cSoil_inter_pos=sum_rh_cSoil_inter_pos,
        rh_cSoil_inter_neg=sum_rh_cSoil_inter_neg,        
        f_v2s_pos=sum_f_v2s_pos,
        f_v2s_neg=sum_f_v2s_neg,        
        percent=percent,        
    )      
    
    
def plot_attribution_2 (
        times,
        delta_nep_pos,
        delta_nep_neg,
        gpp_pos,
        gpp_neg,        
        rh_pos,
        rh_neg,
        ra_pos,
        ra_neg,
        dist_pos,
        dist_neg,     
        cVeg_contrib_pos,
        cVeg_contrib_neg,
        ra_rate_contrib_pos,
        ra_rate_contrib_neg,
        ra_cVeg_inter_pos,
        ra_cVeg_inter_neg,      
        cSoil_contrib_pos,
        cSoil_contrib_neg,   
        rh_rate_contrib_pos,
        rh_rate_contrib_neg,   
        rh_cSoil_inter_pos,
        rh_cSoil_inter_neg, 
        f_v2s_pos,
        f_v2s_neg,
        percent,  
):
        
        gpp_contrib=gpp_pos-gpp_neg 
        rh_contrib=rh_pos-rh_neg
        ra_contrib=ra_pos-ra_neg
        dist_contrib=dist_pos-dist_neg
        cVeg_contrib=cVeg_contrib_pos-cVeg_contrib_neg
        ra_rate_contrib=ra_rate_contrib_pos-ra_rate_contrib_neg
        cSoil_contrib=cSoil_contrib_pos-cSoil_contrib_neg
        rh_rate_contrib=rh_rate_contrib_pos-rh_rate_contrib_neg
        ra_cVeg_inter_contrib=ra_cVeg_inter_pos-ra_cVeg_inter_neg
        rh_cSoil_inter_contrib=rh_cSoil_inter_pos-rh_cSoil_inter_neg
        f_v2s_contrib=f_v2s_pos-f_v2s_neg        
       
        fig5=plt.figure(figsize=(15,5))
        axs=fig5.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]})
        
        # contributions timeline
        ax=axs[0]

        # quadratic trends
            

        z = np.polyfit(times,  gpp_pos, 2)
        p = np.poly1d(z)  
        ax.plot(        
            times, 
            p(times),
            label="$\Delta$ GPP positive",
            color="green",
        ) 
        ax.fill_between(
                times,
                gpp_pos, 
                p(times),
                color="green",
                alpha=0.2                
                )         

        z = np.polyfit(times,  gpp_neg, 2)
        p = np.poly1d(z)  
        ax.plot(        
            times, 
            p(times),
            label="$\Delta$ GPP negative",
            color="green",
            linestyle="dashed",            
        ) 
        ax.fill_between(
                times,
                gpp_neg, 
                p(times),
                color="green",
                #linewidth=0.1,
                alpha=0.2                
                )          
   
        z = np.polyfit(times,  ra_pos, 2)
        p = np.poly1d(z)  
        ax.plot(        
            times, 
            p(times),
            label="$\Delta$ RA positive",
            color="blue",
        ) 
        ax.fill_between(
                times,
                ra_pos, 
                p(times),
                color="blue",
                alpha=0.2                
                )
   
        z = np.polyfit(times,  ra_neg, 2)
        p = np.poly1d(z)  
        ax.plot(        
            times, 
            p(times),
            label="$\Delta$ RA negative",
            color="blue",
            linestyle="dashed",            
        ) 
        ax.fill_between(
                times,
                ra_neg, 
                p(times),
                color="blue",
                alpha=0.2                
                )   
   
        z = np.polyfit(times,  rh_pos, 2)
        p = np.poly1d(z)  
        ax.plot(        
            times, 
            p(times),
            label="$\Delta$ rh positive",
            color="red",
        ) 
        ax.fill_between(
                times,
                rh_pos, 
                p(times),
                color="red",
                alpha=0.2                
                )
   
        z = np.polyfit(times,  rh_neg, 2)
        p = np.poly1d(z)  
        ax.plot(        
            times, 
            p(times),
            label="$\Delta$ rh negative",
            color="red",
            linestyle="dashed",            
        ) 
        ax.fill_between(
                times,
                rh_neg, 
                p(times),
                color="red",
                alpha=0.2                
                )  
                    
        # if abs(np.mean(dist_pos))>0.01:        
            # z = np.polyfit(times,  dist_pos, 2)
            # p = np.poly1d(z)  
            # ax.plot(        
                # times, 
                # p(times),
                # label="$\Delta$ dist positive",
                # color="purple",
            # ) 
            # ax.fill_between(
                    # times,
                    # dist_pos, 
                    # p(times),
                    # color="purple",
                    # alpha=0.2                
                    # )
        # if abs(np.mean(dist_neg))>0.01:        
            # z = np.polyfit(times,  dist_neg, 2)
            # p = np.poly1d(z)  
            # ax.plot(        
                # times, 
                # p(times),
                # label="$\Delta$ dist negative",
                # color="purple",
                # linestyle="dashed",            
            # ) 
            # ax.fill_between(
                    # times,
                    # dist_neg, 
                    # p(times),
                    # color="purple",
                    # alpha=0.2                
                    # )                 

        # bar charts       

        positive_delta_nep=sum(delta_nep_pos)/len(delta_nep_pos)
        negative_delta_nep=sum(delta_nep_neg)/len(delta_nep_neg)
        positive_gpp_contrib=sum(gpp_pos)/len(gpp_pos)
        negative_gpp_contrib=sum(gpp_neg)/len(gpp_neg)        
        positive_ra_contrib=sum(ra_pos)/len(ra_pos)
        negative_ra_contrib=sum(ra_neg)/len(ra_neg)
        positive_rh_contrib=sum(rh_pos)/len(rh_pos)
        negative_rh_contrib=sum(rh_neg)/len(rh_neg)
        positive_dist_contrib=sum(dist_pos)/len(dist_pos)
        negative_dist_contrib=sum(dist_neg)/len(dist_neg)        
        
        ax0=axs[1]  
        ax0.set_title('Temporal average of contributions')       
        ax0.axhline(0, color='black', ls='dashed')
        
        if abs(np.mean(positive_delta_nep+negative_delta_nep))>0.01:
            ax0.bar ('Net $\Delta$ nep', positive_delta_nep+negative_delta_nep, width=0.4, color="black", label='Net $\Delta$ NEP')        
        ax0.bar ('$\Delta$ gpp', positive_gpp_contrib, color="green", label='$\Delta$ gpp')
        ax0.bar ('$\Delta$ gpp', negative_gpp_contrib, color="green")        
        ax0.bar ('$\Delta$ ra', positive_ra_contrib, color="blue", label='$\Delta$ ra')
        ax0.bar ('$\Delta$ ra', negative_ra_contrib, color="blue")
        ax0.bar ('$\Delta$ rh', positive_rh_contrib, color="red", label='$\Delta$ rh')
        ax0.bar ('$\Delta$ rh', negative_rh_contrib, color="red")        
        # ax0.bar ('$\Delta$ dist', positive_dist_contrib, color="purple", label='$\Delta$ dist')
        # ax0.bar ('$\Delta$ dist', negative_dist_contrib, color="purple")   
        ax0.legend(bbox_to_anchor =(1,1))  
        ax0.grid() 

        abs_total=gpp_contrib+rh_contrib+ra_contrib+dist_contrib
                
        percent_gpp=gpp_contrib/abs_total*100
        percent_ra=ra_contrib/abs_total*100
        percent_rh=rh_contrib/abs_total*100
        percent_dist=dist_contrib/abs_total*100
        
        ######### Percents
        if percent==True:   
            fig3=plt.figure(figsize=(15,5))
            axs=fig3.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]})         
            # if delta==False: 
            ax=axs[0]      
            
            # % timeline
       
            # quadratic trends
            z = np.polyfit(times,  percent_gpp, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ gpp",
                color="green",
            ) 
            ax.fill_between(
                    times,
                    percent_gpp, 
                    p(times),
                    color="green",
                    alpha=0.2                
                    )  
                    
            z = np.polyfit(times,  percent_ra, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ NPP",
                color="blue",
            ) 
            ax.fill_between(
                    times,
                    percent_ra, 
                    p(times),
                    color="blue",
                    alpha=0.2                
                    )  

            z = np.polyfit(times,  percent_rh, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ rh",
                color="red",
            )         
            ax.fill_between(
                    times,
                    percent_rh, 
                    p(times),
                    color="red",
                    alpha=0.2                
                    )  

            # z = np.polyfit(times,  percent_dist, 2)
            # p = np.poly1d(z)  
            # ax.plot(        
                # times, 
                # p(times),
                # label="$\Delta$NPP*$\Delta$RT",
                # color="purple",
            # ) 
            # ax.fill_between(
                    # times,
                    # percent_dist, 
                    # p(times),
                    # color="purple",
                    # alpha=0.2                
                    # )  
            
            ax.set_title('% Contributions over time')
            ax.set_ylabel('%')
            ax.grid()
            
            # pie charts
            ax1=axs[1]
            ax1.set_title('Temporal average % contributions')
              

            labels = '$\Delta$ GPP', '$\Delta$ RA', '$\Delta$NPP*$\Delta$ RH'#, '$\Delta$ Dist'
            sizes = [np.mean(percent_gpp), np.mean(percent_ra), 
                np.mean(percent_rh)]#, np.mean(percent_dist)]
            ax1.pie(sizes, autopct='%1.1f%%', 
                startangle=90, counterclock=False, colors=("green", "blue", "red"))#, "purple"))
       
            ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
            ax1.legend(labels, bbox_to_anchor =(1,1))
        
        plt.show()    
        
        
def get_global_mean_uncertainty(dataPath,  
                            experiment_name, # string, e.g. "S2"
                            data_str, # named tuple
                            avd,
                            #lat_var,
                            #lon_var,
                            ):
    # def nc_file_name(nc_var_name, experiment_name):
        # return experiment_name+"{}.nc".format(nc_var_name) if nc_var_name!="npp_nlim" else experiment_name+"npp.nc"

    # def nc_global_mean_file_name(experiment_name):
        # return experiment_name+"gm_all_vars.nc"

    #data_str = msh(model_folder).data_str   
    #names = data_str._fields
    #conf_dict = confDict(model_folder)
    #dataPath=Path(conf_dict["dataPath"])     
    
    # if dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name)).exists():
        # print(""" Found cached global mean files. If you want to recompute the global means
            # remove the following files: """
        # )
        # print( dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name)))

        # def get_cached_global_mean(vn):
            # gm_path=dataPath.joinpath(
                # nc_global_mean_file_name(experiment_name=experiment_name))               
            # return nc.Dataset(str(gm_path)).variables[vn].__array__() 

        # output=data_streams(*map(get_cached_global_mean, data_streams._fields))      
        # return (
            # output
        # )

    # else:
    gcm=globalMask(file_name="common_mask_all_models.nc")
    # # load an example file with mask
    # # special case for YIBS that doesn't have a mask in all files ecept tas
    # if model_folder=="jon_yib":
        # template = nc.Dataset(dataPath.joinpath(msh(model_folder).nc_file_name("tas", 
            # experiment_name=experiment_name))).variables['tas'][0,:,:].mask
    # else:
        # template = nc.Dataset(dataPath.joinpath(msh(model_folder).nc_file_name("cSoil", 
            # experiment_name=experiment_name))).variables['cSoil'][0,:,:].mask
    # # gcm=project_2(
            # # source=gm,
            # # target=CoordMask(
                # # index_mask=np.zeros_like(template),
                # # tr=SymTransformers(
                    # # ctr=msh(model_folder).make_model_coord_transforms(),
                    # # itr=msh(model_folder).make_model_index_transforms()
                # # )
            # # )
    # # )
    # gcm=resample_grid(
        # source_coord_mask=gm, 
        # target_coord_mask=CoordMask(
                # index_mask=np.zeros_like(template),
                # tr=SymTransformers(
                    # ctr=msh(model_folder).make_model_coord_transforms(),
                    # itr=msh(model_folder).make_model_index_transforms(), 
                    # ),
                # ),    
        # var=gm.index_mask, 
        # method="nearest"
        # )

    print("computing means, this may take some minutes...")

    # data_str = namedtuple(
    # 'data_str',
    # ["cVeg", "cSoil_total","cVeg_diff", "cSoil_total_diff","nep", "gpp", "ra", "rh",  "dist", "f_v2s", 
    # "X", "X_c", "X_p","RT", "RT_veg", "RT_soil"]
    # )  
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
        #n_t = var.shape[0]
        res = (weight_mat * var[:, :]).sum() / wms #np.zeros(n_t)
        # for it in tqdm(range(n_t)):
            # el = (weight_mat * var[it, :, :]).sum() / wms
            # res[it] = el
        return res
    
    if avd: suffix = '_avd' 
    else: suffix = '_mean'
    
    def compute_and_cache_global_mean(vn):
        path = Path(dataPath).joinpath(experiment_name+'_'+vn+'_uncertainty.nc')
        #if vn=="npp_nlim": path=dataPath.joinpath(msh(model_folder).nc_file_name("npp", experiment_name=experiment_name))
        print(path)
        ds = nc.Dataset(str(path))
        vs=ds.variables
        lats= vs['lat'].__array__()
        lons= vs['lon'].__array__()
        var=ds.variables[vn+suffix]
        # check if we have a cached version (which is much faster)
        #gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))
        
        #model_mask = msh(model_folder).spatial_mask(dataPath=Path(conf_dict["dataPath"]))
        #combined_mask = combine_masks ([model_mask,gcm])
        gm=global_mean_var(
                lats,
                lons,
                #combined_mask.index_mask,
                gcm.index_mask,
                var      
        )     
        return gm #* 86400 if vn in ["gpp", "npp", "npp_nlim", "rh", "ra"] else gm
        
    #map variables to data
    #print(data_str._fields)
    output=data_str(*map(compute_and_cache_global_mean, data_str._fields)) 
    # cVeg=output.cVeg if output.cVeg.shape[0]<500 else avg_timeline(output.cVeg, 12)
    # if "cLitter" in names:
        # cLitter=output.cLitter if output.cLitter.shape[0]<500 else avg_timeline(output.cLitter, 12)
    # cSoil=output.cSoil if output.cSoil.shape[0]<500 else avg_timeline(output.cSoil, 12)        
    # gpp=output.gpp if output.gpp.shape[0]<500 else avg_timeline(output.gpp, 12)
    # if "npp" in names:
        # npp=output.npp if output.npp.shape[0]<500 else avg_timeline(output.npp, 12)
    # if "npp_nlim" in names:
        # npp=output.npp_nlim if output.npp_nlim.shape[0]<500 else avg_timeline(output.npp_nlim, 12)            
    # if "ra" in names:
        # ra=output.ra if output.ra.shape[0]<500 else avg_timeline(output.ra, 12)    
    # rh=output.rh if output.rh.shape[0]<500 else avg_timeline(output.rh, 12)          
    # # for models like SDGVM where pool data starts earlier than gpp data
    # if cVeg.shape[0]>gpp.shape[0]: cVeg=cVeg[cVeg.shape[0]-gpp.shape[0]:]        
    # if "cLitter" in names and cLitter.shape[0]>gpp.shape[0]: 
        # cLitter=cLitter[cLitter.shape[0]-gpp.shape[0]:]
    # if cSoil.shape[0]>gpp.shape[0]: cSoil=cSoil[cSoil.shape[0]-gpp.shape[0]:]
    
    # output_final=data_streams(
        # cVeg=cVeg,
        # cSoil=cLitter+cSoil if "cLitter" in names else cSoil,
        # gpp=gpp, 
        # npp=npp if ("npp" in names) or ("npp_nlim" in names) else gpp-ra,
        # ra=ra if "ra" in names else gpp-npp,
        # rh=rh,
        # )      
    # gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))
    # write_data_streams_cache(gm_path, output_final)        
    return (
        output
        #output_final
    )         