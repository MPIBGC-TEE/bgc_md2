from cmath import pi, sin
import sys
import json
from pathlib import Path
from collections import namedtuple
import netCDF4 as nc
import numpy as np
from sympy import Symbol
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
    TimeStepIterator2,
)
from copy import copy
from typing import Callable
# from general_helpers import month_2_day_index, monthly_to_yearly
from functools import reduce

sys.path.insert(0, '..')  # necessary to import general_helpers
import general_helpers as gh

# we will use the trendy output names directly in other parts of the output
Observables_annual = namedtuple(
    'Observables_annual',
    ["cVeg", "cLitter", "cSoil"]
)
Observables_monthly = namedtuple(
    'Observables_monthly',
    ["rh", "ra"]
)
Observables = namedtuple(
    "Observables",
    Observables_annual._fields + Observables_monthly._fields
)
# OrgDrivers=namedtuple(
#     "OrgDrivers",
#     ["gpp", "mrso", "tas"]
# )    
Drivers = namedtuple(
    "Drivers",
    ["npp", "mrso", "Ts"]  #
)
# As a safety measure we specify those parameters again as 'namedtuples', which are like a mixture of dictionaries and tuples
# They preserve order as numpy arrays which is great (and fast) for the numeric computations
# and they allow to access values by keys (like dictionaries) which makes it difficult to accidentally mix up values.

Constants = namedtuple(
    "Constants",
    [
        "cLitter_0",
        "cSoil_0",
        "cVeg_0",
        "npp_0",
       # "xi_0",
        "rh_0",
        "ra_0",
        # "r_C_root_litter_2_C_soil_passive",# here  we pretend to know these two rates
        # "r_C_root_litter_2_C_soil_slow",# it would be much better to know more
        "number_of_months"  # necessary to prepare the output in the correct lenght
    ]
)
# Note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated) 
# parameters. In this case we may use only the first entry e.g. to derive startvalues and their length (number of months)
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise 
# It is better to start with only a few

EstimatedParameters = namedtuple(
    "EstimatedParameters", [
        "beta_wood1",
        "beta_wood2",
        "beta_leaf",
        "beta_root",
        "T_0",
        "E",
        "KM",
        "r_C_wood1_2_C_wood3",
        "r_C_wood1_2_C_litter1",
        "r_C_wood2_2_C_wood4",
        "r_C_wood2_2_C_litter2",
        "r_C_wood3_2_C_litter1",
        "r_C_wood4_2_C_litter2",
        "r_C_leaf_2_C_litter3",
        "r_C_leaf_2_C_litter5",
        "r_C_root_2_C_litter4",
        "r_C_root_2_C_litter6",
        "r_C_fruit_2_C_litter3",
        "r_C_fruit_2_C_litter5",
        "r_C_litter1_2_C_som1",
        "r_C_litter1_2_C_som2",
        "r_C_litter2_2_C_som2",
        "r_C_litter2_2_C_som3",
        "r_C_litter3_2_C_som1",
        "r_C_litter3_2_C_som3",
        "r_C_litter4_2_C_som1",
        "r_C_litter4_2_C_som2",
        "r_C_litter5_2_C_som1",
        "r_C_litter6_2_C_som2",
        "r_C_som1_2_C_som3",
        "r_C_som2_2_C_som3",
        "r_C_som2_2_C_som4",
        "r_C_som3_2_C_som2",
        "r_C_som3_2_C_som4",
        "r_C_som4_2_C_som2",
        "r_C_som4_2_C_som3",
        "r_C_litter1_rh",
        "r_C_litter2_rh",
        "r_C_litter3_rh",
        "r_C_litter4_rh",
        "r_C_litter5_rh",
        "r_C_litter6_rh",
        "r_C_som1_rh",
        "r_C_som2_rh",
        "r_C_som3_rh",
        "r_C_som4_rh",
        'C_wood1_0',
        'C_wood2_0',
        'C_wood3_0',
        'C_wood4_0',
        'C_leaf_0',  # for the trendy data also the startvalues have to be estimated but
        'C_root_0',
        # C_root_0 can be inferred as cVeg_0-(C_leaf_0+C_wood_0)
        'C_litter1_0',
        'C_litter2_0',
        'C_litter3_0',
        'C_litter4_0',
        'C_litter5_0',
        # C_root_litter_0 can be inferred
        'C_som1_0',
        'C_som2_0',
        'C_som3_0',
        # C_soil_passive_0 can be inferred
    ]
)


# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated)
# parameters. In this case we may use only the first entry e.g. to derive startvalues. 
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise         
# to guard agains accidentally changed order we use a namedtuple again. Since B_func and u_func rely 
# on the correct ordering of the statevariables we build V dependent on this order 

# create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),  # platform independent path desc. (Windows vs. linux)
        models=['OCN'],
        variables=Observables._fields + OrgDrivers._fields
    )


experiment_name = "OCN_"


def nc_file_name(nc_var_name):
    return experiment_name + "{}.nc".format(nc_var_name)


def nc_global_mean_file_name(nc_var_name):
    return experiment_name + "{}_gm.nc".format(nc_var_name)


def get_global_vars(dataPath):
    # define function to average variables
    def f(tup):
        # define parts of function from nc file
        vn, fn = tup
        path = dataPath.joinpath(fn)
        ds = nc.Dataset(str(path))
        lats = ds.variables["latitude"]
        lons = ds.variables["longitude"]

        # check for npp/gpp/rh/ra to convert from kg/m2/s to kg/m2/day
        if vn in ["npp", "gpp", "rh", "ra"]:
            # for name, variable in ds.variables.items():
            #    for attrname in variable.ncattrs():
            #        print("{} -- {}".format(attrname, getattr(variable, attrname)))
            return (gh.global_mean(lats, lons, ds.variables[vn].__array__()) * 24 * 60 * 60)
        else:
            # for name, variable in ds.variables.items():
            #    for attrname in variable.ncattrs():
            #        print("{} -- {}".format(attrname, getattr(variable, attrname)))
            return (gh.global_mean(lats, lons, ds.variables[vn].__array__()))

    # Link symbols and data:
    # IBIS has annual vs monthly file names so they are linked separately
    # If all your data is similarly named you can do this in one step

    # Create annual file names (single step if files similarly named)
    o_names = [(f, "OCN_{}.nc".format(f)) for f in Observables_annual._fields]

    # Create monthly file names (can remove if done in one step above)
    monthly_names = [(f, "OCN_{}.nc".format(f)) for f in Observables_monthly._fields]
    # Extend name list with monthly names
    o_names.extend(monthly_names)

    # create file names for Drivers
    d_names = [(f, "OCN_{}.nc".format(f)) for f in Drivers._fields]

    # we want to drive with npp and can create it from gpp and ra 
    # observables
    odvs = OrgDrivers(*map(f, d_names))
    obss = Observables(*map(f, o_names))

    dvs = Drivers(
        npp=odvs.gpp - obss.ra,
        mrso=odvs.mrso,
        Ts=odvs.Ts
    )

    # Link symbols and data for Observables/Drivers
    # return (Observables(*map(f, o_names)),Drivers(*map(f,d_names)))
    return (obss, dvs)


def get_example_site_vars(dataPath):
    # pick up 1 site
    s = slice(None, None, None)  # this is the same as :
    t = s, 72, 117  # [t] = [:,49,325]

    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them
        ds = nc.Dataset(str(path))
        return ds.variables[vn][t]

    # note that there use S3 for IBIS model
    o_names = [(f, "OCN_{}.nc".format(f)) for f in Observables._fields]
    d_names = [(f, "OCN_{}.nc".format(f)) for f in OrgDrivers._fields]

    # we want to drive with npp and can create it from gpp and ra 
    # observables
    odvs = OrgDrivers(*map(f, d_names))
    obss = Observables(*map(f, o_names))

    dvs = Drivers(
        npp=odvs.gpp - obss.ra,
        mrso=odvs.mrso,
        Ts=odvs.Ts
    )
    return (obss, dvs)


def get_global_mean_vars(dataPath):
    o_names = Observables._fields
    d_names = Drivers._fields
    names = o_names + d_names

    if all([dataPath.joinpath(nc_global_mean_file_name(vn)).exists() for vn in names]):
        print(""" Found cached global mean files. If you want to recompute the global means
            remove the following files: """
              )
        for vn in names:
            print(dataPath.joinpath(nc_global_mean_file_name(vn)))

        def get_cached_global_mean(vn):
            return gh.get_cached_global_mean(dataPath.joinpath(nc_global_mean_file_name(vn)), vn)

        return (
            Observables(*map(get_cached_global_mean, o_names)),
            Drivers(*map(get_cached_global_mean, d_names))
        )

    else:
        # we now check if any of the arrays has a time lime containing nan values 
        # APART FROM values that are already masked by the fillvalue
        print("computing masks to exclude pixels with nan entries, this may take some minutes...")

        def f(vn):
            path = dataPath.joinpath(nc_file_name(vn))
            ds = nc.Dataset(str(path))
            # scale fluxes vs pools
            var = ds.variables[vn]
            return gh.get_nan_pixel_mask(var)

        masks = [f(name) for name in names]
        # We compute the common mask so that it yields valid pixels for ALL variables 
        combined_mask = reduce(lambda acc, m: np.logical_or(acc, m), masks)
        print("computing means, this may also take some minutes...")

        def compute_and_cache_global_mean(vn):
            path = dataPath.joinpath(nc_file_name(vn))
            ds = nc.Dataset(str(path))
            vs = ds.variables
            print(vs)
            lats = vs["latitude"].__array__()
            lons = vs["longitude"].__array__()
            print(vn)
            var = ds.variables[vn]
            # check if we have a cached version (which is much faster)
            gm_path = dataPath.joinpath(nc_global_mean_file_name(vn))

            gm = gh.global_mean_var(
                lats,
                lons,
                combined_mask,
                var
            )
            gh.write_global_mean_cache(
                gm_path,
                gm,
                vn
            )
            return gm * 86400 if vn in ["gpp", "rh", "ra", "npp"] else gm

        # map variables to data
        return (
            Observables(*map(compute_and_cache_global_mean, o_names)),
            Drivers(*map(compute_and_cache_global_mean, d_names))
        )


def make_npp_func(dvs):
    def func(day):
        month = gh.day_2_month_index(day)
        # kg/m2/s kg/m2/day;
        return (dvs.npp[month])  # * 86400

    return func


def make_xi_func(dvs):
    def func(day):
        return 1.0 # preliminary fake for lack of better data...
    return func


def make_func_dict(mvs, dvs):
    return {
        "NPP": make_npp_func(dvs),
        "xi": make_xi_func(dvs)
    }


# We now build the essential object to run the model forward. We have a
# - startvector $V_0$ and 
# - a function $f$ to compute the next value: $V_{it+1} =f(it,V_{it})$
#   the dependence on the iteration $it$ allows us to represent drivers that
#   are functions of time 
#
# So we can compute all values:
#
# $V_1=f(0,V_0)$
#
# $V_2=f(1,V_1)$
#
# ...
#
# $V_n+1=f(n,V_n)$
#
# Technically this can be implemented as an `iterator` object with a `__next__()` method to move our system one step forward in time. 
#
# What we want to build is an object `it_sym` that can use as follows.
# ```python
# for i in range(10):
#     print(it_sym.__next__())
# ```
# to print out the first 10 values.
#
# If iterators had not been invented yet we would invent them now, because they capture exactly the mathematical concept of an initial value system, 
# without all the nonessential technical details of e.g. how many steps we want to make or which of the results we want to store.
# This is essential if we want to test and use the iterator later in different scenarios but avoid reimplemtation of the basic timestep. 
#
# Remark:
#
# If we were only interested in the timeseries of the pool contents `bgc_md2` could compute the solution automatically without the need to build an iterator ourselves. 
# In our case we are also interested in tracking the autotrophic and heterotrophic respiration and therefore have to recreate and extend the code `bgc_md2` would use in the background.
# We will let `bgc_md2` do part of the work and derive numeric functions for the Compartmental matrix $B$ and the input $u$ and the Outfluxes - from which to compute $ra$ $rh$ - from our symbolic description but we will build our own iterator by combining these functions.    
# We will proceed as follows:
# - create $V_0$ 
# - build the function $f$

def make_iterator_sym(
        mvs,
        V_init,  #: StartVector,
        par_dict,
        func_dict,
        delta_t_val=1  # defaults to 1day timestep
):
    B_func, u_func = gh.make_B_u_funcs_2(mvs, par_dict, func_dict, delta_t_val)

    sv = mvs.get_StateVariableTuple()
    # mass production of output functions

    n = len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct 
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr = np.array(
        [V_init.__getattribute__(str(v)) for v in sv] +
        [V_init.ra, V_init.rh]
    ).reshape(n + 2, 1)  # reshaping is neccessary for matmul (the @ in B @ X)

    # To compute the ra and rh we have to some up the values for autotrophic and heterotrophic respiration we have to sum up outfluxes.
    # We first create numerical functions for all the outfluxes from our symbolic description.
    numOutFluxesBySymbol = {
        k: gh.numfunc(expr_cont, mvs, delta_t_val=delta_t_val, par_dict=par_dict, func_dict=func_dict)
        for k, expr_cont in mvs.get_OutFluxesBySymbol().items()
    }

    def f(it, V):
        X = V[0:n]
        b = u_func(it, X)
        B = B_func(it, X)
        X_new = X + b + B @ X
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep

        ra = np.sum(
            [
                numOutFluxesBySymbol[Symbol(k)](it, *X)
                for k in ["C_leaf", "C_wood1", "C_wood2", "C_fruit", "C_root"]
                if Symbol(k) in numOutFluxesBySymbol.keys()
            ]
        )
        rh = np.sum(
            [
                numOutFluxesBySymbol[Symbol(k)](it, *X)
                for k in
                ["C_litter1","C_litter2","C_litter3","C_litter4","C_litter5","C_litter6", "C_som1","C_som2","C_som3","C_som4"]
                if Symbol(k) in numOutFluxesBySymbol.keys()
            ]
        )
        V_new = np.concatenate((X_new.reshape(n, 1), np.array([ra, rh]).reshape(2, 1)), axis=0)

        return V_new

    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


def make_StartVector(mvs):
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()] +
        ["ra", "rh"]
    )


# +
def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers
) -> Callable[[np.ndarray], np.ndarray]:
    # To compute numeric solutions we will need to build and iterator 
    # as we did before. As it will need numeric values for all the parameters 
    # we will have to create a complete dictionary for all the symbols
    # exept those for the statevariables and time.
    # This set of symbols does not change during the mcmc procedure, since it only
    # depends on the symbolic model.
    # Therefore we create it outside the mcmc loop and bake the result into the 
    # param2res function.
    # The iterator does not care if they are estimated or constant, it only 
    # wants a dictionary containing key: value pairs for all
    # parameters that are not state variables or the symbol for time
    # StartVector=make_StartVector(mvs)

    srm = mvs.get_SmoothReservoirModel()
    model_par_dict_keys = srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))] +
        list(mvs.get_StateVariableTuple())
    )

    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test

    #     seconds_per_day = 86400
    #     def npp_func(day):
    #         month=gh.day_2_month_index(day)
    #         return dvs.npp[month] * seconds_per_day   # kg/m2/s kg/m2/day;

    def param2res(pa):
        epa = EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this
        StartVector = make_StartVector(mvs)
        V_init = StartVector(
            C_wood1=epa.C_wood1_0,
            C_wood2=epa.C_wood2_0,
            C_wood3=epa.C_wood3_0,
            C_wood4= epa.C_wood4_0,
            C_leaf= epa.C_leaf_0,
            C_root= epa.C_root_0,
            C_fruit= cpa.cVeg_0 - (epa.C_wood1_0 + epa.C_wood2_0 + epa.C_wood3_0 + epa.C_wood4_0 + epa.C_leaf_0 + epa.C_root_0),
            C_litter1=epa.C_litter1_0,
            C_litter2= epa.C_litter2_0,
            C_litter3= epa.C_litter3_0,
            C_litter4=epa.C_litter4_0,
            C_litter5= epa.C_litter5_0,
            C_litter6= cpa.cLitter_0 - (epa.C_litter1_0 + epa.C_litter2_0 + epa.C_litter3_0 + epa.C_litter4_0 + epa.C_litter5_0),
            C_som1= epa.C_som1_0,
            C_som2= epa.C_som2_0,
            C_som3= epa.C_som3_0,
            C_som4= cpa.cSoil_0-(epa.C_som1_0 + epa.C_som2_0 + epa.C_som3_0),
            ra = cpa.ra_0,
            rh = cpa.rh_0
        )
        # next we create the parameter dict for the iterator
        # The iterator does not care if they are estimated or not so we look for them
        # in the combination
        apa = {**cpa._asdict(),**epa._asdict()}
        model_par_dict = {
            Symbol(k):v for k,v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }

        # added by cybian for fix some errors
        model_par_dict = {
            Symbol(k): v for k, v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }

        print('model_par_dict:', model_par_dict)
        # from IPython import embed;embed()

        # Beside the par_dict the iterator also needs the python functions to replace the symbolic ones with
        # our fake xi_func could of course be defined outside of param2res but in general it
        # could be build from estimated parameters and would have to live here...
        #         def xi_func(day):
        #             return 1.0 # preliminary fake for lack of better data...

        #         func_dict={
        #             'NPP':npp_func,
        #              'xi':xi_func
        #         }

        func_dict = make_func_dict(mvs, dvs)

        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2) 
        delta_t_val = 1
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
        )

        # Now that we have the iterator we can start to compute.
        # the first thing to notice is that we don't need to store all values (daili yi delta_t_val=1)
        # since the observations are recorded monthly while our iterator possibly has a smaller timestep.
        # - for the pool contents we only need to record the last day of the month
        # - for the respiration(s) ra and rh we want an accumulated value (unit mass) 
        #   have to sum up the products of the values*delta_t over a month
        # 
        # Note: check if TRENDY months are like this...
        # days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        sols = []
        dpm = 30  #
        n = len(V_init)

        # added by cybian 
        rhs = np.zeros(cpa.number_of_months)
        ras = np.zeros(cpa.number_of_months)
        number_of_years = int(cpa.number_of_months / 12)
        print('number_of_years:', number_of_years)

        cVegs = np.zeros(number_of_years)
        cLitters = np.zeros(number_of_years)
        cSoils = np.zeros(number_of_years)
        dpy = 30 * 12
        m_id = 0

        for y in range(number_of_years):
            # print('y:',y)
            cVeg_ave = 0
            cLitter_ave = 0
            cSoil_ave = 0
            for m in range(12):
                # print('y:',y,'m:',m,'m_id:',m_id)
                # dpm = days_per_month[ m % 12]
                mra_ave = 0.0
                mrh_ave = 0.0
                for d in range(dpm):  # int(dpm/delta_t_val)
                    v = it_sym.__next__()  # .reshape(n,)
                    # actually the original datastream seems to be a flux per area (kg m^-2 ^-s )
                    # at the moment the iterator also computes a flux but in kg^-2 ^day

                    V = StartVector(*v)
                    cVeg_ave=np.array(cVeg_ave, dtype=object)+ float(V.C_leaf + V.C_wood1 + V.C_wood2 + V.C_wood3 + V.C_wood4 + V.C_root + V.C_fruit)
                    cLitter_ave =np.array(cLitter_ave, dtype=object) +float(V.C_litter1 + V.C_litter2 + V.C_litter3 + V.C_litter4 + V.C_litter5 + V.C_litter6)
                    cSoil_ave =np.array(cSoil_ave, dtype=object)+ float(V.C_som1 + V.C_som2 + V.C_som3 + V.C_som4)
                    mrh_ave = np.array(mrh_ave) + V.rh
                    mra_ave = np.array(mra_ave) + V.ra

                # print('here:m_id:',m_id)
                rhs[m_id] = mrh_ave / dpm
                ras[m_id] = mrh_ave / dpm
                m_id = m_id + 1

            # print('Here:y:',y)
            cVegs[y] = np.array(cVeg_ave, dtype=object) / dpy
            cLitters[y] = np.array(cLitter_ave, dtype=object) / dpy
            cSoils[y] = np.array(cSoil_ave, dtype=object) / dpy

        return Observables(
            cVeg=cVegs,
            cSoil=cSoils,
            cLitter=cLitters,
            rh=rhs,
            ra=ras
        )

        # comment by cybian

    #                 #from IPython import embed;embed()
    #                 o=Observables(
    #                     cVeg=float(V.C_leaf+V.C_wood+V.C_root),
    #                     cLitter=float(V.C_mll + V.C_mwl + V.C_sll + V.C_swl + V.C_lll + V.C_lwl),
    #                     cSoil=float(V.C_mrl + V.C_srl + V.C_lrl + V.C_mic + V.C_prot + V.C_nonprot + V.C_pass),
    #                     ra=V.ra/seconds_per_day,
    #                     rh=V.rh/seconds_per_day # the data is per second while the time units for the iterator refer to days
    #                 )
    #                 sols.append(o)

    # sol=np.stack(sols)
    # convert to yearly output if necessary (the monthly pool data looks very "yearly")
    # sol_yr=np.zeros(int(cpa.number_of_months/12)*sol.shape[1]).reshape([int(cpa.number_of_months/12),sol.shape[1]])
    # for i in range(sol.shape[1]):
    #   sol_yr[:,i]=monthly_to_yearly(sol[:,i])
    # sol=sol_yr
    # return sol

    return param2res


# -

def make_param_filter_func(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters
) -> Callable[[np.ndarray], bool]:
    # find position of beta_leaf and beta_wood
    beta_leaf_ind = EstimatedParameters._fields.index("beta_leaf")
    beta_wood1_ind = EstimatedParameters._fields.index("beta_wood1")
    beta_wood2_ind = EstimatedParameters._fields.index("beta_wood2")
    beta_root_ind = EstimatedParameters._fields.index("beta_root")

    def isQualified(c):
        beta_leaf_ind
        cond1 = (c >= c_min).all()
        cond2 = (c <= c_max).all()
        cond3 = c[beta_leaf_ind] + c[beta_wood1_ind] + c[beta_wood2_ind] + c[beta_root_ind] <= 1
        return (cond1 and cond2 and cond3)

    return isQualified


def make_weighted_cost_func(
        obs: Observables
) -> Callable[[Observables], np.float64]:
    # first unpack the observation array into its parts
    # cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
    def costfunction(out_simu: np.ndarray) -> np.float64:
        # fixme
        #   as indicated by the fact that the function lives in this
        #   model-specific module it is not apropriate for (all) other models.
        #   There are model specific properties:
        #   1.) The weight for the different observation streams
        #
        number_of_ys = out_simu.cVeg.shape[0]
        number_of_ms = out_simu.rh.shape[0]

        J_obj1 = np.mean((out_simu.cVeg - obs.cVeg) ** 2) / (2 * np.var(obs.cVeg))
        J_obj2 = np.mean((out_simu.cLitter - obs.cLitter) ** 2) / (2 * np.var(obs.cLitter))
        J_obj3 = np.mean((out_simu.cSoil - obs.cSoil) ** 2) / (2 * np.var(obs.cSoil))

        J_obj4 = np.mean((out_simu.rh - obs.rh) ** 2) / (2 * np.var(obs.rh))

        J_new = (J_obj1 + J_obj2 + J_obj3) / 200 + J_obj4 / 4
        # to make this special costfunction comparable (in its effect on the
        # acceptance rate) to the general costfunction proposed by Feng we
        # rescale it by a factor
        return J_new * 400

    return costfunction


def pesudo_yearly_to_monthly(yearly):
    # Added by cybian just for extending the data from yearly to monthly for matching all variables' length
    # not used in InspectModel.py
    yearly_len = len(yearly)
    months_per_year = 12
    monthly_data = np.zeros(yearly_len * months_per_year) + 1

    for iyear in range(0, yearly_len):
        for imonth in range(0, months_per_year):
            monthly_data[iyear * 12:(iyear + 1) * 12] = monthly_data[iyear * 12:(iyear + 1) * 12] * yearly[
                iyear]  # (yearly[iyear]-yearly[iyear]*sin(2*pi/12*imonth - pi/2))

    return monthly_data


def make_traceability_iterator(mvs, dvs, cpa, epa):
    par_dict = {
        Symbol(k): v for k, v in {
            "beta_wood1": epa.beta_wood1,
            "beta_wood2": epa.beta_wood2,
            "beta_leaf": epa.beta_leaf,
            "beta_root": epa.beta_root,
            "T_0": epa.T_0,
            "E": epa.E,
            "KM": epa.KM,

            "r_C_litter1_rh": epa.r_C_litter1_rh,
            "r_C_litter2_rh": epa.r_C_litter2_rh,
            "r_C_litter3_rh": epa.r_C_litter3_rh,
            "r_C_litter4_rh": epa.r_C_litter4_rh,
            "r_C_litter5_rh": epa.r_C_litter5_rh,
            "r_C_litter6_rh": epa.r_C_litter6_rh,
            "r_C_som1_rh": epa.r_C_som1_rh,
            "r_C_som2_rh": epa.r_C_som2_rh,
            "r_C_som3_rh": epa.r_C_som3_rh,
            "r_C_som4_rh": epa.r_C_som4_rh,
            "r_C_wood1_2_C_wood3": epa.r_C_wood1_2_C_wood3,
            "r_C_wood1_2_C_litter1": epa.r_C_wood1_2_C_litter1,
            "r_C_wood2_2_C_wood4": epa.r_C_wood2_2_C_wood4,
            "r_C_wood2_2_C_litter2": epa.r_C_wood2_2_C_litter2,
            "r_C_wood3_2_C_litter1": epa.r_C_wood3_2_C_litter1,
            "r_C_wood4_2_C_litter2": epa.r_C_wood4_2_C_litter2,
            "r_C_leaf_2_C_litter3": epa.r_C_leaf_2_C_litter3,
            "r_C_leaf_2_C_litter5": epa.r_C_leaf_2_C_litter5,
            "r_C_root_2_C_litter4": epa.r_C_root_2_C_litter4,
            "r_C_root_2_C_litter6": epa.r_C_root_2_C_litter6,
            "r_C_fruit_2_C_litter3": epa.r_C_fruit_2_C_litter3,
            "r_C_fruit_2_C_litter5": epa.r_C_fruit_2_C_litter5,
            "r_C_litter1_2_C_som1": epa.r_C_litter1_2_C_som1,
            "r_C_litter1_2_C_som2": epa.r_C_litter1_2_C_som2,
            "r_C_litter2_2_C_som2": epa.r_C_litter2_2_C_som2,
            "r_C_litter2_2_C_som3": epa.r_C_litter2_2_C_som3,
            "r_C_litter3_2_C_som1": epa.r_C_litter3_2_C_som1,
            "r_C_litter3_2_C_som3": epa.r_C_litter3_2_C_som3,
            "r_C_litter4_2_C_som1": epa.r_C_litter4_2_C_som1,
            "r_C_litter4_2_C_som2": epa.r_C_litter4_2_C_som2,
            "r_C_litter5_2_C_som1": epa.r_C_litter5_2_C_som1,
            "r_C_litter6_2_C_som2": epa.r_C_litter6_2_C_som2,
            "r_C_som1_2_C_som3": epa.r_C_som1_2_C_som3,
            "r_C_som2_2_C_som3": epa.r_C_som2_2_C_som3,
            "r_C_som2_2_C_som4": epa.r_C_som2_2_C_som4,
            "r_C_som3_2_C_som2": epa.r_C_som3_2_C_som2,
            "r_C_som3_2_C_som4": epa.r_C_som3_2_C_som4,
            "r_C_som4_2_C_som2": epa.r_C_som4_2_C_som2,
            "r_C_som4_2_C_som3": epa.r_C_som4_2_C_som3,
        }.items()
    }
    X_0_dict = {
        "C_wood1": epa.C_wood1_0,
        "C_wood2": epa.C_wood2_0,
        "C_wood3": epa.C_wood3_0,
        "C_wood4": epa.C_wood4_0,
        "C_leaf": epa.C_leaf_0,
        "C_root": epa.C_root_0,
        "C_fruit": cpa.cVeg_0 - (
                    epa.C_wood1_0 + epa.C_wood2_0 + epa.C_wood3_0 + epa.C_wood4_0 + epa.C_leaf_0 + epa.C_root_0),
        "C_litter1": epa.C_litter1_0,
        "C_litter2": epa.C_litter2_0,
        "C_litter3": epa.C_litter3_0,
        "C_litter4": epa.C_litter4_0,
        "C_litter5": epa.C_litter5_0,
        "C_litter6": cpa.cLitter_0 - (
                    epa.C_litter1_0 + epa.C_litter2_0 + epa.C_litter3_0 + epa.C_litter4_0 + epa.C_litter5_0),
        "C_som1": epa.C_som1_0,
        "C_som2": epa.C_som2_0,
        "C_som3": epa.C_som3_0,
        "C_som4": cpa.cSoil_0 - (epa.C_som1_0 + epa.C_som2_0 + epa.C_som3_0),
    }
    X_0 = np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(17, 1)
    fd = make_func_dict(mvs, dvs)
    V_init = gh.make_InitialStartVectorTrace(
        X_0, mvs,
        par_dict=par_dict,
        func_dict=fd
    )
    it_sym_trace = gh.make_daily_iterator_sym_trace(
        mvs,
        V_init=V_init,
        par_dict=par_dict,
        func_dict=fd
    )
    return it_sym_trace
