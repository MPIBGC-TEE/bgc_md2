import sys
# import json
from pathlib import Path
from collections import namedtuple
import netCDF4 as nc
import numpy as np
from sympy import Symbol
# from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (TimeStepIterator2)
# from copy import copy
from typing import Callable

# from general_helpers import month_2_day_index, monthly_to_yearly
# from functools import reduce

sys.path.insert(0, '..')  # necessary to import general_helpers
import general_helpers as gh

Observables = namedtuple(
    'Observables',
    ["cVeg", "cSoil", "rh", "fVegSoil"]
)
# Driver data streams on TRENDY server
Drivers = namedtuple(
    "Drivers",
    ["npp", "mrso", "tsl", "landCoverFrac"]
)

Constants = namedtuple(
    "Constants",
    [
        'npp_0',  # Initial input/pools
        'rh_0',
        'c_veg_0',
        'c_soil_0',
        'fVegSoil_0',  # Total carbon mass from vegetation directly into the soil
        'nyears'  # Run time (years for my model)
    ]
)

EstimatedParameters = namedtuple(
    'EstimatedParameters',
    [
        'c_leaf_0',  # Names: c_poolname_0
        'c_wood_0',  # Only initial pools that are estimated
        'c_DPM_0',
        'c_RPM_0',
        'c_BIO_0',

        'beta_leaf',
        'beta_wood',
        'Mw',
        'Ms',

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

def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),  # platform independent path desc. (Windows vs. linux)
        models=["JULES-ES"],
        variables=Observables._fields + Drivers._fields
    )

def get_example_site_vars(dataPath):
    # Define single geospatial cell from (3840, 144, 192)
    slayer = slice(None, None, None)  # this is the same as :
    slat = 70
    slon = 160
    t = slayer, slat, slon  # a site in South America

    # Define function to select geospatial cell and scale data
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them
        ds = nc.Dataset(str(path))

        if vn in ["npp_nlim", "gpp", "rh", "ra", "f_veg2soil"]:  # (3840, 144, 192), kg m-2 s-1
            # f_veg2soil: Total carbon mass from vegetation directly into the soil
            print("reading ", vn, ", size is ", ds.variables[vn].shape)
            return ds.variables[vn][t] * 86400  # convert from kg/m2/s to kg/m2/day
        elif vn in ["tsl"]:  # Temperature of Soil - layer, 192x144x4 (depth) x3840, 'K'
            print("reading ", vn, ", size is ", ds.variables[vn].shape)  ## (3840, 4, 144, 192)
            tmp = np.mean(ds.variables[vn][:, :, slat, slon], axis=1)
            return tmp  # average soil temperature at different depth
            print("converted size is ", tmp.shape)
        elif vn in ["landCoverFrac"]:  # Plant Functional Type Grid Fraction, 192x144x17 (vegtype) x3840
            print("reading ", vn, ", size is ", ds.variables[vn].shape)  ## (3840, 17, 144, 192)
            # from IPython import embed;embed()
            var = ds.variables[vn]
            sh = var.shape
            tmp = np.sum(var[:, 0:13, slat, slon], axis=1)
            print("converted size is ", tmp.shape)
            return tmp  # sum the vegetation coverages
        else:
            print("reading ", vn, ", size is ", ds.variables[vn].shape)
            return ds.variables[vn][t]

    # Link symbols and data:

    # Create file names (single step if files similarly named)

    file_name_from_var_name = {
       "npp_nlim":"JULES-ES-1p0_S2_npp.nc",
       **{
            vn: "JULES-ES-1p0_S2_{}.nc".format(vn) 
            for vn in [ "mrso", "tsl","landCoverFrac", "cVeg", "cSoil", "rh","fVegSoil" ]
       }
    }
    #d_name2varname_in_file = {
    #    "npp":'npp_nlim',
    #    **{
    #         vn: vn
    #         for vn in ["fVegSoil", "mrso", "tsl","landCoverFrac", "cVeg", "cSoil", "rh" ]
    #    "npp_nlim": "JULES-ES-1p0_S2_npp.nc",
    #    **{
    #        vn: "JULES-ES-1p0_S2_{}.nc".format(vn)
    #        for vn in ["mrso", "tsl", "landCoverFrac", "cVeg", "cSoil", "rh", "fVegSoil"]
    #    }
    #}
    d_name2varname_in_file = {
        "npp": 'npp_nlim',
        **{
            vn: vn
            for vn in ["fVegSoil", "mrso", "tsl", "landCoverFrac", "cVeg", "cSoil", "rh"]
        }

    }

    o_tuples = [
        (
            d_name2varname_in_file[f],
            file_name_from_var_name[d_name2varname_in_file[f]]
        )
        for f in Observables._fields
    ]

    d_tuples = [
        (
            d_name2varname_in_file[f],
            file_name_from_var_name[d_name2varname_in_file[f]]
        )
        for f in Drivers._fields
    ]

    # Link symbols and data for observables/drivers
    # print(o_tuples)
    return (
        Observables(*map(f, o_tuples)),
        Drivers(*map(f, d_tuples))
    )


def make_StartVector(mvs):
    # fixme mm:
    # add the fraction!
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()] +
        ["rh","fVegSoil"]
    )


def make_xi_func(tsl, Mw, Ms, mrso, landCoverFrac):
    def xi_func(day):
        mi = gh.day_2_month_index(day)
        # alternative FT
            # Q10 function
        FT = 2.0 ** ((tsl[mi] - 298.15) / 10)  # temperature rate modifier
            # RothC temperature function (Jenkinson 1990)
        #FT = 47.9 / (1 + np.exp(106/(tsl[mi] - 254.85)))
        FV = 0.6 + 0.4 * (1 - landCoverFrac[mi] / 100)  # effect of vegetation cover
        # Mw is soil moisture at wilting point as a fraction of saturation
        # Ms is soil moisture content at saturation
        S0 = 0.5 * (1 + Mw)  # optimum soil moisture
        Smin = 1.7 * Mw  # lower threshold soil moisture for soil respiration
        if S0 > Smin:
            FS = 1 - 0.8 * (mrso[mi]/Ms - S0)  # effect of soil moisture
        if (Smin < mrso[mi]/Ms) and (mrso[mi]/Ms <= S0):
            FS = 0.2 + 0.8 * (mrso[mi]/Ms - Smin) / (S0 - Smin)
        if mrso[mi] <= Smin:
            FS = 0.2
        # print("FT,FV,FS", FT, FV, FS)
        rh_factor = FT * FV * FS
        return rh_factor # 1.0     # Set to 1 if no scaler implemented
        # return 1.0

    return xi_func


def make_npp_func(dvs):
    def func(day):
        month = gh.day_2_month_index(day)
        # kg/m2/s kg/m2/day;
        return (dvs.npp[month])

    return func


def make_func_dict(mvs, dvs, cpa, epa):
    return {
        "NPP": make_npp_func(dvs),
        "xi": make_xi_func(dvs.tsl, epa.Mw, epa.Ms, dvs.mrso, dvs.landCoverFrac)
    }


def make_iterator_sym(
        mvs,
        V_init: "StartVector",
        par_dict,
        func_dict,
        delta_t_val=1  # defaults to 1day timestep
):
    B_func, u_func = gh.make_B_u_funcs_2(mvs, par_dict, func_dict, delta_t_val)
    sv = mvs.get_StateVariableTuple()
    n = len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr = np.array(
        [V_init.__getattribute__(str(v)) for v in sv] +
        [V_init.rh, V_init.fVegSoil]
    ).reshape(n + 2, 1)  # reshaping is neccessary for matmul (the @ in B @ X)

    # numOutFluxesBySymbol={
    #    k:numfunc(expr_cont,delta_t_val=delta_t_val)
    #    for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
    # }
    numOutFluxesBySymbol = {
        k: gh.numfunc(
            expr_cont, 
            mvs,
            delta_t_val=delta_t_val,
            par_dict=par_dict,
            func_dict=func_dict
        )
        for k, expr_cont in mvs.get_OutFluxesBySymbol().items()
    }
    v2sfl = {
        k:v for k,v in mvs.get_InternalFluxesBySymbol().items() 
        if k[0] in map(Symbol,["c_leaf","c_wood","c_root"])
    }
    numv2sfl = {
        k: gh.numfunc(
            expr_cont,
            mvs,delta_t_val=delta_t_val,
            par_dict=par_dict,
            func_dict=func_dict
        )
        for k, expr_cont in v2sfl.items()
    }

    def f(it, V):
        X = V[0:n]
        b = u_func(it, X)
        B = B_func(it, X)
        X_new = X + b + B @ X
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep

        rh=np.array([
            f(it, *X.reshape(n, ))
            for k,f in numOutFluxesBySymbol.items()
            if str(k) in ["c_DPM", "c_RPM", "c_BIO", "c_HUM"]
        ]).sum()
        fVegSoil=np.array(
            [
                f(it, *(X.reshape(n, )))
                for f in numv2sfl.values()
            ]
        ).sum()


        V_new = np.concatenate(
            (
                X_new.reshape(n,1),
                np.array([rh]).reshape(1,1),
                np.array([fVegSoil]).reshape(1,1)
            )
            , axis=0
        )
        return V_new
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )

def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers,
) -> Callable[[np.ndarray], np.ndarray]:
    # Build iterator
    # Need dictionary of numeric values for all parameters that are not state variables/time
    srm = mvs.get_SmoothReservoirModel()
    model_par_dict_keys = srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))] +
        list(mvs.get_StateVariableTuple())
    )
    StartVector = make_StartVector(mvs)

    # Time dependent driver function does not change with the estimated parameters
    # Defined once outside param2res function
    # seconds_per_day = 86400

    ########### Ask Markus here
    # Build environmental scaler function  ############### day or monthly, monthly inputs here, Mw and Ms are the to-be-estimated parameters

    # Define actual forward simulation function
    def param2res(pa):
        # Define function dictionary
        # Parameter vector
        epa = EstimatedParameters(*pa)

        func_dict = make_func_dict(mvs, dvs, cpa, epa)
        # Create a startvector for the iterator
        V_init = StartVector(
            c_leaf=epa.c_leaf_0,
            c_wood=epa.c_wood_0,
            c_root=cpa.c_veg_0 - (epa.c_leaf_0 + epa.c_wood_0),

            c_DPM=epa.c_DPM_0,
            c_RPM=epa.c_RPM_0,
            c_BIO=epa.c_BIO_0,
            c_HUM=cpa.c_soil_0 - (epa.c_DPM_0 + epa.c_RPM_0 + epa.c_BIO_0),

            rh=cpa.rh_0,
            fVegSoil=cpa.fVegSoil_0
        )

        # Parameter dictionary for the iterator
        apa = {**cpa._asdict(), **epa._asdict()}
        model_par_dict = {
            Symbol(k): v for k, v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }

        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2)
        delta_t_val = 15
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
        
        # fixme mm:
        # - We want monthly averages instead of the presently implemented last value of the month
        # - We want to compute the flux from Veg 2 Soil.
        sols = []
        dpm = 30  #
        n = len(V_init)
        for m in range(cpa.nyears*12):
            # dpm = days_per_month[ m % 12]
            mrh = 0
            for d in range(int(dpm / delta_t_val)):
                v = it_sym.__next__().reshape(n, )
                # actually the original datastream seems to be a flux per area (kg m^-2 ^-s )
                # at the moment the iterator also computes a flux but in kg^-2 ^day
            V = StartVector(*v)
            # from IPython import embed;embed()
            o = Observables(
                cVeg=float(V.c_leaf + V.c_wood + V.c_root),
                cSoil=float(V.c_DPM + V.c_RPM + V.c_BIO + V.c_HUM),
                rh=V.rh,
                fVegSoil=V.fVegSoil #Total carbon mass from vegetation directly into the soil
            )
            sols.append(o)

        # quick conversion of the list of observables into an  Observables tuple of arrays 
        def arr(field):
            n=len(sols)
            a=np.zeros(n)
            for i in range(n):
                a[i]=sols[i].__getattribute__(field)
            return a

        return Observables(*( arr(field) for field in Observables._fields))
        # sol = np.stack(sols)
        # convert to yearly output if necessary (the monthly pool data looks very "yearly")
        # sol_yr=np.zeros(int(cpa.number_of_months/12)*sol.shape[1]).reshape([int(cpa.number_of_months/12),sol.shape[1]])
        # for i in range(sol.shape[1]):
        #   sol_yr[:,i]=monthly_to_yearly(sol[:,i])
        # sol=sol_yr
        # return sol

    return param2res

def make_weighted_cost_func(
        obs: Observables
) -> Callable[[Observables], np.float64]:
    # first unpack the observation array into its parts
    # cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
    def costfunction(out_simu: Observables) -> np.float64:
        # fixme
        # loop over the fields and move to general_helpers

        J_obj1 = np.mean((out_simu.cVeg - obs.cVeg) ** 2) / (2 * np.var(obs.cVeg))
        J_obj2 = np.mean((out_simu.cSoil - obs.cSoil) ** 2) / (2 * np.var(obs.cSoil))

        J_obj3 = np.mean((out_simu.rh - obs.rh) ** 2) / (2 * np.var(obs.rh))
        J_obj4 = np.mean((out_simu.fVegSoil - obs.fVegSoil) ** 2) / (2 * np.var(obs.fVegSoil))

        J_new = (J_obj1 + J_obj2 + J_obj3 + J_obj4)
        return J_new

    return costfunction
