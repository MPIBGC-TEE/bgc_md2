import pathlib
from typing import Callable
import netCDF4 as nc
import numpy as np
from collections import namedtuple
from functools import reduce
from general_helpers import day_2_month_index, month_2_day_index, months_by_day_arr, TimeStepIterator2, respiration_from_compartmental_matrix

# fixme:
# Your parameters will most likely differ but you can still use the
# destinctions between different sets of parameters. The aim is to make
# the different tasks in the code more obvious. In principal you could
# have a lot of overlapping sets and just have to keep them consistent. 
# 'namedtuples' can be used just as normal tuples by functions
# that are not aware of the names. They can still use the positions like 
# in the original code

# This set is used by the functions that produce the
# specific ingredients (functions) that will be run by
# mcmc alg.
UnEstimatedParameters = namedtuple(
    "UnEstimatedParameters",
    [
        'C_soil_0',
        'C_veg_0',
        'gpp_0',
        'rh_0',
        'npp_0',
        'f_veg2soil_0',
        'npp',
        'number_of_months',
        'mrso',
        'tsl',
        'tas',
        'vegcover'
    ]
)

# This set is used (as argument) by the functions that are called
# inside the mcmc
EstimatedParameters = namedtuple(
    "EstiamatedParameters",
    [
        "beta_leaf",    # 0 (indices uses in original code)
        "beta_wood",    # 1
        "f_leaf2DPM",   # 2, f41
        # "f_leaf2RPM",   #  f51
        "f_wood2DPM",   # 3, f42
        # "f_wood2RPM",   #  f52
        "f_root2DPM",   # 4, f43
        # "f_root2RPM",   #  f53
        "f_DPM2BIO",    # 5, f64
        "f_DPM2HUM",    # 6, f74
        "f_RPM2BIO",    # 7, f65
        "f_RPM2HUM",    # 8, f75
        "f_BIO2HUM",    # 9, f76
        "f_HUM2BIO",    # 10, f67
        "k_leaf",       # 11
        "k_wood",       # 12
        "k_root",       # 13
        "k_DPM",	    # 14
        "k_RPM",	    # 15
        "k_BIO",        # 16
        "k_HUM",        # 17
        "c_leaf0",      # 18
        "c_wood0",      # 19
        "c_RPM0",       # 20
        "c_DPM0",       # 21
        "c_BIO0",       # 22
        'Mw'            # 23
    ]
)

# This is the set off all 
_Parameters = namedtuple(
        "Parameters",
        EstimatedParameters._fields + UnEstimatedParameters._fields
)
class Parameters(_Parameters):
    @classmethod
    def from_EstimatedParametersAndUnEstimatedParameters(
            cls,
            epa :EstimatedParameters,
            cpa :UnEstimatedParameters
        ):
        return cls(*(epa + cpa))


# This set defines the order of the c pools
# The order is crucial for the compatibility
# with the matrices (B and b) If you change ist
# the matrix changes
StateVariables = namedtuple(
    'StateVariables',
    [
        'C_leaf',
        'C_wood',
        'C_root',
        'C_DPM',
        'C_RPM',
        'C_BIO',
        'C_HUM'
    ]
)
Observables = namedtuple(
    'Observables',
    [
        'C_veg',
        'C_soil',
        'rh',
        'f_veg2soil'
    ]
)

# We define another set of parameters which describes
# the parameters of the matrices A,K and the vector b
# and drivers like npp (in form of arrays)
# but does not include start values and hyperparameters like the 'number_of_months'
# This distinction is helpful for the forward simulation where the
# distinction between estimated and constant is irrelevant.
ModelParameters = namedtuple(
    "ModelParameters",
    [
        name for name in Parameters._fields
        if name not in [
            'C_leaf_0',
            'C_wood_0',
            'C_root_0',
            'C_litter_above_0',
            'C_litter_below_0'
            'C_fast_som_0',
            'C_slow_som_0',
            'C_pass_som_0',
            'C_leaf_lit_0',
            'rh_0',
            'f_veg2soil_0',
            'number_of_months'
        ]
    ]
)
# to do
# 1.) make a namedtuple for the yycable data or use xarray to create a multifile dataset
def get_variables_from_files(dataPath):
    ## YZhou: Modified for JULES ##
    # Read NetCDF data. Monthly data, 1700-2019
    # Carbon in Soil (including below-ground litter), 192 (lat) x144 (lon) x3840 (time), kg m-2
    path = dataPath.joinpath("JULES-ES-1p0_S2_cSoil.nc")
    ds = nc.Dataset(str(path))
    var_csoil = ds.variables['cSoil'][:, :, :]
    ds.close()
    # Carbon in Vegetation, 192x144x3840, kg m-2
    path = dataPath.joinpath("JULES-ES-1p0_S2_cVeg.nc")
    ds = nc.Dataset(str(path))
    var_cveg = ds.variables['cVeg'][:, :, :]
    ds.close()
    # Total carbon mass from vegetation directly into the soil, 192x144x3840, kg m-2 s-1
    path = dataPath.joinpath("JULES-ES-1p0_S2_fVegSoil.nc")
    ds = nc.Dataset(str(path))
    var_fvegsoil = ds.variables['fVegSoil'][:, :, :]
    ds.close()
    # Total Soil Moisture Content, 192x144x3840, kg m-2
    path = dataPath.joinpath("JULES-ES-1p0_S2_mrso.nc")
    ds = nc.Dataset(str(path))
    var_mrso= ds.variables['mrso'][:, :, :]
    ds.close()
    # Gross primary production, 192x144x3840, kg m-2 s-1
    path = dataPath.joinpath("JULES-ES-1p0_S2_gpp.nc")
    ds = nc.Dataset(str(path))
    var_gpp = ds.variables['gpp'][:, :, :]
    ds.close()
    # Net Primary Production on Land as Carbon Mass Flux, 192x144x3840, kg m-2 s-1
    path = dataPath.joinpath("JULES-ES-1p0_S2_npp.nc")
    ds = nc.Dataset(str(path))
    var_npp = ds.variables['npp_nlim'][:, :, :]
    # note that NPP is N-limited in S2, there is another file 'JULES-ES-1p0_S2_npp_noNlimitation.nc' with var 'npp'
    ds.close()
    # Heterotrophic Respiration, 192x144x3840, kg m-2 s-1
    path = dataPath.joinpath("JULES-ES-1p0_S2_rh.nc")
    ds = nc.Dataset(str(path))
    var_rh = ds.variables['rh'][:, :, :]
    ds.close()
    # Temperature of Soil - layer, 192x144x4 (depth) x3840, 'K'
    path = dataPath.joinpath("JULES-ES-1p0_S2_tsl.nc")
    ds = nc.Dataset(str(path))
    var_tsl = ds.variables['tsl'][:, :, :, :] # YZhou: check if this works for a 4-D netcdf
    ds.close()
    # Near-Surface Air Temperature (1.5m from JULES metadata), 192x144x3840, 'K'
    path = dataPath.joinpath("JULES-ES-1p0_S2_tas.nc")
    ds = nc.Dataset(str(path))
    var_tas = ds.variables['tas'][:, :, :]
    # Plant Functional Type Grid Fraction, 192x144x17 (vegtype) x3840
    # vegtype= '0.BdlDcd; 1.BdlEvgTrop; 2.BdlEvgTemp; 3.NdlDcd; 4.NdlEvg; 5.c3grass; 6.c3crop;
    # 7.c3pasture; 8.c4grass; 9.c4crop; 10.c4pasture; 11.shrubDcd; 12.shrubEvg; 13.urban; 14.lake; 15.soil; 16.ice'
    path = dataPath.joinpath("JULES-ES-1p0_S2_landCoverFrac.nc")
    ds = nc.Dataset(str(path))
    var_vegcover = ds.variables['landCoverFrac'][:, :, :, :] # YZhou: check if this works for a 4-D netcdf

    lat_data = ds.variables['latitude'][:].data
    lon_data = ds.variables['longitude'][:].data
    ds.close()

    return (var_csoil, var_cveg, var_fvegsoil, var_mrso, var_gpp, var_npp, var_rh, var_tsl, var_tas, var_vegcover, lat_data, lon_data)

def get_example_site_vars(dataPath):
    ( ## YZhou: Modified for JULES ##
        var_csoil, var_cveg, var_fvegsoil, var_mrso, var_gpp, var_npp, var_rh, var_tsl, var_tas, var_vegcover, lat_data, lon_data
    ) = get_variables_from_files(dataPath)
    # pick up 1 site   wombat state forest

    s = slice(None, None, None)  # this is the same as :
    t = s, 50, 33  # [t] = [:,49,325]
    # YZhou: check the factor here. JULES is monthly dataset!!!
    gpp = var_gpp[t] * 86400  # from kg m-2 s-1 to kg m-2 day-1
    npp = var_npp[t] * 86400   # from kg m-2 s-1 to kg m-2 day-1
    rh = var_rh[t]*86400  # from kg m-2 s-1 to kg m-2 day-1
    f_veg2soil = var_fvegsoil[t] * 86400  # from kg m-2 s-1 to kg m-2 day-1
    tsl_mean = np.mean(var_tsl, axis=1)  # average soil temperature at different depth
    vegcover = np.sum(var_vegcover[:, :, 0:13, :], axis = 2)  # sum the vegetation coverages
    (
        C_soil,
        C_veg,
        mrso,
        tsl,
        tas
    ) = map(
        lambda var: var[t],
        (
            var_csoil,
            var_cveg,
            var_mrso,
            tsl_mean,
            var_tas
        )
    )
    return (C_soil, C_veg, f_veg2soil, mrso, gpp, npp, rh, tsl, tas, vegcover)
# get global sum of variables - modification of make_global_average.py from Jon model.


# YZhou: DO NOT USE, haven't tested this function yet!!
# def get_global_sum_vars(dataPath):
#     (  ## YZhou: Modified for JULES ##
#         var_csoil, var_cveg, var_fvegsoil, var_mrso, var_gpp, var_npp, var_rh, var_tsl, var_tas, lat_data, lon_data
#     ) = get_variables_from_files(dataPath)
#
#     # We assume that lat,lon actually specifies the center of the gridcell
#     # and that it extends from lat-0.5 to lat+0.5 and long-0.5
#     for v in ('theta', 'phi', 'delta_theta', 'delta_phi', 'r'):
#         var(v)
#
#     # we compute the are of a delta_phi * delta_theta patch
#     # on the unit ball (which depends also on theta but not
#     # on phi)
#     A_sym = integrate(
#         integrate(
#             sin(theta),
#             (theta, theta - delta_theta / 2, theta + delta_theta / 2)
#         ),
#         (phi, phi - delta_phi / 2, phi + delta_phi / 2)
#     )
#
#     A_fun=lambdify((theta,delta_theta,delta_phi),A_sym)
#
#     def grad2rad(alpha):
#         return np.pi / 180.0 * alpha
#
#     def area(lat):
#         delta_phi = 1
#         delta_theta = 1
#         r = 6378.1370  # km
#         # translate to normal spherical coordinates
#         theta = (90.0 - lat)
#         theta, delta_theta, delta_phi = map(
#             grad2rad,
#             (
#                 theta,
#                 delta_theta,
#                 delta_phi,
#             )
#         )
#         return A_fun(theta, delta_theta, delta_phi) * r
#
#     weights = np.array(list(map(area, var_lat)))
#
#     npp_lon_sum = var_npp.sum(axis=2)*86400 #   kg/m2/s kg/m2/day;
#     npp = (npp_lon_sum*weights).sum(axis=1);
#     print("yo there - the data type is of npp is", type(npp))
#     print("global npp calculated")
#     print("data type is ", type(npp))
#
#     rh_lon_sum = var_rh.sum(axis=2)*86400;   # per s to per day
#     rh = (rh_lon_sum*weights).sum(axis=1);
#     print("global rh calculated")
#
#     clitter_lon_sim = var_clitter.sum(axis=2);
#     clitter = (clitter_lon_sim*weights).sum(axis=1);
#     print("global clitter calculated")
#
#     csoil_lon_sum = var_csoil.sum(axis=2);
#     csoil = (csoil_lon_sum*weights).sum(axis=1);
#     print("global csoil calculated")
#
#     cveg_lon_sum = var_cveg.sum(axis=2);
#     cveg = (cveg_lon_sum * weights).sum(axis=1);
#     print("global cveg calculated")
#
#     cleaf_lon_sum = var_cleaf.sum(axis=2);
#     cleaf= (cleaf_lon_sum * weights).sum(axis=1);
#     print("global cleaf calculated")
#
#     croot_lon_sum = var_croot.sum(axis=2);
#     croot = (croot_lon_sum * weights).sum(axis=1);
#     print("global croot calculated")
#
#     ccwd_lon_sum = var_ccwd.sum(axis=2);
#     ccwd = (ccwd_lon_sum * weights).sum(axis=1);
#     print("global ccwd calculated")
#
#     cwood = cveg - cleaf - croot;
#     print("yo there - the data type is of cwood is", type(cwood))
#
#     return (npp, rh, clitter, csoil, cveg, cleaf, croot, ccwd, cwood)


def make_param_filter_func(
        c_max: np.ndarray,
        c_min: np.ndarray
        ) -> Callable[[np.ndarray], bool]:

    def isQualified(c):
        # fixme
        #   this function is model specific: It discards parameter proposals
        #   where beta1 and beta2 are >0.99
        cond1 =  (c >= c_min).all()
        cond2 =  (c <= c_max).all()
        cond3 = (c[0] + c[1]) < 1
        cond4 = (c[2] + c[3] + c[4]) < 1
        cond5 = (c[5] + c[6] + c[7]) < 1
        cond6 = (c[8] + c[9] + c[10]) < 1
        return (cond1 and cond2 and cond3 and cond4 and cond5 and cond6)


    return isQualified

# def make_weighted_cost_func(
#         obs: np.ndarray
#     ) -> Callable[[np.ndarray],np.float64]:
#     # first unpack the observation array into its parts
#     cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
#     def costfunction(out_simu: np.ndarray) ->np.float64:
#         # fixme
#         #   as indicated by the fact that the function lives in this
#         #   model-specific module it is not apropriate for (all) other models.
#         #   There are model specific properties:
#         #   1.) The weight for the different observation streams
#         #
#         tot_len=out_simu.shape[0]
#         # we assume the model output to be in the same shape and order
#         # as the obeservation
#         # this convention has to be honored by the forwar_simulation as well
#         # which in this instance already compresses the 3 different litter pools
#         # to c_litter and the 3 different soil pools to one
#         c_simu = out_simu[:,0:5]
#
#         # we assume the rh  part to be in the remaining columns again
#         # this convention has to be honored by the forwar_simulation as well
#         rh_simu = out_simu[:,5:]
#         #from IPython import embed; embed()
#
#         J_obj1 = np.mean (( c_simu[:,0] - cleaf[0:tot_len] )**2)/(2*np.var(cleaf[0:tot_len]))
#         J_obj2 = np.mean (( c_simu[:,1] - croot[0:tot_len] )**2)/(2*np.var(croot[0:tot_len]))
#         J_obj3 = np.mean (( c_simu[:,2] - cwood[0:tot_len] )**2)/(2*np.var(cwood[0:tot_len]))
#         J_obj4 = np.mean (( c_simu[:,3]-  clitter[0:tot_len] )**2)/(2*np.var(clitter[0:tot_len]))
#         J_obj5 = np.mean (( c_simu[:,4]-  csoil[0:tot_len] )**2)/(2*np.var(csoil[0:tot_len]))
#
#         J_obj6 = np.mean (( rh_simu[:,0] - rh[0:tot_len] )**2)/(2*np.var(rh[0:tot_len]))
#
#         J_new = (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 )/200+ J_obj6/4
#         return J_new
#     return costfunction


def make_param2res(
        cpa
    ) -> Callable[[np.ndarray], np.ndarray]:
    """The returned function 'param2res' is used by the metropolis hastings algorithm
    to produce a result from a proposed parameter(tuple).
    'param2res' will be called only with the parameters to be estimated 'pe'.
    The constant parameters are the parameters of the present function 'pc' and are automatically
    available in the returned 'param2res' (closure).
    
    In the case of this model 'param2res' has to perform the following 
    tasks.
    -   use both sets of parameters 'pc' and 'pe' to build a forward model
    -   run the forward model to produce output for the times 
        when observations are available.(Here monthly)
    -   project the output of the forward model to the observed variables.
        (here summation of all different soil-pools to compare to the single
        observed soil-pool and the same for the litter pools)
        
    In this version of the function all tasks are performed at once as in 
    the original script.
    See the alternative implementation to see how all the tasks can be 
    delegated to smaller functions.
    """
    # fixme
    # This function is model-specific in several ways:
    # 0. Which are the fixed and which are the estimated variables
    # 1. The matrices for the forward simulation 
    # 2. the driver (here monthly npp)
    # 3. the projection of the simulated pool values to observable variables
    #    (here summation of   
    def param2res(pa):
        # pa is a numpy array when pa comes from the predictor
        # so we transform it to be able to use names instead of positions 
        epa = EstimatedParameters(*pa)
        days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        # Construct b vector
        beta1 = epa.beta_leaf
        beta2 = epa.beta_wood
        beta3 = 1 - beta1 - beta2
        b = np.array([beta1, beta2, beta3, 0, 0, 0, 0]).reshape([7, 1])   # allocation

        # Construct A matrix
        f41 = epa.f_leaf2DPM
        f51 = 1 - epa.f_leaf2DPM  # start with NPP as input, no autotrophic respiration here
        f42 = epa.f_wood2DPM
        f52 = 1 - epa.f_wood2DPM
        f43 = epa.f_root2DPM
        f53 = 1 - epa.f_root2DPM
        f64 = epa.f_DPM2BIO
        f74 = epa.f_DPM2HUM
        f65 = epa.f_RPM2BIO
        f75 = epa.f_RPM2HUM
        f76 = epa.f_BIO2HUM
        f67 = epa.f_HUM2BIO
        A = np.array([  -1,   0,   0,   0,   0,   0,   0,
                         0,  -1,   0,   0,   0,   0,   0,
                         0,   0,  -1,   0,   0,   0,   0,
                       f41, f42, f43,  -1,   0,   0,   0,
                       f51, f52, f53,   0,  -1,   0,   0,
                         0,   0,   0, f64, f65,  -1, f67,
                         0,   0,   0, f74, f75, f76,  -1]).reshape([7, 7])   # transfer

        # Construct K matrix
        # #turnover rate per day of pools:
        temp = [epa.k_leaf, epa.k_wood, epa.k_root, epa.k_DPM, epa.k_RPM, epa.k_BIO, epa.k_HUM]
        # May change to the default values of soil K (s-1):
        # DPM: 3.22 * 10^-7; RPM: 9.65 * 10^-9; BIO: 2.12 * 10^-8; HUM: 6.43 * 10^-10
        K = np.zeros(49).reshape([7, 7])
        for i in range(0, 7):
            K[i][i] = temp[i]

        # initialize X
        # leaf, root , wood,
        # DPM (decomposable plant material), RPM (resistant plant material)
        # BIO (microbial biomass), and HUM (long-lived humified)
        x_init = np.array([
                epa.c_leaf0,
                epa.c_wood0,
                cpa.c_veg - epa.c_leaf0 - epa.c_wood0,
                epa.c_RPM0,
                epa.c_DPM0,
                epa.c_BIO0,
                cpa.c_soil - epa.c_RPM0 - epa.c_DPM0 - epa.c_BIO0]).reshape([7, 1])   # Initial carbon pool size
        # initialize carbon pools 
        X = x_init

        # create output matrices
        x_fin = np.zeros((cpa.number_of_months, 7))
        rh_fin = np.zeros((cpa.number_of_months, 1))
        f_veg2soil_fin = np.zeros((cpa.number_of_months, 1))

        # initialize first respiration value 
        co2_rh = cpa.rh_0
        f_veg2soil = cpa.f_veg2soil_0

        # fixme:
        # slight change to the original
        # I would like to start the solution with the initial values
        # m=0 means after 0 months = in the initial step
        #B=A@K
        #pa=Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)
        #B=make_compartmental_matrix_func(pa)(0,X)

        def Rh_calc (TS, V, Mw, M):
            FT = 2.0 ** ((TS - 298.15) / 10)  # temperature rate modifier
            FV = 0.6 + 0.4 * (1 - V / 100)  # effect of vegetation cover
            # Mw is soil moisture at wilting point as a fraction of saturation
            S0 = 0.5* (1 + Mw)  # optimum soil moisture
            Smin = 1.7 * Mw  # lower threshold soil moisture for soil respiration
            if S0 > Smin:
                FS = 1 - 0.8 * (M - S0)  # effect of soil moisture
            if (Smin < M) and (M <= S0):
                FS = 0.2 + 0.8 * (M - Smin) / (S0 - Smin)
            if M <= Smin:
                FS = 0.2
            rh_factor = FT * FV * FS
            return rh_factor

        for m in range(0, cpa.number_of_months):
            x_fin[m, :] = X.reshape(1, 7)
            npp_in = cpa.npp[m]
            rh_fin[m, 0] = co2_rh
            f_veg2soil_fin[m, 0] = f_veg2soil
            co2_rh = 0
            f_veg2soil = 0
            rh_modifier = Rh_calc(cpa.tsl[m], cpa.vegcover[m], epa.Mw, cpa.mrso[m])
            ksi = np.array(
                [
                    1,
                    1,
                    1,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier
                ]
            ).reshape([7, 1])  # environmental modifiers

            for d in range(0, days[m % 12]):
                K_new = K * ksi
                X = X + b * npp_in + A @ K_new @ X
                co2_rate = [0, 0, 0,
                        (1 - f64 - f74) * K[3, 3] * ksi[3],
                        (1 - f65 - f75) * K[4, 4] * ksi[4],
                        (1 - f76)       * K[5, 5] * ksi[5],
                        (1 - f67)       * K[6, 6] * ksi[6]
                    ]
                co2 = np.sum(co2_rate * X.reshape(1, 9))
                co2_rh = co2_rh + co2 / days[m % 12]   # monthly average rh

                veg2soil_rate = [(f41 + f51) * K[0, 0] * ksi[0],
                                 (f42 + f52) * K[1, 1] * ksi[1],
                                 (f43 + f53) * K[2, 2] * ksi[2],
                                 0, 0, 0, 0]
                veg2soil = np.sum(veg2soil_rate * X.reshape(1, 7))
                f_veg2soil = f_veg2soil + veg2soil / days[m % 12]  # monthly average flux from veg to litter

        # We create an output that has the same shape
        # as the obvervations to make the costfunctions 
        # easier. 
        # To this end we project our 10 output variables of the matrix simulation
        # onto the 6 data streams by summing up the 3 litter pools into one
        # and also the 3 soil pools into one
        #C_litter_above = x_fin[:,3] + x_fin[:,4]
        # C_litter_above = np.sum(x_fin[:, 3:5], axis=1).reshape(cpa.number_of_months,1)
        #C_litter_below = x_fin[:,5]
        #c_litter = np.sum(x_fin[:,3:6],axis=1).reshape(cpa.number_of_months,1)
        #c_soil = np.sum(x_fin[:,6:9],axis=1).reshape(cpa.number_of_months,1)
        #from IPython import embed; embed()
        out_simu = np.concatenate(
            [
                x_fin[:, 0:3],  # Cveg, the first 3 pools are used as they are
                #C_litter_below,
                x_fin[:, 5:7],  # Csoil
                rh_fin,
                f_veg2soil_fin
            ]
            , axis=1
        )
        return out_simu

    return param2res

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# # alternative implementation and helper functions. If you do not want to use it # comment it.
# def make_param2res_2(
#         cpa: UnEstimatedParameters
#     ) -> Callable[[np.ndarray], np.ndarray]:
#     # This function is an alternative implementation with the same results as
#     # make_param2res
#     # Internally it uses a slightly higher level of abstraction and devides the task
#     # into more numerous but smaller testable and reusable units.
#     # Although both properties are not necessary to replicate the original
#     # code they provide a different perspective which is easier to reuse the code
#     # in  a wider range of models e.g. models with daily driver data, time
#     # dependent or even nonlinear matrices and so on.
#     # Thus it allows to test small parts of the code against the respective parts of the
#     # symbolic framework
#     #
#     # It differs in the following respects:
#     # 0.)   We use named tuples for the parameters (just for readability)
#     # 1.)   We use a seperate function to costruct the matrix which will be
#     #       a function B(i,X)  of the iteration it  and the present poolvalues X
#     #       (althouhg in this case unnecessary since B is constant but
#     #       the dependence on 'i' becomes useful if you want to include
#     #       an environmental scalar and the dependence on X makes it possible
#     #       to include nonlinearities.
#     # 2.)   We build a daily advancing model that can provide output for an arbitrary
#     #       selection of days.  To do so we provide all driver data as
#     #       functions of the smalles timestep (here day with index i), which in
#     #       the case of this model means that we provide the same npp value for
#     #       all the days in a given month.
#     # 3.)   We translate the index of a given month to the appropriate day index
#     #       and apply the dayly model of 2.). Again this seems cumbersome for this
#     #       example but allows us to reuse the daily model for all applications.
#     #       This is espacially usefull for testing since we only need some dayly timesteps.
#     #       It makes it also easier to keep an overview over the appropriate
#     #       units: If the smallest timestep is a day, then all time related parameters
#     #       have to be provided in the corresponding  units, regardless of the
#     #       number of available datapoints
#     #
#     def param2res(pa):
#         epa=EstimatedParameters(*pa)
#         # here we want to store only monthly values
#         # although the computation takes place on a daily basis
#         V_init = construct_V0(cpa,epa)
#         # compute the days when we need the results
#         # to be able to compare to the monthly output
#         day_indices = month_2_day_index(range(cpa.number_of_months))
#         print("day_indices=",day_indices)
#         apa = Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)
#
#         mpa = ModelParameters(
#             **{
#                 k:v for k,v in apa._asdict().items()
#                 if k in ModelParameters._fields
#             }
#         )
#         full_res = run_forward_simulation(
#             V_init=V_init,
#             day_indices=day_indices,
#             mpa=mpa
#         )
#
#         # project the litter and soil pools
#         tot_len=cpa.number_of_months
#         c_litter = np.sum(full_res[:,3:6],axis=1).reshape(tot_len,1)
#         c_soil = np.sum(full_res[:,6:9],axis=1).reshape(tot_len,1)
#         out_simu = np.concatenate(
#             [
#                 full_res[:,0:3], # the first 3 pools are used as they are
#                 c_litter,
#                 c_soil,
#                 full_res[:,9:10]
#             ]
#             ,axis=1
#         )
#         return out_simu
#
#     return param2res
#

def run_forward_simulation(
        V_init,
        day_indices,
        mpa
    ):
        tsi = make_daily_iterator(
            V_init,
            mpa=mpa
        )

        def g(acc, i):
            xs,co2s,acc_co2,acc_days = acc
            v = tsi.__next__()
            d_pools = v[0:-1, :]
            d_co2 = v[-1:, :]
            acc_co2 += d_co2
            acc_days += 1
            if i in day_indices:
                xs += [d_pools]
                co2s +=[acc_co2/acc_days]
                acc_co2=np.array([0.0]).reshape(1,1)
                acc_days = 0

            acc = (xs,co2s,acc_co2,acc_days)
            return acc
        xs, co2s, acc_days, _ =  reduce(g, range(max(day_indices)+1), ([], [], 0, 0))

        def h(tup):
            x, co2 = tup
            return np.transpose(np.concatenate([x,co2]))

        values_with_accumulated_co2 = [v for v in map(h, zip(xs, co2s))]
        RES = np.concatenate(values_with_accumulated_co2, axis=0)
        return RES

def make_daily_iterator(
        V_init,
        mpa
    ):

        # Construct npp(day)
        # in general this function can depend on the day i and the state_vector X
        # e.g. typically the size fo X.leaf...
        # In this case it only depends on the day i
        def npp_func(day,X):
            return mpa.npp[day_2_month_index(day)]

        # b (b vector for partial allocation) 
        beta_root = 1- mpa.beta_leaf- mpa.beta_wood
        b = np.array(
            [
                mpa.beta_leaf,
                mpa.beta_wood,
                beta_root,
                0,
                0,
                0,
                0
            ],
        ).reshape(7, 1)
        # Now construct B matrix B=A*K
        B_func  = make_compartmental_matrix_func(
            mpa=mpa,
        )
        # Build the iterator which is the representation of a dynamical system
        # for looping forward over the results of a difference equation
        # X_{i+1}=f(X_{i},i)
        # This is the discrete counterpart to the initial value problem
        # of the continuous ODE
        # dX/dt = f(X,t) and initial values X(t=0)=X_0

        def f(it, V):
            X = V[0:7]
            co2 = V[7]
            npp = npp_func(it, X)
            B = B_func(it, X)
            X_new = X + npp * b + B@X

            # we also compute the respired co2 in every (daily) timestep
            # and use this part of the solution later to sum up the monthly amount
            co2_new = respiration_from_compartmental_matrix(B, X)

            V_new = np.concatenate((X_new, np.array([co2_new]).reshape(1, 1)), axis=0)
            return V_new


        #X_0 = np.array(x_init).reshape(9,1) #transform the tuple to a matrix
        #co2_0 = np.array([0]).reshape(1,1)
        #V_0 = np.concatenate((X_0, co2_0), axis=0)

        return TimeStepIterator2(
                initial_values=V_init,
                f=f,
                #max_it=max(day_indices)+1
        )


def make_compartmental_matrix_func(
        mpa
    ):
    # Now construct A matrix
    # diagonal 
    # make sure to use 1.0 instead of 1 otherwise it will be an interger array
    # and round your values....
    A =np.diag([-1.0 for i in range(7)])
    # because of python indices starting at 0 we have A[i-1,j-1]=fij
    #
    #A = np.array([  -1,   0,   0,   0,   0,   0,   0,   0,   0,
    #                 0,  -1,   0,   0,   0,   0,   0,   0,   0,
    #                 0,   0,  -1,   0,   0,   0,   0,   0,   0,
    #               f41, f42,   0,  -1,   0,   0,   0,   0,   0,
    #               f51, f52,   0,   0,  -1,   0,   0,   0,   0,
    #                 0,   0, f63,   0,   0,  -1,   0,   0,   0,
    #                 0,   0,   0, f74, f75,   0,  -1,   0,   0,
    #                 0,   0,   0,   0, f85, f86, f87,  -1,   0,
    #                 0,   0,   0,   0,   0, f96, f97, f98,  -1 ]).reshape([9,9])   # tranfer

    # note that the matrix description of a model is implicitly related to the
    # order of pools in the state_vector and not independent from it.

    A[3, 0] = mpa.f_leaf2DPM
    A[3, 1] = mpa.f_wood2DPM
    A[3, 2] = mpa.f_root2DPM
    A[4, 0] = 1.0 - mpa.f_leaf2DPM
    A[4, 1] = 1.0 - mpa.f_wood2DPM
    A[4, 2] = 1.0 - mpa.f_root2DPM
    A[5, 3] = mpa.f_DPM2BIO
    A[6, 3] = mpa.f_DPM2HUM
    A[5, 4] = mpa.f_RPM2BIO
    A[6, 4] = mpa.f_RPM2HUM
    A[6, 5] = mpa.f_BIO2HUM
    A[5, 6] = mpa.f_HUM2BIO


    #turnover rate per day of pools:
    K = np.diag([
        mpa.k_leaf,
        mpa.k_wood,
        mpa.k_root,
        mpa.k_DPM,
        mpa.k_RPM,
        mpa.k_BIO,
        mpa.k_HUM
    ])
    # in the general nonautonomous nonlinear case
    # B will depend on an it,X (althouh in this case it does not depend on either
    def B_func(it,X):
        return A @ K

    return B_func
#
# def construct_V0(
#         cpa :UnEstimatedParameters,
#         epa :EstimatedParameters
#     ) -> np.ndarray:
#     """Construct the initial values for the forward simulation
#     from constant and eveluated parameters
#
#     param: cpa : constant parameeters
#     param: epa : estimated parameters
#     """
#     # to make sure that we are in the right order we use the
#     # StateVariables namedtuple
#     X_0 = StateVariables(
#         C_leaf=cpa.C_leaf_0,
#         C_root=cpa.C_root_0,
#         C_wood=cpa.C_wood_0,
#         C_metlit=epa.C_metlit_0,
#         C_strlit=epa.C_strlit_0,
#         C_cwd=C_cwd_0,
#         C_mic=epa.C_mic_0,
#         C_slowsom=cpa.csoil_0- epa.C_mic_0 - epa.C_passom_0,
#         C_passsom=epa.C_passom_0
#     )
#     # add the respiration start value to the tuple
#     V_0 = (*X_0,cpa.rh_0)
#     return np.array(V_0).reshape(10,1)
