import pathlib
from typing import Callable
import netCDF4 as nc
import numpy as np
from collections import namedtuple
from functools import reduce
from general_helpers import day_2_month_index, month_2_day_index, months_by_day_arr, TimeStepIterator2, \
    respiration_from_compartmental_matrix

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
          'C_wood1_0',
          'C_wood2_0',
          'C_wood3_0',
          'C_wood4_0',
          'C_leaf_0',
          'C_root_0',
          'C_fruit_0',
          'C_litter1_0',
          'C_litter2_0',
          'C_litter3_0',
          'C_litter4_0',
          'C_litter5_0',
          'C_litter6_0',
          'C_surface_som_0',
          'C_fast_som_0',
          'C_slow_som_0',
          'C_pass_som_0',
          'rh_0',
          'ra_0',
        #   'gpp',
           'npp',
          'number_of_months',
          'mrso',
        #  'tsl1',
        #  'tsl2',
      #  'beta_fruit' #0.1

        #  'C_wood_0',
        #  'C_leaf_0',
        #  'C_root_0',
        #  'C_litter_0',
        #  'C_soil_0',
        #   'rh_0',
        #  'ra_0',
        #  'npp',
        #  'number_of_months',
        # 'mrso',
          'tsl',
        'beta_fruit'  # 0.1

    ]
)

# This set is used (as argument) by the functions that are called
# inside the mcmc
EstimatedParameters = namedtuple(
    "EstiamatedParameters",
    [

      #  "C_liable_0",
      #  "C_reserve_0",
      #  "C_leaf_0",
      #  "C_root_0",
      #  "C_sapwood1_0",
      #  "C_sapwood2_0",
      #  "C_heartwood1_0",
      #  "C_heartwood2_0",
      #  "C_fruit_0",
      #  "C_litter1_0",
      #  "C_litter2_0",
      #  "C_litter3_0",
      #  "C_litter4_0",
      #  "C_litter5_0",
      #  "C_litter6_0",
      #  "C_surface_som_0",
      #  "C_fast_som_0",
      #  "C_slow_som_0",
      #  "C_pass_som_0",

      "beta_sapwood1",  # 0 ()
      "beta_sapwood2",  # 1
      "beta_leaf",  # 2
      "beta_root",  # 3

      "f_sapwood1_heartwood1",
      "f_sapwood2_heartwood2",


      "f_wood1_liiter1",  # 2
      "f_wood2_liiter2",  # 2
      "f_leaf_liiter3",  # 2
      "f_root_liiter4",  # 2
      "f_fruit_liiter3",  # 2

      "f_litter1_surface_som",  # 2
      "f_litter1_fast_som",  # 2

      "f_litter2_fast_som",  # 2
      "f_litter2_slow_som",  # 2

      "f_litter3_surface_som",  # 2
      "f_litter3_slow_som",  # 2

      "f_litter4_surface_som",  # 2
      "f_litter4_fast_som",  # 2

      "f_litter5_surface_som",  # 2
      "f_litter6_fast_som",  # 2

      "f_surface_som_slow_som",  # 2

      #  "k_labile",  # 11
      "k_leaf",  # 11
      "k_wood1",  # 12
      "k_wood2",  # 12
      "k_wood3",  # 12
      "k_wood4",  # 12
      "k_root",  # 13
      "k_fruit",
      "k_lit1",  # 14
      "k_lit2",  # 15
      "k_lit3",  # 16        "k_leaf_lit",  # 14
      "k_lit4",  # 15
      "k_lit5",  # 16
      "k_lit6",  # 17
      "k_surface_som",  # 18
      "k_fast_som",  # 18
      "k_slow_som",  # 18
      "k_pass_som",  # 19
      "T_0",  # 21
      "E",  # 22
      "KM"  # 23
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
            epa: EstimatedParameters,
            cpa: UnEstimatedParameters
    ):
        return cls(*(epa + cpa))


# This set defines the order of the c pools
# The order is crucial for the compatibility
# with the matrices (B and b) If you change ist
# the matrix changes
StateVariables = namedtuple(
    'StateVariables',
    [
   #     'C_labile',
        'C_sap_wood1',
        'C_sap_wood2',
        'C_heart_wood1',
        'C_heart_wood2',
        'C_leaf',
        'C_root',
        'C_fruit'
        'C_litter1',
        'C_litter2',
        'C_litter3',
        'C_litter4',
        'C_litter5',
        'C_litter6',

        'C_surface_som',
        'C_fast_som',
        'C_slow_som',
        'C_pass_som'
    ]
)
Observables = namedtuple(
    'Observables',
    [
        'C_leaf',
        'C_wood',
        'C_root',
        'C_litter',
        'C_soil',
        'ra',
        'rh',
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
        'C_wood1_0',
        'C_wood2_0',
        'C_wood3_0',
        'C_wood4_0',
        'C_leaf_0',
        'C_root_0',
        'C_fruit_0',
        'C_litter1_0',
        'C_litter2_0',
        'C_litter3_0',
        'C_litter4_0',
        'C_litter5_0',
        'C_litter6_0',
        'C_surface_som_0',
        'C_fast_som_0',
        'C_slow_som_0',
        'C_pass_som_0',
        'number_of_months',
        'beta_fruit'  # 0.1
    ]
    ]
)


# to do
# 1.) make a namedtuple for the yycable data or use xarray to create a multifile dataset
def get_variables_from_files(dataPath):
    # Read NetCDF data  ******************************************************************************************************************************
    names = [
        ('cLeaf', 'cLeaf_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc'),
        ('cLitter', 'cLitter_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc'),
        ('cRoot', 'cRoot_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc'),
        ('cSoil', 'cSoil_Emon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc'),
        ('cVeg', 'cVeg_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc'),
        ('gpp', 'gpp_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc'),
        ('mrso', 'mrso_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc'),
        ('ra', 'ra_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc'),
        ('rh', 'rh_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc'),
        ('tsl', 'tsl_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-020012.nc'),
        ('tsl', 'tsl_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_020101-025012.nc'),
    ]

    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        ds = nc.Dataset(str(path))
        return ds.variables[vn][:, :, :]

    return map(f, names)


#     # Read NetCDF data  ******************************************************************************************************************************

def get_example_site_vars(dataPath):
    (
        C_leaf,
        C_litter,
        C_root,
        C_Soil,
        C_veg,
        gpp,
        mrso,
        ra,
        rh,
        tsl1,
        tsl2
    ) = get_variables_from_files(dataPath)
    # pick up 1 site   wombat state forest
    s = slice(None, None, None)  # this is the same as :
    t = s, 58, 159  # [t] = [:,49,325]
    gpp = gpp[t] * 86400  # kg/m2/s kg/m2/day;
    rh = rh[t] * 86400  # per s to per day
    ra = ra[t] * 86400  # per s to per day
 #   npp = gpp[t] - ra[t]

    tsl_mean1 = np.mean(tsl1, axis=1)  # average soil temperature at different depth

    tsl_mean2 = np.mean(tsl2, axis=1)  # average soil temperature at different depth
  #  print(tsl_mean2.shape)
    tsl_mean = np. append(tsl_mean1,tsl_mean2,axis=0)
  #  print(tsl_mean.shape)
    #print(C_leaf.shape)
    #print(C_leaf[t].shape)

    (
        C_leaf,
        C_litter,
        C_root,
        C_Soil,
        C_veg,
        mrso,
        tsl
    ) = map(
        lambda var: var[t],
        (
            C_leaf,
            C_litter,
            C_root,
            C_Soil,
            C_veg,
            mrso,
            tsl_mean
        )
    )
    C_wood = C_veg - C_leaf - C_root
    npp = gpp - ra
    return (npp, C_leaf, C_wood, C_root, C_litter, C_Soil, rh, ra, mrso, tsl)


def make_param_filter_func(
        c_max: np.ndarray,
        c_min: np.ndarray
) -> Callable[[np.ndarray], bool]:
    def isQualified(c):
        # fixme
        #   this function is model specific: It discards parameter proposals
        #   where beta1 and beta2 are >0.99
        cond1 = (c >= c_min).all()
        cond2 = (c <= c_max).all()
        cond3 = (c[0] + c[1] + c[2] + c[3]) <= 0.9
      #  cond4 = (c[2] + c[3] + c[4]) < 1
      #  cond5 = (c[5] + c[6] + c[7]) < 1
       # cond6 = (c[8] + c[9] + c[10]) < 1
        return (cond1 and cond2 and cond3 )#and cond4 and cond5 and cond6)

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
        beta1 = epa.beta_sapwood1;
        beta2 = epa.beta_sapwood2;
        beta3 = epa.beta_leaf;
        beta4 = 1 - beta1 - beta2 - beta3 - 0.1
        b = np.array([beta1, beta2, 0, 0, beta3, beta4, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).reshape([17, 1])  # allocation
        # Now construct A matrix


        f31 = epa.f_sapwood1_heartwood1;
        f81 = 1-epa.f_sapwood1_heartwood1;
        f42 = epa.f_sapwood2_heartwood2;
        f92 = 1-epa.f_sapwood2_heartwood2;
        f83 = 1;
        f94 = 1;
        f105 = 0.712;
        f125 = 0.288;
        f116 = epa.f_root_liiter4;
        f136 = 1-epa.f_root_liiter4;
        f107 = epa.f_fruit_liiter3;
        f127 = 1-epa.f_fruit_liiter3;
        f158 = epa.f_litter1_surface_som;
        f168 = epa.f_litter1_fast_som;
        f159 = epa.f_litter2_fast_som;
        f169 = epa.f_litter2_slow_som;
        f1410 = epa.f_litter3_surface_som;
        f1610 = epa.f_litter3_slow_som;
        f1511 = epa.f_litter4_fast_som;
        f1411 = epa.f_litter4_surface_som;
        f1412 = 0.4;
        f1513 = 0.45;
        f1614 = epa.f_surface_som_slow_som;
        f1615 = 0.296;
        f1715 = 0.04;
        f1516 = 0.42;
        f1716 = 0.03;
        f1517 = 0.45;


        A = np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      f31, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, f42, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      f81, 0, f83, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, f92, 0, f94, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, f105, 0, f107, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, f116, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, f125, 0, f127, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, f136, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, f1410, 0, f1412, 0, -1, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, f158, f159, 0, f1511, 0, f1513, 0, -1, f1516, f1517,
                      0, 0, 0, 0, 0, 0, 0, f168, f169, f1610, f1411, 0, 0, f1614, f1615, -1,0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, f1715, f1716, -1]).reshape([17, 17])  # transfer

        # turnover rate per day of pools:
        temp = [epa.k_leaf, epa.k_root, epa.k_wood1,epa.k_wood2,epa.k_wood3,epa.k_wood4, epa.k_fruit, epa.k_lit1,epa.k_lit2,epa.k_lit3,epa.k_lit4,epa.k_lit5,epa.k_lit6,
                epa.k_surface_som, epa.k_fast_som, epa.k_slow_som, epa.k_pass_som]
        K = np.zeros(289).reshape([17, 17])
        for i in range(0, 17):
            K[i][i] = temp[i]

        x_fin = np.zeros((cpa.number_of_months, 17))
        rh_fin = np.zeros((cpa.number_of_months, 1))
        ra_fin = np.zeros((cpa.number_of_months, 1))

        # leaf, root , wood, metabolic, structural, CWD, microbial, slow, passive
        x_init = np.array(
            [
                cpa.C_leaf_0,
                cpa.C_root_0,
                cpa.C_wood1_0,
                cpa.C_wood2_0,
                cpa.C_wood3_0,
                cpa.C_wood4_0,
                cpa.C_fruit_0,
                cpa.C_litter1_0,
                cpa.C_litter2_0,
                cpa.C_litter3_0,
                cpa.C_litter4_0,
                cpa.C_litter5_0,
                cpa.C_litter6_0,
                cpa.C_surface_som_0,
                cpa.C_fast_som_0,
                cpa.C_slow_som_0,
                cpa.C_pass_som_0,
            ]
        ).reshape([17, 1])  # Initial carbon pool size
        # initialize carbon pools
        X = x_init

        # initialize first respiration value
        co2_rh = cpa.rh_0
        co2_ra = cpa.ra_0


        # fixme:
        # slight change to the original
        # I would like to start the solution with the initial values
        # m=0 means after 0 months = in the initial step
        # B=A@K
        # pa=Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)
        # B=make_compartmental_matrix_func(pa)(0,X)
        def Rh_calc(TS, M, T0, E, KM):
            TS = TS - 273.15
            if TS > T0:
                rh_out = np.exp(E * (1 / (10 - T0) - 1 / (TS - T0))) * M / (KM + M)
            else:
                rh_out = 0
            return rh_out

        for m in range(0, cpa.number_of_months):
            x_fin[m, :] = X.reshape(1, 17)
            npp_in = cpa.npp[m]
            rh_fin[m, 0] = co2_rh
            ra_fin[m, 0] = co2_ra
            co2_rh = 0;
            co2_ra = 0;

            rh_modifier = Rh_calc(cpa.tsl[m], cpa.mrso[m], epa.T_0, epa.E, epa.KM)
            ksi = np.array(
                [
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier,
                    rh_modifier
                ]
            ).reshape([17, 1])  # environmental modifiers
            for d in range(0, days[m % 12]):
                K_new = K * ksi
                X = X + b * npp_in + A @ K_new @ X
                co2_rate = [0, 0, 0, 0, 0, 0, 0,
                            (1 - f158 - f168) * K[7, 7] * ksi[7],
                            (1 - f159 - f169) * K[8, 8] * ksi[8],
                            (1 - f1410 - f1610) * K[9, 9] * ksi[9],
                            (1 - f1511 - f1411) * K[10, 10] * ksi[10],
                            (1 - f1412) * K[11, 11] * ksi[11],
                            (1 - f1513) * K[12, 12] * ksi[12],
                            K[13, 13] * ksi[13],
                            K[14, 14] * ksi[14],
                            K[15, 15] * ksi[15],
                            K[16, 16] * ksi[16]
                            ]
                co2 = np.sum(co2_rate * X.reshape(1, 17))
                co2_rh = co2_rh + co2 / days[m % 12]  # monthly average rh
                #    litterfall_rate = [(f31 + f81) * K[0, 0] * ksi[0], (f42 + f92) * K[1, 1] * ksi[1], f83 * K[2, 2] * ksi[2],
                #                   f94 * K[3, 3] * ksi[3], (f105 + f125) * K[4, 4] * ksi[4], (f116 + f136) * K[5, 5] * ksi[5],
                 #                  (f107 + f127) * K[6, 6] * ksi[6],
                 #                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                #     litterfall = np.sum(litterfall_rate * X.reshape(1, 17))
           #     f_veg_lit = f_veg_lit + litterfall / days[m % 12]  # monthly average flux from veg to litter
             #   litter_to_soil_rate = [
             #       0, 0, 0, 0, 0, 0, 0,
              #      f158 * K[7, 7] * ksi[7] + f168 * K[7, 7] * ksi[7],
              #      f159 * K[8, 8] * ksi[8] + f169 * K[8, 8] * ksi[8],
               #     f1410 * K[9, 9] * ksi[9] + f1610 * K[9, 9] * ksi[9],
               #     f1511 * K[10, 10] * ksi[10] + f1411 * K[10, 10] * ksi[10],
              #      f1412 * K[11, 11] * ksi[11],
              #      f1513 * K[12, 12] * ksi[12],
              #      0, 0, 0, 0
              #  ]
             #   litter_to_soil = np.sum(litter_to_soil_rate * X.reshape(1, 17))
            #    f_lit_soil = f_lit_soil + litter_to_soil / days[m % 12]  # monthly average flux from litter to soil
        # We create an output that has the same shape
        # as the obvervations to make the costfunctions
        # easier.
        # To this end we project our 10 output variables of the matrix simulation
        # onto the 6 data streams by summing up the 3 litter pools into one
        # and also the 3 soil pools into one
        # C_litter_above = x_fin[:,3] + x_fin[:,4]
        C_wood = np.sum(x_fin[:, 0:3], axis=1).reshape(cpa.number_of_months, 1)
      #  print(np.shape(C_wood))
      #  C_leaf = np.sum(x_fin[:, 4], axis=1)(cpa.number_of_months, 1)
      #  C_root = np.sum(x_fin[:, 5], axis=1).reshape(cpa.number_of_months, 1)
      #  C_fruit = np.sum(x_fin[:, 6], axis=1).reshape(cpa.number_of_months, 1)
        C_litter = np.sum(x_fin[:, 7:12], axis=1).reshape(cpa.number_of_months, 1)
        # C_litter_below = x_fin[:,5]
        # c_litter = np.sum(x_fin[:,3:6],axis=1).reshape(cpa.number_of_months,1)
        C_soil = np.sum(x_fin[:,13:16],axis=1).reshape(cpa.number_of_months,1)
        # from IPython import embed; embed()
        out_simu = np.concatenate(
            [
                C_wood,
            #    C_leaf,
             #   C_root,
             #   C_fruit,
                x_fin[:, 4:6],  # the first 3 pools are used as they are
                C_litter,
                C_soil,
                rh_fin,
                ra_fin
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
        xs, co2s, acc_co2, acc_days = acc
        v = tsi.__next__()
        d_pools = v[0:-1, :]
        d_co2 = v[-1:, :]
        acc_co2 += d_co2
        acc_days += 1
        if i in day_indices:
            xs += [d_pools]
            co2s += [acc_co2 / acc_days]
            acc_co2 = np.array([0.0]).reshape(1, 1)
            acc_days = 0

        acc = (xs, co2s, acc_co2, acc_days)
        return acc

    xs, co2s, acc_days, _ = reduce(g, range(max(day_indices) + 1), ([], [], 0, 0))

    def h(tup):
        x, co2 = tup
        return np.transpose(np.concatenate([x, co2]))

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
    def npp_func(day, X):
        return mpa.npp[day_2_month_index(day)]

        # b (b vector for partial allocation)

    beta_root = 1 - mpa.beta_wood1 - mpa.beta_wood2 - mpa.beta_leaf - 0.1
    b = np.array(
        [
            mpa.beta_wood1,
            mpa.beta_wood2,
            0,
            0,
            mpa.beta_leaf,
            beta_root,
            0.1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0
        ],
    ).reshape(17, 1)
    # Now construct B matrix B=A*K
    B_func = make_compartmental_matrix_func(
        mpa=mpa,
    )

    # Build the iterator which is the representation of a dynamical system
    # for looping forward over the results of a difference equation
    # X_{i+1}=f(X_{i},i)
    # This is the discrete counterpart to the initial value problem
    # of the continuous ODE
    # dX/dt = f(X,t) and initial values X(t=0)=X_0

    def f(it, V):
        X = V[0:17]
        co2 = V[17]
        npp = npp_func(it, X)
        B = B_func(it, X)
        X_new = X + npp * b + B @ X

        # we also compute the respired co2 in every (daily) timestep
        # and use this part of the solution later to sum up the monthly amount
        co2_new = respiration_from_compartmental_matrix(B, X)

        V_new = np.concatenate((X_new, np.array([co2_new]).reshape(1, 1)), axis=0)
        return V_new

    # X_0 = np.array(x_init).reshape(9,1) #transform the tuple to a matrix
    # co2_0 = np.array([0]).reshape(1,1)
    # V_0 = np.concatenate((X_0, co2_0), axis=0)

    return TimeStepIterator2(
        initial_values=V_init,
        f=f,
        # max_it=max(day_indices)+1
    )


def make_compartmental_matrix_func(
        mpa
):
    # Now construct A matrix
    # diagonal
    # make sure to use 1.0 instead of 1 otherwise it will be an interger array
    # and round your values....
    A = np.diag([-1.0 for i in range(17)])
    # because of python indices starting at 0 we have A[i-1,j-1]=fij
    #
    # A = np.array([  -1,   0,   0,   0,   0,   0,   0,   0,   0,
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



    A[2, 0] = mpa.f_sapwood1_heartwood1;
    A[7, 0] = 1 - mpa.f_sapwood1_heartwood1;
    A[3, 1] = mpa.f_sapwood2_heartwood2;
    A[8, 1] = 1 - mpa.f_sapwood2_heartwood2;
    A[7, 2] = 1;
    A[8, 3] = 1;
    A[9, 4] = 0.712;
    A[11, 4] = 0.288;
    A[10, 5] = mpa.f_root_liiter4;
    A[12, 5] = 1 - mpa.f_root_liiter4;
    A[9, 6] = mpa.f_leaf_liiter3;
    A[11, 6] = 1 - mpa.f_leaf_liiter3;
    A[14, 7] = mpa.f_litter1_surface_som;
    A[15, 7] = mpa.f_litter1_fast_som;
    A[14, 8] = mpa.f_litter2_fast_som;
    A[15, 8] = mpa.f_litter2_slow_som;
    A[13, 9] = mpa.f_litter3_surface_som;
    A[15, 9] = mpa.f_litter3_slow_som;
    A[14, 10] = mpa.f_litter4_fast_som;
    A[13, 10] = mpa.f_litter4_surface_som;
    A[13, 11] = 0.4;
    A[14, 12] = 0.45;
    A[15, 13] = mpa.f_surface_som_slow_som;
    A[15, 14] = 0.296;
    A[16, 14] = 0.04;
    A[14, 15] = 0.42;
    A[16, 15] = 0.03;
    A[14, 16] = 0.45;


    # turnover rate per day of pools:
    K = np.diag([
        mpa.k_wood1,
        mpa.k_wood2,
        mpa.k_wood3,
        mpa.k_wood4,
        mpa.k_leaf,
        mpa.k_root,
        mpa.k_fruit,
        mpa.k_lit1,
        mpa.k_lit2,
        mpa.k_lit3,
        mpa.k_lit4,
        mpa.k_lit5,
        mpa.k_lit6,
        mpa.k_surfacesom,
        mpa.k_fastsom,
        mpa.k_slowsom,
        mpa.k_passsom
    ])

    # in the general nonautonomous nonlinear case
    # B will depend on an it,X (althouh in this case it does not depend on either
    def B_func(it, X):
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