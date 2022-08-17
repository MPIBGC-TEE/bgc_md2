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
#cpa should be here
UnEstimatedParameters = namedtuple(
    "UnEstimatedParameters",
    [
        'C_Veg_0',
        'C_root_0',
        'C_litter_0',
        'npp',
        'number_of_months',
        'C_soil_0',
        'rh_o'
    ]
)

# This set is used (as argument) by the functions that are called
# inside the mcmc
#Epa
EstimatedParameters = namedtuple(
    "EstiamatedParameters",
    [
        "beta_leaf",    #  0 (indices uses in original code)
        "beta_wood",    #  1
        "Ca",#  2
        "clay",#  3
        "leachedwtr30",#  4
        "sand",  # 5
        "silt_clay",  # 6
        "k_leaf",  # 7
        "k_wood",  # 8
        "k_root",  # 9
        "C_leaf_0",  # 10
        "C_abvstrlit_0",       #  11
        "C_abvmetlit_0",       #  12
        "C_blwstrlit_0",       #  13
        "C_surfacemic_0",	#  14
        "C_soilmic_0",	#  15
        "C_slow_0"   #  16
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
        'C_abvstructural_lit',
        'C_abvmetabolic_lit',
        'C_belowstructual__lit',
        'C_belowmetabolic_lit',
        'C_surface_microbe',
        'C_soil_microbe',
        'C_slow_soil'
        'C_passive_soil'
    ]
)
Observables = namedtuple(
    'Observables',
    [
        'C_Veg',
        'C_root',
        'C_litter',
        'C_soil',
        'rh'
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
        'C_Veg_0',
        'C_root_0',
        'C_litter_0',
        'number_of_months',
        'C_soil_0'
        'rh_0',
        'C_leaf_0',
        'C_abvstrlit_0',
        'C_abvmetlit_0',
        'C_blwstrlit_0',
        'C_surfacemic_0',
        'C_soilmic_0',
        'C_slow_0'
        ]
    ]
)









# to do
# 1.) make a namedtuple for the SDGVM data or use xarray to create a multifile dataset
def get_variables_from_files(dataPath):
    # Read NetCDF data  ******************************************************************************************************************************
    path = dataPath.joinpath("CABLE-POP_S2_npp.nc")
    ds = nc.Dataset(str(path))
    var_npp = ds.variables['npp'][:, :, :]
    ds.close()

    path = dataPath.joinpath("CABLE-POP_S2_rh.nc")
    ds = nc.Dataset(str(path))
    var_rh = ds.variables['rh'][:, :, :]
    ds.close()

    path = dataPath.joinpath("CABLE-POP_S2_cLeaf.nc")
    ds = nc.Dataset(str(path))
    var_cleaf = ds.variables['cLeaf'][:, :, :]
    ds.close()

    path = dataPath.joinpath("CABLE-POP_S2_cRoot.nc")
    ds = nc.Dataset(str(path))
    var_croot = ds.variables['cRoot'][:, :, :]
    ds.close()

    path = dataPath.joinpath("CABLE-POP_S2_cVeg.nc")
    ds = nc.Dataset(str(path))
    var_cveg = ds.variables['cVeg'][:, :, :]
    ds.close()

    path = dataPath.joinpath("CABLE-POP_S2_cSoil.nc")
    ds = nc.Dataset(str(path))
    var_csoil = ds.variables['cSoil'][:, :, :]
    ds.close()

    path = dataPath.joinpath("CABLE-POP_S2_cLitter.nc")
    ds = nc.Dataset(str(path))
    var_clitter = ds.variables['cLitter'][:, :, :]
    ds.close()

    path = dataPath.joinpath("CABLE-POP_S2_cCwd.nc")
    ds = nc.Dataset(str(path))
    var_ccwd = ds.variables['cCwd'][:, :, :]
    ds.close()

    return (var_npp, var_rh, var_cleaf, var_croot, var_cveg, var_csoil, var_clitter, var_ccwd)


def get_example_site_vars(dataPath):
    var_npp, var_rh, var_cleaf, var_croot, var_cveg, var_csoil, var_clitter, var_ccwd = get_variables_from_files(
        dataPath)
    # pick up 1 site   wombat state forest
    s = slice(None, None, None)  # this is the same as :
    t = s, 49, 325  # [t] = [:,49,325]
    npp = var_npp[t] * 86400  # kg/m2/s kg/m2/day;
    rh = var_rh[t] * 86400;  # per s to per day
    (
        clitter,
        csoil,
        cveg,
        cleaf,
        croot,
        ccwd
    ) = map(
        lambda var: var[t],
        (
            var_clitter,
            var_csoil,
            var_cveg,
            var_cleaf,
            var_croot,
            var_ccwd
        )
    )
    cwood = cveg - cleaf - croot;
    return (npp, rh, clitter, csoil, cveg, cleaf, croot, ccwd, cwood)


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
        cond3 = (c[0] + c[1]) < 1
        return (cond1 and cond2 and cond3)

    return isQualified


def make_weighted_cost_func(
        obs: np.ndarray
) -> Callable[[np.ndarray], np.float64]:
    # first unpack the observation array into its parts
    cleaf, croot, cwood, clitter, csoil, rh = np.split(obs, indices_or_sections=6, axis=1)

    def costfunction(out_simu: np.ndarray) -> np.float64:
        # fixme
        #   as indicated by the fact that the function lives in this
        #   model-specific module it is not apropriate for (all) other models.
        #   There are model specific properties:
        #   1.) The weight for the different observation streams
        #
        tot_len = out_simu.shape[0]
        # we assume the model output to be in the same shape and order
        # as the obeservation
        # this convention has to be honored by the forwar_simulation as well
        # which in this instance already compresses the 3 different litter pools
        # to c_litter and the 3 different soil pools to one
        c_simu = out_simu[:, 0:5]

        # we assume the rh  part to be in the remaining columns again
        # this convention has to be honored by the forwar_simulation as well
        rh_simu = out_simu[:, 5:]
        # from IPython import embed; embed()

        J_obj1 = np.mean((c_simu[:, 0] - cleaf[0:tot_len]) ** 2) / (2 * np.var(cleaf[0:tot_len]))
        J_obj2 = np.mean((c_simu[:, 1] - croot[0:tot_len]) ** 2) / (2 * np.var(croot[0:tot_len]))
        J_obj3 = np.mean((c_simu[:, 2] - cwood[0:tot_len]) ** 2) / (2 * np.var(cwood[0:tot_len]))
        J_obj4 = np.mean((c_simu[:, 3] - clitter[0:tot_len]) ** 2) / (2 * np.var(clitter[0:tot_len]))
        J_obj5 = np.mean((c_simu[:, 4] - csoil[0:tot_len]) ** 2) / (2 * np.var(csoil[0:tot_len]))

        J_obj6 = np.mean((rh_simu[:, 0] - rh[0:tot_len]) ** 2) / (2 * np.var(rh[0:tot_len]))

        J_new = (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5) / 200 + J_obj6 / 4
        return J_new

    return costfunction













def make_param2res(
        cpa
) -> Callable[[np.ndarray], np.ndarray]:
    """The returned function 'param2res' is used by the metropolis hastings algorithm
    to produce a result from a proposed parameter(tuple).
    'param2res' will be called only with the parameters to be estimated 'pe'.
    The constant parameters are the parameters of the present function 'pc' and are automatically available in the returned 'param2res' (closure).

    In the case of this model 'param2res' has to perform the following
    tasks.
    -   use both sets of parameters 'pc' and 'pe' to build a forward model
    -   run the forward model to produce output for the times
        when obeservations are available.(Here monthly)
    -   project the output of the forward model to the obeserved variables.
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
        days = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
        # Construct b vector
        beta1 = epa.beta_leaf;
        beta2 = epa.beta_wood;
        beta3 = 1 - (beta1 + beta2);
        b = np.array([beta1, beta2, beta3, 0, 0, 0, 0, 0, 0, 0, 0]).reshape([11, 11])  # allocation
        # Now construct A matrix
       # lig_leaf = epa.lig_leaf
        ls_aboveground = 0.12
        ls_belowground = 0.35
        f41 = 0.91986 + 0.00324* epa.Ca;
        f51 = 1 - f41;
        f42 = f41;
        f52 = 1 - f42;
        f63 = f41;
        f73 = 1 - f63;
        f84 = (1 - np.exp(-3*ls_aboveground))*0.4;
        f104 = np.exp(-3*ls_aboveground) * 0.7;
        f85 = 0.4;
        f96 = (1 - np.exp(-3*ls_belowground))*0.45;
        f106 = np.exp(-3*ls_aboveground) * 0.7;
        f97 = 0.45;
        f910 = 0.447 + 0.009*epa.clay;
        f911 = 0.45;
        f108 = 0.4;
        f119 = 0.003 + 0.032 * epa.clay;
        f109 = 1 - f119 - (epa.leachedwtr30)/18 * (0.01 + 0.04*epa.sand) - 0.85 - o.68*epa.silt_clay;
        f1110 = 0.003 - 0.009*epa.clay;

        A = np.array([-1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                      0,  -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                      0,   0, -1,  0,  0,  0,  0,  0,  0,  0,  0,
                      f41, f42, 0,-1,  0,  0,  0,  0,  0,  0,  0,
                      f51, f52, 0, 0,  -1, 0,  0,  0,  0,  0,  0,
                      0,  0,  f63, 0,  0,  -1, 0,  0,  0,  0,  0,
                      0,  0,  f73, 0,  0,  0,  -1, 0,  0,  0,  0,
                      0,  0,  0,  f84, f85, 0, 0, -1,  0,  0,  0,
                      0,  0,  0,  0,  0, f96, f97, 0,  -1, f910, f911,
                      0,  0,  0,  f104,  0, f106, 0, f108, f109, -1, 0,
                      0,  0,  0,  0,  0,  0,  0,  0, f119, f1110, -1]).reshape([11, 11])  # tranfer

        # turnover rate per day of pools:
        temp = [epa.k_leaf, epa.k_wood, epa.k_root, 3.9/365, 4.8/365, 7.3/365, 6.0/365, 14.8/365, 18.5/365, 0.2/365, 0.0045/365]
        K = np.zeros(121).reshape([11, 11])
        for i in range(0, 11):
            K[i][i] = temp[i]

        x_fin = np.zeros((cpa.number_of_months, 11))
        rh_fin = np.zeros((cpa.number_of_months, 1))
        # leaf, root , wood, metabolic, structural, CWD, microbial, slow, passive
        x_init = np.array(
            [
                epa.C_leaf_0,
                cpa.C_veg_0 - epa.C_leaf_0 - cpa.C_root_0,
                cpa.C_root_0,
                epa.C_abvstrlit_0,
                epa.C_abvmetlit_0,
                epa.C_blwstrlit_0,
                cpa.C_litter_0 - epa.C_abvstrlit_0 - epa.C_abvmetlit_0 - epa.C_blwstrlit_0,
                epa.C_surfacemic_0,
                epa.C_soilmic_0,
                epa.C_slow_0,
                cpa.C_soil_0 - epa.C_soilmic_0 - epa.C_slow_0
            ]
        ).reshape([11, 1])  # Initial carbon pool size
        # initialize carbon pools
        X = x_init

        # initialize first respiration value
        co2_rh = cpa.rh_0
        # fixme:
        # slight change to the original
        # I would like to start the solution with the initial values
        # m=0 means after 0 moths = in the initial step
        B = A @ K
        # pa=Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)
        # B=make_compartmental_matrix_func(pa)(0,X)

        for m in range(0, cpa.number_of_months):
            x_fin[m, :] = X.reshape(1, 11)
            npp_in = cpa.npp[m]
            rh_fin[m, 0] = co2_rh
            co2_rh = 0
            for d in range(0, days[m % 12]):
                X = X + b * npp_in + B @ X
                co2_rate = [0, 0, 0, (1 - f84) * K[3, 3], (1 - f85) * K[4, 4], (1 - f96 - f106) * K[5, 5],
                            (1 - f97) * K[6, 6], (1 - f108) * K[7, 7], (1 - f109 - f119) * K[8, 8], (1 - f910 - f1110) * K[9, 9], f911* K[10, 10]]
                co2 = np.sum(co2_rate * X.reshape(1, 11))
                co2_rh = co2_rh + co2 / days[m % 12]  # monthly average rh

        # We create an output that has the same shape
        # as the obvervations to make the costfunctions
        # easier.
        # To this end we project our 10 output variables of the matrix simulation
        # onto the 6 data streams by summing up the 3 litter pools into one
        # and also the 3 soil pools into one
        cVeg = np.sum(x_fin[:, 0:3], axis=1).reshape(cpa.number_of_months, 1)
        c_litter = np.sum(x_fin[:, 3:7], axis=1).reshape(cpa.number_of_months, 1)
        c_soil = np.sum(x_fin[:, 8:11], axis=1).reshape(cpa.number_of_months, 1)
        # from IPython import embed; embed()
        out_simu = np.concatenate(
            [
                cVeg,
                x_fin[:, 2],  # root
                c_litter,
                c_soil,
                rh_fin
            ]
            , axis=1
        )
        return out_simu

    return param2res