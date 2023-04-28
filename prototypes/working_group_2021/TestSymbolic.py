# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
# import sys
import shutil
from unittest.case import TestCase, skip
from pathlib import Path
from importlib import import_module
from sympy import (
    Symbol,
    Function,
    sympify,
    simplify,
    lambdify,
    Function,
    Symbol,
    diff,
    exp,
    diag,
)
from plotly.offline import plot
import matplotlib.pyplot as plt
import numpy as np

from CompartmentalSystems.smooth_model_run import SmoothModelRun
from CompartmentalSystems.start_distributions import (
    start_age_moments_from_empty_spinup,
    start_age_moments_from_steady_state,
    start_age_moments_from_zero_initial_content,
    compute_fixedpoint_numerically,
    start_age_distributions_from_steady_state,
    start_age_distributions_from_empty_spinup,
    start_age_distributions_from_zero_initial_content,
)
from CompartmentalSystems.ArrayDictResult import ArrayDictResult
import CompartmentalSystems.helpers_reservoir as hr

from ComputabilityGraphs.CMTVS import CMTVS

from bgc_md2.resolve.mvars import (
    NumericParameterization,
    NumericSimulationTimes,
    NumericStartValueArray,
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.helper as h


import general_helpers as gh
import MIP_output_helpers as moh

model_names = {
    "ab_classic": "CLASSIC",
    "clm5": "CLM5.0",
    "kv_ft_dlem": "DLEM",
    "bian_ibis2": "IBIS",
    "cj_isam": "ISAM",
    "isba-ctrip": "ISBA-CTRIP",
    "jsbach": "JSBACH",
    "yz_jules": "JULES-ES-1p0",
    "lpj-guess": "LPJ-GUESS",
    "lpjwsl": "LPJ",
    "lpx-bern": "LPX-Bern",
    "ORCHIDEE-V2": "OCN",
    "ORCHIDEE": "ORCHIDEE",
    "ORCHIDEE-CNP": "ORCHIDEE-CNP",
    "ORCHIDEEv3": "ORCHIDEEv3",
    "Aneesh_SDGVM": "SDGVM",
    "kv_visit2": "VISIT",
    "jon_yib": "YIBs",
}
experiment_names = {k: v + "_S2_" for k, v in model_names.items()}


class TestSymbolic(TestCase):
    @classmethod
    def clean_dir(cls, testDirPath):
        if testDirPath.exists():
            shutil.rmtree(testDirPath)
        testDirPath.mkdir(parents=True)

    @classmethod
    def data_dir_path(cls, mf):
        confDict = gh.confDict(mf)
        dataPath = Path(confDict["dataPath"])
        return dataPath

    def output_path(self, mf):
        testDirPath = __class__.data_dir_path(mf).joinpath(self.id())
        return testDirPath

    @property
    def model_folders(self):
        return [
            # first tier (best shape)
            "kv_visit2",
            "jon_yib",
            "yz_jules",
            ##
            "Aneesh_SDGVM",  # second tier (not quite ready)
            # "kv_ft_dlem",
            ##
            ##third tier
            ##"cj_isam", # has problems with ODE solution probably due to wrong parameters
            ## msh.numericX0 also yields a negative pool value for the last pool
            # "bian_ibis2",#
            ##"cable-pop",
        ]

    def test_make_model_index_transforms(self):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                msh = gh.msh(mf)
                lats = test_args.lats.data
                lons = test_args.lons.data
                # check that we predict the
                # right latitude values for
                # a given array index
                tr = msh.make_model_index_transforms()
                n_lats = len(lats)
                print("n_lats={}".format(n_lats))
                for i in range(n_lats):
                    self.assertEqual(lats[i], tr.i2lat(i))

                n_lons = len(lons)
                print("n_lons={}".format(n_lons))
                for i in range(n_lons):
                    self.assertEqual(lons[i], tr.i2lon(i))

                ## inverse operation
                ## check that we can predict the index from a given
                ## latitude
                # for i in range(n_lats):
                #    self.assertEqual(tr.lat2i(lats[i]), i)
                ## or longitude
                # for i in range(n_lons):
                #    self.assertEqual(tr.lon2i(lons[i]), i)

                ## check the interpretation of the pixel boundaries
                # for i in range(n_lats - 1):
                #    # minimum belongs to pixel
                #    self.assertEqual(tr.lat2i(tr.i2lat_min_max(i)[0]), i)
                #    # maximum belongs to next pixel (if there is one)
                #    self.assertEqual(tr.lat2i(tr.i2lat_min_max(i)[1]), i + 1)

                # # the last lat pixel contains also the pole
                # last = n_lats - 1
                # lat_min, lat_max = tr.i2lat_min_max(last)
                # print(lat_max)
                # self.assertEqual(tr.lat2i(lat_max), last)

    def test_make_model_coord_transforms(self):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                msh = gh.msh(mf)
                lats = test_args.lats.data
                lons = test_args.lons.data
                # check that we predict the
                # right latitude values for
                # a given array index
                ctr = msh.make_model_coord_transforms()
                # print(lats)
                print([ctr.lat2LAT(lat) for lat in lats])
                print([ctr.lon2LON(lon) for lon in lons])

    def test_age_distributions_and_btt_start_in_ss_3(self):
        # get everything from mvs
        # This is the blueprint for a notebook
        # that does not depend on any data exept those
        # provided in the model folder
        for mf in set(self.model_folders):
        #for mf in set(self.model_folders).intersection(
        #    # ["Aneesh_SDGVM"]
        #    ["kv_visit2"]
        #    # ["jon_yib"]
        #    # ["yz_jules"]
        #    ):
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs

                smr = mvs.get_SmoothModelRun()
                smr.initialize_state_transition_operator_cache(lru_maxsize=None)
                start_mean_age_vec = mvs.get_NumericStartMeanAgeVector()
                sv = mvs.get_StateVariableTuple()
                n_pools = len(sv)
                order = 1
                s_arr, s_func = smr._solve_age_moment_system(
                    order,
                    start_mean_age_vec.reshape(1, n_pools)
                    # method='rk45' does not support kwargs yet
                )

                # the first n colums are the solution
                # solutions = smr.solve()
                solutions = s_arr
                m_a_arr = s_arr[:, n_pools : 2 * n_pools]

                ## plot the continuous solution (by ODE) solver against the  iterator
                ## generated one.
                times = mvs.get_NumericSimulationTimes()
                X_0 = mvs.get_NumericStartValueArray()
                t0 = times[0]
                par_dict = mvs.get_NumericParameterization().par_dict
                func_dict = mvs.get_NumericParameterization().func_dict
                bit = ArrayDictResult(
                    gh.traceability_iterator_internal(
                        mvs, X_0, par_dict, func_dict, delta_t_val=15, t_0=t0
                    )
                )
                fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
                axs = fig1.subplots(n_pools, 2)
                n_steps = 12  # for testing
                vals = bit[0:n_steps]
                for i in range(n_pools):
                    ax = axs[i, 0]
                    ax.plot(times, vals.X[:, i], label="bit")

                    ax.plot(times, solutions[:, i], label="sol")
                    ax.legend()
                    ax = axs[i, 1]
                    ax.plot(times, m_a_arr[:, i])

                testDir = self.output_path(mf)
                self.__class__.clean_dir(testDir)
                fig1.savefig(testDir.joinpath("poolwise.pdf"))

                fig2 = plt.figure(figsize=(10, 10))
                axs2 = fig2.subplots(2, 2)
                mean_btts = smr.backward_transit_time_moment(
                    order=1, start_age_moments=start_mean_age_vec.reshape(1, n_pools)
                )
                ax = axs2[0, 0]
                ax.plot(times, mean_btts, label="mean backward transit time")
                ax.plot(
                    times, vals.system_RT_sum, label="$\sum_i (RT)_i$"
                )  # steady state transit times
                ax.plot(
                    times, vals.rt, label="rt of surrogate one pool system"
                )  # steady state transit times
                ax.legend()
                fig2.savefig(testDir.joinpath("system.pdf"))
                # construct a function p that takes an age array "ages" as argument
                # and gives back a three-dimensional ndarray (ages x times x pools)
                # from the a array-valued function of a single age a_dens_function
                srm = mvs.get_SmoothReservoirModel()
                a_dens_function, X_fix = start_age_distributions_from_steady_state(
                    srm, t0=t0, parameter_dict=par_dict, func_set=func_dict, x0=X_0
                )
                p = smr.pool_age_densities_func(a_dens_function)
                ages = np.linspace(
                    0,
                    (np.array(start_mean_age_vec, dtype=float).reshape(-1)).max() * 2,
                    2,
                )  # 21)
                age_densities = p(ages)

                for n in range(srm.nr_pools):
                    max_ind = np.argmin(ages < start_mean_age_vec[n] * 2)
                    fig = smr.plot_3d_density_plotly(
                        "age distribution pool {0}".format(sv[n]),
                        age_densities[0:max_ind, :, n],
                        ages[0:max_ind],
                    )
                    # plot the computed start age density for t0 on top
                    fig.add_scatter3d(
                        x=np.array([-t0 for a in ages]),
                        y=np.array([a for a in ages]),
                        z=np.array([a_dens_function(a)[n] for a in ages]),
                        mode="lines",
                        line=dict(color="#FF0000", width=15),
                    )
                    smr.add_line_to_density_plot_plotly(
                        fig,
                        data=m_a_arr[:, n],
                        color="#FF0000",
                        name="mean age",
                        time_stride=1,
                        on_surface=True,
                        bottom=True,
                        legend_on_surface=True,
                        legend_bottom=False,
                    )

                    plot(
                        fig,
                        filename=str(
                            testDir.joinpath("age_distribution_{0}.html".format(sv[n]))
                        ),
                        auto_open=False,
                    )

                btt_dens = smr.backward_transit_time_density(age_densities)
                fig_btt = smr.plot_3d_density_plotly(
                    "backward_transit_time_density_steady_state",
                    btt_dens,
                    ages,
                    y_label="transit time",
                )
                smr.add_line_to_density_plot_plotly(
                    fig_btt,
                    data=mean_btts,
                    color="#FF0000",
                    name="mean age",
                    time_stride=1,
                    on_surface=True,
                    bottom=True,
                    legend_on_surface=True,
                    legend_bottom=False,
                )
                plot(
                    fig_btt,
                    filename=str(testDir.joinpath("btt_distribution.html")),
                    auto_open=False,
                )

    def test_get_parameterization_from_data_assimilation(self):
        for mf in set(self.model_folders):
        #for mf in set(self.model_folders).intersection(
        #    # ["Aneesh_SDGVM"]
        #    # ["kv_visit2"]
        #    ["jon_yib"]
        #    # ["yz_jules"]
        #):
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs
                conf_dict = gh.confDict(mf)
                data_path = Path(conf_dict["dataPath"])
                tr_path = Path(mf).joinpath(
                    "data_assimilation_parameters_from_test_args"
                )
                cpa = msh.Constants(
                    **h.load_dict_from_json_path(tr_path.joinpath("cpa.json"))
                )

                epa_min, epa_max, epa_0 = tuple(
                    map(
                        lambda p: msh.EstimatedParameters(
                            **h.load_dict_from_json_path(p)
                        ),
                        [
                            tr_path.joinpath(f"{s}.json")
                            for s in ["epa_min", "epa_max", "epa_0"]
                        ],
                    )
                )
                target_path = Path(mf).joinpath("global_means")
                # when it runs the first time it will create the global mean
                # cache files there
                svs, dvs = msh.get_global_mean_vars_2(conf_dict, target_path)
                # run different data assimilation procedures
                for func in [
                    gh.get_parameterization_from_data_1,
                    # gh.get_parameterization_from_data_2,
                ]:
                    with self.subTest(func=func):
                        cpfd = gh.make_cached_data_assimilation_func_1(
                            func,
                            msh,
                        )
                        cp, _, _, _ = cpfd(
                            data_path,
                            msh,
                            mvs,
                            svs,
                            dvs,
                            cpa,
                            epa_min,
                            epa_max,
                            epa_0,
                            nsimu=15,
                        )

    def test_minimal_iterator_internal_args(self):
        # this is a smaller test and a stepping stone for da_iterator and
        # ultimately param2res
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs
                delta_t_val = 5

                CP = import_module(
                    f"{msh.model_mod}.CachedParameterization"
                ).CachedParameterization
                dir_path = Path(mf).joinpath("parameterization_from_test_args")
                cp = CP.from_path(dir_path)
                X_0 = np.array(
                    [cp.X_0_dict[str(s)] for s in mvs.get_StateVariableTuple()]
                )
                func_dict = cp.func_dict
                par_dict = cp.parameter_dict
                bit = gh.minimal_iterator_internal(
                    mvs,
                    X_0,
                    par_dict,
                    func_dict,
                    delta_t_val=delta_t_val,
                )
                start = 0
                stop = 10
                step = 2
                adr = ArrayDictResult(bit)
                vals = adr[start:stop]

                parts = gh.partitions(start, stop, step)
                res1 = adr.averaged_values(parts)  # using the result class
                res2 = vals.averaged_values(
                    parts
                )  # using the already computed valued d

                # if the number (stop-start) is a multiple of step then
                # we can create an averaging iterator and wrap it with the
                # ArrayDictResult class which provides the index notation

                res3 = ArrayDictResult(hr.average_iterator(bit, 2))[
                    start : int(stop / step)
                ]
                self.assertTrue(np.allclose(res1.X, res2.X))

    def test_make_da_iterator(self):
        # this is a stepping stone for ultimately param2res
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs
                CP = import_module(
                    f"{msh.model_mod}.CachedParameterization"
                ).CachedParameterization
                dir_path = Path(mf).joinpath("parameterization_from_test_args")
                cp = CP.from_path(dir_path)
                X_0 = np.array(
                    [cp.X_0_dict[str(s)] for s in mvs.get_StateVariableTuple()]
                )
                func_dict = cp.func_dict
                par_dict = cp.parameter_dict
                delta_t_val = 5
                #par_dict = gh.make_param_dict(mvs, cpa, epa)
                bit = msh.make_da_iterator(
                    mvs,
                    X_0,
                    par_dict,
                    func_dict,
                    delta_t_val=delta_t_val,
                )
                start = 0
                stop = 10
                step = 2
                adr = ArrayDictResult(bit)
                vals = adr[start:stop]

    def test_aggregate_surrogate_systems(self):
        for mf in set(self.model_folders):
            # for mf in set(self.model_folders).intersection(
            #    #["jon_yib"]
            #    ["kv_visit2"]
            # ):
            with self.subTest(mf=mf):
                # test_args = gh.test_args_2(mf)
                # cpa = test_args.cpa
                # epa = test_args.epa_opt
                # dvs = test_args.dvs
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs
                state_vector = mvs.get_StateVariableTuple()
                n_pools = len(state_vector)
                CP = import_module(
                    f"{msh.model_mod}.CachedParameterization"
                ).CachedParameterization
                dir_path = Path(mf).joinpath("parameterization_from_test_args")
                cp = CP.from_path(dir_path)
                X_0 = np.array(
                    [cp.X_0_dict[str(s)] for s in mvs.get_StateVariableTuple()]
                )
                func_dict = cp.func_dict
                par_dict = cp.parameter_dict
                testDir = self.output_path(mf)
                self.__class__.clean_dir(testDir)
                conf_dict = gh.confDict(mf)
                t_min = 0
                t_max = 1400
                delta_t_val = 5  # this affects the precision of the iterator
                stride = 6  # this does not affects the precision of the iterator
                # but makes it more effiecient (the values in between the strides
                # are computed but not stored)
                # the presision of the continuous solution is not affected by either
                # number since the solver decides where to compute the next value
                # the times argument just tells it where WE want to know the values...

                n_steps = int((t_max - t_min) / delta_t_val)
                # times = np.linspace(t_min, t_max, n_steps)
                times = np.arange(t_min, t_max, delta_t_val * stride)
                # in our case we want to compute a steady state start value
                a_dens_function, X_fix = start_age_distributions_from_steady_state(
                    t0=t_min,  # time of freezing (determines B(t_0) and u(t_0) from which to compute the equilibrium)
                    srm=mvs.get_SmoothReservoirModel(),
                    parameter_dict=par_dict,
                    func_set=func_dict,
                )
                mvs = mvs.update(
                    {
                        NumericParameterization(par_dict=par_dict, func_dict=func_dict),
                        NumericStartValueArray(X_fix),
                        NumericSimulationTimes(times),
                    }
                )
                t0 = times[0]
                sv = mvs.get_StateVariableTuple()
                n_pools = len(sv)
                smr = mvs.get_SmoothModelRun()
                # now create two submodels
                start_mean_age_vec = start_age_moments_from_steady_state(
                    smr.model,
                    t0=0,
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
                m_a_arr = s_arr[:, n_pools : 2 * n_pools]
                t = mvs.get_TimeSymbol()

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
                        lambda d: {
                            k: sympify(flux_ex).subs(subs_dict)
                            for k, flux_ex in d.items()
                        },
                        (sub_in_fluxes, sub_internal_fluxes, sub_out_fluxes),
                    )

                    # now prepare the extension of the func_dict by the solutions
                    # for the variables
                    def s_func_maker(sym):
                        return lambda t: s_func(t)[svl.index(sym)]

                    outer_pool_funcs = {
                        # note that we create the function in another
                        # function (s_func_maker) since the normal dictionary
                        # comprehension
                        # fails due to pythons lazy evaluation not
                        # protecting the function specific variable
                        # (in this case sym) so that all flux functions
                        # would be the same
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
                    # for the start mean ages we can not just take the
                    # respective parts of the system mean ages since they refer
                    # to the times stuff has spent in the SYSTEM and not in the
                    # SUB system.  Although the pools in the subsystem are
                    # pools in the
                    sub_start_mean_age_vec = start_age_moments_from_steady_state(
                        sub_smr.model,
                        t0=0,
                        parameter_dict=par_dict,
                        func_set=sub_func_dict,
                        max_order=1,
                    )
                    return sub_mvs, sub_smr, sub_start_mean_age_vec

                veg_sv = mvs.get_VegetationCarbonStateVariableTuple()
                veg_mvs, veg_smr, veg_smav = sub_mr_smav(mvs, veg_sv)
                veg_btt = veg_smr.backward_transit_time_moment(
                    order=1, start_age_moments=veg_smav
                )
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
                # def btt(mvs, svt):
                #    sub_smr, sub_start_mean_age_vec = sub_mr_smav(mvs, svt)
                #    return sub_smr.backward_transit_time_moment(
                #        order=1, start_age_moments=sub_start_mean_age_vec
                #    )

                # compute the start age distribution

                # smr.initialize_state_transition_operator_cache(lru_maxsize=None)
                order = 1
                # From a_dens_function we should be able to extract the mean age
                # but the numeric integration does not converge for this
                # numerically computed example.  It usually has no problem
                # with an age distribution given as an explicit function but we can
                # compute the mean age directly as (1,..,1)B^{-1}
                print(start_mean_age_vec)
                ##ages = np.linspace(-1,10,12)
                # ages = np.linspace(0,500,3)
                # construct a function p that takes an age array "ages" as argument
                # and gives back a three-dimensional ndarray (ages x times x pools)
                # from the a array-valued function representing the start age density
                # p = smr.pool_age_densities_func(
                #    a_dens_function
                #    #start_age_distributions_from_zero_initial_content(srm)
                #    #lambda a: np.exp(-a)*X0
                # )
                # age_densities = p(ages)
                # btt_dens = smr.backward_transit_time_density(age_densities)
                # for n in range(srm.nr_pools):
                #    fig = smr.plot_3d_density_plotly(
                #        "pool {0}".format(n), age_densities[:, :, n], ages
                #    )
                #    # plot the computed start age density for t0 on top
                #    fig.add_scatter3d(
                #        x=np.array([-t0 for a in ages]),
                #        y=np.array([a for a in ages]),
                #  14      z=np.array([a_dens_function(a)[n] for a in ages]),
                nn#        mode='lines',
                #        line=dict(
                #            color='#FF0000',
                #            width=15
                #            )
                #    )
                #    plot(fig, filename="test_{0}.html".format(n), auto_open=False)
                # LuoRT=s1.mvs.get_LuoRT()
                # Ot=s1.mvs.get_OutputTuple()

                # system_rt=LuoRT.dot(np.ones_like(X_fix))
                it = gh.traceability_iterator_internal(
                    mvs, X_fix, par_dict, func_dict, t_min, delta_t_val
                )
                bit = ArrayDictResult(it)
                vals = bit[0 : (n_steps + 1) : stride]
                disc_times = vals.t

                fig2 = plt.figure(figsize=(20, 20))
                axs2 = fig2.subplots(4, 1)
                mean_btts = smr.backward_transit_time_moment(
                    order=1, start_age_moments=start_mean_age_vec.reshape(1, n_pools)
                )
                color_dict = {"veg": "green", "soil": "brown", "system": "black"}
                marker_dict = {"RT": "*", "btt": "+", "tot": "o"}
                #############################################################
                ax = axs2[0]
                ax.set_title("turnover times vs $\sum_i (RT)_i$")
                ax.plot(
                    disc_times,
                    # vals.RT.sum(axis=1),
                    vals.system_RT_sum,
                    label="$\sum_i (RT)_i$",
                    color=color_dict["system"],
                    marker=marker_dict["RT"],
                )
                ax.plot(
                    disc_times,
                    vals.rt,
                    color=color_dict["system"],
                    marker=marker_dict["tot"],
                    label="tot = rt of surrogate one pool system",
                )
                ax.plot(
                    disc_times,
                    vals.veg_tot,
                    label="tot of veg pool",
                    color=color_dict["veg"],
                    marker=marker_dict["tot"],
                )
                key = "soil_tot"
                ax.plot(
                    disc_times,
                    vals[key],
                    label=key,
                    color=color_dict["soil"],
                    marker=marker_dict["tot"],
                )
                ax.legend()

                #############################################################
                ax = axs2[1]
                ax.set_title("mean transit times vs $\sum_i (RT)_i$")
                ax.plot(
                    times,
                    mean_btts,
                    color=color_dict["system"],
                    marker=marker_dict["btt"],
                    label="mean system backward transit time",
                )
                for i, r in enumerate((veg_btt,)):
                    ax.plot(
                        times,
                        r,
                        color=color_dict["veg"],
                        marker=marker_dict["btt"],
                        label=f"mean veg {0} subsystem backward transit time".format(i),
                    )
                for i, r in enumerate((soil_btt,)):
                    ax.plot(
                        times,
                        r,
                        color=color_dict["soil"],
                        marker=marker_dict["btt"],
                        label=f"mean soil {0} sub system backward transit time".format(
                            i
                        ),
                    )
                ax.plot(
                    times,
                    vals.system_RT_sum,
                    color=color_dict["system"],
                    marker=marker_dict["RT"],
                    label="$\sum_i (RT)_i$",
                )  # steady state transit times
                ax.legend()

                #############################################################
                ax = axs2[2]
                ax.set_title("Fluxes")
                for l in ["out_2_veg", "veg_2_out", "veg_2_soil", "soil_2_out"]:
                    ax.plot(times, vals[l], label="{0}".format(l))
                ax.legend()

                #############################################################
                ax = axs2[3]
                ax.set_title("Stocks")
                ax.plot(times, vals.x, color=color_dict["system"], label="$x$")
                ax.plot(times, vals.veg_x, color=color_dict["veg"], label="$x_veg$")
                ax.plot(times, vals.soil_x, color=color_dict["soil"], label="$x_soil$")
                ax.legend()
                fig2.savefig(testDir.joinpath("system.pdf"))

                fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
                axs = fig1.subplots(n_pools, 2)
                vnp = veg_smr.nr_pools
                snp = soil_smr.nr_pools
                for i, sym in enumerate(mvs.get_StateVariableTuple()):
                    print(sym)
                    ax = axs[i, 0]
                    ax.set_title("solutions {0}".format(sym))
                    ax.plot(
                        vals.t, vals.X[:, i], color=color_dict["system"], label="bit"
                    )
                    ax.plot(
                        times, solutions[:, i], color=color_dict["system"], label="sol"
                    )
                    if sym in set(veg_sv):
                        i_veg = list(veg_sv).index(sym)
                        ax.plot(
                            times,
                            veg_arr[:, i_veg],
                            color=color_dict["veg"],
                            label="veg_sol",
                        )
                    if sym in set(soil_sv):
                        i_soil = list(soil_sv).index(sym)
                        ax.plot(
                            times,
                            soil_arr[:, i_soil],
                            color=color_dict["soil"],
                            label="soil_sol",
                        )
                    ax.legend()

                    ax = axs[i, 1]
                    ax.set_title("mean age {0}".format(sym))
                    ax.plot(
                        times, m_a_arr[:, i], color=color_dict["system"], label="sol"
                    )
                    if sym in set(veg_sv):
                        i_veg = list(veg_sv).index(sym)
                        ax.plot(
                            times,
                            veg_arr[:, i_veg + vnp],
                            color=color_dict["veg"],
                            label="veg_ma",
                        )
                    if sym in set(soil_sv):
                        i_soil = list(soil_sv).index(sym)
                        ax.plot(
                            times,
                            soil_arr[:, i_soil + snp],
                            color=color_dict["soil"],
                            label="soil_ma",
                        )
                    ax.legend()

                fig1.savefig(testDir.joinpath("poolwise.pdf"))
