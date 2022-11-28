# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
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
import CompartmentalSystems.helpers_reservoir as hr
import unittest
import pathlib
import inspect
import shutil
import matplotlib.pyplot as plt
import numpy as np
from unittest.case import TestCase, skip
from pathlib import Path
from importlib import import_module
from collections import OrderedDict
import json
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
import matplotlib.pyplot as plt
from plotly.offline import plot

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
import bgc_md2.resolve.computers as bgc_c


import general_helpers as gh
import MIP_output_helpers as moh


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
            "kv_visit2",
            # "jon_yib",
            # "Aneesh_SDGVM",
            ##"cable-pop",
            ##"cj_isam", # has problems with ODE solution probably due to wrong parameters
            ## msh.numericX0 also yields a negative pool value for the last pool
            # "yz_jules",
            # "kv_ft_dlem",
            # "bian_ibis2",
        ]

    def test_aggregate_surrogate_systems(self):
        # at the moment only for visit
        # mf = "kv_visit2"
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                # subs_1 = import_module("{}.subs_1".format(mf))
                # s1 = import_module("{}.source_1".format(mf))
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                state_vector = mvs.get_StateVariableTuple()
                n_pools = len(state_vector)
                cpa = test_args.cpa
                epa = test_args.epa_opt
                dvs = test_args.dvs
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
                par_dict = gh.make_param_dict(mvs, test_args.cpa, epa)
                func_dict = msh.make_func_dict(mvs, test_args.dvs, test_args.cpa, epa)
                X0 = msh.numeric_X_0(mvs, test_args.dvs, test_args.cpa, epa).reshape(-1)
                # in our case we want to compute a steady state start value
                a_dens_function, X_fix = start_age_distributions_from_steady_state(
                    t0=0,
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
                # from IPython import embed; embed()
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
                        t0=0,
                        parameter_dict=par_dict,
                        func_set=sub_func_dict,
                        max_order=1,
                    )
                    from IPython import embed; embed()
                    return sub_mvs, sub_smr, sub_start_mean_age_vec
                        
                veg_sv = mvs.get_VegetationCarbonStateVariableTuple()
                veg_mvs, veg_smr, veg_smav = sub_mr_smav(mvs, veg_sv) 
                veg_btt = veg_smr.backward_transit_time_moment(
                    order=1,
                    start_age_moments=veg_smav
                )
                veg_arr, veg_func = veg_smr._solve_age_moment_system(
                    order,
                    start_age_moments=veg_smav
                )
                #from IPython import embed;embed()
                soil_sv = mvs.get_SoilCarbonStateVariableTuple()
                soil_mvs, soil_smr, soil_smav = sub_mr_smav(mvs, soil_sv) 
                #soil_smr.initialize_state_transition_operator_cache(lru_maxsize=None)
                soil_arr, soil_func = soil_smr._solve_age_moment_system(
                    order,
                    start_age_moments=soil_smav
                )
                soil_btt = soil_smr.backward_transit_time_moment(
                    order=1, start_age_moments=soil_smav
                )
                #def btt(mvs, svt):
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
                #        mode='lines',
                #        line=dict(
                #            color='#FF0000',
                #            width=15
                #            )
                #    )
                #    plot(fig, filename="test_{0}.html".format(n), auto_open=False)
                # LuoRT=s1.mvs.get_LuoRT()
                # Ot=s1.mvs.get_OutputTuple()

                # system_rt=LuoRT.dot(np.ones_like(X_fix))
                bit = gh.traceability_iterator(
                    X_fix,
                    func_dict,
                    mvs,
                    dvs,
                    test_args.cpa,
                    epa,
                    delta_t_val=delta_t_val,
                )
                vals = bit[0 : (n_steps + 1) : stride]
                disc_times = vals.t
                # from IPython import embed; embed()

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
                    vals.RT.sum(axis=1),
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
                    vals.tot_veg,
                    label="tot of veg pool",
                    color=color_dict["veg"],
                    marker=marker_dict["tot"],
                )
                key = "tot_soil"
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
                    vals.RT.sum(axis=1),
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
                ax.plot(times, vals.x_veg, color=color_dict["veg"], label="$x_veg$")
                ax.plot(times, vals.x_soil, color=color_dict["soil"], label="$x_soil$")
                ax.legend()
                fig2.savefig(testDir.joinpath("system.pdf"))

                fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
                axs = fig1.subplots(n_pools, 2)
                vnp = veg_smr.nr_pools
                snp = soil_smr.nr_pools
                for i,sym in enumerate(mvs.get_StateVariableTuple()):
                    print(sym)
                    ax = axs[i, 0]
                    ax.set_title(
                        'solutions {0}'.format(sym)
                    )
                    ax.plot(
                        vals.t,
                        vals.X[:,i],
                        color=color_dict["system"],
                        label="bit"
                    )
                    ax.plot(
                        times,
                        solutions[:, i],
                        color=color_dict["system"],
                        label="sol"
                    )
                    if sym in set(veg_sv):
                        i_veg=list(veg_sv).index(sym)
                        ax.plot(
                            times,
                            veg_arr[:, i_veg],
                            color=color_dict["veg"],
                            label="veg_sol"
                        )
                    if sym in set(soil_sv):
                        i_soil=list(soil_sv).index(sym)
                        ax.plot(
                            times,
                            soil_arr[:, i_soil],
                            color=color_dict["soil"],
                            label="soil_sol"
                        )
                    ax.legend()

                    ax = axs[i, 1]
                    ax.set_title("mean age {0}".format(sym))
                    ax.plot(
                        times,
                        m_a_arr[:, i],
                        color=color_dict["system"],
                        label="sol"
                    )
                    if sym in set(veg_sv):
                        i_veg=list(veg_sv).index(sym)
                        ax.plot(
                            times,
                            veg_arr[:, i_veg + vnp],
                            color=color_dict["veg"],
                            label="veg_ma"
                        )
                    if sym in set(soil_sv):
                        i_soil=list(soil_sv).index(sym)
                        ax.plot(
                            times,
                            soil_arr[:, i_soil + snp],
                            color=color_dict["soil"],
                            label="soil_ma"
                        )
                    ax.legend()

                fig1.savefig(testDir.joinpath("poolwise.pdf"))
                # svs, dvs = msh.get_example_site_vars(Path(conf_dict["dataPath"]))
                ## create the aggregated veg soil model
                # C__Veg,C__Soil=map(Symbol,["C__Veg","C__Soil"])
                # mvs_vs=CMTVS(
                #    {
                #        mvs.get_TimeSymbol(),
                #        StateVariableTuple([C__Veg,C__Soil]),
                #        InFluxesBySymbol({C__Veg: mvs.get_AggregatedVegetationCarbonInFlux()}),
                #        OutFluxesBySymbol(
                #            {
                #                C__Veg: mvs.get_AggregatedVegetationCarbonOutFlux(),
                #                C__Soil: mvs.get_AggregatedSoilCarbonOutFlux()
                #            }
                #        ),
                #        InternalFluxesBySymbol({
                #            (C__Veg,C__Soil): mvs.get_AggregatedVegetation2SoilCarbonFlux(),
                #            (C__Soil, C__Veg): mvs.get_AggregatedSoil2VegetationCarbonFlux()
                #        })
                #    },
                #    computers=h.module_computers(bgc_c)
                # )

                # h.compartmental_graph(mvs_vs)
                # func_dict=msh.make_func_dict(mvs,dvs, test_args.cpa, test_args.epa_0)
                # get the continuous system to check
