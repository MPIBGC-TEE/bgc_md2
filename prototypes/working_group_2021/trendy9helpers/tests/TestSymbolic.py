# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
# import sys
import shutil
import re
from unittest.case import TestCase, skip
from pathlib import Path
from importlib import import_module
from importlib.resources import files as mod_files
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
from testinfrastructure.InDirTest import InDirTest

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


from trendy9helpers import general_helpers as gh

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


class TestSymbolic(InDirTest):

    @property
    def model_folders(self):
        return [
            # first tier (best shape)
            "kv_visit2",
            "jon_yib",
            "yz_jules",
            ###
            "Aneesh_SDGVM",  # second tier (not quite ready)
            ### "kv_ft_dlem",
            ####
            ####third tier
            ####"cj_isam", # has problems with ODE solution probably due to wrong parameters
            #### msh.numericX0 also yields a negative pool value for the last pool
            ### "bian_ibis2",#
            ####"cable-pop",
        ]


    def test_start_dt(self):
        # get everything from mvs
        # This is the blueprint for a notebook
        # that does not depend on any data exept those
        # provided in the model folder
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                print(gh.msh(mf).start_dt())

    def test_age_distributions_and_btt_start_in_ss_3(self):
        # get everything from mvs
        # This is the blueprint for a notebook
        # that does not depend on any data exept those
        # provided in the model folder
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                Path(mf).mkdir()
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs

                smr = mvs.get_SmoothModelRun()
                smr.initialize_state_transition_operator_cache(lru_maxsize=None)
                start_mean_age_vec = mvs.get_NumericStartMeanAgeTuple()
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
                stride = 5
                delta_t_val = (times[1]-times[0])/stride
                par_dict = mvs.get_NumericParameterization().par_dict
                func_dict = mvs.get_NumericParameterization().func_dict
                bit = ArrayDictResult(
                    gh.traceability_iterator_internal(
                        mvs,
                        X_0,
                        par_dict,
                        func_dict,
                        delta_t_val=delta_t_val,
                        t_0=t0
                    )
                )
                fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
                axs = fig1.subplots(n_pools, 2)
                times_max_index=min(len(times),100)# for testing
                iterator_max_step_index = times_max_index*stride  # for testing
                # produces vals.t = times[0:times_max_index]
                vals = bit[0:iterator_max_step_index:stride]
                for i in range(n_pools):
                    ax = axs[i, 0]

                    ax.plot(
                        times[: times_max_index],
                        solutions[: times_max_index, i],
                        label="sol"
                    )
                    ax.plot(vals.t[:], vals.X[:, i], label="bit")
                    ax.legend()
                    ax = axs[i, 1]
                    ax.plot(times, m_a_arr[:, i])

                fig1.savefig(Path(mf).joinpath("poolwise.pdf"))

                self.assertTrue(
                    np.allclose(
                        times[:times_max_index],
                        vals.t
                   )
                ) 
                #from IPython import embed; embed()
                self.assertTrue(
                    np.allclose(
                        solutions[:times_max_index,0:n_pools],
                        vals.X,
                        rtol=2e-2 #2%
                   )
                )   

                fig2 = plt.figure(figsize=(10, 10))
                axs2 = fig2.subplots(2, 2)
                mean_btts = smr.backward_transit_time_moment(
                    order=1, start_age_moments=start_mean_age_vec.reshape(1, n_pools)
                )
                ax = axs2[0, 0]
                ax.plot(times, mean_btts, label="mean backward transit time")
                ax.plot(
                    vals.t, vals.system_RT_sum, label="$\sum_i (RT)_i$"
                )  # steady state transit times
                ax.plot(
                    vals.t, vals.rt, label="rt of surrogate one pool system"
                )  # steady state transit times
                ax.legend()
                fig2.savefig(Path(mf).joinpath("system.pdf"))
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
                            "age_distribution_{0}.html".format(sv[n])
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
                    filename="btt_distribution.html",
                    auto_open=False,
                )

    def test_param2res(self):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                msh = import_module(f"trendy9helpers.{mf}.model_specific_helpers_2")
                mvs = import_module(f"{msh.model_mod}.source").mvs
                mdp = mod_files(f"trendy9helpers.{mf}")
                # where to look for global_means
                gm_cache_path = mdp.joinpath("global_means")
                svs, dvs = msh.get_global_mean_vars(
                    None, 
                    gm_cache_path, 
                    flash_cache=False
                )
                # we look for data assimilation submodules (directories) with
                # names da_0, da_1 ... ,da_9 and test all of them
                ex = re.compile("da_[0-9]")
                da_mod_names = [
                    f.stem 
                    for f in mod_files(f"trendy9helpers.{mf}").iterdir()
                    if bool(ex.match(f.stem)) & f.is_dir()
                ]
                #da_mod_names = ["da_0"]#,"da_1"] 
                for name in da_mod_names:
                    with self.subTest(name):
                        outp=Path(mf).joinpath(name)
                        outp.mkdir(parents=True)
                        da_mod = import_module(
                            f"trendy9helpers.{mf}.{name}.mod"
                        )
                        # read the constans and ingredients 
                        da_dir_p = mod_files(f"trendy9helpers.{mf}.{name}")
                        da_param_path = da_dir_p.joinpath( "par_1","in") 
                        fcpa = da_mod.FreeConstants(
                            **h.load_dict_from_json_path(da_param_path.joinpath("cpa.json"))
                        )
                        fepa_0 = da_mod.FreeEstimatedParameters(
                            **h.load_dict_from_json_path(da_param_path.joinpath("epa_0.json"))
                        )
                        mdp = mod_files(f"trendy9helpers.{mf}")
                        # where to look for global_means
                        gm_cache_path = mdp.joinpath("global_means")
                        svs, dvs = msh.get_global_mean_vars(
                            None, 
                            gm_cache_path, 
                            flash_cache=False
                        )
                        epa_0 = da_mod.epa_0(fepa_0, dvs, svs)
                        param2res = da_mod.make_param2res_sym(
                            mvs,
                            fcpa,
                            dvs,
                            svs 
                        )    
                        sim_0 = param2res(np.array(epa_0))
                        fig = plt.figure(figsize=(10,50))
                        axs=fig.subplots(len(svs._fields),1)
                        
                        
                        for ind,f in enumerate(svs._fields):
                            val_sim_0=sim_0.__getattribute__(f)
                            val_obs=svs.__getattribute__(f)
                            axs[ind].plot(range(len(val_sim_0)),val_sim_0,label=f+"_sim_0")
                            axs[ind].plot(range(len(val_obs)),val_obs,label=f+"_obs")
                            axs[ind].legend()
                            
                        fig.savefig(outp.joinpath('param2res.pdf'))
                        # check that the length of the returned arrays
                        # matches the length of the observations
                        for f in sim_0._fields:
                            self.assertEqual(
                                sim_0.__getattribute__(f).shape,
                                svs.__getattribute__(f).shape
                            )

    def test_da_res(self):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                Path(mf).mkdir()
                msh = import_module(f"trendy9helpers.{mf}.model_specific_helpers_2")
                mvs = import_module(f"{msh.model_mod}.source").mvs
                mdp = mod_files(f"trendy9helpers.{mf}")
                # where to look for global_means
                gm_cache_path = mdp.joinpath("global_means")
                svs, dvs = msh.get_global_mean_vars(
                    None, 
                    gm_cache_path, 
                    flash_cache=False
                )
                # we look for data assimilation submodules (directories) with
                # names da_0, da_1 ... ,da_9 and test all of them
                ex = re.compile("da_[0-9]")
                da_mod_names = [
                    f.stem 
                    for f in mod_files(f"trendy9helpers.{mf}").iterdir()
                    if bool(ex.match(f.stem)) & f.is_dir()
                ]
                #da_mod_names = ["da_1"]
                for name in da_mod_names:
                    with self.subTest(name):
                        da_mod = import_module(
                            f"trendy9helpers.{mf}.{name}.mod"
                        )
                        
                        func=gh.cached_da_res_1_maker(
                            da_mod.make_proposer,
                            da_mod.make_param_filter_func,
                            da_mod.make_param2res_sym,
                            msh.make_weighted_cost_func,
                            da_mod.numeric_X_0,
                            #msh.CachedParameterization,
                            da_mod.EstimatedParameters,
                        )    
                        # specify output dir (inside our the testdir)
                        output_cache_path=Path(mf).joinpath(name,"out")
                        # read the constans and ingredients 
                        da_dir_p = mod_files(f"trendy9helpers.{mf}.{name}")
                        da_param_path = da_dir_p.joinpath( "par_1","in") 
                        fcpa = da_mod.FreeConstants(
                            **h.load_dict_from_json_path(da_param_path.joinpath("cpa.json"))
                        )
                        def ep(p):
                            d = h.load_dict_from_json_path(p)
                            ep = da_mod.FreeEstimatedParameters(**d)
                            return ep

                        fepa_min, fepa_max, fepa_0 = tuple(
                            map(
                                ep,
                                [
                                    da_param_path.joinpath(f"{s}.json")
                                    for s in ["epa_min", "epa_max", "epa_0"]
                                ],
                            )
                        )
                        # check that we can read the hyperparameters
                        # change the name to hyper.json for the real
                        # hyperparameters her we use values for
                        # quick testing
                        hyper_dict = h.load_dict_from_json_path(
                            da_param_path.joinpath("hyper_test.json")
                        )
                        cpa = da_mod.cpa(fcpa, dvs, svs)
                        Cs, Js, epa_opt = func(
                            output_cache_path,
                            mvs,
                            svs,
                            dvs,
                            fcpa,
                            da_mod.epa_min(fepa_min,dvs,svs),
                            da_mod.epa_max(fepa_max,dvs,svs),
                            da_mod.epa_0(fepa_0,dvs,svs),
                            nsimu=hyper_dict['nsimu'],
                            acceptance_rate=hyper_dict['acceptance_rate'],
                            chunk_size=hyper_dict['chunk_size'],
                            D_init=hyper_dict['D_init'],
                            K=hyper_dict['K'],
                        )
                        
                        #param_dict=gh.make_param_dict(mvs,cpa,epa_opt) 
                        apa = {**cpa._asdict(), **epa_opt._asdict()}
                        param_dict=gh.make_param_dict2(mvs,apa) 
                        X_0=da_mod.numeric_X_0(mvs,dvs,apa)
                        X_0_dict={
                            sym: X_0[i,0] 
                            for i,sym in enumerate(
                                mvs.get_StateVariableTuple()
                            )
                        }
                        # some models (e.g. yz_jules) need extra parameters to build
                        # the func_dict for the parameterization
                        func_dict_param_dict = { 
                            str(k): v 
                            for k, v in apa.items() 
                            if str(k) in msh.CachedParameterization.func_dict_param_keys 
                        }
                        cp = msh.CachedParameterization(
                            param_dict,
                            dvs,
                            X_0_dict,
                            func_dict_param_dict
                        )
                        cp.write(output_cache_path)
                        # check if we can recreate cp from the cache directory
                        # check that we can run the model 
                        cp_fp = msh.CachedParameterization.from_path(output_cache_path)

                        X_0 = np.array(
                            [cp.X_0_dict[s] for s in mvs.get_StateVariableTuple()]
                        )
                        func_dict = cp.func_dict
                        par_dict = cp.parameter_dict
                        bit = gh.minimal_iterator_internal(
                            mvs,
                            X_0,
                            par_dict,
                            func_dict,
                            delta_t_val=2,
                        )
                        start = 0
                        stop = 10
                        step = 2
                        adr = ArrayDictResult(bit)
                        vals = adr[start:stop]

    def test_minimal_iterator_internal_args(self):
        # this is a smaller test and a stepping stone for da_iterator and
        # ultimately param2res
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                Path(mf).mkdir()
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs
                delta_t_val = 5

                CP = import_module(
                    f"{msh.model_mod}.CachedParameterization"
                ).CachedParameterization
                dir_path = mod_files(f"trendy9helpers.{mf}").joinpath("parameterization_from_test_args")
                cp = CP.from_path(dir_path)
                X_0 = np.array(
                    [cp.X_0_dict[s] for s in mvs.get_StateVariableTuple()]
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
                Path(mf).mkdir()
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs
                CP = import_module(
                    f"{msh.model_mod}.CachedParameterization"
                ).CachedParameterization
                dir_path = mod_files(f"trendy9helpers.{mf}").joinpath("parameterization_from_test_args")
                cp = CP.from_path(dir_path)
                X_0 = np.array(
                    [cp.X_0_dict[s] for s in mvs.get_StateVariableTuple()]
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
    
    def test_synthetic_observables(self):
        # this is a stepping stone for ultimately param2res
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                outp=Path(mf)
                outp.mkdir()
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs
                CP = import_module(
                    f"{msh.model_mod}.CachedParameterization"
                ).CachedParameterization
                dir_path = mod_files(f"trendy9helpers.{mf}").joinpath("parameterization_from_test_args")
                mdp = mod_files(f"trendy9helpers.{mf}")
                gm_cache_path = mdp.joinpath("global_means")
                svs, dvs = msh.get_global_mean_vars(
                    None, 
                    gm_cache_path, 
                    flash_cache=False
                )
                cp=CP(
                    parameter_dict=CP.parameter_dict_from_path(dir_path),
                    drivers=dvs,
                    X_0_dict=CP.X_0_dict_from_path(dir_path),
                    func_dict_param_dict=CP.func_dict_param_dict_from_path(dir_path)
                )
                #cp = CP.from_path(dir_path)
                X_0 = np.array(
                    [cp.X_0_dict[s] for s in mvs.get_StateVariableTuple()]
                )
                func_dict = cp.func_dict
                par_dict = cp.parameter_dict
                synth_obs=msh.synthetic_observables(
                    mvs,
                    X_0,
                    par_dict=par_dict,
                    func_dict=func_dict,
                    dvs=dvs
                )
                fig = plt.figure(figsize=(10,50))
                axs=fig.subplots(len(svs._fields),1)
                
                
                for ind,f in enumerate(svs._fields):
                    val_sim_0=synth_obs.__getattribute__(f)
                    val_obs=svs.__getattribute__(f)
                    axs[ind].plot(range(len(val_sim_0)),val_sim_0,label=f+"_synth")
                    axs[ind].plot(range(len(val_obs)),val_obs,label=f+"_obs")
                    axs[ind].legend()
                    
                fig.savefig(outp.joinpath('plot.pdf'))
                self.assertTrue(
                    all(
                        [
                            svs.__getattribute__(f).shape == synth_obs.__getattribute__(f).shape
                            for f in  svs._fields
                        ]
                    )
                )

    def test_aggregate_surrogate_systems(self):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                Path(mf).mkdir()
                msh = gh.msh(mf)
                mvs = import_module(f"{msh.model_mod}.source").mvs
                state_vector = mvs.get_StateVariableTuple()
                n_pools = len(state_vector)
                CP = import_module(
                    f"{msh.model_mod}.CachedParameterization"
                ).CachedParameterization
                dir_path = mod_files(f"trendy9helpers.{mf}").joinpath("parameterization_from_test_args")
                cp = CP.from_path(dir_path)
                X_0 = np.array(
                    [cp.X_0_dict[s] for s in mvs.get_StateVariableTuple()]
                )
                func_dict = cp.func_dict
                par_dict = cp.parameter_dict
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
                fig2.savefig(Path(mf).joinpath("system.pdf"))

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

                fig1.savefig(Path(mf).joinpath("poolwise.pdf"))
# If you want to run the script from this location (where the file lives)
