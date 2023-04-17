# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
import unittest
import pathlib
import inspect
import shutil
import matplotlib.pyplot as plt
import numpy as np
import json
import matplotlib.pyplot as plt
from unittest.case import TestCase, skip
from pathlib import Path
from importlib import import_module
from collections import OrderedDict
from sympy import  (
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
from functools import partial

from CompartmentalSystems.smooth_model_run import SmoothModelRun
from CompartmentalSystems.start_distributions import (
    start_age_moments_from_empty_spinup,
    start_age_moments_from_steady_state,
    start_age_moments_from_zero_initial_content,
    compute_fixedpoint_numerically,
    start_age_distributions_from_steady_state,
    start_age_distributions_from_empty_spinup,
    start_age_distributions_from_zero_initial_content
)
from CompartmentalSystems.ArrayDictResult import  ArrayDictResult
import CompartmentalSystems.helpers_reservoir as  hr

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
    "jon_yib": "YIBs"    
}
experiment_names = {k:v+"_S2_" for k,v in model_names.items() }

class TestSymbolic(TestCase):

    @classmethod
    def clean_dir(cls, testDirPath):
        if testDirPath.exists():
            shutil.rmtree(testDirPath)
        testDirPath.mkdir(parents=True)

    @classmethod
    def data_dir_path(cls, mf):
        confDict = gh.confDict(mf)
        dataPath = Path(confDict['dataPath'])
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
            "Aneesh_SDGVM", #second tier (not quite ready)
            #"kv_ft_dlem",
            ##
            ##third tier
            ##"cj_isam", # has problems with ODE solution probably due to wrong parameters
            ## msh.numericX0 also yields a negative pool value for the last pool
            #"bian_ibis2",# 
            ##"cable-pop", 
        ]

    def test_age_distributions_and_btt_start_in_ss_4(self):
        # this test checks that we can reproduce our 
        # mean age and age distribution from cached data 
        # produced by one of the dataassimilation procedures
        # expressed as a function # (here get_parameterization_from_data_1 
        # but it could be another one of the same kind 
        # the point is that we reconstruct a model run 
        # from cached data
        for mf in set(self.model_folders).intersection(["kv_visit2"]):
            with self.subTest(mf=mf):
                th = gh.th(mf)
                msh = gh.msh(mf)
                mvs=import_module(f"{msh.model_mod}.source").mvs
                conf_dict=gh.confDict(mf)
                data_path=Path(conf_dict['dataPath'])
                # obviously has to be here
                tr_path=Path(mf).joinpath("data_assimilation_parameters_from_test_args")
                cpa= msh.Constants(
                    **h.load_dict_from_json_path(
                        tr_path.joinpath("cpa.json") 
                    )
                )

                epa_min,epa_max,epa_0=tuple(
                    map(
                        lambda p:msh.EstimatedParameters(**h.load_dict_from_json_path(p)),
                        [
                            tr_path.joinpath(f"{s}.json") 
                            for s in ['epa_min','epa_max','epa_0']
                        ]
                    )
                )
                svs, dvs = msh.get_global_mean_vars_2(conf_dict)
                cpfd_1 = gh.make_cached_data_assimilation_func(
                    msh.get_parameterization_from_data,
                    msh
                )
                cp,_,_,_= cpfd_1(
                    data_path.joinpath(),
                    mvs,
                    svs,
                    dvs,
                    cpa,
                    epa_min,
                    epa_max,
                    epa_0,
                )
                
    def test_age_distributions_and_btt_start_in_ss_3(self):
        # get everything from mvs
        for mf in set(self.model_folders).intersection(["kv_visit2"]):
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs=import_module(f"{msh.model_mod}.source").mvs

                smr = mvs.get_SmoothModelRun()
                smr.initialize_state_transition_operator_cache(lru_maxsize=None)
                start_mean_age_vec=mvs.get_NumericStartMeanAgeVector()
                sv = mvs.get_StateVariableTuple()
                n_pools=len(sv)
                order = 1
                s_arr, s_func = smr._solve_age_moment_system(
                    order, 
                    start_mean_age_vec.reshape(1,n_pools)
                    # method='rk45' does not support kwargs yet
                )
                
                # the first n colums are the solution
                #solutions = smr.solve()
                solutions = s_arr
                m_a_arr = s_arr[:, n_pools:2 * n_pools]

                ## plot the continuous solution (by ODE) solver against the  iterator 
                ## generated one.
                times = mvs.get_NumericSimulationTimes()
                X_0 = mvs.get_NumericStartValueArray()
                t0 = times[0]
                par_dict = mvs.get_NumericParameterization().par_dict
                func_dict = mvs.get_NumericParameterization().func_dict
                bit = gh.traceability_iterator_result(
                    mvs,
                    X_0,
                    par_dict,
                    func_dict, 
                    delta_t_val=15,
                    t_0=t0
                )
                fig1 = plt.figure(figsize=(2*10, n_pools*10))
                axs = fig1.subplots(n_pools, 2)
                n_steps=12 # for testing
                vals=bit[0:n_steps]
                for i in range(n_pools):
                    ax = axs[i,0]
                    ax.plot(times,vals.X[:,i],label="bit")

                    ax.plot(times,solutions[:,i],label="sol")
                    ax.legend()
                    ax = axs[i,1]
                    ax.plot(times,m_a_arr[:,i])

                testDir = self.output_path(mf)
                self.__class__.clean_dir(testDir)
                fig1.savefig(
                    testDir.joinpath(
                        "poolwise.pdf"
                    )
                )

                fig2 = plt.figure(figsize=(10, 10))
                axs2 = fig2.subplots(2, 2)
                mean_btts= smr.backward_transit_time_moment(
                    order=1,
                    start_age_moments=start_mean_age_vec.reshape(1,n_pools)
                )
                ax = axs2[0,0]
                ax.plot(times, mean_btts,label="mean backward transit time")
                ax.plot(times, vals.system_RT_sum,label="$\sum_i (RT)_i$") #steady state transit times
                ax.plot(times, vals.rt,label="rt of surrogate one pool system") #steady state transit times
                ax.legend()
                fig2.savefig(
                            testDir.joinpath(
                                "system.pdf"
                            )
                )
                # construct a function p that takes an age array "ages" as argument
                # and gives back a three-dimensional ndarray (ages x times x pools)
                # from the a array-valued function of a single age a_dens_function
                srm = mvs.get_SmoothReservoirModel()
                a_dens_function, X_fix = start_age_distributions_from_steady_state(
                    srm,
                    t0=t0,
                    parameter_dict=par_dict,
                    func_set=func_dict,
                    x0=X_0
                )
                p = smr.pool_age_densities_func(a_dens_function)
                #from IPython import embed; embed()
                ages = np.linspace(0, (np.array(start_mean_age_vec,dtype=float).reshape(-1)).max()*2, 2)# 21)
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
                        mode='lines',
                        line=dict(
                            color='#FF0000',
                            width=15
                            )
                    )
                    smr.add_line_to_density_plot_plotly(
                        fig,
                        data=m_a_arr[:, n],
                        color='#FF0000',
                        name="mean age",
                        time_stride       = 1,
                        on_surface        = True,
                        bottom            = True,
                        legend_on_surface = True,
                        legend_bottom     = False,
                    )

                    plot(
                        fig, 
                        filename=str(
                            testDir.joinpath(
                                "age_distribution_{0}.html".format(sv[n])
                            )
                        ),
                        auto_open=False
                    )
                
                btt_dens = smr.backward_transit_time_density(age_densities)
                fig_btt = smr.plot_3d_density_plotly(
                    "backward_transit_time_density_steady_state",
                    btt_dens,
                    ages,
                    y_label="transit time"
                )
                smr.add_line_to_density_plot_plotly(
                    fig_btt,
                    data=mean_btts,
                    color='#FF0000',
                    name="mean age",
                    time_stride       = 1,
                    on_surface        = True,
                    bottom            = True,
                    legend_on_surface = True,
                    legend_bottom     = False
                )
                plot(
                    fig_btt, 
                    filename=str(testDir.joinpath("btt_distribution.html")), 
                    auto_open=False
                )
                #svs, dvs = msh.get_example_site_vars(Path(conf_dict["dataPath"]))
                #from IPython import embed; embed()
    
    def test_get_parameterization_from_data(self):
        #for mf in set(self.model_folders):
        for mf in set(self.model_folders).intersection(["yz_jules"]):
        #for mf in set(self.model_folders).intersection(["Aneesh_SDGVM"]):
        #for mf in set(self.model_folders).intersection(["kv_visit2"]):
        #for mf in set(self.model_folders).intersection(["jon_yib"]):
            with self.subTest(mf=mf):
                th = gh.th(mf)
                msh = gh.msh(mf)
                #CP= import_module(f"{msh.model_mod}.CachedParameterization").CachedParameterization
                mvs=import_module(f"{msh.model_mod}.source").mvs
                conf_dict=gh.confDict(mf)
                data_path=Path(conf_dict['dataPath'])
                tr_path=Path(mf).joinpath("data_assimilation_parameters_from_test_args")
                cpa= msh.Constants(
                    **h.load_dict_from_json_path(
                        tr_path.joinpath("cpa.json") 
                    )
                )

                epa_min,epa_max,epa_0=tuple(
                    map(
                        lambda p:msh.EstimatedParameters(**h.load_dict_from_json_path(p)),
                        [
                            tr_path.joinpath(f"{s}.json") 
                            for s in ['epa_min','epa_max','epa_0']
                        ]
                    )
                )
                svs, dvs = msh.get_global_mean_vars_2(conf_dict)
                # run different data assimilation procedures 
                for func in [
                    gh.get_parameterization_from_data,
                    #gh.get_parameterization_from_data,
                ]:
                    with self.subTest(func=func):
                        cpfd = gh.make_cached_data_assimilation_func(
                            func,   
                            msh,
                        )
                        cp,_,_,_= cpfd(
                            data_path,
                            msh,
                            mvs,
                            svs,
                            dvs,
                            cpa,
                            epa_min,
                            epa_max,
                            epa_0,
                            nsimu=15
                        )
    



    def test_parameterization_from_test_args(self):
        # temporary function to write the files to replace the production run information
        # stored in testargs up to now
        for mf in set(self.model_folders): #.intersection(["yz_jules"]):
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs=import_module(f"{msh.model_mod}.source").mvs
                test_args = gh.test_args_2(mf)
                CP= import_module(f"{msh.model_mod}.CachedParameterization").CachedParameterization
                #extract the additional parameters for make_func
                apa = {**test_args.cpa._asdict(), **test_args.epa_opt._asdict()}
                
                func_dict_param_dict = {
                    k: v 
                    for k, v in apa.items() 
                    if k in CP.func_dict_param_keys
                } if CP.func_dict_param_keys is not None else None
                name = "parameterization_from_test_args"
                p = Path(mf).joinpath(name)
                gh.write_parameterization_from_test_args(
                    test_args,
                    func_dict_param_dict,
                    CP,
                    p 
                )
                cp = CP.from_path(p)
                for k in test_args.dvs._fields:
                    self.assertTrue(
                        np.allclose(
                            cp.drivers.__getattribute__(k),
                            test_args.dvs.__getattribute__(k)
                        )
                    )

    # temporary function to write the files to replace the production run information
    # stored in the model testargs up to now
    def test_data_assimilation_parameters_from_test_args(self):
        for mf in set(self.model_folders).intersection(
            #["yz_jules","kv_visit2", "jon_yib"]#, "Aneesh_SDGVM"]
            ["yz_jules"]
        ):
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs=import_module(f"{msh.model_mod}.source").mvs
                f = gh.data_assimilation_parameters_from_test_args
                CachedParameterization = import_module(f"{msh.model_mod}.CachedParameterization").CachedParameterization
                dir_path = Path(mf).joinpath(f.__name__)
                try:
                    cp=CachedParameterization.from_path(dir_path)
                except:    
                    test_args = gh.test_args_2(mf)

                    f(
                        test_args,
                        CachedParameterization,
                        dir_path
                    )
   

    def test_minimal_iterator_args(self):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                cpa = test_args.cpa
                epa = test_args.epa_0
                dvs = test_args.dvs
                testDir = self.output_path(mf)
                self.__class__.clean_dir(testDir)
                sv = mvs.get_StateVariableTuple()
                func_dict = msh.make_func_dict(dvs , cpa=cpa, epa=epa)
                X0 = msh.numeric_X_0(mvs,test_args.dvs,test_args.cpa,test_args.epa_0)
                delta_t_val=5
                par_dict = gh.make_param_dict(mvs, cpa, epa)
                bit = gh.minimal_iterator(
                    X0,
                    func_dict, 
                    mvs,
                    dvs,
                    cpa,
                    epa,
                    delta_t_val=delta_t_val,
                )
                start = 0
                stop = 10 
                step = 2
                adr = ArrayDictResult(bit)
                vals = adr[start: stop]

                parts = gh.partitions(start, stop, step)
                res1 = adr.averaged_values(parts) # using the result class 
                res2 = vals.averaged_values(parts) # using the already computed valued d
                

                # if the number (stop-start) is a multiple of step then 
                # we can create an averaging iterator and wrap it with the
                # ArrayDictResult class which provides the index notation
                
                res3 = ArrayDictResult(hr.average_iterator(bit,2))[start:int(stop/step)] 
                #from IPython import embed; embed() 
                self.assertTrue(np.allclose(res1.X,res2.X))

                for k,v in cp.X_0_dict.items():
                    self.assertTrue(
                        np.allclose( v,test_args.V_init.__getattribute__(k))
                    )
                
                for k,v in cp.parameter_dict.items():
                    self.assertTrue(
                        np.allclose( v,test_args.par_dict[Symbol(k)])
                    )

    def test_make_da_iterator(self):
        #for mf in set(
        #        self.model_folders
        #    ).intersection(["kv_visit2","Aneesh_SDGVM","cj_isam","yz_jules","kv_ft_dlem","bian_ibis2"]):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                test_args = gh.test_args(mf)
                cpa = test_args.cpa
                epa = test_args.epa_0
                dvs = test_args.dvs
                testDir = self.output_path(mf)
                self.__class__.clean_dir(testDir)
                sv = mvs.get_StateVariableTuple()
                func_dict = msh.make_func_dict(dvs , cpa=cpa, epa=epa)
                X0 = msh.numeric_X_0(mvs,test_args.dvs,test_args.cpa,test_args.epa_0)
                delta_t_val=5
                par_dict = gh.make_param_dict(mvs, cpa, epa)
                bit = msh.make_da_iterator(
                    mvs,
                    X0,
                    par_dict,
                    func_dict, 
                    delta_t_val=delta_t_val,
                )
                start = 0
                stop = 10 
                step = 2
                adr = ArrayDictResult(bit)
                vals = adr[start: stop]
                #from IPython import embed; embed()

