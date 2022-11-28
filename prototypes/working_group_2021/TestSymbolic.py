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
    start_age_distributions_from_zero_initial_content
)
import CompartmentalSystems.helpers_reservoir as  hr
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
        dataPath = Path(confDict['dataPath'])
        return dataPath

    def output_path(self, mf):
        testDirPath = __class__.data_dir_path(mf).joinpath(self.id())
        return testDirPath

    @property
    def model_folders(self):
        return [
            "kv_visit2",
            "jon_yib",
            "Aneesh_SDGVM",
            #"cable-pop",
            #"cj_isam", # has problems with ODE solution probably due to wrong parameters
            # msh.numericX0 also yields a negative pool value for the last pool
            "yz_jules",
            "kv_ft_dlem",
            "bian_ibis2",
        ]

    def test_symobolic_description(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                mvs = gh.mvs(mf)
                # we assert that some variables are present
                mvs.get_SmoothReservoirModel()
                mvs.get_BibInfo()

    # @skip
    def test_download_data(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                conf_dict = gh.confDict(mf)
                msh = gh.msh(mf)
                msh.download_my_TRENDY_output(conf_dict)

    def test_get_example_site_vars(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                conf_dict = gh.confDict(mf)
                msh = gh.msh(mf)
                svs, dvs = msh.get_example_site_vars(Path(conf_dict["dataPath"]))
                # print(svs)

    def test_get_global_mean_vars(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                conf_dict = gh.confDict(mf)
                msh = gh.msh(mf)
                svs, dvs = msh.get_global_mean_vars(Path(conf_dict["dataPath"]))
                # print(svs)

    def test_make_func_dict_interface(self):
        # Purpose: this test checks that your make_func_dict accepts the
        # maximum set of parameters mvs,dvs,cpa,epa because some of the models
        # need some of them and we want to call the function in the same way
        # from the comparison notebooks
        #
        # How to fix your make_func_dict: If the test fails and tells you that
        # you called the function with too many parameters just add it to the
        # parameter list of your make_func_dict and ignore it in the rest of
        # the function...
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                conf_dict = gh.confDict(mf)
                test_args = gh.test_args(mf)
                svs, dvs = msh.get_example_site_vars(Path(conf_dict["dataPath"]))
                msh.make_func_dict(mvs, dvs, test_args.cpa, test_args.epa_0)

    def test_make_iterator_sym(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):

                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                conf_dict = gh.confDict(mf)
                test_args = gh.test_args(mf)
                t_min = 0
                t_max = 1400
                delta_t_val = 15 # this affects the precision of the iterator
                stride = 2 # this does not affects the precision of the iterator
                # but makes it more effiecient (the values in between the strides
                # are computed but not stored)
                # the presision of the continuous solution is not affected by either
                # number since the solver decides where to compute the next value
                # the times argument just tells it where WE want to know the values...

                n_steps = int((t_max - t_min)/delta_t_val)
                #times = np.linspace(t_min, t_max, n_steps)
                times = np.arange(t_min, t_max, delta_t_val*stride)
                # create a smooth model run to check the results against 
                mvs = mvs.update({
                    NumericParameterization(
                        par_dict=test_args.par_dict,
                        func_dict=test_args.func_dict,
                    ),
                    NumericStartValueArray(
                        msh.numeric_X_0(
                            mvs,
                            test_args.dvs,
                            test_args.cpa,
                            test_args.epa_0,
                        )
                    ),
                    NumericSimulationTimes(times)
                })
                
                smr=mvs.get_SmoothModelRun()
                solution=smr.solve()

                V_init = test_args.V_init
                it_sym_2 = msh.make_iterator_sym(
                    mvs=mvs,
                    V_init=V_init,
                    par_dict=test_args.par_dict,
                    func_dict=test_args.func_dict,
                    delta_t_val=delta_t_val,
                )
                ns = delta_t_val * 3
                times_2 = np.arange(0, ns, delta_t_val)
                res_2 = np.zeros((len(times_2), len(V_init)))
                res_2[0, :] = V_init
                for i in range(1, len(times_2) - 1):
                    res_2[i, :] = it_sym_2.__next__().reshape(
                        len(V_init),
                    )
                # create a smooth model run to check the results against 
                mvs = mvs.update({
                    NumericParameterization(
                        par_dict=test_args.par_dict,
                        func_dict=test_args.func_dict,
                    ),
                    NumericStartValueArray(
                        msh.numeric_X_0(
                            mvs,
                            test_args.dvs,
                            test_args.cpa,
                            test_args.epa_0,
                        )
                    ),
                    NumericSimulationTimes(times)
                })

    def test_param2res_sym(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                mvs = gh.mvs(mf)
                msh = gh.msh(mf)
                conf_dict = gh.confDict(mf)
                test_args = gh.test_args(mf)

                cpa = test_args.cpa
                dvs = test_args.dvs
                svs = test_args.svs
                epa_0 = test_args.epa_0
                param2res_sym = msh.make_param2res_sym(mvs, cpa, dvs)
                xs = param2res_sym(epa_0)

    # @skip
    # def test_autostep_mcmc_array_cost_func(self):
    #    # this test is only performed for certain models which have (or have created) monthly data
    #    # for all observed variables an can therefore use the simpler general costfunctions in general
    #    # helpers.
    #    # Most other models implement their own costfunctions in model_specific_helpers_2 and are
    #    # are tested with different arguments to the mcmc
    #    for mf in set(self.model_folders).intersection(['cj_isam']):
    #        #print("############################  {}  ############################".format(mf))
    #        with self.subTest(mf=mf):
    #            mvs = gh.mvs(mf)
    #            msh = gh.msh(mf)
    #            conf_dict=gh.confDict(mf)
    #            test_args = gh.test_args(mf)
    #            cpa = test_args.cpa
    #            dvs = test_args.dvs
    #            svs = test_args.svs
    #            epa_min = test_args.epa_min
    #            epa_max = test_args.epa_max
    #            epa_0 = test_args.epa_0

    #            isQualified = gh.make_param_filter_func(epa_max, epa_min)
    #            param2res = msh.make_param2res_sym( mvs, cpa, dvs)

    #            obs=test_args.obs_arr
    #            #obs=np.column_stack([ np.array(v) for v in svs])
    #            obs=obs[0:cpa.number_of_months,:] #cut
    #            # Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
    #            C_autostep, J_autostep = gh.autostep_mcmc(
    #                initial_parameters=epa_0,
    #                filter_func=isQualified,
    #                param2res=param2res,
    #                costfunction=gh.make_feng_cost_func(obs),
    #                nsimu=20, # for testing and tuning mcmc
    #                #nsimu=20000,
    #                c_max=np.array(epa_max),
    #                c_min=np.array(epa_min),
    #                acceptance_rate=15,   # default value | target acceptance rate in %
    #                chunk_size=10,  # default value | number of iterations to calculate current acceptance ratio and update step size
    #                D_init=1,   # default value | increase value to reduce initial step size
    #                K=2 # default value | increase value to reduce acceptance of higher cost functions
    #            )

    def test_autostep_mcmc_tupel_cost_func(self):
        for mf in set(self.model_folders).intersection(["kv_visit2"]):
            with self.subTest(mf=mf):
                mvs = gh.mvs(mf)
                msh = gh.msh(mf)
                conf_dict = gh.confDict(mf)
                test_args = gh.test_args(mf)
                cpa = test_args.cpa
                dvs = test_args.dvs
                svs = test_args.svs
                epa_min = test_args.epa_min
                epa_max = test_args.epa_max
                epa_0 = test_args.epa_0

                isQualified = gh.make_param_filter_func(epa_max, epa_min)
                param2res = msh.make_param2res_sym(mvs, cpa, dvs)
                # obs=np.column_stack([ np.array(v) for v in svs])
                # obs=obs[0:cpa.number_of_months,:] #cut
                # Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
                C_autostep, J_autostep = gh.autostep_mcmc(
                    initial_parameters=epa_0,
                    filter_func=isQualified,
                    param2res=param2res,
                    # costfunction=msh.make_feng_cost_func_2(svs),
                    costfunction=gh.make_feng_cost_func_2(svs),
                    nsimu=20,  # for testing and tuning mcmc
                    # nsimu=20000,
                    c_max=np.array(epa_max),
                    c_min=np.array(epa_min),
                    acceptance_rate=15,  # default value | target acceptance rate in %
                    chunk_size=10,  # default value | number of iterations to calculate current acceptance ratio and update step size
                    D_init=1,  # default value | increase value to reduce initial step size
                    K=1,  # default value | increase value to reduce acceptance of higher cost functions
                )

    # def test_autostep_mcmc_with_model_specific_costfunction(self):
    #
    #    for mf in set(self.model_folders).intersection(['Aneesh_SDGVM']):
    #        with self.subTest(mf=mf):
    #            mvs = gh.mvs(mf)
    #            msh = gh.msh(mf)
    #            conf_dict=gh.confDict(mf)
    #            test_args = gh.test_args(mf)
    #            cpa = test_args.cpa
    #            dvs = test_args.dvs
    #            svs = test_args.svs
    #            epa_min = test_args.epa_min
    #            epa_max = test_args.epa_max
    #            epa_0 = test_args.epa_0

    #            isQualified = gh.make_param_filter_func(epa_max, epa_min)
    #            param2res = msh.make_param2res_sym( mvs, cpa, dvs)
    #            # Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
    #            C_autostep, J_autostep = gh.autostep_mcmc(
    #                initial_parameters=epa_0,
    #                filter_func=isQualified,
    #                param2res=param2res,
    #                costfunction=msh.make_weighted_cost_func(svs),
    #                nsimu=20, # for testing and tuning mcmc
    #                #nsimu=20000,
    #                c_max=np.array(epa_max),
    #                c_min=np.array(epa_min),
    #                acceptance_rate=15,   # default value | target acceptance rate in %
    #                chunk_size=10,  # default value | number of iterations to calculate current acceptance ratio and update step size
    #                D_init=1,   # default value | increase value to reduce initial step size
    #                K=2 # default value | increase value to reduce acceptance of higher cost functions
    #            )

    def test_epa_opt_presence(self):
        # Purpose:
        # For the model comparison it's nice to have your best and shiniest parameters ;-)
        # available without having to go through the dataassimilation.
        #
        # How to fix this test:
        # 1.)   print out your optimal parameters after the data assimilation (e.g from your inspectModel.py)
        # 2.)   Add the epa_opt field to your TestArgs tupel definition in your test_helpers
        # 3.)   To avoid  your make_test_args to throw an error at you you have to set the epa_opt value...
        #       Just paste your painstakenly obtained  parameters there and change them if you have
        #       new even shinier ones
        #       you can look at Kostia's "kv_visit2/test_helpers.py
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                mvs = gh.mvs(mf)
                msh = gh.msh(mf)
                test_args = gh.test_args(mf)
                test_args.epa_opt

    def test_start_date_presence(self):
        # Purpose:
        # For the model comparison plots it's necessary to have a common timeline
        #
        # How to fix this test:
        # 1.)   implement a function start_date in your model_specific_helpers_2 file
        #       you can look at Kostia's "kv_visit2/model_specific_helpers_2.py as an
        #       example
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                mvs = gh.mvs(mf)
                msh = gh.msh(mf)
                msh.start_date()

    # def test_make_param_filter_func_presence(self):
    #    #Purpose:
    #    # Although we do not have to

    def test_numeric_X0(self):
        # Purpose:
        # This function assembles the startvector for several iterators
        # especially the one that computes the variables for the tracebility analysis.
        # but could also be used in your param2res (for the pools part)
        #
        # How to make it work:
        # look at kv_visit2/model_specific_helpers_2.py
        # or yz_jules/model_specific_helpers_2.py
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                mvs = gh.mvs(mf)
                msh = gh.msh(mf)
                ta = gh.test_args(mf)
                mvs_t = gh.mvs(mf)
                dvs_t = ta.dvs
                cpa_t = ta.cpa
                epa_t = ta.epa_0
                X_0 = gh.msh(mf).numeric_X_0(mvs_t, dvs_t, cpa_t, epa_t)
                # it is important that  X_0 and the result of I(t,X) have 
                # the same dimension
                # fixme mm-11-17-2022 
                # the older iterators rely on a n,1 matrix
                # where as CompartmentalSystems expects a n, vector
                # until the old code is completely removed we check for n,1
                self.assertEqual(X_0.ndim,2)

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
                #for i in range(n_lats):
                #    self.assertEqual(tr.lat2i(lats[i]), i)
                ## or longitude
                #for i in range(n_lons):
                #    self.assertEqual(tr.lon2i(lons[i]), i)

                ## check the interpretation of the pixel boundaries
                #for i in range(n_lats - 1):
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


    def test_mask(self):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                # test_args = gh.test_args(mf)
                msh = gh.msh(mf)

                c_mask = msh.spatial_mask(Path(gh.confDict(mf)["dataPath"]))
                # c_mask.write_netCDF4(mf+".nc")
                f = plt.figure()
                ax = f.add_subplot(1, 1, 1)
                c_mask.plot_dots(ax)
                f.savefig(str(mf) + ".pdf")


    def test_func_dict_old_vs_interp(self):
        # at the moment only for visit
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                conf_dict = gh.confDict(mf)

                fd = msh.make_func_dict(mvs,test_args.dvs, test_args.cpa, test_args.epa_0)
                fd_old=msh.make_func_dict_old(mvs,test_args.dvs, test_args.cpa, test_args.epa_0)
                months=range(0,len(test_args.dvs.npp))
                days=[m*30 for m in months]
                for fn in fd.keys():
                    arr=np.array(list(map(fd[fn],days)))
                    arr_old=np.array(list(map(fd_old[fn],days)))
                    
                    self.assertTrue(
                        np.allclose(
                            arr,
                            arr_old,
                            rtol=0.01
                        )
                    )

    def test_iterator_instance(self):
        # this test just makes sure that Kostias plot functions still get the data
        # at the moment only for visit
        delta_t_val = 15
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                itr = moh.traceability_iterator_instance(
                    mf,
                    test_args=gh.test_args(mf),
                    delta_t_val=delta_t_val
                )
                start, stop = 0, 5
                vals = itr[start:stop]
                #check the availablility of certain properties 
                target_props=[
                    "x",
                    "x_c",
                    "x_p",
                    "u",
                    "rt",
                ]
                #from IPython import embed; embed()
                for p in target_props:
                    #v1=vals[p]
                    # check the . interface
                    v2=vals.__getattribute__(p)


    def test_age_distributions_and_btt_start_in_ss_2(self):
        # mf = "kv_visit2"
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                #subs_1= import_module("{}.subs_1".format(mf)) 
                #s1 = import_module("{}.source_1".format(mf)) 
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                cpa = test_args.cpa
                epa = test_args.epa_0
                dvs = test_args.dvs
                testDir = self.output_path(mf)
                self.__class__.clean_dir(testDir)
                sv = mvs.get_StateVariableTuple()
                n_pools = len(sv)
                func_dict = msh.make_func_dict(mvs,dvs , cpa, epa)
                # X0 = msh.numeric_X_0(mvs,test_args.dvs,test_args.cpa,test_args.epa)
                par_dict = gh.make_param_dict(mvs, cpa, epa)
                # compute the start age distribution
                srm = mvs.get_SmoothReservoirModel()
                t0 = 0 #3/2*np.pi
                # we assume that the system was in steady state 
                # at t_0 with X_fix given by X_fix = M^{-1} I
                # This might not be a good fit for the data since the data
                # assimilation does not necessaryly make use of the steady state assumption
                # (If it does the startvector does not contain any estimable components)
                a_dens_function, X_fix = start_age_distributions_from_steady_state(
                    srm,
                    t0=t0,
                    parameter_dict=par_dict,
                    func_set=func_dict
                )
                X_fix=X_fix.reshape(-1)
                # create a model run that starts at X_fix and t0
                n_steps = 12  # 2881
                t_max = 144 
                times = np.linspace(t0, t_max, n_steps)
                delta_t_val = (t_max - t0)/n_steps
                mvs = mvs.update({
                    NumericParameterization(
                        par_dict=gh.make_param_dict(mvs, cpa, epa),
                        func_dict=func_dict
                    ),
                    NumericStartValueArray(X_fix),
                    NumericSimulationTimes(times)
                })
                # smr = mvs.get_SmoothModelRun()
                smr = SmoothModelRun(
                    srm,
                    parameter_dict=par_dict,
                    start_values=X_fix,
                    times=times,
                    func_set=func_dict
                )
                smr.initialize_state_transition_operator_cache(lru_maxsize=None)
                # We cant compute the mean by integration since for this example the numeric integration does not
                # converge) It usually has no problem with an age distribution
                # given as an explicit function but we can compute the mean age
                # directly as (1,..,1)B^{-1} 
                M0 = mvs.get_NumericCompartmentalMatrixFunc()(t0,X_fix)
                M0_inv = np.linalg.inv(M0)
                start_mean_age_vec = -1*(np.ones_like(X_fix) @ M0_inv)
                print(start_mean_age_vec)
                
                #compute solutions for the mean age system starting at X_fix and star
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

                # alternatively we can compute the symbolic inverse
                #tau = s1.mvs.get_LuoTau()
                #tau_num=gh.numfunc(tau,mvs,delta_t_val,par_dict,func_dict)
                #start_mean_age_vec = np.ones_like(X_fix)@tau_num(0,*X_fix)
                #system_rt = LuoRT.dot(np.ones_like(X_fix))
                # plot the continuous solution (by ODE) solver against the  iterator 
                # generated one.
                bit = gh.traceability_iterator(
                    X_fix,
                    func_dict, 
                    mvs,
                    dvs,
                    cpa,
                    epa,
                    delta_t_val=delta_t_val,
                )
                fig1 = plt.figure(figsize=(2*10, n_pools*10))
                axs = fig1.subplots(n_pools, 2)
                vals=bit[0:n_steps]
                for i in range(n_pools):
                    ax = axs[i,0]
                    ax.plot(times,vals.X[:,i],label="bit")

                    ax.plot(times,solutions[:,i],label="sol")
                    ax.legend()
                    ax = axs[i,1]
                    ax.plot(times,m_a_arr[:,i])

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
                ax.plot(times,vals.RT.sum(axis=1),label="$\sum_i (RT)_i$") #steady state transit times
                ax.plot(times,vals.rt,label="rt of surrogate one pool system") #steady state transit times
                ax.legend()
                fig2.savefig(
                            testDir.joinpath(
                                "system.pdf"
                            )
                )
                # construct a function p that takes an age array "ages" as argument
                # and gives back a three-dimensional ndarray (ages x times x pools)
                # from the a array-valued function of a single age a_dens_function
                p = smr.pool_age_densities_func(a_dens_function)
                ages = np.linspace(0, np.max(start_mean_age_vec)*2,2)# 21)
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
                

    @skip 
    def test_age_distributions_and_btt_from_zero_initial_content(self):
        pass
        #a_dens_function, X_fix = start_age_distributions_from_zero_initial_content(srm)
            
    
    @skip 
    def test_age_distributions_and_btt_from_empty_spinup(self):
        pass
        #a_dens_function, X_fix = start_age_distributions_from_empty_spinup(srm)

    def test_aggregate_surrogate_systems(self):
        # at the moment only for visit
        #mf = "kv_visit2"
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                #subs_1 = import_module("{}.subs_1".format(mf)) 
                #s1 = import_module("{}.source_1".format(mf)) 
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                state_vector=mvs.get_StateVariableTuple()
                n_pools=len(state_vector)
                cpa = test_args.cpa
                epa = test_args.epa_opt
                dvs = test_args.dvs
                testDir = self.output_path(mf)
                self.__class__.clean_dir(testDir)
                conf_dict = gh.confDict(mf)
                t_min = 0
                t_max = 1400
                delta_t_val = 5 # this affects the precision of the iterator
                stride = 6 # this does not affects the precision of the iterator
                # but makes it more effiecient (the values in between the strides
                # are computed but not stored)
                # the presision of the continuous solution is not affected by either
                # number since the solver decides where to compute the next value
                # the times argument just tells it where WE want to know the values...

                n_steps = int((t_max - t_min)/delta_t_val)
                #times = np.linspace(t_min, t_max, n_steps)
                times = np.arange(t_min, t_max, delta_t_val*stride)
                par_dict=gh.make_param_dict(mvs, test_args.cpa, epa)
                func_dict=msh.make_func_dict(mvs,test_args.dvs, test_args.cpa, epa)
                X0 = msh.numeric_X_0(mvs,test_args.dvs,test_args.cpa,epa).reshape(-1)
                # in our case we want to compute a steady state start value
                #from IPython import embed;embed()
                a_dens_function, X_fix = start_age_distributions_from_steady_state(
                    t0=0,
                    srm=mvs.get_SmoothReservoirModel(),
                    parameter_dict=par_dict,
                    func_set=func_dict
                )
                mvs=mvs.update({
                    NumericParameterization(
                        par_dict=par_dict,
                        func_dict=func_dict
                    ),
                    NumericStartValueArray(X_fix),
                    NumericSimulationTimes(times)
                })
                t0 = times[0]
                sv = mvs.get_StateVariableTuple()
                n_pools = len(sv)
                smr=mvs.get_SmoothModelRun()
                from IPython import embed; embed()
                # now create two submodels
                M0 = mvs.get_NumericCompartmentalMatrixFunc()(t0,X_fix)
                print(M0)
                M0_inv = np.linalg.inv(M0)
                start_mean_age_vec = start_age_moments_from_steady_state(
                    srm,
                    t0=0,
                    parameter_dict=par_dict,
                    func_set=func_dict,
                    max_order=1
                ).reshape(-1)
                #compute solutions for the mean age system starting at X_fix and star
                order = 1
                s_arr, s_func = smr._solve_age_moment_system(
                    order, 
                    start_mean_age_vec.reshape(1,n_pools),
                )
                # the first n colums are the solution
                #solutions = smr.solve()
                solutions = s_arr
                m_a_arr = s_arr[:, n_pools:2 * n_pools]
                t=mvs.get_TimeSymbol()
                def sub_mr_smav(smr, svt): 
                    combined = (
                        set(mvs.get_StateVariableTuple()),
                        mvs.get_InFluxesBySymbol(),
                        mvs.get_OutFluxesBySymbol(),
                        mvs.get_InternalFluxesBySymbol()
                    )

                    svl=list(sv)
                    sv_set, in_fluxes, out_fluxes, internal_fluxes = hr.extract(
                        combined,
                        set(svt)
                    )    

                    # we have to provide solutions (functions of time)  for the
                    # state variables of the outer system ( In our case soil
                    # influxes depend on the poolvalues of the veg system. In
                    # general every subsystem could depend on pools of every other
                    # subsystem
                    outer_pools=set(sv).difference(sv_set)
                    # first replace the symbols by functions
                    subs_dict={sym: Function(str(sym))(t) for sym in outer_pools}
                    in_fluxes_f, internal_fluxes_f, out_fluxes_f = map(
                        lambda d: {
                            k: sympify(flux_ex).subs(subs_dict) 
                            for k,flux_ex in d.items()
                        },
                        (in_fluxes, internal_fluxes, out_fluxes)
                    )    
                    # now prepare the extension of the func_dict by the solutions
                    # for the variables
                    outer_pool_funcs={
                        Function(str(sym)): lambda t: s_func(t)[svl.index(sym)]
                        for sym in outer_pools
                    }
                    sub_func_dict = { **func_dict, **outer_pool_funcs}
                    sub_X_fix = X_fix[[svl.index(w) for w in svt]]
                    sub_mvs = CMTVS(
                        {
                            t,
                            StateVariableTuple(svt),
                            InFluxesBySymbol(in_fluxes_f),
                            InternalFluxesBySymbol(internal_fluxes_f),
                            OutFluxesBySymbol(out_fluxes_f),
                            NumericParameterization(
                                par_dict=par_dict,
                                func_dict=sub_func_dict,
                            ),
                            NumericStartValueArray(sub_X_fix),
                            NumericSimulationTimes(times)
                        },
                        mvs.computers
                    )
                    return sub_mvs.get_SmoothModelRun()
                
                def sub_btt(sub_svt):
                    sub_smr = sub_mr(
                        mvs.get_StateVariableTuple(),
                        sub_svt,
                        combined,
                        s_func,
                        X_fix,
                    )    
                    # for the start mean ages we can not just take the respective parts of the system mean ages
                    # since they refer to the times stuff has spent in the SYSTEM and not in the SUB system.
                    # Although the pools in the subsystem are pools in the 
                    sub_start_mean_age_vec = start_age_moments_from_steady_state(
                        sub_smr.model,
                        t0=0,
                        parameter_dict=par_dict,
                        func_set=sub_func_dict,
                        max_order=1
                    )
                    return  sub_smr.backward_transit_time_moment(
                        order=1,
                        start_age_moments=sub_start_mean_age_vec
                    )
                veg_res = sub_btt(mvs.get_VegetationCarbonStateVariableTuple())
                sub_smr.backward_transit_time_moment(
                    order=1,
                    start_age_moments=sub_smav
                )
                soil_res= sub_system_results(mvs.get_SoilCarbonStateVariableTuple())
                # compute the start age distribution

                #smr.initialize_state_transition_operator_cache(lru_maxsize=None)
                order = 1
                # From a_dens_function we should be able to extract the mean age
                # but the numeric integration does not converge for this
                # numerically computed example.  It usually has no problem
                # with an age distribution given as an explicit function but we can
                # compute the mean age directly as (1,..,1)B^{-1} 
                print(start_mean_age_vec)
                ##ages = np.linspace(-1,10,12)
                #ages = np.linspace(0,500,3)
                # construct a function p that takes an age array "ages" as argument
                # and gives back a three-dimensional ndarray (ages x times x pools)
                # from the a array-valued function representing the start age density
                #p = smr.pool_age_densities_func(
                #    a_dens_function
                #    #start_age_distributions_from_zero_initial_content(srm)
                #    #lambda a: np.exp(-a)*X0
                #)
                #age_densities = p(ages)
                #btt_dens = smr.backward_transit_time_density(age_densities)
                #for n in range(srm.nr_pools):
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
                #LuoRT=s1.mvs.get_LuoRT()
                #Ot=s1.mvs.get_OutputTuple()

                #system_rt=LuoRT.dot(np.ones_like(X_fix))
                bit = gh.traceability_iterator(
                    X_fix,
                    func_dict, 
                    mvs,
                    dvs,
                    test_args.cpa,
                    epa,
                    delta_t_val=delta_t_val,
                )
                vals = bit[0:(n_steps+1):stride] 
                disc_times = vals.t
                # from IPython import embed; embed()

                fig2 = plt.figure(figsize=(20, 20))
                axs2 = fig2.subplots(4, 1)
                mean_btts= smr.backward_transit_time_moment(
                    order=1,
                    start_age_moments=start_mean_age_vec.reshape(1,n_pools)
                )
                color_dict={"veg":"green", "soil": "brown", "system": "black"}
                marker_dict={"RT":"*", "btt": "+", "tot": "o"}
                #############################################################
                ax = axs2[0]
                ax.set_title("turnover times vs $\sum_i (RT)_i$")
                ax.plot(
                    disc_times, 
                    vals.RT.sum(axis=1),
                    label="$\sum_i (RT)_i$",
                    color=color_dict["system"],
                    marker=marker_dict["RT"]
                )
                ax.plot(
                    disc_times,
		            vals.rt,
                    color=color_dict["system"],
                    marker=marker_dict["tot"],
		            label="tot = rt of surrogate one pool system"
                )
                ax.plot(
                    disc_times,
		            vals.tot_veg,
		            label="tot of veg pool",
                    color=color_dict["veg"],
                    marker=marker_dict["tot"]
                )
                key="tot_soil"
                ax.plot(
                    disc_times,
		    		vals[key],
		    		label=key,
                    color=color_dict["soil"],
                    marker=marker_dict["tot"]
                )
                ax.legend()
                
                #############################################################
                ax = axs2[1]
                ax.set_title("mean transit times vs $\sum_i (RT)_i$")
                ax.plot(
                    times,
		    soil_res,
                    color=color_dict["soil"],
                    marker=marker_dict["btt"],
		    label="mean soil sub system backward transit time"
                )
                ax.plot(
                    times,
		    mean_btts,
                    color=color_dict["system"],
                    marker=marker_dict["btt"],
		    label="mean system backward transit time"
                )
                for i,r in enumerate((veg_res,)):
                    ax.plot(
                        times,
		        r,
                        color=color_dict["veg"],
                        marker=marker_dict["btt"],
		        label=f"mean veg {0} subsystem backward transit time".format(i)
                    )
                for i,r in enumerate((soil_res,)):
                    ax.plot(
                        times,
		        r,
                        color=color_dict["soil"],
                        marker=marker_dict["btt"],
		        label=f"mean soil {0} sub system backward transit time".format(i)
                    )
                ax.plot(
                    times,
		    vals.RT.sum(axis=1),
                    color=color_dict["system"],
                    marker=marker_dict["RT"],
		    label="$\sum_i (RT)_i$"
                ) #steady state transit times
                ax.legend()
                
                #############################################################
                ax = axs2[2]
                ax.set_title("Fluxes")
                for l in ["out_2_veg","veg_2_out", "veg_2_soil", "soil_2_out"]:
                    ax.plot(
                        times,
		    	vals[l],
		    	label='{0}'.format(l)
                    )
                ax.legend()

                
                #############################################################
                ax = axs2[3]
                ax.set_title("Stocks")
                ax.plot(
                    times,
		    vals.x,
                    color=color_dict["system"],
		    label="$x$"
                ) 
                ax.plot(
                    times,
		    vals.x_veg,
                    color=color_dict["veg"],
		    label="$x_veg$"
                ) 
                ax.plot(times,
		    vals.x_soil,
                    color=color_dict["soil"],
		    label="$x_soil$") 
                ax.legend()
                fig2.savefig(testDir.joinpath( "system.pdf"))
                
                fig1 = plt.figure(figsize=(2*10, n_pools*10))
                axs = fig1.subplots(n_pools, 2)
                vals=bit[0:n_steps]
                for i in range(n_pools):
                    ax = axs[i,0]
                    #ax.plot(times,vals.X[:,i],label="bit")

                    ax.plot(times,solutions[:,i],label="sol")
                    ax.legend()
                    ax = axs[i,1]
                    ax.plot(times,m_a_arr[:,i])

                fig1.savefig(
                    testDir.joinpath(
                        "poolwise.pdf"
                    )
                )
                #svs, dvs = msh.get_example_site_vars(Path(conf_dict["dataPath"]))
                ## create the aggregated veg soil model
                #C__Veg,C__Soil=map(Symbol,["C__Veg","C__Soil"])
                #mvs_vs=CMTVS(
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
                #)

                #h.compartmental_graph(mvs_vs)
                #func_dict=msh.make_func_dict(mvs,dvs, test_args.cpa, test_args.epa_0)
                # get the continuous system to check


    def test_aggregate_surrogate_systems_2(self):
        # at the moment only for visit
        #mf = "kv_visit2"
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                #subs_1 = import_module("{}.subs_1".format(mf)) 
                #s1 = import_module("{}.source_1".format(mf)) 
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                state_vector=mvs.get_StateVariableTuple()
                n_pools=len(state_vector)
                cpa = test_args.cpa
                epa = test_args.epa_opt
                dvs = test_args.dvs
                testDir = self.output_path(mf)
                self.__class__.clean_dir(testDir)
                conf_dict = gh.confDict(mf)
                t_min = 0
                t_max = 1400*10
                delta_t_val = 5 # this affects the precision of the iterator
                stride = 6 # this does not affects the precision of the iterator
                # but makes it more effiecient (the values in between the strides
                # are computed but not stored)
                # the presision of the continuous solution is not affected by either
                # number since the solver decides where to compute the next value
                # the times argument just tells it where WE want to know the values...

                n_steps = int((t_max - t_min)/delta_t_val)
                #times = np.linspace(t_min, t_max, n_steps)
                times = np.arange(t_min, t_max, delta_t_val*stride)
                par_dict=gh.make_param_dict(mvs, test_args.cpa, epa)
                func_dict=msh.make_func_dict(mvs,test_args.dvs, test_args.cpa, epa)
                #X0 = msh.numeric_X_0(mvs,test_args.dvs,test_args.cpa,epa).reshape(-1)
                # in our case we want to compute a steady state start value
                #from IPython import embed;embed()
                a_dens_function, X_fix = start_age_distributions_from_steady_state(
                    srm=mvs.get_SmoothReservoirModel(),
                    t0=times[0],
                    parameter_dict=par_dict,
                    func_set=func_dict
                )
                mvs=mvs.update({
                    NumericParameterization(
                        par_dict=par_dict,
                        func_dict=func_dict
                    ),
                    NumericStartValueArray(X_fix),
                    NumericSimulationTimes(times)
                })
                t0 = times[0]
                sv = mvs.get_StateVariableTuple()
                n_pools = len(sv)
                smr=mvs.get_SmoothModelRun()
                srm = smr.model
                # now create two submodels
                combined = (
                    set(mvs.get_StateVariableTuple()),
                    mvs.get_InFluxesBySymbol(),
                    mvs.get_OutFluxesBySymbol(),
                    mvs.get_InternalFluxesBySymbol()
                )
                start_mean_age_vec = start_age_moments_from_steady_state(
                    srm,
                    t0=0,
                    parameter_dict=par_dict,
                    func_set=func_dict,
                    max_order=1
                ).reshape(-1)

                #Ot=s1.mvs.get_OutputTuple()

                #system_rt=LuoRT.dot(np.ones_like(X_fix))
                order = 1
                s_arr, s_func = smr._solve_age_moment_system(
                    order, 
                    start_mean_age_vec.reshape(1,n_pools),
                )
                # the first n colums are the solution
                #solutions = smr.solve()
                solutions = s_arr
                m_a_arr = s_arr[:, n_pools:2 * n_pools]
                t=mvs.get_TimeSymbol()

                    

                #smr.initialize_state_transition_operator_cache(lru_maxsize=None)
                print(start_mean_age_vec)
                bit = gh.traceability_iterator(
                    X_fix,
                    func_dict, 
                    mvs,
                    dvs,
                    test_args.cpa,
                    epa,
                    delta_t_val=delta_t_val,
                )
                vals = bit[0:(n_steps+1):stride] 
                disc_times = vals.t

                fig1 = plt.figure(figsize=(2*10, n_pools*10))
                axs = fig1.subplots(n_pools, 2)
                for i in range(n_pools):
                    ax = axs[i,0]
                    ax.plot(disc_times,vals.X[:,i],label="bit")
                    ax.plot(times,solutions[:,i],label="sol")
                    ax.legend()
                    #++++++++++++++++++++
                    ax = axs[i,1]
                    ax.plot(times,m_a_arr[:,i])

                fig1.savefig(
                    testDir.joinpath(
                        "poolwise.pdf"
                    )
                )
                #############################################################
                #############################################################

                fig2 = plt.figure(figsize=(20, 20))
                axs2 = fig2.subplots(1, 1)
                mean_btts= smr.backward_transit_time_moment(
                    order=1,
                    start_age_moments=start_mean_age_vec.reshape(1,n_pools)
                )
                color_dict={"veg":"green", "soil": "brown", "system": "black"}
                marker_dict={"RT":"*", "btt": "+", "tot": "o"}
                #############################################################
                ax = axs2#[0]
                ax.set_title("turnover times vs $\sum_i (RT)_i$ vs$ vs. mean transit times")
                ax.plot(
                    disc_times, 
                    vals.RT.sum(axis=1),
                    label="$\sum_i (RT)_i$",
                    color=color_dict["system"],
                    marker=marker_dict["RT"]
                )
                ax.plot(
                    disc_times,
		            vals.rt,
                    color=color_dict["system"],
                    marker=marker_dict["tot"],
		            label="tot = rt of surrogate one pool system"
                )
                
                
                ax.plot(
                    times,
		    		mean_btts,
                    color=color_dict["system"],
                    marker=marker_dict["btt"],
		    		label="mean system backward transit time"
                )
                ax.legend()
                
                fig2.savefig(testDir.joinpath( "system.pdf"))
                

