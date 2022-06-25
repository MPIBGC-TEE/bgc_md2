# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
import matplotlib.pyplot as plt
import numpy as np
from unittest.case import TestCase, skip
from pathlib import Path
from importlib import import_module
import json 

import general_helpers as gh


class TestSymbolic(TestCase):
    
    @property
    def model_folders(self):
        return [
            'kv_visit2',
	     'jon_yib',
	     'Aneesh_SDGVM',
	     'cable-pop',
	     'cj_isam',
	     'yz_jules',
	     'kv_ft_dlem'
        ]

    def test_symobolic_description(self):
        for mf in self.model_folders: 
            with self.subTest(mf=mf):
                mvs = gh.mvs(mf)
                # we assert that some variables are present
                mvs.get_SmoothReservoirModel()
                mvs.get_BibInfo()
    
    #@skip
    def test_download_data(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                conf_dict=gh.confDict(mf)
                msh = gh.msh(mf)
                msh.download_my_TRENDY_output(conf_dict)
    
    def test_get_example_site_vars(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                conf_dict=gh.confDict(mf)
                msh = gh.msh(mf)
                svs, dvs = msh.get_example_site_vars(Path(conf_dict['dataPath']))
                #print(svs)
    

    def test_get_global_mean_vars(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                conf_dict=gh.confDict(mf)
                msh = gh.msh(mf)
                svs, dvs = msh.get_global_mean_vars(Path(conf_dict['dataPath']))
                #print(svs)
    
    def test_make_func_dict(self):
        # Purpose:
        # this test checks that your make_func_dict accepts the maximum set of parameters
        # mvs,dvs,cpa,epa because some of the models need some of them 
        # and we want to call the function in the same way from the comparison notebooks
        # 
        # How to fix your make_func_dict:
        # If the test fails and tells you that you called the function with too many parameters
        # just add it to the parameter list of your make_func_dict and ignore it in the rest 
        # of the function...
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                conf_dict=gh.confDict(mf)
                test_args = gh.test_args(mf)
                svs, dvs = msh.get_example_site_vars(Path(conf_dict['dataPath']))
                msh.make_func_dict(mvs, dvs, test_args.cpa, test_args.epa_0)

    def test_make_iterator_sym(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                
                msh = gh.msh(mf)
                mvs = gh.mvs(mf)
                conf_dict=gh.confDict(mf)
                test_args = gh.test_args(mf)
                delta_t_val=30 #n_day iterator
                V_init=test_args.V_init
                it_sym_2 = msh.make_iterator_sym(
                    mvs=test_args.mvs,
                    V_init=V_init,
                    par_dict=test_args.par_dict,
                    func_dict=test_args.func_dict,
                    delta_t_val=delta_t_val
                )
                ns=delta_t_val*3
                times_2= np.arange(0,ns,delta_t_val)
                res_2= np.zeros((len(times_2),len(V_init)))
                res_2[0,:]=V_init
                for i in range(1,len(times_2)-1):
                    res_2[i,:]=it_sym_2.__next__().reshape(len(V_init),)

    def test_param2res_sym(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                mvs = gh.mvs(mf)
                msh = gh.msh(mf)
                conf_dict=gh.confDict(mf)
                test_args = gh.test_args(mf)

                cpa = test_args.cpa
                dvs = test_args.dvs
                svs = test_args.svs
                epa_0 = test_args.epa_0
                param2res_sym = msh.make_param2res_sym( mvs, cpa, dvs)
                xs= param2res_sym(epa_0)


    #@skip
    #def test_autostep_mcmc_array_cost_func(self):
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
        for mf in set(self.model_folders).intersection(['kv_visit2']):
            with self.subTest(mf=mf):
                mvs = gh.mvs(mf)
                msh = gh.msh(mf)
                conf_dict=gh.confDict(mf)
                test_args = gh.test_args(mf)
                cpa = test_args.cpa
                dvs = test_args.dvs
                svs = test_args.svs
                epa_min = test_args.epa_min
                epa_max = test_args.epa_max
                epa_0 = test_args.epa_0

                isQualified = gh.make_param_filter_func(epa_max, epa_min)
                param2res = msh.make_param2res_sym( mvs, cpa, dvs)
                #obs=np.column_stack([ np.array(v) for v in svs])
                #obs=obs[0:cpa.number_of_months,:] #cut
                # Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
                C_autostep, J_autostep = gh.autostep_mcmc(
                    initial_parameters=epa_0,
                    filter_func=isQualified,
                    param2res=param2res,
                    #costfunction=msh.make_feng_cost_func_2(svs),
                    costfunction=gh.make_feng_cost_func_2(svs),
                    nsimu=20, # for testing and tuning mcmc
                    #nsimu=20000,
                    c_max=np.array(epa_max),
                    c_min=np.array(epa_min),
                    acceptance_rate=15,   # default value | target acceptance rate in %
                    chunk_size=10,  # default value | number of iterations to calculate current acceptance ratio and update step size
                    D_init=1,   # default value | increase value to reduce initial step size
                    K=1 # default value | increase value to reduce acceptance of higher cost functions
                )
    
    #def test_autostep_mcmc_with_model_specific_costfunction(self):
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



    #def test_make_param_filter_func_presence(self):
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
                ta=gh.test_args(mf)
                mvs_t=gh.mvs(mf)
                dvs_t=ta.dvs
                cpa_t=ta.cpa
                epa_t=ta.epa_0
                X_0=gh.msh(mf).numeric_X_0(mvs_t,dvs_t,cpa_t,epa_t)



    def test_make_model_index_transforms(self):
        model_folders=['kv_visit2','kv_ft_dlem','Aneesh_SDGVM','cj_isam','jon_yib']#,'kv_visit2', 'jon_yib''Aneesh_SDGVM','cable-pop','yz_jules',]
        #for mf in set(self.model_folders):
        for mf in set(model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                msh = gh.msh(mf)
                lats=test_args.lats.data
                lons=test_args.lons.data
                # check that we predict the
                # right latitude values for
                # a given array index 
                tr = msh.make_model_index_transforms()
                n_lats=len(lats)
                print("n_lats={}".format(n_lats))
                for i in range(n_lats):
                    self.assertEqual(lats[i],tr.i2lat(i))
                
                n_lons=len(lons)
                print("n_lons={}".format(n_lons))
                for i in range(n_lons):
                    self.assertEqual(lons[i],tr.i2lon(i))

                # inverse operation
                # check that we can predict the index from a given
                # latitude
                for i in range(n_lats):
                    self.assertEqual(tr.lat2i(lats[i]),i)
                # or longitude
                for i in range(n_lons):
                    self.assertEqual(tr.lon2i(lons[i]),i)

                #check the interpretation of the pixel boundaries
                for i in range(n_lats-1):
                    #minimum belongs to pixel
                    self.assertEqual(tr.lat2i(tr.i2lat_min_max(i)[0]),i)
                    #maximum belongs to next pixel (if there is one)
                    self.assertEqual(tr.lat2i(tr.i2lat_min_max(i)[1]),i+1)

               
                # the last lat pixel contains also the pole
                last=n_lats-1
                lat_min, lat_max=tr.i2lat_min_max(last)
                print(lat_max)
                self.assertEqual(tr.lat2i(lat_max),last)
                

    def test_make_model_coord_transforms(self):
        model_folders=['kv_visit2','kv_ft_dlem','Aneesh_SDGVM','cj_isam','jon_yib','yz_jules']#,'cable-pop']
        #for mf in set(self.model_folders):
        for mf in set(model_folders):
            with self.subTest(mf=mf):
                test_args = gh.test_args(mf)
                msh = gh.msh(mf)
                lats=test_args.lats.data
                lons=test_args.lons.data
                # check that we predict the
                # right latitude values for
                # a given array index 
                ctr = msh.make_model_coord_transforms()
                #print(lats)
                print(ctr.lat2LAT(lats))
                print(ctr.lon2LON(lons))

    def test_mask(self):
        #model_folders=['cj_isam','kv_visit2']#,'kv_ft_dlem','Aneesh_SDGVM','cj_isam','jon_yib','yz_jules']#,'cable-pop']
        #for mf in set(model_folders):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                #test_args = gh.test_args(mf)
                msh = gh.msh(mf)
                
                c_mask=msh.spatial_mask(Path(gh.confDict(mf)['dataPath']))
                #c_mask.write_netCDF4(mf+".nc")
                f=plt.figure()
                ax=f.add_subplot(1,1,1)
                c_mask.plot_dots(ax)
                f.savefig(str(mf)+".pdf")



                    
                    

