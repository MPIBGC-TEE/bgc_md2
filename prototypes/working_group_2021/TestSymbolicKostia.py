# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
from unittest.case import TestCase, skip
from pathlib import Path
import json 

import matplotlib.pyplot as plt
import numpy as np
from importlib import import_module
from general_helpers import (
    plot_solutions,
    autostep_mcmc,
    make_param_filter_func,
    make_feng_cost_func
)

class TestSymbolic(TestCase):
    #define a class_variable that will be overloaded by subclasses
    model_folders=['kv_visit2', 'jon_yib','Aneesh_SDGVM','cable-pop','cj_isam','yz_jules','kv_ft_dlem']

    @property
    def model_folders(self):
        return self.__class__.model_folders

    def test_symobolic_description(self):
        for mf in self.model_folders: 
            with self.subTest(mf=mf):
                mvs= import_module('{}.source'.format(mf)).mvs
                print(mvs.get_SmoothReservoirModel())
                print(mvs.get_BibInfo())
    
    #@skip
    def test_download_data(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                #sys.path.insert(0,mf)
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                msh= import_module('{}.model_specific_helpers_2'.format(mf))
                msh.download_my_TRENDY_output(conf_dict)
    
    def test_get_example_site_vars(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                msh= import_module('{}.model_specific_helpers_2'.format(mf))
                svs, dvs = msh.get_example_site_vars(Path(conf_dict['dataPath']))
                #print(svs)
    

    def test_get_global_mean_vars(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                msh= import_module('{}.model_specific_helpers_2'.format(mf))
                svs, dvs = msh.get_global_mean_vars(Path(conf_dict['dataPath']))
                #print(svs)
    

    def test_make_func_dict(self):
        #model_folders=['kv_visit2', 'Aneesh_SDGVM', 'cable-pop']
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f)                 
                #sys.path.insert(0,mf)
                msh= import_module('{}.model_specific_helpers_2'.format(mf))
                mvs = import_module('{}.source'.format(mf)).mvs
                th= import_module('{}.test_helpers'.format(mf))
                test_args=th.make_test_args(conf_dict,msh,mvs)
                svs, dvs = msh.get_example_site_vars(Path(conf_dict['dataPath']))
                msh.make_func_dict(mvs,dvs,test_args.epa_0)
    #@skip
    def test_make_iterator_sym(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                
                #sys.path.insert(0,mf)
                msh= import_module('{}.model_specific_helpers_2'.format(mf))
                th= import_module('{}.test_helpers'.format(mf))
                mvs = import_module('{}.source'.format(mf)).mvs
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                test_args=th.make_test_args(conf_dict,msh,mvs)
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
                #sys.path.insert(0,mf)
                mvs = import_module('{}.source'.format(mf)).mvs
                msh= import_module('{}.model_specific_helpers_2'.format(mf))
                th= import_module('{}.test_helpers'.format(mf))
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                test_args=th.make_test_args(conf_dict,msh,mvs)

                cpa = test_args.cpa
                dvs = test_args.dvs
                svs = test_args.svs
                epa_0 = test_args.epa_0
                param2res_sym = msh.make_param2res_sym( mvs, cpa, dvs)
                xs= param2res_sym(epa_0)


    #@skip
    def test_autostep_mcmc(self):
        for mf in self.model_folders:
            with self.subTest(mf=mf):
                #sys.path.insert(0,mf)
                mvs = import_module('{}.source'.format(mf)).mvs
                msh= import_module('{}.model_specific_helpers_2'.format(mf))
                th= import_module('{}.test_helpers'.format(mf))
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                test_args=th.make_test_args(conf_dict,msh,mvs)
                cpa = test_args.cpa
                dvs = test_args.dvs
                svs = test_args.svs
                epa_min = test_args.epa_min
                epa_max = test_args.epa_max
                epa_0 = test_args.epa_0

                isQualified = make_param_filter_func(epa_max, epa_min)
                param2res = msh.make_param2res_sym( mvs, cpa, dvs)
                #obs=np.column_stack([ np.array(v) for v in svs])
                #obs=obs[0:cpa.number_of_months,:] #cut
                # Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
                C_autostep, J_autostep = autostep_mcmc(
                    initial_parameters=epa_0,
                    filter_func=isQualified,
                    param2res=param2res,
                    costfunction=msh.make_feng_cost_func_2(svs),
                    nsimu=20, # for testing and tuning mcmc
                    #nsimu=20000,
                    c_max=np.array(epa_max),
                    c_min=np.array(epa_min),
                    acceptance_rate=15,   # default value | target acceptance rate in %
                    chunk_size=10,  # default value | number of iterations to calculate current acceptance ratio and update step size
                    D_init=1,   # default value | increase value to reduce initial step size
                    K=1 # default value | increase value to reduce acceptance of higher cost functions
                )

