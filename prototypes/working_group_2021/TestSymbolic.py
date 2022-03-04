# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
from unittest.case import TestCase, skip
from pathlib import Path
import json 

import matplotlib.pyplot as plt
import numpy as np
from general_helpers import (
    month_2_day_index,
    plot_solutions,
    autostep_mcmc,
    make_param_filter_func,
    make_feng_cost_func
)

class TestSymbolic(TestCase):

    def test_autostep_mcmc(self):
        model_folders=['kv_visit2', 'Aneesh_SDGVM']
        for mf in model_folders:
            with self.subTest(mf=mf):
                sys.path.insert(0,mf)
                from source import mvs
                import model_specific_helpers as msh
                import testhelpers as th
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                test_args=th.make_test_args(conf_dict)
                cpa = test_args.cpa
                dvs = test_args.dvs
                svs = test_args.svs
                epa_min = test_args.epa_min
                epa_max = test_args.epa_max
                epa_0 = test_args.epa_0

                isQualified = make_param_filter_func(epa_max, epa_min)
                param2res = msh.make_param2res_sym(cpa,dvs,svs)
                obs=np.column_stack([ np.array(v) for v in svs])
                obs=obs[0:cpa.number_of_months,:] #cut
                # Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
                C_autostep, J_autostep = autostep_mcmc(
                    initial_parameters=epa_0,
                    filter_func=isQualified,
                    param2res=param2res,
                    costfunction=make_feng_cost_func(obs),
                    nsimu=200, # for testing and tuning mcmc
                    #nsimu=20000,
                    c_max=np.array(epa_max),
                    c_min=np.array(epa_min),
                    acceptance_rate=15,   # default value | target acceptance rate in %
                    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
                    D_init=1,   # default value | increase value to reduce initial step size
                    K=2 # default value | increase value to reduce acceptance of higher cost functions
                )


    def test_param2res_sym(self):
        model_folders=['kv_visit2']
        for mf in model_folders:
            with self.subTest(mf=mf):
                sys.path.insert(0,mf)
                from source import mvs
                import model_specific_helpers as msh
                import testhelpers as th
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                test_args=th.make_test_args(conf_dict)

                cpa = test_args.cpa
                dvs = test_args.dvs
                svs = test_args.svs
                param2res_sym = msh.make_param2res_sym(
                    cpa,
                    dvs,
                    svs
                )
                xs= param2res_sym(test_args.epa_0)
                obs=np.column_stack([ np.array(v) for v in svs])
                
                obs=obs[0:cpa.number_of_months,:] #cut 
                obs[:,3:4]=obs[:,3:4]
                n=cpa.number_of_months
                
                
                print("Forward run with initial parameters (blue) vs TRENDY output (orange)")
                fig = plt.figure(figsize=(12, 4), dpi=80)
                plot_solutions(
                        fig,
                        times=range(n),
                        var_names=msh.Observables._fields,
                        tup=(xs,obs)
                        #tup=(obs,)
                )
                fig.savefig('solutions.pdf')

    #@skip
    def test_make_iterator_sym(self):
        model_folders=['kv_visit2']
        for mf in model_folders:
            with self.subTest(mf=mf):
                
                sys.path.insert(0,mf)
                from source import mvs
                import testhelpers as th
                import model_specific_helpers as msh
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                test_args=th.make_test_args(conf_dict)
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
                
    def test_symobolic_description(self):
        model_folders=['kv_visit2']
        for mf in model_folders:
            with self.subTest(mf=mf):
                sys.path.insert(0,mf)
                from source import mvs
                mvs.get_SmoothReservoirModel()
    
    @skip
    def test_donwload_data(self):
        model_folders=['kv_visit2']
        for mf in model_folders:
            with self.subTest(mf=mf):
                sys.path.insert(0,mf)
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                import model_specific_helpers as msh
                msh.download_my_TRENDY_output(conf_dict)

    def test_get_example_site_vars(self):
        model_folders=['kv_visit2']
        for mf in model_folders:
            with self.subTest(mf=mf):
                sys.path.insert(0,mf)
                with Path(mf).joinpath('config.json').open(mode='r') as f:
                    conf_dict=json.load(f) 
                import model_specific_helpers as msh
                msh.get_example_site_vars(Path(conf_dict['dataPath']))
