# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
from unittest.case import TestCase, skip
from pathlib import Path
import json 

# +
import matplotlib.pyplot as plt
import numpy as np
from importlib import import_module
from general_helpers import (
    plot_solutions,
    autostep_mcmc,
    make_param_filter_func,
    make_feng_cost_func
)

#model_folders=['kv_visit2', 'jon_yib','Aneesh_SDGVM','cable-pop','cj_isam','yz_jules','kv_ft_dlem']
model_folders=['Aneesh_SDGVM']
#model_folders=['kv_ft_dlem']
# -

class TestSymbolic(TestCase):

    def test_autostep_mcmc(self):
        #model_folders=['kv_visit2']#, 'Aneesh_SDGVM']
        model_folders=['Aneesh_SDGVM']
        for mf in model_folders:
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
                    #costfunction=make_feng_cost_func(obs),
                    costfunction=msh.make_weighted_cost_func(svs),
                    nsimu=20, # for testing and tuning mcmc
                    #nsimu=20000,
                    c_max=np.array(epa_max),
                    c_min=np.array(epa_min),
                    acceptance_rate=15,   # default value | target acceptance rate in %
                    chunk_size=10,  # default value | number of iterations to calculate current acceptance ratio and update step size
                    D_init=1,   # default value | increase value to reduce initial step size
                    K=2 # default value | increase value to reduce acceptance of higher cost functions
                )
