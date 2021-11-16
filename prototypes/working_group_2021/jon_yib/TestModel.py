# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
sys.path.insert(0,'..')

import unittest
import json
from sympy import var, Symbol
from time import time
import numpy as np
from testinfrastructure.InDirTest import InDirTest
#from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
from pathlib import Path
import matplotlib.pyplot as plt
from collections import OrderedDict

from model_specific_helpers import (
    EstimatedParameters, 
    UnEstimatedParameters, 
    Parameters, 
    StateVariables,
    ModelParameters,
    Observables,
    monthly_to_yearly,
    pseudo_daily_to_yearly,
    get_example_site_vars,
    get_variables_from_files,
    make_param2res,
)
#from bgc_md2.models.cable_yuanyuan.helpers import (
#    construct_matrix_func_sym
#)
from general_helpers import (
    respiration_from_compartmental_matrix,
    month_2_day_index,
    plot_solutions
)

# read the config file and create the config dictionary
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

class TestModel(InDirTest):
    def test_get_variables_from_files(self):
        var_npp, var_rh, var_ra, var_cveg, var_csoil =  get_variables_from_files(Path(conf_dict['dataPath']))

    def test_monthly_to_yearly(self):
        ex_monthly=np.ones(shape=(3840,181,360))
        ref=np.ones(shape=(320,181,360))
        res=monthly_to_yearly(ex_monthly)
        self.assertTrue((res==ref).all())

        ex_monthly=np.ones(shape=(3840,))
        ref=np.ones(shape=(320,))
        res=monthly_to_yearly(ex_monthly)
        self.assertTrue((res==ref).all())
    
    def test_pseudo_daily_to_yearly(self):
        ex_daily=np.ones(shape=(720,11))
        ref=np.ones(shape=(2,11))
        res=pseudo_daily_to_yearly(ex_daily)
        #from IPython import embed; embed()
        self.assertTrue((res==ref).all())


    def test_forward_simulation(self):
        # compare stored monthly timesteps (although the computation happens in daily steps)
        t0=time()
        npp, rh, ra, csoil, cveg = get_example_site_vars(Path(conf_dict['dataPath']))
        print("data_loaded after",time()-t0)
        print(list(map(lambda a:a.shape,(npp,rh,ra,csoil,cveg))))
        nyears = 320
        #nyears = 2
        obs_tup=Observables(
            c_veg=cveg,
            c_soil=csoil,
            a_respiration=monthly_to_yearly(ra),
            h_respiration=monthly_to_yearly(rh)
        )

        obs = np.stack(obs_tup, axis=1)[0:nyears,:]
        
        epa0 = EstimatedParameters(
            beta_leaf=0.15,             # 0         
            beta_root=0.2,              # 1      
            k_leaf=1/365,               # 2      
            k_root=1/(365*5),           # 3         
            k_wood=1/(365*40),          # 4
            k_cwd=1/(365*5),            # 5      
            k_samet=0.5/(365*0.1),      # 6      
            k_sastr=0.5/(365*0.1),      # 7      
            k_samic=0.3/(365*0.137),    # 8      
            k_slmet=0.3/(365),          # 9      
            k_slstr=0.3/(365),          # 10      
            k_slmic=0.3/(365),          # 11      
            k_slow=0.3/(365*5),         # 12      
            k_arm=0.3/(222*365),        # 13      
            f_samet_leaf=0.3,           # 14      
            f_slmet_root=0.3,           # 15      
            f_samic_cwd=0.3,            # 16     
            C_leaf_0=cveg[0]/5,         # 17      
            C_root_0=cveg[0]/5,         # 18      
            C_cwd_0=cveg[0]/50,         # 19      
            C_samet_0=cveg[0]/300,      # 20      
            C_sastr_0=cveg[0]/300,      # 21      
            C_samic_0=cveg[0]/500,      # 22      
            C_slmet_0=csoil[0]/10,      # 23      
            C_slstr_0=csoil[0]/10,      # 24      
            C_slmic_0=csoil[0]/10,      # 25      
            C_slow_0=csoil[0]/10        # 26 
        )

        cpa = UnEstimatedParameters(
            C_soil_0=csoil[0],
            C_veg_0=cveg[0],
            rh_0 = rh[0],
            ra_0 = ra[0],
            npp=npp,
            clay=0.2028,
            silt=0.2808,
            nyears=320
        )
        param2res = make_param2res(cpa)
        t1=time()
        res = param2res(epa0)
        print(time()-t1) 

        #from IPython import embed; embed()
        # the times have to be computed in days
        fig = plt.figure()
        plot_solutions(
            fig,
            times=np.array(range(nyears)), 
            var_names = Observables._fields,
            tup=(
                res ,
                obs 
            ),
            names=["solution with initial params","observations"]
        )
        fig.savefig('solutions.pdf')
