# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
sys.path.insert(0,'..')

import unittest
import json
from sympy import var, Symbol
import numpy as np
from testinfrastructure.InDirTest import InDirTest
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
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
    run_forward_simulation,
    construct_V0,
    make_compartmental_matrix_func,
    get_example_site_vars,
    month_2_day_index,
    day_2_month_index,
    make_param2res,
    make_param2res_2
)
from bgc_md2.models.cable_yuanyuan.helpers import (
    construct_matrix_func_sym
)
from general_helpers import (
    respiration_from_compartmental_matrix,
    month_2_day_index,
    plot_solutions
)

# read the config file and create the config dictionary
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

class TestModel(InDirTest):
    @unittest.skip
    def test_forward_simulation(self):
        # compare stored monthly timesteps (although the computation happens in daily steps)
        npp, rh, clitter, ccwd, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(Path(conf_dict['dataPath']))
        nyears = 10
        tot_len = 12*nyears
        obs_tup=Observables(
            C_leaf=cleaf,
            C_root=croot,
            C_wood=cwood,
            c_litter=clitter,
            c_soil=csoil,
            respiration=rh
        )
        obs = np.stack(obs_tup, axis=1)[0:tot_len,:]
        
        epa0 = EstimatedParameters(
            beta_leaf=0.15,
	    beta_root=0.2,
	    lig_leaf=0.15,
	    f_leaf2metlit=0.28,
	    f_root2metlit=0.6,
	    k_leaf=1/365,
	    k_root=1/(365*5),
	    k_wood=1/(365*40),
	    k_metlit=0.5/(365*0.1),
	    k_mic=0.3/(365*0.137),
	    k_slowsom=0.3/(365*5),
	    k_passsom=0.3/(222.22*365),
	    C_metlit_0=0.05,
	    C_strlit_0=0.1,
	    C_mic_0=1,
	    C_passom_0=5,
        )

        cpa = UnEstimatedParameters(
            C_leaf_0=cleaf[0],
            C_root_0=croot[0],
            C_wood_0=cwood[0],
            C_cwd_0=ccwd[0],
            c_litter_0=clitter[0],
            c_soil_0=csoil[0],
            rh_0 = rh[0],
            npp=npp,
            number_of_months=tot_len,
            clay=0.2028,
            silt=0.2808,
            lig_wood=0.4,
            f_wood2CWD=1, 
            f_metlit2mic=0.45
        )
        param2res = make_param2res(cpa)
        
        res = param2res(epa0)
        # the times have to be computed in days
        fig = plt.figure()
        plot_solutions(
            fig,
            times=np.array(range(tot_len)), 
            var_names = Observables._fields,
            tup=(
                res ,
                obs 
            ),
            names=["solution with initial params","observations"]
        )
        fig.savefig('solutions.pdf')
