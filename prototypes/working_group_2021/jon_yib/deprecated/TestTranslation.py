# Intermediate Tests for the translation from the original code# that asserted that the results did not change. These Tests will be obsolete after the Transition to the symbolice code


# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
sys.path.insert(0,'..')

from sympy import var, Symbol
import numpy as np
import pandas as pd
from pathlib import Path
import json 
import matplotlib.pyplot as plt
from unittest import TestCase, skip
from testinfrastructure.InDirTest import InDirTest
from model_specific_helpers import (
    EstimatedParameters, 
    UnEstimatedParameters, 
    Parameters, 
    StateVariables,
    ModelParameters,
    Observables,
    run_forward_simulation,
    run_forward_simulation_sym,
    construct_V0,
    make_compartmental_matrix_func,
    get_example_site_vars,
    month_2_day_index,
    day_2_month_index,
    make_param2res,
    make_param2res_2,
    make_param_filter_func,
    make_weighted_cost_func,
    construct_matrix_func_sym
)

from general_helpers import (
    respiration_from_compartmental_matrix,
    month_2_day_index,
    year_2_day_index,
    plot_solutions
)
from TestCommon import TestCommon

class TestTranslation(TestCommon):
    #@skip
    def test_daily_forward_simulation(self):
        # This test does not have a fixture (a result to test agains)
        # but it proves that we can compute arbitrary time steps
        #
        cpa=self.cpa  # make UnEstimatedParametersparameters available from setup
        epa=self.epa0 # EstimatedParameters
        mpa=self.mpa  # ModelParameters
        # Initial carbon pool size
        V_init = jon_yib.construct_V0(cpa,epa)
        jon_yib.run_forward_simulation(
                V_init=V_init,
                day_indices=year_2_day_index(range(self.pa.nyears)),
                mpa=mpa
        )
    @skip   
    def test_param2res_variants(self):
        param2res = jon_yib.make_param2res(self.cpa)
        res = param2res(self.epa0)

        #from IPython import embed; embed()
        # the times have to be computed in days
        fig = plt.figure()
        plot_solutions(
            fig,
            times=np.array(range(nyears)), 
            var_names = jon_yib.Observables._fields,
            tup=(
                res ,
                self.obs 
            ),
            names=["solution with initial params","observations"]
        )
        fig.savefig('solutions.pdf')
        
        param2res_2 = jon_yib.make_param2res_2(self.cpa)
        res = param2res(self.epa0)
        res_2 = param2res_2(self.epa0)
        
        day_indices=month_2_day_index(range(self.pa.number_of_months)),

        fig = plt.figure()
        plot_solutions(
                fig,
                times=day_indices,
                var_names=jon_yib.Observables._fields,
                tup=(res, res_2)
        )
        fig.savefig('solutions.pdf')
        self.assertTrue(
                np.allclose(
                    res,
                    res_2,
                    rtol=1e-2
                ),
        )
