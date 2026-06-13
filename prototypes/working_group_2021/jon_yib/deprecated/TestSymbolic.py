# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
sys.path.insert(0,'..')

from unittest import skip
from sympy import Symbol

import matplotlib.pyplot as plt
from model_specific_helpers import (
    Observables,
    month_2_day_index,
    make_param2res_sym,
    make_daily_iterator_sym,
    monthly_to_yearly,
    construct_V0
)

from general_helpers import (
    day_2_year_index,
    month_2_day_index,
    plot_solutions
)
from TestCommon import TestCommon

class TestSymbolic(TestCommon):
    def test_construct_V0(self):
        V_init = construct_V0(
            self.cpa,
            self.epa0
        )


    @skip   
    def test_daily_iterator_sym(self):
        mpa = self.mpa
        # Construct npp(day)
        # in general this function can depend on the day i and the state_vector X
        # e.g. typically the size fo X.leaf...
        # In this case it only depends on the day i 
        def npp_func(day,X):
            #return mpa.npp[day_2_month_index(day)] 
            return mpa.npp[day_2_year_index(day)] 

        func_dict = {Symbol('npp'):npp_func}
        res = make_daily_iterator_sym(
                V_init=construct_V0(
                    self.cpa,
                    self.epa0
                ),
                mpa=self.mpa,
                func_dict={}
        )
    @skip   
    def test_param2res_sym(self):
        param2res = jon_yib.make_param2res_sym(self.cpa)
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
    
