# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
sys.path.insert(0,'..')

import matplotlib.pyplot as plt
from model_specific_helpers import (
    Observables,
    month_2_day_index,
    make_param2res_sym,
)

from general_helpers import (
    month_2_day_index,
    plot_solutions
)
from TestCommon import TestCommon

class TestSymbolic(TestCommon):

    def test_param2res_sym(self):
        const_params = self.cpa
        
        param2res_sym = make_param2res_sym(const_params)
        res_sym = param2res_sym(self.epa0)
        
        day_indices=month_2_day_index(range(self.pa.number_of_months)),

        fig = plt.figure()
        plot_solutions(
                fig,
                times=day_indices,
                var_names=Observables._fields,
                tup=(res_sym,)
        )
        fig.savefig('solutions.pdf')

