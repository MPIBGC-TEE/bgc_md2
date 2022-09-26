# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
sys.path.insert(0,'..')

import unittest
import numpy as np
import matplotlib.pyplot as plt

from model_specific_helpers import (
    Observables,
    month_2_day_index,
    make_param2res,
    make_param2res_2,
)
from general_helpers import (
    month_2_day_index,
    plot_solutions
)

from TestCommon import TestCommon


class TestModel(TestCommon):

    def test_param2res(self):
        const_params = self.cpa
        param2res = make_param2res(const_params)
        param2res_2 = make_param2res_2(const_params)
        res = param2res(self.epa0)
        res_2 = param2res_2(self.epa0)
        
        day_indices=month_2_day_index(range(self.pa.number_of_months)),

        fig = plt.figure()
        plot_solutions(
                fig,
                times=day_indices,
                var_names=Observables._fields,
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

