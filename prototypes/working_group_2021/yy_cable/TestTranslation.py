# Intermediate Tests for the translation from the original code# that asserted that the results did not change


# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
sys.path.insert(0,'..')

from sympy import var, Symbol
import numpy as np
import matplotlib.pyplot as plt
from unittest import skip
from copy import copy
from model_specific_helpers import (
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
    make_param2res_sym,
    make_daily_iterator,
    make_daily_iterator_sym,
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
    def test_param2res_vs_sym(self):
        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(self.dataPath)
        const_params = self.cpa
        
        param2res = make_param2res(const_params)
        param2res_sym = make_param2res_sym(const_params)
        res = param2res(self.epa0)
        res_sym = param2res_sym(self.epa0)
        
        day_indices=month_2_day_index(range(self.pa.number_of_months)),

        fig = plt.figure()
        plot_solutions(
                fig,
                times=day_indices,
                var_names=Observables._fields,
                tup=(res, res_sym)
        )
        fig.savefig('solutions.pdf')
        
        self.assertTrue(
                np.allclose(
                    res,
                    res_sym,
                    rtol=1e-2
                ),
        )

    #######################################################################
    ## All tests from here on are only required to track down the cause of failure 
    ## if test_param2res_versions fails

    def test_param2res_versions(self):
        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(self.dataPath)
        const_params = self.cpa
        
        param2res = make_param2res(const_params)
        param2res_2 = make_param2res_2(const_params)
        param2res_sym = make_param2res_sym(const_params)
        res = param2res(self.epa0)
        res_2 = param2res_2(self.epa0)
        res_sym = param2res_sym(self.epa0)
        
        day_indices=month_2_day_index(range(self.pa.number_of_months)),

        fig = plt.figure()
        plot_solutions(
                fig,
                times=day_indices,
                var_names=Observables._fields,
                tup=(res, res_2, res_sym)
        )
        fig.savefig('solutions.pdf')
        self.assertTrue(
                np.allclose(
                    res,
                    res_2,
                    rtol=1e-2
                ),
        )
        self.assertTrue(
                np.allclose(
                    res,
                    res_sym,
                    rtol=1e-2
                ),
        )

    #@skip   
    def test_daily_iterator_versions(self):
        mpa = self.mpa
        # Construct npp(day)
        # in general this function can depend on the day i and the state_vector X
        # e.g. typically the size fo X.leaf...
        # In this case it only depends on the day i 
        def npp_func(day):
            #xs = tup[1:]
            return mpa.npp[day_2_month_index(day)] 
        it= make_daily_iterator(
            V_init=construct_V0(
                self.cpa,
                self.epa0
            ),
            mpa=mpa,
        )

        it_sym = make_daily_iterator_sym(
            V_init=construct_V0(
                self.cpa,
                self.epa0
            ),
            mpa=mpa,
            func_dict={
                'NPP':npp_func
            },
        )
        n=5
        res= np.zeros((n,10))
        res_sym = copy(res)
        for i in range(n):
            res[i,:]=it.__next__().reshape(10,)
            res_sym[i,:]=it_sym.__next__().reshape(10,)
        self.assertTrue(
            np.allclose(
                res,
                res_sym
            )
        )






    def test_npp_func(self):
        pa = self.pa
        months = range(pa.number_of_months)
        npp = pa.npp
        def npp_func(day):
            return npp[day_2_month_index(day)] 
        
        day_indices=month_2_day_index(months)
        res = npp[months]

        res_2 = np.array([npp_func(d) for d in day_indices])
        
        fig = plt.figure()
        plot_solutions(
                fig,
                times=day_indices,
                var_names=['nnp'],
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

    def test_respiration(self):
        # The respiration is a function of the Matrix B = A*K
        # here we test that this leads to the same result as yuanyuan's original version
        # which computes it directly from the parameters
        def respiration_from_params(pa,X):
            lig_leaf = pa.lig_leaf
            f41 = pa.f_leaf2metlit
            f42 = pa.f_root2metlit
            f51 = 1 - f41
            f52 = 1 - f42
            f63 = 1
            f74 = 0.45
            f75 = 0.45 * (1 - lig_leaf)
            f85 = 0.7 * lig_leaf
            f86 = 0.4 * (1 - pa.lig_wood)
            f96 = 0.7 * pa.lig_wood
            f87=(0.85 - 0.68 * (pa.clay+pa.silt)) * (0.997 - 0.032 * pa.clay)
            f97=(0.85 - 0.68 * (pa.clay+pa.silt)) * (0.003 + 0.032 * pa.clay)
            f98=0.45 * (0.003 + 0.009 * pa.clay)
            temp = [
                pa.k_leaf,
                pa.k_root,
                pa.k_wood,
                pa.k_metlit,
                pa.k_metlit/(5.75*np.exp(-3*pa.lig_leaf)),
                pa.k_metlit/20.6,
                pa.k_mic,
                pa.k_slowsom,
                pa.k_passsom
            ]
            K = np.zeros(81).reshape([9, 9])
            for i in range(0, 9):
                K[i][i] = temp[i]
        
            co2_rate = np.array(
                    [0,0,0, (1-f74)*K[3,3],(1-f75-f85)*K[4,4],(1-f86-f96)*K[5,5], (1- f87 - f97)*K[6,6], (1-f98)*K[7,7], K[8,8] ])#.reshape(9,1)
        
            return np.sum(co2_rate*X.reshape(1,9))

        Xnt=self.x_init
        X=np.array(Xnt).reshape(9,1)
        orh = respiration_from_params(self.mpa,X)
        
        B_func  = make_compartmental_matrix_func(self.mpa)
        B = B_func(0, self.x_init)
        nrh = respiration_from_compartmental_matrix(B,X)
        self.assertTrue(np.allclose(orh-nrh,np.zeros((9,1))))

    def test_construct_matrix(self):
        B_func  = make_compartmental_matrix_func(
            self.mpa
        )
        B_func_sym = construct_matrix_func_sym(self.mpa)
        B = B_func(0, self.x_init)
        B_sym = B_func_sym(0, self.x_init)
        b=np.array(np.abs(B-B_sym)>1e-3)
        l=[ (i,j) for i in range(9) for j in range(9) if b[i,j]]
        print(l)
        self.assertTrue(np.allclose(B, B_sym))

    def test_daily_forward_simulation(self):
        cpa=self.cpa  # make UnEstimatedParametersparameters available from setup
        epa=self.epa0 # EstimatedParameters
        mpa=self.mpa  # ModelParameters
        # Initial carbon pool size
        V_init = construct_V0(cpa,epa)
        res1=run_forward_simulation(
                V_init=V_init,
                day_indices=month_2_day_index(range(self.pa.number_of_months)),
                mpa=mpa
        )
        res2=run_forward_simulation_sym(
                V_init=V_init,
                day_indices=month_2_day_index(range(self.pa.number_of_months)),
                mpa=mpa,
        )
        self.assertTrue(np.allclose(res1,res2))





