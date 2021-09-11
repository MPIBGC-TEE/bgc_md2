import unittest
import numpy as np
from testinfrastructure.InDirTest import InDirTest
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
from pathlib import Path
import matplotlib.pyplot as plt
import mm_helpers as h
import yy_cable_specific_helpers as yyh
def plot_solutions(times, var_names,tup):
    n_vars = tup[0].shape[1]
    fig, axs = plt.subplots(n_vars,1,figsize=(10,n_vars*10))
    colors =('red','blue','green','organge')
    #from IPython import embed;embed()
    for j in range(n_vars):
        for i,sol in enumerate(tup):
            axs[j].plot(times, sol[:, j], color = colors[i])
            axs[j].set_title(var_names[i])
    fig.savefig('solutions.pdf')

class TestMatrixSimul(InDirTest):
    def setUp(self):
        self.dataPath = Path('/home/data/yuanyuan/')
        ####################################################
        # The estimated parameters include start values
        # for some of the pools
        epa0 = h.EstimatedParameters(
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
	    cmet_init=0.05,
	    cstr_init=0.1,
	    cmic_init=1,
	    cpassive_init=5,
        )

        mpa = h.ModelParameters(
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
	    k_passsom=0.3/(222.22*365)
        )

        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = yyh.get_example_site_vars(self.dataPath)

        self.x_init = h.InitialValues(
            C_leaf=cleaf[0],
            C_root=croot[0],
            C_wood=cwood[0],
            C_metlit=epa0[12],
            C_stlit=epa0[13],
            CWD=clitter[0]-epa0[12]-epa0[13],
            C_mic=epa0[14],
            C_slowsom=csoil[0]- epa0[14] - epa0[15],
            C_passsom=epa0[15]
        )
       
        
        self.npp= npp
        self.epa0 = epa0 
        self.mpa = mpa 

    @unittest.skip
    def test_one_step(self):
        x, r = h.one_day_matrix_simu(
                pa=self.mpa,
                Xnt=self.x_init,
                npp_in=self.npp
        )
    def test_construct_matrix(self):
        B_func = h.construct_matrix_func(self.mpa)
        B_func_sym = h.construct_matrix_func_sym(self.mpa)
        B = B_func(0, self.x_init),
        B_sym = B_func_sym(0, self.x_init)
        self.assertTrue(np.allclose(B, B_sym))


    def test_construct_allocation_vector(self):
        #only checks that the function works
        b_func = h.construct_allocation_vector_func(self.mpa)


    
    def test_respiration(self):
        # The respiration is a function of the Matrix B = A*K
        # here we test that this leads to the same result as yuanyuan's original version
        # which computes it directly from the parameters
        Xnt=self.x_init,
        X=np.array(Xnt).reshape(9,1)
        orh = h.respiration_from_params(self.mpa,X)
        
        B_func = h.construct_matrix_func(self.mpa)
        B = B_func(0, self.x_init)
        nrh =h.respiration_from_compartmental_matrix(B,X)
        self.assertTrue(np.allclose(orh-nrh,np.zeros((9,1))))

    
    def test_daily_forward_simulation(self):
        # compute dayly timesteps with constant npp input
        # this is a preliminary step since npp changes monthly
        n = 20
        n_pools = len(self.x_init)

        # Construct b vector for allocation
        b_func  = h.construct_allocation_vector_func(self.mpa)
        # Now construct B matrix B=A*K
        B_func  = h.construct_matrix_func(self.mpa)

        # construct a function that returns the 
        # npp value for day it
        npp_func= lambda d: self.npp[h.day_2_month_index(d)]
        
        # V_i+1 = f(it,V_i) 
        def f(it,V):
            X = V[0:9]
            co2 = V[9]
            b = b_func(it,X)
            B = B_func(it,X)
            X_new = X + b * npp_func(it) + B@X
            co2_new = h.respiration_from_compartmental_matrix(B,X) #accumulate respiration
            V_new = np.concatenate((X_new,np.array([co2_new]).reshape(1,1)), axis=0)
            return V_new


        X_0 = np.array(self.x_init).reshape(9,1) #transform the tuple to a matrix
   
        co2_0 = np.array([0]).reshape(1,1)
        tsi = TimeStepIterator2(
                initial_values=np.concatenate((X_0, co2_0), axis=0),
                f=f,
                max_it=n
        )
        res = np.concatenate([ np.transpose(val) for val in tsi], axis=0)  
        
        times = np.arange(0,n,step=1)
        res_sym = h.matrix_simul_from_symbolic(
            pa=self.mpa,
            X0=self.x_init,
            npp_in=self.npp,
            times=times
        ) 

        plot_solutions(
            times, 
            var_names = list(self.x_init._asdict().keys())
                        + ['respiration'],
            tup=(
                res,
                res_sym
            )
        )
        self.assertTrue(np.allclose(res, res_sym, rtol=1e-2))

    def test_monthly_forward_simulation(self):
        # compare stored monthly timesteps (although the computation happens in daily steps)
        ns = [n for n in range(20)] 
        n_pools = len(self.x_init)

        res = h.monthly_matrix_simu(
                pa=self.mpa,
                Xnt=self.x_init,
                npp_in=self.npp,
                ns=ns
        )
        # the times have to be computed in days
        times = h.month_2_day_index(ns)
        res_sym = h.matrix_simul_from_symbolic(
            pa=self.mpa,
            X0=self.x_init,
            npp_in=self.npp,
            times=times
        ) 
        plot_solutions(
            times, 
            var_names = list(self.x_init._asdict().keys())
                        + ['respiration'],
            tup=(
                res,
                res_sym
            )
        )
        self.assertTrue(np.allclose(res, res_sym,rtol=1e-2))

    @unittest.skip
    def test_mcmc(self):

        data = h.get_example_site_vars(self.dataPath)
        df, df_j = h.mcmc(data, start_pa=self.epa0, nyears=140)
        df.to_csv(str(self.dataPath.joinpath('cable_demo_da_aa.csv')),sep=',')
        df_j.to_csv(str(self.dataPath.joinpath('cable_demo_da_j_aa.csv')),sep=',')


    def test_month_2_day_index(self):
        self.assertEqual(
                h.month_2_day_index([0]),
                [0]
        ) 
        self.assertEqual(
                h.month_2_day_index([1]),
                [31]
        ) 
        self.assertEqual(
                h.month_2_day_index([2]),
                [59]
        ) 
        self.assertEqual(
                h.month_2_day_index([3]),
                [90]
        ) 
        self.assertEqual(
                h.month_2_day_index([1,3]),
                [31,90]
        ) 

    def test_day_2_month_index(self):
        # note that days are counted from zero so day 30 is January 31.
        self.assertEqual(h.day_2_month_index( 0), 0) 
        self.assertEqual(h.day_2_month_index(30), 0) 
        self.assertEqual(h.day_2_month_index(31), 1) 
        self.assertEqual(h.day_2_month_index(60), 2) 

    #@skip
    #def test_results(self):
    #    # This test checks that we can produce the data assimilation for all models in the working group.
    #    param_sets = [ expected_parameters(mdc) for mdc in cmip_5 model_data_combis ]  
    #    # with
    #    def expected_parameters(model_data_combi):
    #        pass
    
    def test_matrix_simu_maker_org(self):
        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = h.get_example_site_vars(self.dataPath)
        const_params= h.UnEstimatedParameters(
            cleaf_0 = cleaf[0],
            croot_0 = croot[0],
            cwood_0 = cwood[0],
            clitter_0 = clitter[0],
            csoil_0 = csoil[0],
            npp = npp,
            number_of_months = 3
        )
        matrix_simu_org =  h.matrix_simu_maker_org(const_params)
        x,rh = matrix_simu_org(self.epa0)
        

