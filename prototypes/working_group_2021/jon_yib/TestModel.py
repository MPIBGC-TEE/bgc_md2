# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
sys.path.insert(0,'..')

import unittest
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



class TestModel(InDirTest):
    def setUp(self):
        self.dataPath = Path('/home/data/yuanyuan/')
        ####################################################
        # The estimated parameters include start values
        # for some of the pools
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
	    CWD_0=0.1,
	    C_mic_0=1,
	    C_passom_0=5,
        )


        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(self.dataPath)

        self.x_init = StateVariables(
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
       
        cpa= UnEstimatedParameters(
            C_leaf_0 = cleaf[0],
            C_root_0 = croot[0],
            C_wood_0 = cwood[0],
            clitter_0 = clitter[0],
            csoil_0 = csoil[0],
            rh_0 = rh[0],
            clay=0.2028,
            silt=0.2808,
            lig_wood=0.4,
            f_wood2CWD=1,
            f_metlit2mic=0.45,
            npp=npp,
            number_of_months=10
        )
        pa = Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa0,cpa)

        mpa = ModelParameters(**{k:v for k,v in pa._asdict().items() if k in ModelParameters._fields})
        
        self.npp= npp
        self.epa0 = epa0 
        self.cpa = cpa
        self.pa = pa 
        self.mpa = mpa 

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

    
    def test_daily_forward_simulation(self):
        # This test does not have a fixture (a result to test again)
        # but it proves that we can compute arbitrary time steps
        #
        cpa=self.cpa  # make UnEstimatedParametersparameters available from setup
        epa=self.epa0 # EstimatedParameters
        mpa=self.mpa  # ModelParameters
        # Initial carbon pool size
        V_init = construct_V0(cpa,epa)
        run_forward_simulation(
                V_init=V_init,
                day_indices=month_2_day_index(range(self.pa.number_of_months)),
                mpa=mpa
        )
        #def make_smr(**kwargs):
        #    #    x_init,
        #    #    npp,
        #    #    clay,
        #    #    silt,
        #    #    lig_wood,
        #    #    beta_leaf,
        #    #    beta_root,
        #    #    lig_leaf,
        #    #    f_leaf2metlit,
        #    #    f_root2metlit,
        #    #    f_wood2CWD,
        #    #    f_metlit2mic,
        #    #    k_leaf,
        #    #    k_root,
        #    #    k_wood,
        #    #    k_metlit,
        #    #    k_mic,
        #    #    k_slowsom,
        #    #    k_passsom
        #    #):
        #    from bgc_md2.models.cable_yuanyuan.source import mvs 
        #    
        #    # we now use all the different parameters of this function
        #    # to build the constituents of a SmoothModelRun 
        #    # 1.) parDict
        #    # 2.) funcDict
        #    # 3.) initial values 
        #    
        #    # 1.)
        #    symbol_names = mvs.get_BibInfo().sym_dict.keys()   
        #    for name in symbol_names:
        #        var(name)

        #    
        #    par_dict={
        #        # here we assume that the symbols have the same name
        #        # as the parameters of this function
        #        Symbol(k): kwargs[k]  
        #        for k in [
        #            'clay',
        #            'silt',
        #            'lig_wood',
        #            'beta_leaf',
        #            'beta_root',
        #            'lig_leaf',
        #            'f_leaf2metlit',
        #            'f_root2metlit',
        #            'f_wood2CWD',
        #            'f_metlit2mic',
        #            'k_leaf',
        #            'k_root',
        #            'k_wood',
        #            'k_metlit',
        #            'k_mic',
        #            'k_slowsom',
        #            'k_passsom'
        #        ]
        #    }
        #    
        #    # 2.) 
        #    def npp_func(t):
        #        # this is a piecewise (monthwise) constant function that
        #        # expects its argument t in unit days which in turn requires
        #        # that all other time related parameters are compatible with
        #        # this assumption
        #        return npp[day_2_month_index(int(t))]

        #    funcDict={NPP :npp_func}

        #    # 3.)
        #    X0=kwargs['x_init']
        #    start_values=np.array(
        #        [
        #           X0._asdict() for k in mvs.get_StateVariableTuple()
        #        ]
        #    )
        #    # 4.) 
        #    times = day_indices

        #    # 5.)
        #    srm = mvs.get_SmoothReservoirModel()
        #    smr = SmoothModelRun(
        #            srm,
        #            parameter_dict=parDict,
        #            start_values=start_values,
        #            times=times,
        #            func_set={}
        #    )
        #    return smr

        #smr = make_smr(**args)
        #sol = smr.solve()
        #RESP_vec = smr.acc_gross_external_output_vector()
        #RESP = np.sum(RESP_vec,axis=1)
        ## in order to attach it to the solution we 
        ## have to have equal length. (in timesteps)
        ## since the accumulated respiration vector
        ## does not return the zeros for the t=0
        ## but the solution contains the start values
        ## we add the zero at the first time step 
        #RESP_w0 = np.concatenate([np.array([0]),RESP]).reshape(sol.shape[0],1)
        ## have to add zeros at the start (
        #res_sym = np.concatenate([sol,RESP_w0], axis=1)
        #res_sym = h.matrix_simul_from_symbolic(
        #    pa=self.mpa,
        #    X0=self.x_init,
        #    npp_in=self.npp,
        #    times=day_indices
        #) 
        ##from IPython import embed;embed()

        #plot_solutions(
        #    day_indices, 
        #    var_names = list(self.x_init._asdict().keys())
        #                + ['respiration'],
        #    tup=(
        #        res,
        #        #res_sym
        #    )
        #)
        #self.assertTrue(np.allclose(res, res_sym, rtol=1e-2))

#    def test_monthly_forward_simulation(self):
#        # compare stored monthly timesteps (although the computation happens in daily steps)
#        ns = [n for n in range(20)] 
#        n_pools = len(self.x_init)
#
#        res = h.monthly_matrix_simu(
#                pa=self.mpa,
#                Xnt=self.x_init,
#                npp_in=self.npp,
#                ns=ns
#        )
#        # the times have to be computed in days
#        times = month_2_day_index(ns)
#        res_sym = h.matrix_simul_from_symbolic(
#            pa=self.mpa,
#            X0=self.x_init,
#            npp_in=self.npp,
#            times=times
#        ) 
#        plot_solutions(
#            times, 
#            var_names = list(self.x_init._asdict().keys())
#                        + ['respiration'],
#            tup=(
#                res,
#                res_sym
#            )
#        )
#        self.assertTrue(np.allclose(res, res_sym,rtol=1e-2))

    @unittest.skip
    def test_mcmc(self):

        data = h.get_example_site_vars(self.dataPath)
        df, df_j = h.mcmc(data, start_pa=self.epa0, nyears=140)
        df.to_csv(str(self.dataPath.joinpath('cable_demo_da_aa.csv')),sep=',')
        df_j.to_csv(str(self.dataPath.joinpath('cable_demo_da_j_aa.csv')),sep=',')


    def test_month_2_day_index(self):
        self.assertEqual(
                month_2_day_index([0]),
                [0]
        ) 
        self.assertEqual(
                month_2_day_index([1]),
                [31]
        ) 
        self.assertEqual(
                month_2_day_index([2]),
                [59]
        ) 
        self.assertEqual(
                month_2_day_index([3]),
                [90]
        ) 
        self.assertEqual(
                month_2_day_index([1,3]),
                [31,90]
        ) 

    def test_day_2_month_index(self):
        # note that days are counted from zero so day 30 is January 31.
        self.assertEqual(day_2_month_index( 0), 0) 
        self.assertEqual(day_2_month_index(30), 0) 
        self.assertEqual(day_2_month_index(31), 1) 
        self.assertEqual(day_2_month_index(60), 2) 

    #@skip
    #def test_results(self):
    #    # This test checks that we can produce the data assimilation for all models in the working group.
    #    param_sets = [ expected_parameters(mdc) for mdc in cmip_5 model_data_combis ]  
    #    # with
    #    def expected_parameters(model_data_combi):
    #        pass
    
    def test_param2res(self):
        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(self.dataPath)
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

