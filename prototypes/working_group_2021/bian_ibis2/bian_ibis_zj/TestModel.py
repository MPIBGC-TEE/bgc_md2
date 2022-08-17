# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
import sys
sys.path.insert(0,'..')

import unittest
from sympy import var, Symbol
import numpy as np
import pandas as pd
from testinfrastructure.InDirTest import InDirTest
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
from pathlib import Path
import json 
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
    make_param2res_2,
    make_param_filter_func,
    make_weighted_cost_func
)
from bgc_md2.models.cable_yuanyuan.helpers import (
    construct_matrix_func_sym
)
from general_helpers import (
    respiration_from_compartmental_matrix,
    month_2_day_index,
    make_uniform_proposer,
    make_multivariate_normal_proposer,
    mcmc,
    make_feng_cost_func,
    plot_solutions
)


with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 
dataPath = Path(conf_dict['dataPath'])

class TestModel(InDirTest):
    def setUp(self):
        self.dataPath = dataPath
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
        # This test does not have a fixture (a result to test agains)
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

    def test_mcmc(self):
        
        dataPath = self.dataPath
        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(dataPath)
        
        nyears = 10
        tot_len = 12*nyears
        obs = np.stack([cleaf, croot, cwood, clitter, csoil, rh], axis=1)[0:tot_len,:]
        c_min=np.array([0.09, 0.09,0.09,0.01,0.01,  1/(2*365), 1/(365*10), 1/(60*365), 0.1/(0.1*365), 0.06/(0.137*365), 0.06/(5*365), 0.06/(222.22*365),clitter[0]/100,clitter[0]/100,csoil[0]/100,csoil[0]/2])
        c_max=np.array([1   ,    1,0.21,   1,   1,1/(0.3*365),1/(0.8*365),      1/365,   1/(365*0.1),  0.6/(365*0.137),  0.6/(365*5),  0.6/(222.22*365),    clitter[0],    clitter[0],  csoil[0]/3,csoil[0]])
        
        isQualified = make_param_filter_func(c_max,c_min)
        uniform_prop = make_uniform_proposer(
            c_min,
            c_max,
            #D=10.0,
            D=20.0,
            filter_func=isQualified
        )
        
        cpa = UnEstimatedParameters(
            C_leaf_0=cleaf[0],
            C_root_0=croot[0],
            C_wood_0=cwood[0],
            clitter_0=clitter[0],
            csoil_0=csoil[0],
            rh_0 = rh[0],
            npp=npp,
            number_of_months=tot_len,
            clay=0.2028,
            silt=0.2808,
            lig_wood=0.4,
            f_wood2CWD=1, 
            f_metlit2mic=0.45
        )
        param2res = make_param2res(cpa) #pa=[beta1,beta2, lig_leaf, f41,f42, kleaf,kroot,kwood,kmet,kmic, kslow,kpass, cmet_init, cstr_init, cmic_init, cpassive_init ]
#        pa=            [0.15,  0.2,0.15,0.28, 0.6,      1/365,  1/(365*5), 1/(365*40), 0.5/(365*0.1),  0.3/(365*0.137),  0.3/(365*5),  0.3/(222.22*365),          0.05,           0.1,           1,         5]
        epa_0 = EstimatedParameters(
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
        # save the parameters and costfunctionvalues for postprocessing 
        demo_aa_path = Path('cable_demo_da_aa.csv')
        demo_aa_j_path = Path('cable_demo_da_j_aa.csv')
        if not demo_aa_path.exists():

            print("did not find demo run results. Will perform  demo run")
            nsimu_demo = 200    
            C_demo, J_demo = mcmc(
                    initial_parameters=epa_0,
                    proposer=uniform_prop,
                    param2res=param2res,
                    #costfunction=make_weighted_cost_func(obs)
                    costfunction=make_feng_cost_func(obs),
                    nsimu=nsimu_demo
            )
            # save the parameters and costfunctionvalues for postprocessing 
            pd.DataFrame(C_demo).to_csv(demo_aa_path,sep=',')
            pd.DataFrame(J_demo).to_csv(demo_aa_j_path,sep=',')
        else:
            print("""Found {p} from a previous demo run. 
            If you also want to recreate the demo output move the file!
            """.format(p = demo_aa_path)) 
            C_demo = pd.read_csv(demo_aa_path).to_numpy()
            J_demo = pd.read_csv(demo_aa_j_path).to_numpy() 
        
        # build a new proposer based on a multivariate_normal distribution using the estimated covariance of the previous run if available
        # parameter values of the previous run
        # first we check how many accepted parameters we got 
        n_accept=C_demo.shape[1]
        # and then use part of them to compute a covariance matrix for the 
        # formal run
        covv = np.cov(C_demo[:, int(n_accept/10):]) 

        
        normal_prop = make_multivariate_normal_proposer(
            covv = covv,
            filter_func=isQualified
        )
        C_formal, J_formal = mcmc(
                initial_parameters=epa_0,
                proposer=normal_prop,
                param2res=param2res,
                #costfunction=make_weighted_cost_func(obs)
                costfunction=make_feng_cost_func(obs),
                nsimu=200
        )
        # save the parameters and costfunctionvalues for postprocessing 
        formal_aa_path = Path('cable_formal_da_aa.csv')
        formal_aa_j_path = Path('cable_formal_da_j_aa.csv')
        pd.DataFrame(C_formal).to_csv(formal_aa_path,sep=',')
        pd.DataFrame(J_formal).to_csv(formal_aa_j_path,sep=',')




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

