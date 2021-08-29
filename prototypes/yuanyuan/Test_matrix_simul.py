from mm_helpers import matrix_simul_from_symbolic
import unittest
import numpy as np
from testinfrastructure.InDirTest import InDirTest
from pathlib import Path
from collections import namedtuple
from mm_helpers import (
    get_example_site_vars,
    one_step_matrix_simu,
    mcmc
)

class TestMatrixSimul(InDirTest):
    def setUp(self):
        self.dataPath = Path('/home/data/yuanyuan/')
        ####################################################
        # The estimated parameters include start values
        # for some of the pools
        EstimatedParameters = namedtuple(
            "EstiamatedParameters",
            [
                "beta_leaf",
                "beta_root",
                "lig_leaf",
	        "f_leaf2metlit",
                "f_root2metlit",
	        "k_leaf",
	        "k_root",
	        "k_wood",
	        "k_metlit",
	        "k_mic",
	        "k_slowsom",
                "k_passsom",
	        "cmet_init",
	        "cstr_init",
	        "cmic_init",
	        "cpassive_init"
            ]
        )

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
	    cmet_init=0.05,
	    cstr_init=0.1,
	    cmic_init=1,
	    cpassive_init=5,
        )

        # the model parameters are parameters are just
        # constants in the model
        ModelParameters = namedtuple(
            "ModelParameters",
            [
                "beta_leaf",
                "beta_root",
                "lig_leaf",
	        "f_leaf2metlit",
                "f_root2metlit",
	        "k_leaf",
	        "k_root",
	        "k_wood",
	        "k_metlit",
	        "k_mic",
	        "k_slowsom",
                "k_passsom",
            ]
        )
        mpa = ModelParameters(
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

        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(self.dataPath)
        InitialValues = namedtuple(
            'IntitialValues',
            [
                'C_leaf',
                'C_root',
                'C_wood',
                'C_metlit',
                'C_stlit',
                'CWD',
                'C_mic',
                'C_slowsom',
                'C_passsom',
            ]
        )

        self.x_init = InitialValues(
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
       
        
        self.npp_init = npp[0]
        self.epa0 = epa0 
        self.mpa = mpa 

    def test_one_step(self):
        x, r = one_step_matrix_simu(
                pa=self.mpa,
                X=self.x_init,
                npp_in=self.npp_init
        )
    
    def test_from_symbolic_forward_simulation(self):
        matrix_simul_from_symbolic(
            pa=self.mpa,
            X0=self.x_init,
            npp_in=self.npp_init
        ) 

    def test_mcmc(self):

        data = get_example_site_vars(self.dataPath)
        df, df_j = mcmc(data, start_pa=self.pa, nyears=140)
        df.to_csv(str(self.dataPath.joinpath('cable_demo_da_aa.csv')),sep=',')
        df_j.to_csv(str(self.dataPath.joinpath('cable_demo_da_j_aa.csv')),sep=',')
