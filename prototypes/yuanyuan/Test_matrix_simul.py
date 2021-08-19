import unittest
import numpy as np
from testinfrastructure.InDirTest import InDirTest
from pathlib import Path
from mm_helpers import (
    get_example_site_vars,
    one_step_matrix_simu,
    mcmc,
    pa1
)

class TestMatrixSimul(InDirTest):
    def setUp(self):
        self.dataPath = Path('/home/data/yuanyuan/')
    def test_one_step(self):
        npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = get_example_site_vars(self.dataPath)
        pa = pa1 
        
        x_init = np.array([cleaf[0],croot[0],cwood[0],pa[12],pa[13],clitter[0]-pa[12]-pa[13],pa[14],csoil[0]- pa[14] - pa[15], pa[15]]).reshape([9,1])   # Initial carbon pool size
        npp_in = npp[0]
        x, r = one_step_matrix_simu(pa=pa,X=x_init,npp_in=npp_in)

    def test_mcmc(self):

        data = get_example_site_vars(self.dataPath)
        df, df_j = mcmc(data, start_pa=pa1, nyears=140)
        df.to_csv(str(self.dataPath.joinpath('cable_demo_da_aa.csv')),sep=',')
        df_j.to_csv(str(self.dataPath.joinpath('cable_demo_da_j_aa.csv')),sep=',')
