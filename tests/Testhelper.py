import unittest
import numpy as np
from functools import reduce
from bgc_md2.helper import batchSlices




class TestCableUtils(unittest.TestCase):

    def test_batchSlices(self):
        nland = 30

        #  1D array
        testarr = np.arange(0, nland)
        #  make sure that the combined batches combine to the original array
        for nproc in [5, 7]:  # uniform and nonuniform batch sizes
            with self.subTest(nproc=nproc): 
                batches = [testarr[s] for s in batchSlices(nland, nproc)]
                combined = reduce(lambda cum,el:np.append(cum,el),batches)
                self.assertTrue(
                        np.alltrue(
                            combined==testarr
                        )
                )
        # 2D arrays
        a1 = np.arange(0, nland)
        a2 = 10*a1
        testarr = np.stack([a1, a2], axis=0)

        # make sure that the combined batches combine to the original array
        for nproc in [5]:#,7]: # uniform and nonuniform batch sizes
            with self.subTest(nproc=nproc): 
                batches = [testarr[:,s] for s in batchSlices(nland, nproc)]
                combined = reduce(lambda cum,el:np.append(cum,el,axis=1),batches)
                self.assertTrue(
                        np.alltrue(
                            combined==testarr
                        )
                )

        
