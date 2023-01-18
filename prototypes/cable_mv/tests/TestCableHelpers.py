import unittest

import dask.array
from testinfrastructure.InDirTest import InDirTest

import bgc_md2.models.cable_all.cableHelpers as cH

class TestCableHelpers(InDirTest):

    def test_valid_combis_parallel(self):
        n = 10
        A = dask.array.stack(
            [
                dask.array.from_array(
                    [
                       [3, 0, 0],
                       [0, 4, 0],
                       [5, 6, 0]
                    ]
                )
                for i in range(n)
            ],
            axis=0,
        )
        cond = A[0,...] != 0
        nz_tuple = dask.array.nonzero(cond) # will be of length cond.ndim
        for el in nz_tuple:
            el.compute_chunk_sizes()
        res = cH.valid_combies_parallel(nz_tuple, A)
        self.assertEqual(res.shape, (10,4))

    def test_unmasked_part(self):
        n = 10
        A = dask.array.stack(
            [
                dask.array.from_array(
                    [
                       [3, 0, 0],
                       [0, 4, 0],
                       [5, 6, 0]
                    ]
                )
                for i in range(n)
            ],
            axis=0,
        )
        cond = A[0, ...] != 0
        mask = dask.array.broadcast_to(~cond, A.shape)
        A_ma = dask.array.ma.masked_array(A, mask=mask)
        res = A_ma[~mask]
        res.compute_chunk_sizes()
        final = res.reshape(n,-1)
        print(final)

if __name__ == '__main__':
    ## single test
    s = unittest.TestSuite()
    #s.addTest(TestCableHelpers("test_valid_combis_parallel"))
    s.addTest(TestCableHelpers("test_unmasked_part"))
    runner = unittest.TextTestRunner()
    result = runner.run(s)

