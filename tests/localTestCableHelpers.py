# import unittest
# import zarr as zr
import numpy as np
import dask.array

# from dask.distributed import Client, LocalCluster
import bgc_md2.models.cable_all.cableHelpers as cH
from bgc_md2.helper import batchSlices
from testinfrastructure.InDirTest import InDirTest


class TestCableHelpers(InDirTest):
    def test_test(self):
        a = dask.array.from_array(np.array([1, 2]))
        b = dask.array.from_array(np.array([1, 0]))
        c = dask.array.from_array(np.array([1, 0]))

        self.assertFalse(dask.array.all(a == b).compute())
        self.assertTrue(dask.array.all(b == c).compute())

    def test_valid_combies(self):
        ifv = -99999
        nland = 23
        npatch = 2
        ntime = 30
        npool = 9

        # fake an iveg array where half the patches are valid
        iveg = dask.array.stack(
            [
                np.concatenate(
                    [
                        np.full((int(npatch / 2),), ifv, dtype=np.int32),
                        np.full((int(npatch / 2),), 1, dtype=np.int32)
                        # np.random(1, 3, size=((int(npatch/2), ), 1, dtype=np.int32)
                    ]
                )
                for i in range(nland)
            ],
            axis=-1,
        ).rechunk((npatch, 1))

        # fake B and u arrays
        B = dask.array.stack(
            [
                # np.zeros((ntime, npool, npool, npatch))
                np.random.uniform(1, 3, size=(ntime, npool, npool, npatch))
                for i in range(nland)
            ],
            axis=-1,
        )

        u = dask.array.stack(
            [
                # np.zeros((ntime, npool, npool, npatch))
                np.random.uniform(1, 3, size=(ntime, npool, npatch))
                for i in range(nland)
            ],
            axis=-1,
        )

        x0 = dask.array.stack(
            [
                # np.zeros((ntime, npool, npool, npatch))
                np.random.uniform(1, 3, size=(npool, npatch))
                for i in range(nland)
            ],
            axis=-1,
        )

        nz = dask.array.nonzero((iveg == ifv))
        n_val = int(npatch * nland / 2)
        # print(nz)
        B_val = cH.valid_combies(nz, B)
        u_val = cH.valid_combies(nz, u)
        x0_val = cH.valid_combies(nz, x0)
        self.assertTrue(B_val.shape == (ntime, npool, npool, n_val))
        self.assertTrue(u_val.shape == (ntime, npool, n_val))
        self.assertTrue(x0_val.shape == (npool, n_val))

        # write them to disk
        cH.batchwise_to_zarr(B_val, "B.zarr")
        cH.batchwise_to_zarr(u_val, "u.zarr")
        cH.batchwise_to_zarr(x0_val, "x0.zarr")

        # and read them again
        self.assertTrue(
            dask.array.all(dask.array.from_zarr("B.zarr") == B_val).compute()
        )
        self.assertTrue(
            dask.array.all(dask.array.from_zarr("u.zarr") == u_val).compute()
        )
        self.assertTrue(
            dask.array.all(dask.array.from_zarr("x0.zarr") == x0_val).compute()
        )

    def test_valid_combies_parallel(self):
        ifv = -99999
        nland = 56
        npatch = 4
        ntime = 30
        npool = 9

        # fake an iveg array where half the patches are valid
        iveg = dask.array.stack(
            [
                np.concatenate(
                    [
                        np.full((int(npatch / 2),), ifv, dtype=np.int32),
                        np.full((int(npatch / 2),), 1, dtype=np.int32)
                        # np.random(1, 3, size=((int(npatch/2), ), 1, dtype=np.int32)
                    ]
                )
                for i in range(nland)
            ],
            axis=-1,
        ).rechunk((npatch, 1))

        # fake B and u arrays
        B = dask.array.stack(
            [
                # np.zeros((ntime, npool, npool, npatch))
                np.random.uniform(1, 3, size=(ntime, npool, npool, npatch))
                for i in range(nland)
            ],
            axis=-1,
        )

        u = dask.array.stack(
            [
                # np.zeros((ntime, npool, npool, npatch))
                np.random.uniform(1, 3, size=(ntime, npool, npatch))
                for i in range(nland)
            ],
            axis=-1,
        )

        x0 = dask.array.stack(
            [
                # np.zeros((ntime, npool, npool, npatch))
                np.random.uniform(1, 3, size=(npool, npatch))
                for i in range(nland)
            ],
            axis=-1,
        )

        nz = dask.array.nonzero((iveg == ifv))
        n_val = int(npatch * nland / 2)
        # print(nz)
        B_val = cH.valid_combies_parallel(nz, B)
        u_val = cH.valid_combies_parallel(nz, u)
        x0_val = cH.valid_combies_parallel(nz, x0)
        self.assertTrue(B_val.shape == (ntime, npool, npool, n_val))
        self.assertTrue(u_val.shape == (ntime, npool, n_val))
        self.assertTrue(x0_val.shape == (npool, n_val))

        # write them to disk
        cH.batchwise_to_zarr(B_val, "B.zarr")
        cH.batchwise_to_zarr(u_val, "u.zarr")
        cH.batchwise_to_zarr(x0_val, "x0.zarr")

        # and read them again
        self.assertTrue(
            dask.array.all(dask.array.from_zarr("B.zarr") == B_val).compute()
        )
        self.assertTrue(
            dask.array.all(dask.array.from_zarr("u.zarr") == u_val).compute()
        )
        self.assertTrue(
            dask.array.all(dask.array.from_zarr("x0.zarr") == x0_val).compute()
        )
