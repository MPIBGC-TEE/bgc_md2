import dask.array
from pathlib import Path

from testinfrastructure.InDirTest import InDirTest
from bgc_md2.models.cable_all.cableHelpers import cacheWrapper, fromCachedArgsMaker


def A_from_raw(cable_out_path, zarr_cache_path):
    # Example for a variable that does not depend on other cached vars but
    # only on the original dataset which is hardcoded like 'iveg' or 'Cplant'
    raw_data = dask.array.from_zarr('raw_data')
    A = dask.array.add(raw_data,raw_data)
    return A

def B_from_raw(cable_out_path, zarr_cache_path):
    # Example for a variable that does not depend on other cached vars but
    # only on the original dataset which is hardcoded like 'iveg' or 'Cplant'
    raw_data = dask.array.from_zarr('raw_data')
    B = dask.array.add(raw_data,raw_data)
    return B

def C_from_cached_args(cable_out_path, zarr_cache_path):
    # example for a dependent function
    A = cacheWrapper(A_from_raw,cable_out_path, zarr_cache_path)
    B = cacheWrapper(B_from_raw,cable_out_path, zarr_cache_path)
    return A+B


def D_from_cached_args(cable_out_path, zarr_cache_path):
    # example for an indirectly dependent function
    C = cacheWrapper(C_from_cached_args,cable_out_path,zarr_cache_path)
    return 5*C


class TestRecursiveCache(InDirTest):
    def setUp(self):
        # simulate a cable data_set that we do not touch
        # this is a placeholder for the fortran generated cable output dir

        cable_out_path = Path(".").joinpath('org')
        # create a file
        cable_out_path.mkdir()
        raw_data = dask.array.ones((3, 3))
        dask.array.to_zarr(arr=raw_data, url='raw_data')

        # create a directory for the zarr files as cache for the computed variables
        zarr_cache_path = Path(".").joinpath('cache')
        zarr_cache_path.mkdir()

        self.cable_out_path = cable_out_path
        self.zarr_cache_path = zarr_cache_path
        self.raw_data = raw_data

    def test_root(self):

        # compute variables directly from the original data
        res_A = cacheWrapper(A_from_raw, self.cable_out_path, self.zarr_cache_path)

        self.assertTrue((res_A == 2 * self.raw_data).all().compute())
        zarr_dirs = [f.name for f in self.zarr_cache_path.iterdir()]
        self.assertTrue('A_from_raw' in zarr_dirs)


        res_B = cacheWrapper(B_from_raw, self.cable_out_path, self.zarr_cache_path)

        zarr_dirs = [f.name for f in self.zarr_cache_path.iterdir()]
        self.assertTrue((res_B == 2 * self.raw_data).all().compute())
        self.assertTrue('B_from_raw' in zarr_dirs)

    def test_directly_dependent(self):

        # compute something that depends on cached variables
        res_C = cacheWrapper(C_from_cached_args, self.cable_out_path, self.zarr_cache_path)
        zarr_dirs = [f.name for f in self.zarr_cache_path.iterdir()]

        # all the variables have been cached
        for name in ['A_from_raw', 'B_from_raw', 'C_from_cached_args']:
            self.assertTrue(name in zarr_dirs)

    def test_indirectly_dependent(self):

        # compute something that only
        res_C = cacheWrapper(D_from_cached_args, self.cable_out_path, self.zarr_cache_path)
        zarr_dirs = [f.name for f in self.zarr_cache_path.iterdir()]
        for name in ['A_from_raw', 'B_from_raw', 'C_from_cached_args', 'D_from_cached_args']:
            self.assertTrue(name in zarr_dirs)

    def test_auto_created_root(self):
        import testFuncs  # the module with the functions doing the computing
        cached_A = fromCachedArgsMaker(
                testFuncs.A, 
                testFuncs
        )
        res_A = cacheWrapper(cached_A, self.cable_out_path, self.zarr_cache_path)

    def test_auto_created_directly_dependent(self):
        import testFuncs  # the module with the functions doing the computing
        C_from_cached = fromCachedArgsMaker(
                testFuncs.C, 
                testFuncs
        )
        res_C = cacheWrapper(C_from_cached, self.cable_out_path, self.zarr_cache_path)

    def test_auto_created_indirectly_dependent(self):
        import testFuncs  # the module with the functions doing the computing
        D_from_cached = fromCachedArgsMaker(
                testFuncs.D, 
                testFuncs
        )
        res_D = cacheWrapper(D_from_cached, self.cable_out_path, self.zarr_cache_path)
