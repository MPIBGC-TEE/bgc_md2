import dask.array
from unittest import skip
from pathlib import Path

from bgc_md2.models.cable_all.cableHelpers import (
        cacheWrapper,
        fromCachedArgsMaker,
        cachedFromCachedArgsMaker
)
from TestCache import TestCache, time_dict

########################################################################
# now we create the complete recursively caching functions automatically.
class TestCacheFromCachedArgsMaker(TestCache):

    def test_root(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        import testFuncs  # the module with the functions doing the computing
        cached_A = cachedFromCachedArgsMaker(
                testFuncs.A,
                testFuncs
        )
        res_A = cached_A(cable_out_path, zarr_cache_path)
        self.assertTrue((res_A == 2 * raw_data).all().compute())

        cached_B = cachedFromCachedArgsMaker(
                testFuncs.B,
                testFuncs
        )
        res_B = cached_B(
                cable_out_path,
                zarr_cache_path,
                batch_size=32
        )
        self.assertTrue((res_B == 2 * raw_data).all().compute())

        # check if  zarr dirs have been created
        zarr_dirs = [f.name for f in zarr_cache_path.iterdir()]
        for name in [
            'A',
            'B'
        ]:
            self.assertTrue(name in zarr_dirs)

    def test_root_slice(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        import testFuncs  # the module with the functions doing the computing
        landpoint_slice = slice(0,10)
        raw_data_sliced = raw_data[...,landpoint_slice,:]

        cached_A = cachedFromCachedArgsMaker(
                testFuncs.A,
                testFuncs
        )
        res_A = cached_A(
                cable_out_path,
                zarr_cache_path,
                landpoint_slice=landpoint_slice
        )
        self.assertTrue((res_A == 2 * raw_data_sliced).all().compute())

        cached_B = cachedFromCachedArgsMaker(
                testFuncs.B,
                testFuncs
        )
        res_B = cached_B(
                cable_out_path,
                zarr_cache_path,
                landpoint_slice=landpoint_slice,
                batch_size=32
        )
        self.assertTrue((res_B == 2 * raw_data_sliced).all().compute())

        # check if  zarr dirs have been created
        zarr_dirs = [f.name for f in zarr_cache_path.iterdir()]
        for name in [
            'A',
            'B'
        ]:
            self.assertTrue(name in zarr_dirs)

    def test_root_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        import testFuncs  # the module with the functions doing the computing

        cached_B = cachedFromCachedArgsMaker(
                testFuncs.B,
                testFuncs
        )
        res_B = cached_B(cable_out_path, zarr_cache_path)
        # check the times the cache files have been touched
        # first to make sure that they are not touched 
        td_before = time_dict(zarr_cache_path)

        #now with rm flag
        res_B = cached_B(
            cable_out_path,
            zarr_cache_path,
            rm=True
        )
        td_after_rm = time_dict(zarr_cache_path)
        key = 'B'
        self.assertTrue(td_before[key] < td_after_rm[key])

    def test_root_rec_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        import testFuncs  # the module with the functions doing the computing

        cached_B = cachedFromCachedArgsMaker(
                testFuncs.B,
                testFuncs
        )
        res_B = cached_B(cable_out_path, zarr_cache_path)
        # check the times the cache files have been touched
        td_before = time_dict(zarr_cache_path)

        #now with rec_rm flag
        res_B = cached_B(
            cable_out_path,
            zarr_cache_path,
            rec_rm=True
        )
        td_after_rm = time_dict(zarr_cache_path)
        key = 'B'
        # for directly computable vars rec_rm has the same
        # effect as rm and only affects B 
        self.assertTrue(td_before[key] < td_after_rm[key])


    def test_indirectly_dependent(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        import testFuncs  # the module with the functions doing the computing
        cached_D = cachedFromCachedArgsMaker(
                testFuncs.D,
                testFuncs
        )
        res_D = cached_D(cable_out_path, zarr_cache_path)
        self.assertTrue((res_D == 20 * raw_data).all().compute())

        zarr_dirs = [f.name for f in zarr_cache_path.iterdir()]
        for name in [
            'A',
            'B',
            'C',
            'D'
        ]:
            self.assertTrue(name in zarr_dirs)

    def test_indirectly_dependent_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        import testFuncs  # the module with the functions doing the computing
        cached_D = cachedFromCachedArgsMaker(
                testFuncs.D,
                testFuncs
        )
        res_D = cached_D(cable_out_path, zarr_cache_path)
        self.assertTrue((res_D == 20 * raw_data).all().compute())

        zarr_dirs = [f.name for f in zarr_cache_path.iterdir()]
        for name in [
            'A',
            'B',
            'C',
            'D'
        ]:
            self.assertTrue(name in zarr_dirs)

        # check the times the cache files have been touched
        # first to make sure that they are not touched
        td_before = time_dict(zarr_cache_path)
        
        res_D = cached_D(
            cable_out_path,
            zarr_cache_path,
            rm=True
        )

        td_after_rm = time_dict(zarr_cache_path)
        # check that the newly created result is newer
        key = 'D'
        self.assertTrue(td_before[key] < td_after_rm[key])

        # check that the dependencies are untouched
        rest = frozenset.difference(
                frozenset(td_after_rm.keys()),
                frozenset([key])
        )
        for k in rest:
            self.assertEqual(td_before[k], td_after_rm[k])

    def test_indirectly_dependent_rec_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        import testFuncs  # the module with the functions doing the computing
        cached_D = cachedFromCachedArgsMaker(
                testFuncs.D,
                testFuncs
        )
        res_D = cached_D(cable_out_path, zarr_cache_path)
        self.assertTrue((res_D == 20 * raw_data).all().compute())

        zarr_dirs = [f.name for f in zarr_cache_path.iterdir()]
        for name in [
            'A',
            'B',
            'C',
            'D'
        ]:
            self.assertTrue(name in zarr_dirs)

        # check the times the cache files have been touched
        # first to make sure that they are not touched
        td_before = time_dict(zarr_cache_path)
        # now with rec_rm flag
        res_D = cached_D(
            cable_out_path,
            zarr_cache_path,
            rec_rm=True
        )
        td_after_rec_rm = time_dict(zarr_cache_path)
        # check that all dependencies have been recreated
        self.assertTrue(
            all(
                [
                    td_before[key] < td_after_rec_rm[key]
                    for key in td_before.keys()
                ]
            )
        )

