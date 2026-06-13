import dask.array
from unittest import skip
from pathlib import Path
from testinfrastructure.InDirTest import InDirTest

from helpers import (
        cacheWrapper,
        fromCachedArgsMaker2,
        cachedFromCachedArgsMaker2
)

# helpers
def time_dict(zarr_cache_path):
    subPaths = zarr_cache_path.iterdir()
    return {
        p.name: (p.stat()).st_mtime
        for p in subPaths
    }

class TestCache(InDirTest):
    """helper class for the setUp method shared by all the childred classes"""
    def setUp(self):
        # simulate a cable data_set that we do not touch
        # this is a placeholder for the fortran generated cable output dir

        cable_out_path = Path(".").joinpath('org')
        # create a file
        cable_out_path.mkdir()
        raw_data = dask.array.ones((3, 3, 100, 9))
        dask.array.to_zarr(arr=raw_data, url=str(cable_out_path.joinpath('raw_data')))

        # create a directory for the zarr files as cache for the computed
        # variables
        zarr_cache_path = Path(".").joinpath('cache')

        self.cable_out_path = cable_out_path
        self.zarr_cache_path = zarr_cache_path
        self.raw_data = raw_data

#from TestCache import TestCache, time_dict

###################################################################
# Example fucntions
# These could be created automatically later
# from simpler functions that do not contain
# cache specific arguments but this is how they look like
# if created manually.


def A_from_raw(
        cable_out_path,
        landpoint_slice=slice(None,None,None),
        **kwargs
    ):
    # Example for a cached variable that does not depend on other cached vars but
    # only on the original dataset which is hardcoded like 'iveg' or 'Cplant'
    url = str(cable_out_path.joinpath('raw_data'))
    raw_data = dask.array.from_zarr(url)[..., landpoint_slice, :]
    A = dask.array.add(raw_data, raw_data)
    return A


def B_from_raw(
        cable_out_path, 
        landpoint_slice=slice(None,None,None),
        **kwargs
    ):
    # Example for a variable that does not depend on other cached vars but
    # only on the original dataset which is hardcoded like 'iveg' or 'Cplant'
    url = str(cable_out_path.joinpath('raw_data'))
    raw_data = dask.array.from_zarr(url )[...,landpoint_slice,:]
    B = dask.array.add(raw_data, raw_data)
    return B


def C_from_cached_args(**kwargs):
    # example for a dependent function
    A = cacheWrapper(A_from_raw, **kwargs)
    B = cacheWrapper(B_from_raw, **kwargs)
    return A+B


def D_from_cached_args(**kwargs):
    C = cacheWrapper(C_from_cached_args, **kwargs)
    return 5*C

# helpers
def time_dict(zarr_cache_path):
    subPaths = zarr_cache_path.iterdir()
    return {
        p.name: (p.stat()).st_mtime
        for p in subPaths
    }

#class TestCache(InDirTest):
#    """helper class for the setUp method shared by all the childred classes"""
#    def setUp(self):
#        # simulate a cable data_set that we do not touch
#        # this is a placeholder for the fortran generated cable output dir
#
#        cable_out_path = Path(".").joinpath('org')
#        # create a file
#        cable_out_path.mkdir()
#        raw_data = dask.array.ones((3, 3, 100, 9))
#        dask.array.to_zarr(arr=raw_data, url=str(cable_out_path.joinpath('raw_data')))
#
#        # create a directory for the zarr files as cache for the computed
#        # variables
#        zarr_cache_path = Path(".").joinpath('cache')
#
#        self.cable_out_path = cable_out_path
#        self.zarr_cache_path = zarr_cache_path
#        self.raw_data = raw_data


class TestCacheWrapper(TestCache):

    def test_root(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data

        # compute variables directly from the original data
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
            'batch_size': 32
        }
        res_A = cacheWrapper(
            A_from_raw,
            **kwargs
        )
        self.assertTrue((res_A == 2 * raw_data).all().compute())

        # compute variables directly from the original data
        # with batchsize specified
        res_B = cacheWrapper(
            B_from_raw,
            **kwargs
        )
        self.assertTrue((res_B == 2 * raw_data).all().compute())

        # check if a zarr dir has been created
        zarr_dirs = [f.name for f in self.zarr_cache_path.iterdir()]
        for name in [
            'A_from_raw',
            'B_from_raw'
        ]:
            self.assertTrue(name in zarr_dirs)
    
    def test_root_sliced(self):
        raw_data = self.raw_data

        # compute variables directly from the original data
        landpoint_slice = slice(0,10)
        kwargs={
            'cable_out_path': self.cable_out_path,
            'zarr_cache_path': self.zarr_cache_path,
            'landpoint_slice': landpoint_slice
        }
        raw_data_sliced = raw_data[...,landpoint_slice,:]
        res_A = cacheWrapper(
            A_from_raw,
            **kwargs
        )
        self.assertTrue((res_A == 2 * raw_data_sliced).all().compute())

        # compute variables directly from the original data
        # with batchsize specified
        kwargs={
            'cable_out_path': self.cable_out_path,
            'zarr_cache_path': self.zarr_cache_path,
            'batch_size': 32,
            'landpoint_slice': landpoint_slice,
        }
        res_B = cacheWrapper(
            B_from_raw,
            **kwargs
        )
        self.assertTrue((res_B == 2 * raw_data_sliced).all().compute())
 
        # check if a zarr dir has been created
        zarr_dirs = [f.name for f in self.zarr_cache_path.iterdir()]
        for name in [
            'A_from_raw',
            'B_from_raw'
        ]:
            self.assertTrue(name in zarr_dirs)
 
    def test_root_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
            'batch_size': 32
        }

        # compute variables directly from the original data
        res_B = cacheWrapper(
            B_from_raw,
            **kwargs
        )

        # check the times the cache files have been touched
        # first to make sure that they are not touched 
        td_before = time_dict(zarr_cache_path)

        #now with rm flag
        kwargs.update({'rm':True})
        res_B = cacheWrapper(
            B_from_raw,
            **kwargs
        )
        td_after_rm = time_dict(zarr_cache_path)
        key = 'B_from_raw'

        self.assertTrue(td_before[key] < td_after_rm[key])

    def test_root_rec_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path
        }

        # compute variables directly from the original data
        res_B = cacheWrapper(
            B_from_raw,
            **kwargs
        )

        # check the times the cache files have been touched
        td_before = time_dict(zarr_cache_path)

        # now with rec_rm flag which in this case just implies rm
        # but has no other impact since there is no recursion
        kwargs.update({'rec_rm':True})
        res_B = cacheWrapper(
            B_from_raw,
            **kwargs
        )
        td_after_rm = time_dict(zarr_cache_path)
        key = 'B_from_raw'
        # for directly computable vars rec_rm has the same
        # effect as rm and only affects B 
        self.assertTrue(td_before[key] < td_after_rm[key])

    def test_indirectly_dependent(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        kwargs = {
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path
        }

        res_D = cacheWrapper(
            D_from_cached_args,
            **kwargs
        )
        self.assertTrue((res_D == 20 * raw_data).all().compute())

        zarr_dirs = [f.name for f in zarr_cache_path.iterdir()]
        for name in [
            'A_from_raw',
            'B_from_raw',
            'C_from_cached_args',
            'D_from_cached_args'
        ]:
            self.assertTrue(name in zarr_dirs)

        # check the times the cache files have been touched
        # first show that they are untouched by the recent computation
        td_before = time_dict(zarr_cache_path)
        res_D = cacheWrapper(
            D_from_cached_args,
            **kwargs
        )

        td_after = time_dict(zarr_cache_path)
        self.assertEqual(td_before, td_after)

    def test_indirectly_dependent_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data

        kwargs = {
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path
        }
        # compute something that only
        res_D = cacheWrapper(
            D_from_cached_args,
            **kwargs
        )
        self.assertTrue((res_D == 20 * raw_data).all().compute())

        # check the times the cache files have been touched
        # first show that they are untouched by the recent computation
        td_before = time_dict(zarr_cache_path)
        # now with rm flag
        kwargs.update({'rm':True})
        res_D = cacheWrapper(
            D_from_cached_args,
            **kwargs
        )
        td_after_rm = time_dict(zarr_cache_path)
        # check that the newly created result is newer
        key = 'D_from_cached_args'
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
        kwargs = {
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path
        }

        # compute something that only
        res_D = cacheWrapper(
            D_from_cached_args,
            **kwargs
        )
        self.assertTrue((res_D == 20 * raw_data).all().compute())

        # check the times the cache files have been touched
        # first show that they are untouched by the recent computation
        td_before = time_dict(zarr_cache_path)

        # now with rec_rm flag
        kwargs.update({'rec_rm':True})
        res_D = cacheWrapper(
            D_from_cached_args,
            **kwargs
        )
        # check that all dependencies have been recreated
        td_after_rec_rm = time_dict(zarr_cache_path)
        self.assertTrue(
            all(
                [
                    td_before[key] < td_after_rec_rm[key]
                    for key in td_before.keys()
                ]
            )
        )


########################################################################
# now we create the complete recursively caching functions automatically.
class TestCacheFromCachedArgsMaker2(TestCache):

    def test_root_scalar(self):
        # assume we have a function that only needs the 
        # directory without a slice like the fillvalue detectors
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data

        
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
            'batch_size': 32
        }

        import testFuncs  # the module with the functions doing the computing
        cached_F = cachedFromCachedArgsMaker2(
                testFuncs.F,
                testFuncs
        )
        res_F = cached_F(**kwargs)

    def test_root(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
            'batch_size': 32
        }
        import testFuncs  # the module with the functions doing the computing
        cached_A = cachedFromCachedArgsMaker2(
                testFuncs.A,
                testFuncs
        )
        res_A = cached_A(**kwargs)
        self.assertTrue((res_A == 2 * raw_data).all().compute())

        cached_B = cachedFromCachedArgsMaker2(
                testFuncs.B,
                testFuncs
        )
        res_B = cached_B(**kwargs)
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
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
            'landpoint_slice': landpoint_slice
        }
        raw_data_sliced = raw_data[...,landpoint_slice,:]

        cached_A = cachedFromCachedArgsMaker2(
                testFuncs.A,
                testFuncs
        )
        res_A = cached_A(**kwargs)
        self.assertTrue((res_A == 2 * raw_data_sliced).all().compute())

        cached_B = cachedFromCachedArgsMaker2(
                testFuncs.B,
                testFuncs
        )
        res_B = cached_B(**kwargs)
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
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
        }
        import testFuncs  # the module with the functions doing the computing

        cached_B = cachedFromCachedArgsMaker2(
                testFuncs.B,
                testFuncs
        )
        res_B = cached_B(**kwargs)
        # check the times the cache files have been touched
        # first to make sure that they are not touched 
        td_before = time_dict(zarr_cache_path)

        #now with rm flag
        kwargs.update({'rm':True})
        res_B = cached_B(**kwargs)
        td_after_rm = time_dict(zarr_cache_path)
        key = 'B'
        self.assertTrue(td_before[key] < td_after_rm[key])

    def test_root_rec_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
        }
        import testFuncs  # the module with the functions doing the computing

        cached_B = cachedFromCachedArgsMaker2(
                testFuncs.B,
                testFuncs
        )
        res_B = cached_B(**kwargs)
        # check the times the cache files have been touched
        td_before = time_dict(zarr_cache_path)

        #now with rec_rm flag
        kwargs.update({'rec_rm':True})
        res_B = cached_B(**kwargs)
        td_after_rm = time_dict(zarr_cache_path)
        key = 'B'
        # for directly computable vars rec_rm has the same
        # effect as rm and only affects B 
        self.assertTrue(td_before[key] < td_after_rm[key])


    def test_indirectly_dependent(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
        }
        import testFuncs  # the module with the functions doing the computing
        cached_D = cachedFromCachedArgsMaker2(
                testFuncs.D,
                testFuncs
        )
        res_D = cached_D(**kwargs)
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
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
        }
        import testFuncs  # the module with the functions doing the computing
        cached_D = cachedFromCachedArgsMaker2(
                testFuncs.D,
                testFuncs
        )
        res_D = cached_D(**kwargs)
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
        
        kwargs.update({'rm':True})
        res_D = cached_D(**kwargs)

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
        kwargs={
            'cable_out_path': cable_out_path,
            'zarr_cache_path': zarr_cache_path,
        }
        import testFuncs  # the module with the functions doing the computing
        cached_D = cachedFromCachedArgsMaker2(
                testFuncs.D,
                testFuncs
        )
        res_D = cached_D(**kwargs)
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
        kwargs.update({'rec_rm':True})
        res_D = cached_D(**kwargs)
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

