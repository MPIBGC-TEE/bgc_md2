import dask.array
from pathlib import Path

from bgc_md2.models.cable_all.cableHelpers import (
        cacheWrapper,
        fromCachedArgsMaker,
        cachedFromCachedArgsMaker
)
from TestCache import TestCache, time_dict

###################################################################
# Example fucntions
# These will be created automatically later
# from simpler functions that do not contain
# cache specific arguments but this is how they look like
# if created manually.
# It becomes also obvious that the necessity to pass through
# arguments is error prone and obscures the actual computation

def A_from_raw(
        cable_out_path,
        zarr_cache_path, # ignored, but has to be present for the caller
        rec_rm: bool,   # ignored, but has to be present for the caller
        landpoint_slice: slice = slice(None,None,None), 
        batch_size: int = 1  # ignored, but has to be present for the caller
):
    # Example for a variable that does not depend on other cached vars but
    # only on the original dataset which is hardcoded like 'iveg' or 'Cplant'
    raw_data = dask.array.from_zarr('raw_data')[...,landpoint_slice,:]
    A = dask.array.add(raw_data, raw_data)
    return A


def B_from_raw(
        cable_out_path,
        zarr_cache_path, # ignored, but has to be present for the caller
        rec_rm: bool,   # ignored, but has to be present for the caller
        landpoint_slice: slice = slice(None, None, None), 
        batch_size: int = 1 # ignored, but has to be present for the caller
):

    # Example for a variable that does not depend on other cached vars but
    # only on the original dataset which is hardcoded like 'iveg' or 'Cplant'
    raw_data = dask.array.from_zarr('raw_data')[...,landpoint_slice,:]
    B = dask.array.add(raw_data, raw_data)
    return B


def C_from_cached_args(
        cable_out_path,
        zarr_cache_path,
        rec_rm: bool = False,
        landpoint_slice: slice = slice(None,None,None),
        batch_size: int = 1  # has to be passed through
):
    # example for a dependent function
    A = cacheWrapper(
            A_from_raw,
            cable_out_path,
            zarr_cache_path,
            rec_rm=rec_rm,
            landpoint_slice=landpoint_slice,
            batch_size=batch_size
    )
    B = cacheWrapper(
            B_from_raw,
            cable_out_path,
            zarr_cache_path,
            rec_rm=rec_rm,
            landpoint_slice=landpoint_slice,
            batch_size=batch_size
        )
    return A+B


def D_from_cached_args(
    cable_out_path,
    zarr_cache_path,
    rec_rm: bool = False,  # has to be passed on
    landpoint_slice: slice = slice(None, None, None),  # has to be passed on
    batch_size: int = 1 # has to be passed on
):
    # example for an indirectly dependent function
    C = cacheWrapper(
            C_from_cached_args,
            cable_out_path,
            zarr_cache_path,
            rec_rm=rec_rm,
            landpoint_slice=landpoint_slice,
            batch_size=batch_size
    )
    return 5*C


class TestCacheWrapper(TestCache):

    def test_root(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data

        # compute variables directly from the original data
        res_A = cacheWrapper(
            A_from_raw,
            cable_out_path,
            zarr_cache_path,
            batch_size=32
        )
        self.assertTrue((res_A == 2 * raw_data).all().compute())

        # compute variables directly from the original data
        # with batchsize specified
        res_B = cacheWrapper(
            B_from_raw,
            cable_out_path,
            zarr_cache_path,
            batch_size=32
        )
        self.assertTrue((res_B == 2 * raw_data).all().compute())

        # check if a zarr dir has been created
        zarr_dirs = [f.name for f in zarr_cache_path.iterdir()]
        for name in [
            'A_from_raw',
            'B_from_raw'
        ]:
            self.assertTrue(name in zarr_dirs)
    
    def test_root_sliced(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data

        # compute variables directly from the original data
        landpoint_slice = slice(0,10)
        raw_data_sliced = raw_data[...,landpoint_slice,:]
        res_A = cacheWrapper(
            A_from_raw,
            cable_out_path,
            zarr_cache_path,
            landpoint_slice=landpoint_slice,
            batch_size=32
        )
        self.assertTrue((res_A == 2 * raw_data_sliced).all().compute())

        # compute variables directly from the original data
        # with batchsize specified
        res_B = cacheWrapper(
            B_from_raw,
            cable_out_path,
            zarr_cache_path,
            landpoint_slice=landpoint_slice,
            batch_size=32
        )
        self.assertTrue((res_B == 2 * raw_data_sliced).all().compute())

        # check if a zarr dir has been created
        zarr_dirs = [f.name for f in zarr_cache_path.iterdir()]
        for name in [
            'A_from_raw',
            'B_from_raw'
        ]:
            self.assertTrue(name in zarr_dirs)

    def test_root_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data

        # compute variables directly from the original data
        res_B = cacheWrapper(
            B_from_raw,
            cable_out_path,
            zarr_cache_path
        )

        # check the times the cache files have been touched
        # first to make sure that they are not touched 
        td_before = time_dict(zarr_cache_path)

        #now with rm flag
        res_B = cacheWrapper(
            B_from_raw,
            cable_out_path,
            zarr_cache_path,
            rm=True
        )
        td_after_rm = time_dict(zarr_cache_path)
        key = 'B_from_raw'

        self.assertTrue(td_before[key] < td_after_rm[key])

    def test_root_rec_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data

        # compute variables directly from the original data
        res_B = cacheWrapper(
            B_from_raw,
            cable_out_path,
            zarr_cache_path
        )

        # check the times the cache files have been touched
        td_before = time_dict(zarr_cache_path)

        # now with rec_rm flag which in this case just implies rm
        # but has no other impact since there is no recursion
        res_B = cacheWrapper(
            B_from_raw,
            cable_out_path,
            zarr_cache_path,
            rec_rm=True
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

        # compute something that only
        res_D = cacheWrapper(
            D_from_cached_args,
            cable_out_path,
            zarr_cache_path
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
            cable_out_path,
            zarr_cache_path
        )

        td_after = time_dict(zarr_cache_path)
        self.assertEqual(td_before, td_after)

    def test_indirectly_dependent_rm(self):
        zarr_cache_path = self.zarr_cache_path
        cable_out_path = self.cable_out_path
        raw_data = self.raw_data

        # compute something that only
        res_D = cacheWrapper(
            D_from_cached_args,
            cable_out_path,
            zarr_cache_path
        )
        self.assertTrue((res_D == 20 * raw_data).all().compute())

        # check the times the cache files have been touched
        # first show that they are untouched by the recent computation
        td_before = time_dict(zarr_cache_path)
        # now with rm flag
        res_D = cacheWrapper(
            D_from_cached_args,
            cable_out_path,
            zarr_cache_path,
            rm=True
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

        # compute something that only
        res_D = cacheWrapper(
            D_from_cached_args,
            cable_out_path,
            zarr_cache_path
        )
        self.assertTrue((res_D == 20 * raw_data).all().compute())

        # check the times the cache files have been touched
        # first show that they are untouched by the recent computation
        td_before = time_dict(zarr_cache_path)

        # now with rec_rm flag
        res_D = cacheWrapper(
            D_from_cached_args,
            cable_out_path,
            zarr_cache_path,
            rec_rm=True
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

