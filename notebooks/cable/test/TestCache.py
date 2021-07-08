from pathlib import Path
from testinfrastructure.InDirTest import InDirTest
import dask.array

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
        dask.array.to_zarr(arr=raw_data, url='raw_data')

        # create a directory for the zarr files as cache for the computed
        # variables
        zarr_cache_path = Path(".").joinpath('cache')

        self.cable_out_path = cable_out_path
        self.zarr_cache_path = zarr_cache_path
        self.raw_data = raw_data

