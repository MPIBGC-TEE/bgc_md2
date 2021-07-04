from pathlib import Path

"""This module serves to implement a policy for data storage locations relative to 
a result data set produced by a cable run.
The purpose is 
- to make using the same files from different notebooks a bit easier.
- to localize changes to the policy here and avoid having to change 
  all the notebooks.

At the moment the policy consists in the following assumptions:
    1. The original cable output directory  (containing the NetCdf) files
       is not used to write files. 
       Instead directories are created under it. 
       1. The directory containing zarr directories for all computed/converted 
       variables
       1. Directories for testing containing slices of the above
    1. The names of the cached zarr arrays are determined from the name
       of the function that computes them and have to be unique

    
"""


valid_combi_dir_name = "zarr_valid_landpoint_patch_combis"

def zarr_path(netcdf_dir_path: Path):
    return netcdf_dir_path.joinpath("zarr")

def slice_dir_path(
    """Returns the default path to the directory with the zarr files
    for the valid combinations of Landpoints and patches"""
    
        netcdf_dir_path: 'Path',
        sl: slice=slice(None,None,None)
    ) -> Path:
    
    from_str = "" if sl.start is None else "_from_" + str(sl.start)
    to_str = "" if sl.stop is None else "_to_" + str(sl.stop)
    step_str = "" if sl.step is None else "_step_" + str(sl.step)
    suffix = from_str + to_str + step_str 
    return netcdf_dir_path.joinpath(valid_combi_dir_name+suffix)

def slice_suffix(
        sl :slice):
    from_str = "" if sl.start is None else "_from_" + str(sl.start)
    to_str = "" if sl.stop is None else "_to_" + str(sl.stop)
    step_str = "" if sl.step is None else "_step_" + str(sl.step)
    suffix = from_str + to_str + step_str 
    return suffix
