from pathlib import Path

"""This module serves to implement a policy for data storage locations
relative to a result data set produced by a cable run.
The purpose is:
- to make using the same cache files from different notebooks a bit easier.
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

def zarr_path(
        cable_out_path: Path,
        sub_dir_trunk: str = 'zarr'
):
    return cable_out_path.joinpath(sub_dir_trunk)

def slice_dir_path(
        cable_out_path,
        sub_dir_trunk,
        time_slice: slice = slice(None, None, None),
        landpoint_slice: slice = slice(None, None, None)
) -> Path:
    """Returns the default path to the directory with the zarr files
    as constructed from the time and landpoint slices to avoid
    confusion with the unsliced data.

    The sliced arrays could of cause be derived from the complete 
    arrays and need not be stored. 
    But the sliced arrays are much quicker to compute, especially if
    there are changes in one of the base variables which require 
    the recursive regeneration of all the derived arrays.

    It is just a help for experimenting 
    """
    time_suff = "" if time_slice == slice(None, None, None) else "_time_" + slice_suffix(time_slice)
    landpoint_suff = "" if landpoint_slice == slice(None, None, None) else "_landpoints_" + slice_suffix(landpoint_slice)
    return cable_out_path.joinpath(sub_dir_trunk+time_suff+landpoint_suff)


def slice_suffix(
    sl: slice
):
    from_str = "" if sl.start is None else "_from_" + str(sl.start)
    to_str = "" if sl.stop is None else "_to_" + str(sl.stop)
    step_str = "" if sl.step is None else "_step_" + str(sl.step)
    suffix = from_str + to_str + step_str
    return suffix
