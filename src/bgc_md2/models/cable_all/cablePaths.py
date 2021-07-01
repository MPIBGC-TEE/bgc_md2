from pathlib import Path
valid_combi_dir_name = "zarr_valid_landpoint_patch_combis"
def zarr_path(netcdf_dir_path: Path):
    return netcdf_dir_path.joinpath("zarr")

def zarr_valid_landpoint_patch_combis_slice_path(
        netcdf_dir_path: 'Path',
        sl: slice=slice(None,None,None)
    ) -> Path:
    
    from_str = "" if sl.start is None else "_from_" + str(sl.start)
    to_str = "" if sl.stop is None else "_to_" + str(sl.start)
    step_str = "" if sl.step is None else "_step_" + str(sl.step)
    suffix = from_str + to_str + step_str 
    return netcdf_dir_path.joinpath(valid_combi_dir_name+suffix)


