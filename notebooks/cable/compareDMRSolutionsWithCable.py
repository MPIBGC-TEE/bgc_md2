import xarray as xr
from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
import bgc_md2.models.cable_all.cableCache as cC
import bgc_md2.models.cable_all.cablePaths as cP
import bgc_md2.models.cable_all.cableHelpers as cH

if __name__ == '__main__':
    import dask.array
    from dask.distributed import LocalCluster,Client 
    if 'cluster' not in dir():
        cluster = LocalCluster()

    client = Client(cluster)

# +
try:
    from ports.server_helpers import print_commands
    print_commands(cluster, local_port=8880)
    
except ImportError as e:
    pass  # module doesn't exist,dont make a fuss 
# -

# chose the cable output directory you want to work with
cable_out_path = Path('/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4')
time_slice=slice(0,1000)
landpoint_slice = slice(0,1024)

zarr_cache_path = cP.slice_dir_path(
    cable_out_path,
    sub_dir_trunk='zarr',
    time_slice=time_slice,
    landpoint_slice=landpoint_slice
)

cable_data_set = cH.cable_ds(cable_out_path)
args={
    'cable_data_set': cable_data_set,
    'zarr_cache_path': zarr_cache_path,
    'landpoint_slice': landpoint_slice,
    'time_slice': time_slice,
    'batch_size': 128,
    'rec_rm': True
}

#Clitter= cH.cacheWrapper(
#    cC.Clitter,
#    **args
#)
#sol_val= cH.cacheWrapper(
#    cC.sol_val,
#    **args
#)
x_val= cH.cacheWrapper(
    cC.x_val,
    **args
)
