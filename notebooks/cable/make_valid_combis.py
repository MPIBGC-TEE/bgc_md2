
import sys
from getpass import getuser
import xarray as xr
import numpy as np
import dask.array
import zarr as zr
from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
import shutil
import matplotlib.pyplot as plt
import bgc_md2.models.cable_all.cableHelpers as cH
from bgc_md2.helper import batchSlices
from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
from time import time
# %load_ext autoreload 
# %autoreload 2

if __name__ == '__main__':
    import dask.array
    from dask.distributed import LocalCluster,Client 
    if 'cluster' not in dir():
        cluster = LocalCluster()

    client = Client(cluster)

try:
    from ports.server_helpers import print_commands
    print_commands(cluster,local_port=8880)
    
except ImportError as e:
    pass  # module doesn't exist,dont make a fuss 
    

out_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4'
out_path = Path(out_dir)
zarr_sub_dir_name = 'zarr'
zarr_dir_path= out_path.joinpath(zarr_sub_dir_name)
dad=cH.load_or_make_cable_dask_array_dict(
        out_path,
        zarr_dir_path,
        rm=False,
        batch_size=64
)

sl=slice(0,32)
B_val, u_val ,x0_val = cH.load_or_make_valid_B_u_x0_slice(
    out_path,
    zarr_sub_dir_name,
    sl,
    rm=False
)
