from bgc_md2.models.cable_all.cableHelpers import cable_ds
from dask.distributed import Client
from pathlib import Path
import bgc_md2.models.cable_all.cableHelpers as cH

#if __name__ == '__main__':
from bgc_md2.sitespecificHelpers import get_client
client = get_client()
#client =Client(scheduler_file='/home/mm/scheduler.json')

example_run_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4'
ds = cable_ds(example_run_dir)

dir_name= example_run_dir+'/zarr_vars'
cH.write_vars_as_zarr(ds,dir_name)

# +

# This code never finishes
# dir_name2= example_run_dir+'xarray_ds_zarr'
# ds.to_zarr(dir_name2) #serial
# -


