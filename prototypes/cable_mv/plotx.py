import matplotlib.pyplot as plt
import dask.array
from sympy import latex
from pathlib import Path

from dask.distributed import Client
from dask.distributed import LocalCluster
import bgc_md2.models.cable_all.cableCache as cC
import bgc_md2.models.cable_all.cablePaths as cP
import bgc_md2.models.cable_all.cableHelpers as cH

if __name__ == "__main__":
    if "cluster" not in dir():
        # cluster = LocalCluster(n_workers = 1)
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
cable_out_path = Path(
    "/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4"
)
time_slice = slice(0, 2000)
landpoint_slice = slice(0, 2028)
# landpoint_slice = slice(None,None)
# time_slice=slice(None,None,None)

zarr_cache_path = cP.slice_dir_path(
    cable_out_path,
    sub_dir_trunk="zarr_mm11",
    # sub_dir_trunk='zarr_mm',
    time_slice=time_slice,
    landpoint_slice=landpoint_slice,
)
if "cable_data_set" not in dir():
    cable_data_set = cH.cable_ds(cable_out_path)

args = {
    "cable_data_set": cable_data_set,
    "zarr_cache_path": zarr_cache_path,
    "landpoint_slice": landpoint_slice,
    "time_slice": time_slice,
    #'batch_size': 128,
    "batch_size": 8,
    #'rm': True
}
x_org_iveg = cH.cacheWrapper(cC.x_org_iveg, **args)
times = cH.cacheWrapper(cC.time, **args)
patches, landpoints = cC.all_pools_vary_cond_nz(**args)
#lets find the first combinations that fulfill our requirements
my_sl= slice(0,96)
pcs = patches[my_sl].compute()
lpcs = landpoints[my_sl].compute()
pcs,lpcs
# p = Path('plots')
# p.mkdir(exist_ok=True)
# for lp in range(landpoint_slice.start, landpoint_slice.stop):
#    for patch in range(10):
#        print(x_org_iveg[0,:,patch,lp])
#        fig=plt.figure()
#        for pool in range(n):
#            ax = fig.add_subplot(n+1, 1, 2+pool)
#            title =latex(pool)
#            ax.plot(times, x_org_iveg[:,pool,patch,lp])
#            fontsize = 10
#            ax.set_title(title, fontsize=fontsize)
#        fig.savefig(p.joinpath('sol_'+str(patch)+'_'+str(lp)+'.pdf'))
