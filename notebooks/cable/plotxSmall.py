import matplotlib.pyplot as plt
from pathlib import Path

from dask.distributed import Client
from dask.distributed import LocalCluster
import bgc_md2.models.cable_all.cableCache as cC
import bgc_md2.models.cable_all.cablePaths as cP
import bgc_md2.models.cable_all.cableHelpers as cH

if __name__ == "__main__":
    if "cluster" not in dir():
        cluster = LocalCluster(
            n_workers=1, 
            memory_limit="14GB"
        )
        # cluster = LocalCluster()

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
time_slice = slice(0, None)
landpoint_slice = slice(1590, 1637) # cheated
# landpoint_slice = slice(None,None)
# time_slice=slice(None,None,None)

zarr_cache_path = cP.slice_dir_path(
    cable_out_path,
    sub_dir_trunk="zarr_mm11",
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
    "batch_size": 12,
    #'rm': True
}
#iveg_mask = cH.cacheWrapper(cC.iveg_mask, **args)
#x_org= cH.cacheWrapper(cC.x_org, **args)
x_org_iveg = cH.cacheWrapper(cC.x_org_iveg, **args)
sol_org_iveg = cH.cacheWrapper(cC.sol_org_iveg, **args)
time = cH.cacheWrapper(cC.time, **args)
patches, landpoints = cC.all_pools_vary_cond_nz(**args)
pcs = patches.compute()
lpcs = landpoints.compute()
pcs,lpcs
p = Path('plots')
p.mkdir(exist_ok=True)
ind = 0 # first pair that has a nonconstant solution for 
lp = lpcs[ind]
patch = lpcs[ind]
print(x_org_iveg[0,:,patch,lp])
n=9
fig=plt.figure(figsize=(15,25))
my_time_sl = slice(0,2000)
for pool in range(n):
    ax = fig.add_subplot(n+1, 1, 1+pool)
    #title =latex(pool)
    ax.plot(time[my_time_sl], x_org_iveg[my_time_sl ,pool, patch,lp], color='b')
    ax.plot(time[my_time_sl], sol_org_iveg[my_time_sl ,pool, patch,lp], color='r')
    #fontsize = 10
    #ax.set_title(title) #, fontsize=fontsize)
fig.savefig(p.joinpath('sol_'+str(patch)+'_'+str(lp)+'.pdf'))
