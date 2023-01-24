from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
import bgc_md2.models.cable_all.cableCache as cC
import bgc_md2.models.cable_all.cablePaths as cP
import bgc_md2.models.cable_all.cableHelpers as cH

if __name__ == "__main__":
    if "cluster" not in dir():
        cluster = LocalCluster(n_workers=1)
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
    "batch_size": 64,
    #'batch_size': 8,
    #'rm': True
}
sol_org_iveg = cH.cacheWrapper(cC.sol_org_iveg, **args)

sol_org_iveg

# x_val = cH.cacheWrapper(
#     cC.x_val,
#     **args
# )
# # check that the solution is non negative
# assert(((x_val>=0).all()).compute())
#
# # now check that the solution from the reconstructed Bs is non negative
# sol_val = cH.cacheWrapper(
#     cC.sol_val,
#     **args
# )
# this is unfortunately not true
# assert(((sol_val>=0).all()).compute())

# # assess absolute deviation
# diff = x_val-sol_val
# print("maximum absolute deviation",(dask.array.absolute(diff).max()).compute())
#
# #
# # assess relative deviation
# x_val_non_zero = dask.array.ma.masked_array(x_val,mask=x_val==0)
#
# print("maximum relative deviation",((dask.array.absolute(diff)/x_val_non_zero).max()).compute())
#
# the objective is to check if the errors in the reproduction have a spatial pattern, In the best case scenario the problem can be traced back
# to some grid points
# rel_err = dask.array.absolute(diff)/x_val_non_zero
# max_rel_err_per_trajectory = dask.array.max(rel_err,(0,1))
# cond=max_rel_err_per_trajectory < 1e-1
# inds, = dask.array.nonzero(cond)


# x0_val = cH.cacheWrapper(
#    cC.x0_val,
#    **args
# )
# u_val = cH.cacheWrapper(
#    cC.u_val,
#    **args
# )
# B_val = cH.cacheWrapper(
#    cC.B_val,
#    **args
# )
# x_org= cH.cacheWrapper(
#    cC.x_org,
#    **args
# )

# x0_org= cH.cacheWrapper(
#     cC.x0_org_iveg,
#     **args
# )
# B_org= cH.cacheWrapper(
#     cC.B_org_iveg,
#     **args
# )
# u_org= cH.cacheWrapper(
#     cC.u_org_iveg,
#     **args
# )
# #sol_org= cH.cacheWrapper(
# #    cC.sol_org,
# #    **args
# #
# # sol_org
# # osl = slice(0,1)
# X0 = x0_org.rechunk((x0_org.shape[0],1,1))
# U = u_org.rechunk((u_org.shape[0],u_org.shape[1],1,1))
# B = B_org.rechunk((B_org.shape[0],B_org.shape[1],B_org.shape[2],1,1))
# times = cH.cacheWrapper(cC.time,**args)
# sol_org = dask.array.blockwise(
#              cH.trajectory_org, 'ijkl',
#              X0, 'jkl',
#              times, 'i',
#              B, 'ijjkl',
#              U, 'ijkl',
#              dtype='f8'
#     )
# sol_org
# sol_org.compute()
#
# print(cH.trajectory_org(X0[...,osl,osl], times, B[...,osl,osl], U[...,osl,osl]))
# nz = cH.cacheWrapper(
#    cC.nz,
#    **args
# )
# Csoil= cH.cacheWrapper(
#    cC.Csoil,
#    **args
# )
# Csoil0= cH.cacheWrapper(
#    cC.Csoil0,
#    **args
# )
