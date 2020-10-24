# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from getpass import getuser
import xarray as xr
import numpy as np
import dask
from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
from sympy import Symbol,var
import bgc_md2.models.cable_all.cableHelpers as cH

# To be able to work independently we
# set up independent dask clusters for different users 
# with different dashboards on different ports 
# You can use the same cluster from different notebooks though ()
#  
# Note that only the dashboard port (in addition of your jupyterlab or jupyter  port) 
# has to be forwarded to your local maschine since the scheduler port  
# will be used for communication between jupyter and the dask cluster who run on the same machine (matagorda or antakya)
# It is only specified in order to reuse it from different notebooks running simultaneously on the same machine
port_dict = {
    'mm':(8689,8789),       # first is the port of the actual scheduler, second the port for the dashboard
    'hmetzler':(8690,8790), # change at will to a port you forward via ssh to your local machine
    'cs':(8691,8791)        # change at will
}
my_user_name = getuser()
addr = 'localhost:'+str(port_dict[my_user_name][0])
try:
    Client(addr) 
    # my cluster is listening already on my port (started by another running notebook)
    # so I can use it and its dashboard
except IOError:
    # I have to instanciate one, but do so on a specific port that I can use in another notebook
    my_cluster = LocalCluster(
        scheduler_port= port_dict[my_user_name][0],   
        dashboard_address='localhost:'+str(port_dict[my_user_name][1])
    )
    Client(my_cluster)#same as Client(addr)
print(my_user_name)



# cableDataDir = os.environ["cableDataDir"] # within nix-shell
cableDataDir = '/home/data/cable-data/example_runs'
runId = "parallel_1901_2004_with_spinup"
outDir = "output/new4"
first_yr = 1901
last_yr = 2004
fns = ["out_ncar_" + str(yr) + "_ndep.nc" for yr in range(first_yr, last_yr + 1)]
outpath = Path(cableDataDir).joinpath(runId, outDir)
ps = [outpath.joinpath(fn) for fn in fns]
# just have a peek at the 
dat0 = xr.open_dataset(ps[0])
# have a peek at the first file
dat0

# assemble all files into one dataset
ds = cH.cable_ds(outpath,first_yr,last_yr)
for s in ('leaf','wood','fine_root','metabolic_lit','structural_lit','cwd','fast_soil','slow_soil','passive_soil'):
    var(s)
stateVariableTuple=(leaf,fine_root,wood,metabolic_lit,structural_lit,cwd,fast_soil,slow_soil,passive_soil)
npool=len(stateVariableTuple)
npatch=ds.dims['patch']
nland=ds.dims['land']
ntime=ds.dims['time']

ds.fromLeaftoL.sel(litter_casa_pools=[2],land=100)[0:-1:5000].plot(hue='patch')
ds.fromLeaftoL.sel(litter_casa_pools=2,land=100).mean('patch')[0:-1:5000]

post_processing_dir= "../src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_org"
dfac = xr.open_dataset(Path(post_processing_dir).joinpath('outAC.nc'))
dfac


for p in range(npatch):
   print(
       dfac.A.sel(patch=p).mean('land').compute()
   )

# +
# this is a reconstruction of 
# `bgc_md2/src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_org/GetFluxFromCABLEoutput.txt`
# Looking at the "fluxes" computed here we should be able to infer the  compartmental  matrix B
# Note: 
# 1.) B!=A*diag(C)  since this would be constant over a year as is apparent from the 
#     way A and C are constructed.
# 2.) We do not know yet if the variables here called Fluxes are actually accumulated over a Timestep or momentary
#     values taken at the time


C=cH.reconstruct_C_diag(ps[0])
A=cH.reconstruct_A(ps[0])
leaf_ind = 0
wood_ind = 1
root_ind = 2

metabolic_ind =0
structural_ind=1
cwd_ind=2

fast_ind=0
slow_ind=1
passive_ind=2

# Input fluxes
InputFluxes={
    leaf      : (ds.NPP * ds.fracCalloc.sel(plant_casa_pools=leaf_ind)),
    wood      : (ds.NPP * ds.fracCalloc.sel(plant_casa_pools=wood_ind)),
    fine_root : (ds.NPP * ds.fracCalloc.sel(plant_casa_pools=root_ind))
}

#Internal Fluxes

#helper
Flux_from_Leaf = ds.kplant.sel(plant_casa_pools=leaf_ind) * ds.Cplant.sel(plant_casa_pools=leaf_ind)
Flux_from_Root = ds.kplant.sel(plant_casa_pools=root_ind) * ds.Cplant.sel(plant_casa_pools=root_ind)
Flux_from_Wood = ds.kplant.sel(plant_casa_pools=wood_ind) * ds.Cplant.sel(plant_casa_pools=wood_ind)
xk_t_w         = ds.xktemp*ds.xkwater

Flux_from_fast_soil =  ds.Csoil.sel(soil_casa_pool=fast_ind)*xk_t_w
Flux_from_slow_soil =  ds.Csoil.sel(soil_casa_pool=slow_ind)*xk_t_w

xk_t_w_Nl      = xk_t_w*ds.xkNlimiting
Flux_from_cwd  = ds.Clitter.sel(litter_casa_pools=cwd_ind)*xk_t_w_Nl

Flux_from_structural_lit = ds.Clitter.sel(litter_casa_pools=structural_ind)*xk_t_w_Nl
Flux_from_metabolic_lit = ds.Clitter.sel(litter_casa_pools=metabolic_ind)*xk_t_w_Nl

InternalFluxes= {
    # 4. Leaf to Metoblic litter
    # fin->fromLeaftoL(:,0,:,:)*fin->Cplant(:,0,:,:)*fin->kplant(:,0,:,:)
    (leaf,metabolic_lit)       : ds.fromLeaftoL.sel(litter_casa_pools=metabolic_ind) * Flux_from_Leaf, 
    #
    # 5. Leaf to Structural litter
    # fin->fromLeaftoL(:,1,:,:)*fin->Cplant(:,0,:,:)*fin->kplant(:,0,:,:)
    (leaf,structural_lit)      : ds.fromLeaftoL.sel(litter_casa_pools=structural_ind)* Flux_from_Leaf, 
    #
    # 6. Root to Metoblic litter
    # fin->fromRoottoL(0,0,:,:)*fin->Cplant(:,2,:,:)*fin->kplant(:,2,:,:)
    (fine_root,metabolic_lit)  : ds.fromRoottoL.sel(litter_casa_pools=metabolic_ind) * Flux_from_Root, 
    #
    #7. Root to Structural litter
    # fin->fromRoottoL(0,1,:,:)*fin->Cplant(:,2,:,:)*fin->kplant(:,2,:,:)
    (fine_root,structural_lit) : ds.fromRoottoL.sel(litter_casa_pools=structural_ind)* Flux_from_Root, 
    #
    # 8. Wood to CWD
    # fin->fromWoodtoL(:,2,:,:)*fin->Cplant(:,1,:,:)*fin->kplant(:,1,:,:)
    (wood,cwd)                 : ds.fromWoodtoL.sel(litter_casa_pools=cwd_ind)       * Flux_from_Wood, 
    #
    # 9. Metabolic litter to Fast soil
    # fAC->A(6,3,:,:)*fAC->C(3,:,:)*fin->Clitter(:,0,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    (metabolic_lit,fast_soil)  : A[6,3,:,:] * C[3,:,:] * Flux_from_metabolic_lit,
    #
    # 11. Structural Litter to Fast soil: 
    #fAC->A(6,4,:,:)*fAC->C(4,:,:) *fin->Clitter(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    (structural_lit,fast_soil) : A[6,4,:,:] * C[4,:,:] * Flux_from_structural_lit,
    #
    # 12. Structural Litter to Slow soil: 
    # fAC->A(7,4,:,:)*fAC->C(4,:,:) *fin->Clitter(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    (structural_lit,slow_soil) : A[7,4,:,:] * C[4,:,:] * Flux_from_structural_lit,
    #
    # 14. CWD to fast soil: 
    # fAC->A(6,5,:,:)*fAC->C(5,:,:) *fin->Clitter(:,2,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    (cwd,fast_soil)            : A[6,5,:,:] * C[5,:,:] * Flux_from_cwd,
    #
    # 15. CWD to slow soil: 
    # fAC->A(7,5,:,:)*fAC->C(5,:,:) missing *fin->Clitter(:,2,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    (cwd,slow_soil)            : A[7,5,:,:] * C[5,:,:] * Flux_from_cwd,
    #
    # 17. fast soil to slow soil 
    # fAC->A(7,6,:,:)*fAC->C(6,:,:) *fin->Csoil(:,0,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    (fast_soil,slow_soil)      : A[7,6,:,:] * C[6,:,:] * Flux_from_fast_soil,
    #
    # 18. fast soil to passive soil 
    # fAC->A(8,6,:,:)*fAC->C(6,:,:) *fin->Csoil(:,0,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    (fast_soil,passive_soil)   : A[8,6,:,:] * C[6,:,:] * Flux_from_fast_soil,
    #
    # 20. slow soil to passive soil 
    # fAC->A(8,7,:,:)*fAC->C(7,:,:) *fin->Csoil(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    (slow_soil,passive_soil)   : A[8,7,:,:] * C[7,:,:] * Flux_from_slow_soil,


    
 }
OutFluxes = {
    # 10. Metabolic litter to atmosphere
    # (1-fAC->A(6,3,:,:))*fAC->C(3,:,:) *fin->Clitter(:,0,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    metabolic_lit  : A[6,3,:,:] * C[3,:,:] * Flux_from_metabolic_lit,
    #
    # 13. structural Litter to atmosphere 
    # (1-fAC->A(6,4,:,:)-fAC->A(7,4,:,:))*fAC->C(4,:,:) *fin->Clitter(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    structural_lit : (1-A[6,4,:,:]-A[7,4,:,:])* C[4,:,:] * Flux_from_structural_lit,
    #
    # 16. CWD to atmosphere 
    # (1-fAC->A(6,5,:,:)-fAC->A(7,5,:,:))*fAC->C(5,:,:) 
    # *fin->Clitter(:,2,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    cwd            : (1-A[6,5,:,:]-A[7,5,:,:])*C[5,:,:] * Flux_from_cwd,
    #
    # 19. fast soil to atmosphere (1-fAC->A(7,6,:,:)-fAC->A(8,6,:,:))*fAC->C(6,:,:) 
    # *fin->Csoil(:,0,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    fast_soil      : (1-A[7,6,:,:]-A[8,6,:,:])*C[6,:,:] * Flux_from_fast_soil,
    #
    # 21. slow soil to atmosphere (1-fAC->A(8,7,:,:))*fAC->C(7,:,:) 
    # *fin->Csoil(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    slow_soil      : (1-A[8,7,:,:])*C[7,:,:]* Flux_from_slow_soil
}
# -

# We now reconstruct the matrix B by factoring out the pool contents from The Flux_from_... terms
#
# The first task is to identify the state variables in the output file
#  A@pool_name  = (/"leaf,root,wood,metabolic,structure,CWD,fast,slow,passive"/)
#  C@pool_name  = (/"leaf,root,wood,metabolic,structure,CWD,fast,slow,passive"/)

#time dependent B matrix (The first index indexes time):
B_chunk = (ntime,npool,npool,npatch)
B_shape = B_chunk+(nland,)
#B=dask.array.full(B_shape,np.nan,dtype='float64',chunksize=B_chunk)
B=dask.array.where(
    ds.iveg==ds.iveg.attrs['_FillValue'],
    dask.array.full(B_shape,np.nan,chunks=B_chunk+(1,)),
    dask.array.zeros(B_shape,chunks=B_chunk+(1,)),
)
B
