# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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
from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
from sympy import Symbol,var

# +
# To be able to work independently.
# set up independent dask clusters for different users 
# with different dashboards on different ports 
#  
# The resulting port number has to be forwarded via ssh to look at the dashboard 
#(Alternativly we probably could use a common cluster though) 

port_dict = {
    'mm':8789,
    'hmetzler':8790, # change at will
    'cs':8791        # change at will
}
my_user_name = getuser()
print(my_user_name)

my_port = port_dict[my_user_name]
print(my_port)

my_cluster = LocalCluster(dashboard_address='localhost:'+str(my_port))

# -

Client(my_cluster)

# +
# cableDataDir = os.environ["cableDataDir"] # within nix-shell
cableDataDir = '/home/data/cable-data/example_runs'
runId = "parallel_1901_2004_with_spinup"
outDir = "output/new4"
first_yr = 1901
last_yr = 2004
fns = ["out_ncar_" + str(yr) + "_ndep.nc" for yr in range(first_yr, last_yr + 1)]
outpath = Path(cableDataDir).joinpath(runId, outDir)
ps = [outpath.joinpath(fn) for fn in fns]

dat0 = xr.open_dataset(ps[0])
dat0

# +
xr.open_dataset(Path(cableDataDir).joinpath(runId,'restart_in.nc'))


# -

# assemble all files into one dataset
ds = xr.open_mfdataset(
    paths=ps, 
    combine="by_coords",
    parallel=True # use dask
);ds

# now start postprocessing the raw data to get the fluxes
# following ../src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_org/GetFluxFromCABLEoutput.txt
# We get results for 10 patches per pixel and have to sum over the results to get the total fluxes per pixes 
# To align the results with the symbolic flux dict descritption we create dictionaries 
# index by the target and source pools according to the naming used in the work with Chris
# ~/bgc-md/prototypes/ModelsAsScriptsWithSpecialVarNames/models/pseudoCable
for s in ('leaf','wood','fine_root','metabolic_lit','structural_lit','cwd','fast_soil','slow_soil','passive_soil'):
    var(s)

# +
#ds=data

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
    (leaf,metabolic_lit)       : ds.fromLeaftoL.sel(litter_casa_pools=metabolic_ind) * Flux_from_Leaf, #4
    (leaf,structural_lit)      : ds.fromLeaftoL.sel(litter_casa_pools=structural_ind)* Flux_from_Leaf, #5
    (fine_root,metabolic_lit)  : ds.fromRoottoL.sel(litter_casa_pools=metabolic_ind) * Flux_from_Root, #6
    (fine_root,structural_lit) : ds.fromRoottoL.sel(litter_casa_pools=structural_ind)* Flux_from_Root, #7
    (wood,cwd)                 : ds.fromWoodtoL.sel(litter_casa_pools=cwd_ind)       * Flux_from_Wood, #8
    #
    # 9. Metabolic litter to Fast soil: fAC->A(6,3,:,:)*fAC->C(3,:,:) missing 
    (metabolic_lit,fast_soil)  : Flux_from_metabolic_lit,
    #
    # 11. Structural Litter to Fast soil: fAC->A(6,4,:,:)*fAC->C(4,:,:) missing
    # *fin->Clitter(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    (structural_lit,fast_soil) : Flux_from_structural_lit,
    #
    # 12. Structural Litter to Slow soil: fAC->A(7,4,:,:)*fAC->C(4,:,:) missing
    # *fin->Clitter(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    (structural_lit,slow_soil) : Flux_from_structural_lit,
    #
    # 14. CWD to fast soil: fAC->A(6,5,:,:)*fAC->C(5,:,:) missing
    # *fin->Clitter(:,2,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    (cwd,fast_soil)            : Flux_from_cwd,
    #
    # 15. CWD to slow soil: fAC->A(7,5,:,:)*fAC->C(5,:,:) missing
    # *fin->Clitter(:,2,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    (cwd,slow_soil)            : Flux_from_cwd,
    #
    # 17. fast soil to slow soil fAC->A(7,6,:,:)*fAC->C(6,:,:) missing
    # *fin->Csoil(:,0,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    (fast_soil,slow_soil)      : Flux_from_fast_soil,
    #
    # 18. fast soil to passive soil fAC->A(8,6,:,:)*fAC->C(6,:,:) missing
    # *fin->Csoil(:,0,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    (fast_soil,passive_soil)   : Flux_from_fast_soil,
    #
    # 20. slow soil to passive soil fAC->A(8,7,:,:)*fAC->C(7,:,:) missing
    # *fin->Csoil(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    (slow_soil,passive_soil)   : Flux_from_slow_soil,


    
 }
OutFluxes = {
    # 10. Metabolic litter to atmosphere
    # (1-fAC->A(6,3,:,:))*fAC->C(3,:,:) #missing
    # *fin->Clitter(:,0,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    metabolic_lit  : Flux_from_metabolic_lit,
    #
    # 13. structural Litter to atmosphere (1-fAC->A(6,4,:,:)-fAC->A(7,4,:,:))*fAC->C(4,:,:) missing
    # *fin->Clitter(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    structural_lit : Flux_from_structural_lit,
    #
    # 16. CWD to atmosphere (1-fAC->A(6,5,:,:)-fAC->A(7,5,:,:))*fAC->C(5,:,:) missing
    # *fin->Clitter(:,2,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)*fin->xkNlimiting(:,:,:)
    cwd            : Flux_from_cwd,
    #
    # 19. fast soil to atmosphere (1-fAC->A(7,6,:,:)-fAC->A(8,6,:,:))*fAC->C(6,:,:) missing
    # *fin->Csoil(:,0,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    fast_soil      : Flux_from_fast_soil,
    #
    # 21. slow soil to atmosphere (1-fAC->A(8,7,:,:))*fAC->C(7,:,:) missing
    # *fin->Csoil(:,1,:,:)*fin->xktemp(:,:,:)*fin->xkwater(:,:,:)
    slow_soil      : Flux_from_slow_soil
}

# -





