import xarray as xr
import numpy as np
import dask.array as da
from bgc_md2.sitespecificHelpers import get_client
from bgc_md2.models.cable_all.cableHelpers import cable_ds
from dask.distributed import Client
#client = Client(scheduler_file='/home/mm/scheduler.json')
client = get_client()
print(client)
ds = cable_ds()
dd =ds.dims

npool = dd['plant_casa_pools']+dd['litter_casa_pools']+dd['soil_casa_pool']
npatch = dd['patch']
nland = dd['land']
ntime = dd['time']
# npool = 9
# npatch = 10
# nland = 500
# ntime = 300
A_chunk_shape = (npool,npool,npatch ,1)
A_d = da.stack(
    [
        np.stack(
            [np.eye(npool) for i in range(npatch)],

            axis=2
        )
        for j in range(nland)
    ],
    axis=3
).rechunk(A_chunk_shape)
fLtoL_chunk_shape = (3, npatch, 1)
def append_land_dim(tup):
    return tup[:-1]+(nland,)

fLtoL = da.from_array(
    np.ones(append_land_dim(fLtoL_chunk_shape)),
    chunks=fLtoL_chunk_shape
)

print(A_d)
#narr=np.ones(chunk_shape)
#arr1=da.from_array( narr)
#garr1=da.stack([ arr1 for i in range(nland)]) 


def initA(
    # note that the chunks have different sizes
    # Ac.shape =(9,9,10)
    # LeaftoL.shape =(3,10)
    Ac,
    LeaftoL,
    RoottoL,
    WoodtoL,
    MettoS,
    StrtoS,
    CWDtoS
):
    Ac[3:6, 0, :] = LeaftoL
    #Ac[3:6, 1, :] = RoottoL
    #Ac[3:6, 2, :] = WoodtoL
    #Ac[6:9, 3, :] = MettoS
    #Ac[6:9, 4, :] = StrtoS
    #Ac[6:9, 5, :] = CWDtoS

    #Ac[7  ,6,:,:] = fromSOMtoSOM[1,0,:,:]    
    #Ac[8  ,6,:,:] = fromSOMtoSOM[1,1,:,:]    
    #Ac[8  ,7,:,:] = fromSOMtoSOM[1,2,:,:]   
    return Ac
#
ds_first=ds.isel(time=0)
res = da.map_blocks(
    initA,
    A_d,
    ds_first.fromLeaftoL,
    ds_first.fromRoottoL,
    ds_first.fromWoodtoL,
    ds_first.fromMettoS,
    ds_first.fromStrtoS,
    ds_first.fromCWDtoS,
    dtype='f8'
)
print(res)

