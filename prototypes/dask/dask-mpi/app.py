import dask.array as da
from dask.distributed import Client
import numpy as np
client = Client(scheduler_file='/home/mm/scheduler.json')
client

ntime=200
npool=9
npatch=10
nland=1000
chunk_shape=(ntime,npool,npool,npatch ,1)
narr=np.ones(chunk_shape)
#print(narr)
arr1=da.from_array( narr)
#print(arr1)
garr1=da.stack([ arr1 for i in range(nland)],axis=4) 
#print(garr1)
def myfunc(chunk1,chunk2):
    return chunk1+chunk2

res=da.map_blocks(
   myfunc,
   garr1,
   garr1
   )
res

# %time res.compute()


