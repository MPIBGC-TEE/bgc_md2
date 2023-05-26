import netCDF4 as nc 
import numpy as np
import dask.array as da
from tqdm import tqdm
from functools import reduce
tsl=nc.Dataset("/home/data/yu/JULES-ES-1p0_S2_tsl.nc").variables['tsl']
#tsl_da=da.from_array(tsl)
#def reduce_is_infinite_over_first_dim(arr):
#    return  reduce(
#        lambda acc,el: np.logical_or(
#            acc,
#            np.logical_not(
#                np.isfinite(el)
#            )
#        )
#        ,
#        map(lambda i:arr[i], range(arr.shape[0]))
#    )
#
#def chunk_func(chunk):
#    final=reduce_is_infinite_over_first_dim(
#        reduce_is_infinite_over_first_dim(
#            chunk
#        )
#    )
#    #print("g final.shape",final.shape)
#    return final
#
#def j(chunk_list_list): 
#    print("len(chunk_list_list)",len(chunk_list_list))
#    def f(chunk_list):
#        #print(len(chunk_list))
#        res_f=reduce(
#            lambda acc,el:np.logical_or(acc,el),
#            map(
#                chunk_func,
#                chunk_list
#            )
#        )
#        #print("f res.shape",res_f.shape)
#        return res_f
#    
#    res_j = reduce(
#        lambda acc,el:np.logical_or(acc,el),
#        tqdm(map(f,chunk_list_list))
#    )
#    #print("j res.shape",res_j.shape)
#    return res_j
#    
#res=da.blockwise(j,'kl',tsl_da,'ijkl',adjust_chunks={"i": lambda i:1,"j": lambda j:1},dtype=np.float64).compute()
#res.shape

# NOW WE DO THE SAME WITHOUT dask
partitions=gh.partitions(start=0,stop=tsl.shape[0],nr_acc=tsl.chunking()[0])
# no list comprehension but a lazy map object to save memory
subarrs=map(lambda t: tsl[t[0]:t[1]],partitions)

[sa.shape for sa in subarrs]
