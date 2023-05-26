import dask.array as da, numpy as np
from functools import reduce


# This is an application of blockwise where the result has one dimension less than the
# argument. 
# There are 2 explicit reduce calls involved
# 1.) sum over the squares of the colums inside the 2-d chunk (inner)
# 2.) sum over the chunk colums sums (outer)

# There are also implicit array reassembly activities performed by DASK in the background
# the function g below gets called 1+6 times by da.blockwise
# once for testing and once for every chunk in the first dimensions (2/1)*(3/1)=6
# 
# The first reduction happens definitely in parallel
# The second one gathers the results of the first computation
# The result array is 
# 
# Note:
# There is a dask.array function for reduction which would do the gathering recursively in 
# parallel

# Task:
# For every row of a 3 dimensional array a_{i,j}
# we want to compute the sum of the squares of the columns values.
# res_i=sum_j a_{i,j}**2
# For a numpy array we could do it like this:

def f(acc,col):
    return acc+col
line=np.array([1,2,3,4,5,6])
arr=np.array(
    [
        np.array(
            [
                line, 
                2*line,
                3*line,
            ]
        ),
        np.array(
            [
                line, 
                2*line,
                3*line,
            ]
        ),
    ]
)
print("arr.shape",arr.shape)

def sq_col_sum(arr):
    return reduce(
        f,
        # we use a generator comprehension rather than a list to save memory
        (arr[:,:,i]**2 for i in range(arr.shape[2])), 
    )

print("res serial",sq_col_sum(arr))

## Now we want to parallalize this with blockwise
## We can apply our original technique on the blocks but
#
#
x=da.from_array(arr,chunks=(1,1,2))
def g(chunks):
    print("chunks",chunks)

        
    chunk_results=[sq_col_sum(c)  for c in chunks]
    print(
        "chunk_results = ",
        list(chunk_results)
    )
    res=reduce(
            lambda acc,el: acc+el,
            chunk_results
    )
    print("res",res)
    return res

z = da.blockwise(
    g,'jk',
    x,'jkl',
    #adjust_chunks=('l',lambda l:1),
    dtype=np.float64
)
res=z.compute()
print(res)

# alternatively we can use map_block recursively
# the result blocks can be smaller than the ones the function is mapped to
# but have to have the same number of dimensions
# as the original block
result_chunk_shape=(1,1,1)
def h(chunk,block_info=None):
    return sq_col_sum(chunk).reshape(result_chunk_shape)
inter=da.map_blocks(h,x,dtype=np.float64,chunks=result_chunk_shape)

#the chunksize of the result is (1,1,1) to aggregate further we have to drop an axis 
def i(chunk):
    #print(chunk.shape)
    return chunk.sum(axis=2)
res2=da.map_blocks(i,inter,dtype=np.float64,drop_axis=2).compute()
print(res2)

#We could have done this also by
res3=da.map_blocks(h,x,dtype=np.float64,chunks=result_chunk_shape).sum(axis=2).compute()
print(res3)
