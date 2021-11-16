import dask.array as da
from dask.distributed import Client
from dask.distributed import LocalCluster
import numpy as np

# 1.)
# we check how map_blocks works
# assume a dask array of inputs that contains the startvectors for different chains
# for testing we just stack the same vector over and over again.
n_chains=10
l_chain=100
iv=np.array([1.0,2.0,3.0])
npar=len(iv)
IV=da.stack([iv for i in range(n_chains)]) 
# the shape of the array is now (10,3) and it is already distributed over 10 workers
# we could have enforced the latter by saying
IV.rechunk(1,3)

# now we will assume that the function working on the chunk will receive 
# one sample (start values iv) but 
# return an array of l_chain samples, where l_chain is fixed e.g.=100
# (however it comes up with them) 
# so the result_array of f_chunk is of shape (3,100)
def f_chunk(iv):
    # this would be a function like mcmc but with a fixed number of 
    # succesfull samples
    return np.ones((len(iv),l_chain))*2

# this means that the result of the blockwise operation of f_chunk on IV
# leads to (len(iv),l_chain) = (30,100)  chunks
RES = IV.map_blocks(f_chunk,chunks=(len(iv),l_chain))
RES.shape

# up to now no computation has happened
RES.compute()

# Interestingly if the inquire for the chunksize of the result we
# get (3,100) so every processor has a complete chain
RES.chunksize

# 2.)
# We now proceed to replace the f_chunk with something closer to a real mcmc
# sampler.
# A real mcmc sampler gets a startValue and produces an unpredictable number
# of accepted samples C and their associated costfunction values J.
# This # unfortunate for a parallel computation for two reasons
# 1.)  
# dask.map_blocks assignes the part of the result array to 
# be handled by one of the many parallel workers 
# BEFORE the computation starts (the chunks argument)
# 
# 2.)
# It is also unfortunate that each processor might take a different
# amount of time to reach a certain fixed number of accepted values,
# which makes the loadbalancing unpredictable.
# To avoid the situation that all processors wait for the one with the slowest
# converging chain we want to achieve an 'equal effort' parallelisation.
#
# One way to achieve this is to use arrays of fixed size 
# (But the results of the parallel computation will contain zero vectors
# in unknown positions)

# p.s. the actual implementation is still fake at the moment but
# simulates the same behavior w.r.t. the parallelisation

def f_chunk(iv):
    #print(iv.shape)
    g = np.random.default_rng()
    npar =iv.shape[1]
    res = np.zeros(shape=(npar,l_chain))
    # choose an acceptance rate randomly to model the worst case scenario.
    # that the parallel chains behave very differently
    acr = g.random()
    # if acr is smaller than nt/nss we will have 'blank' spots in the output
    # it acr is bigger we discard some succesfull samples (not too bad)
    i = 0
    for i in range(l_chain): 
        v = g.random()
        if v < acr:
            res[:,i]=np.ones(npar)*v

    return res

RES = IV.map_blocks(f_chunk,chunks=(len(iv),l_chain),meta=np.array(()))
RES.shape
RES.compute()
    

# Now we have a result that is computed very efficiently in parallel but contains zeros in
# unknown locations. We have several ways to deal with this.
# 1.) 
# Filtering 
# In our case this is not particularly expensive compared to the generation
# of the chains but it is not completely parallel operation. Fortunately 
# there is a dask function to deal with this effiently "dask.nonzero"
# but its input is a logical array to mark the values we want to keep
# we will do this in parallel

def filter_chunk(chunk):
    npar,l_chain = chunk.shape
    return np.array([(chunk[:,i] != np.zeros(npar)).all() for i in range(l_chain)]).reshape(1,l_chain)

b=RES.map_blocks(filter_chunk,chunks=(1,l_chain),dtype=np.int64)
cb=b.compute()
# b has no a value 0 or 1 for every 3 element paramter value 
# so the total shape is (10,100)
cr,cc=cb.nonzero()
#rows.compute_chunk_sizes()
#cols.compute_chunk_sizes()
#cr=rows.compute()
#cc=cols.compute()
# we can now collect all the nonzero results in a combined chain.
def extract_vec(i):
    b_row=cr[i]
    RES_row_start=b_row*npar
    RES_row_end=(b_row+1)*npar
    col = cc[i]
    return RES[RES_row_start:RES_row_end,col]

vec_list = [extract_vec(i) for i in range(len(cr))]
valid_params=da.stack(vec_list)
valid_params.compute()
