import dask.array

# Note:
# There are two kind of functions 
# 1.) Functions that directly deal with the cable output
#     and produce dask arrays 
#     They accept the arguments
#     - cable_out_path
#     - land_point_slice


def A(
        cable_out_path,
        landpoint_slice: slice = slice(None ,None, None)
):
    url = str(cable_out_path.joinpath('raw_data'))
    raw_data = dask.array.from_zarr(url)[..., landpoint_slice, :]
    A = dask.array.add(raw_data, raw_data)
    return A


def B(
        cable_out_path,
        landpoint_slice: slice = slice(None, None, None)
):
    url = str(cable_out_path.joinpath('raw_data'))
    raw_data = dask.array.from_zarr(url )[...,landpoint_slice,:]
    B = dask.array.add(raw_data, raw_data)
    return B

def F(cable_out_path):
    return 99999

# 2 .) Functions that depend on other variables produced 
#      by other functions whose size will also depend on the
#      size of the source variables.
#      They follow the convention that parameter names in the signature must
#      be equal to function names in this module. 
#      E.g. for C(A,B) there must be
#      two functions A, B.  
#      This is how the recursive link is established.
def C(A,B):
    return A+B

def D(C): 
    return 5*C

