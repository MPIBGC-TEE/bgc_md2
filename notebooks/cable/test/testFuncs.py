import dask.array

def D(C): # the parameter name must be equal to a function in the same module
    return 3*C

def C(A,B):
    return A+B

def A(cable_out_path):
    raw_data = dask.array.from_zarr('raw_data')
    A = dask.array.add(raw_data,raw_data)
    return A
def B(cable_out_path):
    raw_data = dask.array.from_zarr('raw_data')
    B = dask.array.add(raw_data,raw_data)
    return B
