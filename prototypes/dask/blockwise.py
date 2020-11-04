# ---
# jupyter:
#   jupytext:
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

# +
import xarray as xr
import numpy as np
import numpy.ma as ma
import dask
from pathlib import Path
from sympy import Symbol,var
import bgc_md2.models.cable_all.cableHelpers as cH
import netCDF4

from functools import reduce

# -

from bgc_md2.sitespecificHelpers import get_client,getCluster
from dask.distributed import Client

cluster=getCluster()

client=Client(cluster)
client

# +
#cluster.close()

# +
x=dask.array.ones((10000,1000,96),chunks=(10000,1000,1))

# an example without contraction
def f(x):
    print(x.shape)
    # just do something expensive
    return reduce(lambda acc,el:acc+el,[x for i in range(5)])

res=dask.array.blockwise(f,'ijk',x,'ijk',dtype='f8').compute()
x


# -

# if the result has  one  index less than the argument (here the i)
# this implies a contracion
# then f has to accept an interable for x 
def f(xs):
    y=xs[0].sum(0)
    print(
        '##################',
        "len",len(xs),
        "xs[0].shape",
        xs[0].shape,
        'y.shape=',
        y.shape
    )
    return y
print('#######################################################################################')
res=dask.array.blockwise(f,'jk',x,'ijk',dtype='f8')
res.compute()


