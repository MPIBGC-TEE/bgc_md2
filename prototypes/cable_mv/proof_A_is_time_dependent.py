# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# The notebook reimplements the ncl scripts to generate the matrix C in python.
# There are the following differences:
#
# - In the process of rewriting the creation of C we change the supposedly unintentionally produced `inf` 
#   values to `NaN` (or `_FillValue`) which we think were actually intended to mark the landpoint patch combinations
#   where no data is available.
#   (The infs occure because the fillValue is sometimes multiplied causing an overflow of the int32 datatype.
#   Instead we avoid processing the `_FillValue` altogeher and thus can also use it in the resulting C to mark
#   unavailable data points)
#

import xarray as xr
import numpy as np
import numpy.ma as ma
import dask
from pathlib import Path
from sympy import Symbol,var
import bgc_md2.models.cable_all.cableHelpers as cH
import netCDF4

from bgc_md2.sitespecificHelpers import get_client,getCluster
from dask.distributed import Client


# +
#cluster=getCluster()
#client=Client(cluster)

# +

#alternatively
client = get_client()
# -

# now lets repeat this exercise for the whole multifile dataset
example_run_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/'
outDir = "output/new4"
outpath = Path(example_run_dir).joinpath(outDir)
#ds=cH.cable_ds(outpath)
#dad = cH.da_dict_from_dataset(ds)
zarr_dir_path= outpath.joinpath('zarr_vars')
dad=cH.cable_dask_array_dict(zarr_dir_path)
#dad['Cplant']
dad['fromLeaftoL']
dad.keys()
dad['ksoil']


# +
def max_diff(x):
    return dask.array.max(x)-dask.array.min(x)
    
def fnum(chs):
    ch=chs[0]
    s=ch.shape
    bool_ind = np.logical_not(np.isnan(ch[0,:,:,:]))
    y=np.zeros_like(bool_ind,dtype=np.float64)
    v_inds,p_inds,l_inds = np.nonzero(bool_ind)
    for i in range(len(p_inds)):
        v_ind = v_inds[i]
        p_ind = p_inds[i]
        l_ind = l_inds[i]
        #print(ch[0,v_ind,p_ind,l_ind])
        y[v_ind,p_ind,l_ind] =   max_diff(ch[:,v_ind,p_ind,l_ind])
    #
    #print(len(chs),s,y.shape)
    return y   
res = dask.array.blockwise(fnum,'jkl',dad['fromLeaftoL'],'ijkl',dtype=np.float64)
res
# -

#resc=res.compute()
resc=res[:,0:10].compute()
ind=np.unravel_index(np.argmax(resc,axis=None),resc.shape)
ind ,resc[ind]

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
timeline= dad['fromLeaftoL'][5318:5332,ind[0],ind[1],ind[2]]
ax.plot(timeline) 

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
timeline= dad['fromRoottoL'][:,ind[0],ind[1],ind[2]]
ax.plot(timeline) 
