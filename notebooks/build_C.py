# ---
# jupyter:
#   jupytext:
#     formats: py:light
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

# The notebook reimplements the ncl scripts to generate the matrices A and C in python.
# There are the following differences:
# - the original ncl script `mm_mkinitialmatrices_A_C_1901_to_2004.ncl`
#   only uses data of the first day of each year to create only one A per year.
#   We produce a completely time dependent A that changes daily.
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
ds=cH.cable_ds(outpath)

dd=ds.dims
npool=sum([dd[k] for k in ['plant_casa_pools','soil_casa_pool','litter_casa_pools']])
npatch=dd['patch']
nland=dd['land']
npool,npatch,nland
dd

# +
# build C
# 1.) similar to the ncl script
iveg_m1=np.where(
    ds.iveg.data==ifv,
    18,
    ds.iveg.data
)-1
flat_m1=iveg_m1.flatten()

tau = dask.array.from_array(
    [
        [3.72837    ,10,20,0.04,0.23,0.824,0.137,5,222],
        [1.65467    ,10,30,0.04,0.23,0.824,0.137,5,222],
        [0.52343    ,10,20,0.04,0.23,0.824,0.137,5,222],
        [0.50679    ,10,10,0.04,0.23,0.824,0.137,5,222],
        [1.44000    , 2, 4,0.04,0.23,0.824,0.137,5,222],
        [0.2910     , 0.28918, 1,0.04,0.23,0.824,0.137,5,222],
        [0.21420    , 0.21404, 1,0.04,0.23,0.824,0.137,5,222],
        [0.54065    , 0.54030, 1,0.04,0.23,0.824,0.137,5,222],
        [0.28935    , 0.28935, 1,0.04,0.23,0.824,0.137,5,222],
        [       0.37, 0.37000, 1,0.04,0.23,0.824,0.137,5,222],
        [          1, 1, 1,0.04,0.23,0.824,0.137,5,222],
        [          1, 1, 1,0.04,0.23,0.824,0.137,5,222],
        [          1, 1, 1,0.04,0.23,0.824,0.137,5,222],
        [ 0.43293   , 2, 5,0.04,0.23,0.824,0.137,5,222],
        [          1, 1, 1,0.04,0.23,0.824,0.137,5,222],
        [          1, 1, 1,0.04,0.23,0.824,0.137,5,222],
        [          1, 1, 1,0.04,0.23,0.824,0.137,5,222],
        [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999]
    ]
)
    
    
xkoptlitter = dask.array.from_array(
    [
        [1,1,1,0.4,0.4,0.4,1,1,1],
        [1,1,1,0.4,0.4,0.4,1,1,1],
        [1,1,1,0.3,0.3,0.3,1,1,1],
        [1,1,1,0.6,0.6,0.6,1,1,1],
        [1,1,1,0.3,0.3,0.3,1,1,1],
        [1,1,1,0.3,0.3,0.3,1,1,1],
        [1,1,1,0.3,0.3,0.3,1,1,1],
        [1,1,1,0.2,0.2,0.2,1,1,1],
        [1,1,1,0.2,0.2,0.2,1,1,1],
        [1,1,1,0.4,0.4,0.4,1,1,1],
        [1,1,1,0.4,0.4,0.4,1,1,1],
        [1,1,1,0.4,0.4,0.4,1,1,1],
        [1,1,1,0.4,0.4,0.4,1,1,1],
        [1,1,1,2.0,2.0,2.0,1,1,1],
        [1,1,1,0.4,0.4,0.4,1,1,1],
        [1,1,1,0.4,0.4,0.4,1,1,1],
        [1,1,1,0.4,0.4,0.4,1,1,1],
        [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999]
    ]
)

xkoptsoil = dask.array.from_array(
    [
        [1,1,1,1,1,1,0.40,0.40,0.40],
        [1,1,1,1,1,1,0.40,0.40,0.40],
        [1,1,1,1,1,1,0.30,0.30,0.30],
        [1,1,1,1,1,1,0.60,0.60,0.60],
        [1,1,1,1,1,1,0.30,0.30,0.30],
        [1,1,1,1,1,1, 0.3, 0.3, 0.3],
        [1,1,1,1,1,1, 0.3, 0.3, 0.3],
        [1,1,1,1,1,1, 0.2, 0.2, 0.2],
        [1,1,1,1,1,1,0.25, 0.3, 0.3],  #crop *1.25;1.5;1.5 of original number
        [1,1,1,1,1,1,0.25,0.25,0.25],
        [1,1,1,1,1,1,   1,   1,   1],
        [1,1,1,1,1,1,0.65,0.65,0.65],
        [1,1,1,1,1,1, 0.5, 0.5, 0.5],
        [1,1,1,1,1,1,   2,   2,   2],
        [1,1,1,1,1,1, 0.5, 0.5, 0.5],
        [1,1,1,1,1,1, 1.0, 1.0, 1.0],
        [1,1,1,1,1,1, 1.0, 1.0, 1.0],
        [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999]
    ]
)

fracLigninplant = dask.array.from_array(
    [
        [0,0,0,0,0.25,0,0,0,0],
        [0,0,0,0,0.2 ,0,0,0,0],
        [0,0,0,0,0.2 ,0,0,0,0],
        [0,0,0,0,0.2 ,0,0,0,0],
        [0,0,0,0,0.2 ,0,0,0,0],
        [0,0,0,0,0.1 ,0,0,0,0],
        [0,0,0,0,0.1 ,0,0,0,0],
        [0,0,0,0,0.1 ,0,0,0,0],
        [0,0,0,0,0.1 ,0,0,0,0],
        [0,0,0,0,0.1 ,0,0,0,0],
        [0,0,0,0,0.15,0,0,0,0],
        [0,0,0,0,0.15,0,0,0,0],
        [0,0,0,0,0.15,0,0,0,0],
        [0,0,0,0,0.15,0,0,0,0],
        [0,0,0,0,0.15,0,0,0,0],
        [0,0,0,0,0.25,0,0,0,0],
        [0,0,0,0,0.1 ,0,0,0,0],
        [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999]
    ]
)

silt = np.array([0,0,0,0,0,0,0.33,0,0])
clay = np.array([0,0,0,0,0,0,0.3 ,0,0])



# to avoid confusion by the multiplication of fillvalues
# we replace  the fillvalue by np.nan which has the advantage that 
# nan*nan =nan instead of some fill value that changes into something else 
# if multiplied. This is the real problem of the ncl script
# in case iveg is undefined (fillvalue) the fillvalue gets multiplied by every where clause.
# which leads to an overflow and causes the inf values
C1= np.zeros(mat_shape,dtype='float32')
for ipool in range(npool):
    C1[ipool,ipool,:,:] = (
        np.where(
            tau[flat_m1,ipool] == -9999,
            np.nan,
            1.0/tau[flat_m1,ipool]/365.0 
        )
        *np.where(
            xkoptlitter[flat_m1,ipool] == -9999,
            np.nan,
            xkoptlitter[flat_m1,ipool]
        ) 
        *np.where(
            xkoptsoil[flat_m1,ipool] == -9999,
            np.nan,
            xkoptsoil[flat_m1,ipool]
        ) 
        *np.where(
            fracLigninplant[flat_m1,ipool] == -9999,
            np.nan,
            np.exp(-3.0*fracLigninplant[flat_m1,ipool])
        ) 
        *(1-0.75*(silt[ipool]+clay[ipool]))
    ).reshape((npatch,nland))
# we compare it to the one read from the ncl produced file  
C_inf=np.asanyarray(C)
mask=C_inf==np.inf
import numpy.ma as ma
np.max(ma.masked_array(C_inf,mask) -ma.masked_array(C1,np.isnan(C1)))
# -


# The actual aim is to deliver the expression of looked up values where iveg is defined and 
# nan everywhere else. This can be expressed much more clearly in the following code
# which is finally implemented
C2= np.zeros(mat_shape,dtype='float32')
for ipool in range(npool):
    C2[ipool,ipool,:,:] = (
        np.where(
            ds.iveg.data.flatten() == ifv,
            np.nan,
            1.0/tau[flat_m1,ipool]/365.0 
            *xkoptlitter[flat_m1,ipool]
            *xkoptsoil[flat_m1,ipool]
            *np.exp(-3.0*fracLigninplant[flat_m1,ipool])
            *(1-0.75*(silt[ipool]+clay[ipool]))
        ) 
    ).reshape((npatch,nland))
np.max(np.where(np.isnan(C1),0,C1-C2))


# +
# even simpler without artifical enlargement of tau ...
# to demonstrate that we do not need it we cut off the last line [-9990,....] from the arrays
tau_short,xkoptlitter_short,xkoptsoil_short,fracLigninplant_short =  ( 
    arr[0:17,:] for arr in [tau,xkoptlitter,xkoptsoil, fracLigninplant])

# instead of mapping the _fillvalue of iveg to 18  (and later using it to pick the 18 th  line of tau)
# we map it to 1 (the first line of tau) since we mask the result anyway
# so we can as well use one that is easy to deal with
data=ds.iveg.data
mask= data==ifv
flat_0=(np.where(mask,1,data)-1).flatten()
# a
C3= np.zeros(mat_shape,dtype='float32')

for ipool in range(npool):
    C3[ipool,ipool,:,:] = np.where(
        mask,
        np.nan,
        (
            1.0/tau_short[flat_0,ipool]/365.0 
            *xkoptlitter_short[flat_0,ipool]
            *xkoptsoil_short[flat_0,ipool]
            *np.exp(-3.0*fracLigninplant_short[flat_0,ipool])
            *(1-0.75*(silt[ipool]+clay[ipool]))
        ).reshape((npatch,nland)),
    ) 
np.max(np.where(np.isnan(C1),0,C1-C3))
# -

# actually up to now the script starts with a Zero array even if no 
# iveg value is available. This treats the offdiagonal entries differenlty than the diagonal entries
# It would be more appropriate if we either have the whole matrix or nan
C4=np.where(
    mask,
    np.nan,
    np.zeros(mat_shape,dtype='float32')
)
for ipool in range(npool):
    C4[ipool,ipool,:,:] = np.where(
        mask,
        np.nan,
        (
            1.0/tau_short[flat_0,ipool]/365.0 
            *xkoptlitter_short[flat_0,ipool]
            *xkoptsoil_short[flat_0,ipool]
            *np.exp(-3.0*fracLigninplant_short[flat_0,ipool])
            *(1-0.75*(silt[ipool]+clay[ipool]))
        ).reshape((npatch,nland)),
    ) 
# this ensures finally that we only set values different from nan where
# we have iveg data (and can actually compute the matrix)
for i in range(npool):
    for j in range(npool):
        print(np.all(np.isnan(C4[j,i,:,:])==(ds.iveg.data==ifv)))

# This is identical to the version implemented in the library
C5=cH.reconstruct_C(fns[0])
np.all(np.isnan(C5)==np.isnan(C4))
np.max(np.where(np.isnan(C5),0,C4-C5))


# +
# 6.) with function 
# This translation is probably the most easy to read 
# iveg is an array of shape (npatch,nland) 
# one can thus see it as a integer valued 
# function iveg_func(ipatch,iland)
def iveg_func(ipatch,iland):
    return ds.iveg.data[ipatch,iland]

# tau is a function of iveg and ipool
# (and therefore implicitly dependent on (ipatch,iland,ipool))
# same for the other arrays
def func_maker(arr):
    def lookup(iveg,ipool):
        return (np.nan if iveg == ifv else arr[iveg-1,ipool]) 
    return lookup


tau_func,xkoptlitter_func,xkoptsoil_func,fracLigninplant_func =  ( 
    func_maker(arr) for arr in [tau,xkoptlitter,xkoptsoil, fracLigninplant])

# Unfortunately this is very slow and therefore not used in the end 
# Array masking is orders of magnitude faster 
# although computing in vain all the values that are later masked 
#C_map=np.where(
#    mask,
#    np.nan,
#    np.zeros(mat_shape,dtype='float32')
#)
#for iland in range(nland): 
#    for ipatch in range(npatch):
#        iveg=ds.iveg.data[ipatch,iland]
#        for ipool in range(npool):
#            C_map[ipool,ipool,:,:] = 1.0/365.0/tau_func(iveg,ipool)\
#            *xkoptlitter_func(iveg,ipool)\
#            *xkoptsoil_func(iveg,ipool)\
#            *fracLigninplant_func(iveg,ipool)\
#            *(1-0.75*(silt[ipool]+clay[ipool]))
# -

