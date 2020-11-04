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

# The goals are the following:
# - check if  the matrices A and C  produced by 
#   `bgc_md2/src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_orgmm_mkinitialmatrices_A_C_1901_to_2004.ncl`
#   change from year to year.
#   The orirginal script from Chris was`/mkinitialmatrix.ncl`
#   and has been adapted to our paths and changed to loop over all the years of the simulation and write 
#   files of the form `mm_${year}_outAC.nc`  
#
# - check if the original 'outAC.nc' is reproduced by the script
#
# - reimplement it in python.  
#   In the process of rewriting also change the supposedly unintentionally produced `inf` 
#   values where actually  `NaN` are intended.
#   (This happens by multiplicating the fillValue and thereby causing an overflow)
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

# +
# We first look at the data produced by the ncl script
# ddir  = "/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4/"
fns_AC =[ddir+"mm_"+str(yr)+"_outAC.nc" for yr in range(1901,2004)]
C_fs=[np.array(netCDF4.Dataset(fn)['C']) for fn in fns_AC] 
A_fs=[np.array(netCDF4.Dataset(fn)['A']) for fn in fns_AC] 
#netCDF4.Dataset(fns_AC[0])

#for our later reconstruction we remember the fill value
C0= netCDF4.Dataset(fns_AC[0])['C']
C_fv = C0.getncattr('_FillValue')
# -


# for reference use the file sent by Chris an check if it is reproduced by the script
nds_ref=netCDF4.Dataset('../src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_org/outAC.nc')
print(nds_ref['A'].getncattr('_FillValue')==C_fv)
A_ref=np.array(nds_ref['A'])
C_ref=np.array(nds_ref['C'])
print(np.all([np.all(C_ref==M) for M in C_fs])) #all are equal to
print(np.any([np.all(A_ref==M) for M in A_fs])) #not any one is equal

# So while the C matrices produced by the script are all completely identical to the one in the "outAC.nc" file, 
# none of the A matrices produced by the script matches the A of "outAC.nc"

# +
# We now check how A_fs[0] differs from A of the reference file
#print(A_ref.size)
#print(np.sum(A_fs[0]==A_ref))
#print(A_ref.size-np.sum(A_fs[0]==A_ref))

diff_poss = np.nonzero(A_ref!=A_fs[0]) 
# since A is fourdimensional this yields four arrays with the indeces in every dimension

# Analysing which parts of the matrix A_ref and A_fs[0] differ
# we see that the differences occure only in rows {3,4} and {}  
print([set(index_arr) for index_arr in diff_poss ])

#n_diff=len(diff_poss[0])
#print(n_diff)
#for i_diff in range(n_diff):
#    i_row   = diff_poss[0][i_diff]
#    i_col   = diff_poss[1][i_diff]
#    i_patch = diff_poss[2][i_diff]
#    i_land  = diff_poss[3][i_diff]
#    print(i_row,i_col,i_patch,i_land)
#    #print(A_ref[i_row,i_col,i_patch,i_land])
#    #print(A_fs[0][i_row,i_col,i_patch,i_land])
# -


# The As produced by the script differ between them (they have yearly differences)
print(np.all([np.all(A_fs[0]==M) for M in A_fs]))
print([np.max(A_fs[0]-M) for M in A_fs])


# Now we analyse how the matrices differ among them selfes
A0=A_fs[0]
for M in A_fs[0:10]:
    diff_poss = np.nonzero(A0!=M)
    print([set(index_arr) for index_arr in diff_poss ][0:2])

# But at least the structure is the same (non zero entries in the same places)
print(np.all([(A_fs[0]!=0) == (M !=0) for M in A_fs]))


# +
fns =[ddir+"out_ncar_"+str(yr)+"_ndep.nc" for yr in (1901,2004)]

rootgrp = netCDF4.Dataset(fns[0])
np.array(rootgrp.variables['iveg'])
# -


# this is confirmed if we open it by xarray but we have to give 
# mask_and_scale=False to prevent it from silently converting iveg from int32 to float64
# thereby making it useless as an index
dsts=[ 
    xr.open_dataset(
        fn
        ,mask_and_scale=False
    ) for fn in fns
]
ds=dsts[0]
ifv= ds.iveg.attrs['_FillValue']
imv = ds.iveg.attrs['missing_value']
ds

# +
#ds.fromLeaftoL.sel({'litter_casa_pools':0,'patch':ipatch,'land':0})
ltl=ds.fromLeaftoL
fv=ds.fromLeaftoL.attrs['_FillValue']
mv=ds.fromLeaftoL.attrs['missing_value']
ltl0=ltl.sel({'litter_casa_pools':0})
ntime, npatch, nland = ltl0.shape 
bool_ind=(ltl0.data!=fv) * (ltl0.data!=mv)#neither 
times,patches,landpoints=np.nonzero(bool_ind)
#set(times),set(patches),set(landpoints)

def is_timeline_constant(timeline):
    return np.all(timeline==timeline[0])

# we now iterate over some patch landpoint combinations that are nonzero
for i in range(10):
    ip=patches[i]
    ilp=landpoints[i]
    timeline=ltl0.sel({'patch':ip,'land':ilp}).data
    print(is_timeline_constant(timeline))
    
# -

# however there are quite a few
len(times),len(patches)

ltl.data.shape

# to do this for all combinations we better do it in parallel 
#cs1 = (ntime,npatch,1)
cs1 = (ntime,3,npatch,int(nland/96))
ltl_chunked = dask.array.from_array(ltl.data,chunks=(cs1))
#ltl0.data.shape
ltl_chunked


# +
def f(chs):
    ch=chs[0]
    s=ch.shape
    bool_ind = (ch[0,:,:,:]!=fv)*(ch[0,:,:,:]!=mv)
    y=np.ones(shape=s[1:],dtype=np.int32)*ifv
    v_inds,p_inds,l_inds = np.nonzero(bool_ind)
    for i in range(len(p_inds)):
        v_ind = v_inds[i]
        p_ind = p_inds[i]
        l_ind = l_inds[i]
        print(p_ind,l_ind)
        y[v_ind,p_ind,l_ind] =   is_timeline_constant(ch[:,v_ind,p_ind,l_ind])
    #
    print(len(chs),s,y.shape)
    return y   

res = dask.array.blockwise(f,'jkl',ltl_chunked,'ijkl',dtype=np.int32)
# We now check all available entries 
#dask.array.nonzero(res[res!=fv]).compute()
np.unique(res.compute())
# -

veg_inds,patch_inds,land_inds=[arr.compute() for arr in np.nonzero(res==0)]
#i_l=land_inds.compute()
patch_inds


# +
def max_diff(timeline):
    tl_min=np.min(timeline)
    tl_max=np.max(timeline)
    res=(tl_max  - tl_min)#/tl_max
    return res

def fnum(chs):
    ch=chs[0]
    s=ch.shape
    bool_ind = (ch[0,:,:,:]!=fv)*(ch[0,:,:,:]!=mv)
    #y=np.zeros(shape=s[1:],dtype=np.float64)#*ifv
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

res = dask.array.blockwise(fnum,'jkl',ltl_chunked,'ijkl',dtype=np.float64)
#res.compute()
resc=res.compute()
#resc=res[:,:,0:56].compute()
# -

# now lets find the maximum of the non nan.part
#validpartresc[np.logical_not(np.isnan(resc))]
ind=np.unravel_index(np.argmax(resc,axis=None),resc.shape)
ind ,resc[ind]

# now plot the timeline with the gratest change
from matplotlib import pyplot as plt
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
timeline= ltl_chunked[:,ind[0],ind[1],ind[2]]
ax.plot(timeline) 


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


def fnum(chs):
    ch=chs[0]
    s=ch.shape
    bool_ind = np.logical_not(np.isnan(ch[0,:,:,:]))
    #bool_ind = (ch[0,:,:,:]!=fv)*(ch[0,:,:,:]!=mv)
    #y=np.zeros(shape=s[1:],dtype=np.float64)#*ifv
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

resc=res.compute()
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

dd=ds.dims
npool=sum([dd[k] for k in ['plant_casa_pools','soil_casa_pool','litter_casa_pools']])
npatch=dd['patch']
nland=dd['land']
npool,npatch,nland
dd

# +
mat_dims = ('poolx','pooly','patch','land')

chunk_shape =(npool,npool,npatch,1)
mat_shape = chunk_shape[0:3]+(nland,)

#Ax = xr.DataArray(
#    np.stack([ np.stack( [ np.eye(npool) for i in range(npatch) ] ,axis=2 )  for j in range(nland) ], axis=3 ),
#
#    dims=mat_dims
#)
#Ax= xr.DataArray(
#    np.zeros(mat_shape),
#    dims=mat_dims
#)
#
#Ax.load() # this is necessarry if we want to use the assignment later but usually expensive but here A is small

Ax = np.zeros(mat_shape)
for i in range(npool):
    Ax[i,i,:,:]=-1

Ax[3:6,0,:,:] = ds.fromLeaftoL [1,:,:,:]
Ax[3:6,1,:,:] = ds.fromRoottoL [1,:,:,:]
Ax[3:6,2,:,:] = ds.fromWoodtoL [1,:,:,:]
Ax[6:9,3,:,:] = ds.fromMettoS  [1,:,:,:]    
Ax[6:9,4,:,:] = ds.fromStrtoS  [1,:,:,:]    
Ax[6:9,5,:,:] = ds.fromCWDtoS  [1,:,:,:]    
Ax[7  ,6,:,:] = ds.fromSOMtoSOM[1,0,:,:]    
Ax[8  ,6,:,:] = ds.fromSOMtoSOM[1,1,:,:]    
Ax[8  ,7,:,:] = ds.fromSOMtoSOM[1,2,:,:]   
#Ax
# -

type(Ax)


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



