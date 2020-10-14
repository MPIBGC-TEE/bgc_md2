import dask
import xarray as xr
from pathlib import Path
import numpy as np
# fixme mm 10-14-2020
# This module contains hardcoded paths which will ultimately have to go.

def cable_ds(out_dir_path,first_yr = 1901,last_yr = 2004):
    cableDataDir = '/home/data/cable-data/example_runs'
    fns = ["out_ncar_" + str(yr) + "_ndep.nc" for yr in range(first_yr, last_yr + 1)]
    ps = [out_dir_path.joinpath(fn) for fn in fns]

    # data_vars minimal 
    # is critical to avoid 
    # appending a time dimension to variables 
    # that do not have it in a single file DataSet 
    # iveg 
    ds = xr.open_mfdataset(
        paths=ps, 
        concat_dim='time',
        data_vars='minimal', 
        mask_and_scale=False, # prevent the int32 iveg part to be automatically 
        # translated to float32
        parallel=True # use dask delayed a
    )
    dd=ds.dims
    chunk_dict = {'time':dd['time'],'patch':dd['patch'],'land':1}
    # We can rechunk the whole dataset (or could have given a chunk argument to  xr.open_mfdataset)
    ds=ds.chunk(chunk_dict)
    return ds

#def cable_dask_array_dict():
#    dirpath_str='/home/data/cable-data/zarr/parallel_1901_2004_with_spinup.zarr'
#    paths=[p for p in Path(dirpath_str).iterdir() if p.is_dir()]
#    return { p.name:dask.array.from_zarr(str(p)) for p in paths}

def reconstruct_A(fn:str) -> np.ndarray: 
    ds=xr.open_dataset(
        fn,
        mask_and_scale=False #avoids the unintended conversion of iveg from int32 to float32
    )
    dd = ds.dims
    npool=sum([dd[k] for k in ['plant_casa_pools','soil_casa_pool','litter_casa_pools']])
    npatch = dd['patch']
    nland = dd['land']
    mat_shape = (npool, npool, npatch, nland)
    A = np.zeros(mat_shape)
    for i in range(npool):
        A[i, i, :, :] = -1

    A[3:6,0,:,:] = ds.fromLeaftoL [1,:,:,:]
    A[3:6,1,:,:] = ds.fromRoottoL [1,:,:,:]
    A[3:6,2,:,:] = ds.fromWoodtoL [1,:,:,:]
    A[6:9,3,:,:] = ds.fromMettoS  [1,:,:,:]    
    A[6:9,4,:,:] = ds.fromStrtoS  [1,:,:,:]    
    A[6:9,5,:,:] = ds.fromCWDtoS  [1,:,:,:]    
    A[7  ,6,:,:] = ds.fromSOMtoSOM[1,0,:,:]    
    A[8  ,6,:,:] = ds.fromSOMtoSOM[1,1,:,:]    
    A[8  ,7,:,:] = ds.fromSOMtoSOM[1,2,:,:]
    return A

def single_year_ds(fn:str)->xr.core.dataset.Dataset:
    return xr.open_dataset(
        fn, 
        mask_and_scale=False, # prevent the int32 iveg part to be automatically 
        # translated to float32
    )


def reconstruct_C(fn:str) -> np.ndarray: 
    ds= single_year_ds(fn)
    ifv = ds.iveg.attrs['_FillValue']
    C_diag = reconstruct_C_diag(fn)
    arr_shape = C_diag.shape
    npool = arr_shape[0]
    mat_shape = (npool,)+(arr_shape)
    C = np.where(
            ds.iveg.data ==ifv,
            np.nan,
            np.zeros(mat_shape,dtype='float32')
    )
    for ipool in range(npool):
        C[ipool,ipool,:,:] = C_diag[ipool,:,:]
    return C

def reconstruct_C_diag(fn:str) -> np.ndarray: 
    '''Reconstruct the output of the original ncl script
    `bgc_md2/src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_org/mkinitialmatrix.ncl`
    The original algorithm supposedly unintentionally produces `inf` values
    where actually `NaN` are appropriate because no `iveg` value is available.
    This happens by multiplicating the fillValue and causing an overflow.  This
    function implements the expected behaviour, where lookups with a non
    availabel index just returns a nan for the looked up value.
    '''
    ds= single_year_ds(fn)
    dd=ds.dims
    npool=sum([dd[k] for k in ['plant_casa_pools','soil_casa_pool','litter_casa_pools']])
    npatch=dd['patch']
    nland=dd['land']
    arr_shape =(npool,npatch,nland)

    nans=   [np.nan for i in range(npool)]
    tau = np.array(
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
            nans 
            # original:
            # [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999] 
        ]
    )

    xkoptlitter = np.array(
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
            nans
            # original:
            # [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999] 
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
            nans 
            # original:
            # [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999] 
        ]
    )

    fracLigninplant = np.array(
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
            nans
            # original:
            # [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999] 
            # although the results computed from this line are masked
            # in the computation that used the array 
            # the masking happens AFTER they have been processed
            # which leads to errors and warnings in the computation of
            # values that will be thrown away afterwards.
            # In other words np.where(cond(A),exp1(A,..),exp2(A,..) always computes exp2(A) 
            # even if it is later discarded because cond(A) was true for this part of the array.
            # 
            # We could actually leave out this line of the array completely and map the ivge _FillValue to 
            # one any of the regular values (0,1,...16) 
            # So we would always find a tau ,.. fracLigninpant and could compute the result 
            # mapping the _FillValue to the 17 th line of the arrays and
            # putting a value of -9999 there can be seen as some safety strategy since the results will
            # look wired if we forget to mask them.
            # unfortunately this (original) approch leads
            # to the expression  np.exp(-9999) to be evaluated which leads to a value too big for float32
            # and thus to inf and warnings
            # We use np.nan instead
        ]
    )

    silt = np.array([0,0,0,0,0,0,0.33,0,0])
    clay = np.array([0,0,0,0,0,0,0.3 ,0,0])

    # since we want to use iveg as an index array with values between 0 and 17 
    # to look up value in tau and the other arrays we have to convert the non available parts attrs[_FillValue]
    # to an integer (in this case 17 where 17 refers to the line with values that will never be used)
    # this is necessarry since numpy tries to lookup the values BEFORE masking them...
    ifv = ds.iveg.attrs['_FillValue']
    iveg_m1=np.where(
        ds.iveg.data==ifv,
        18,
        ds.iveg.data
    )-1
    flat_m1=iveg_m1.flatten() 

    C = np.zeros(arr_shape,dtype='float32')
    for ipool in range(npool):
        C[ipool,:,:] = (
            np.where(
                ds.iveg.data.flatten == ifv,
                np.nan,
                1.0/tau[flat_m1, ipool]/365.0 
                *xkoptlitter[flat_m1, ipool]
                *xkoptsoil[flat_m1, ipool]
                *np.exp(-3.0*fracLigninplant[flat_m1, ipool ])
                *(1-0.75*(silt[ipool]+clay[ipool]))
            )
        ).reshape((npatch,nland))
    return C
