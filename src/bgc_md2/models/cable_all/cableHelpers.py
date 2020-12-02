import dask
import shutil
import xarray as xr
import zarr as zr
from pathlib import Path
import numpy as np
from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
from bgc_md2.helper import batchSlices

# fixme mm 10-14-2020
# This module contains hardcoded paths which will ultimately have to go.


def cable_ds(out_dir_path, first_yr=1901, last_yr=2004):
    """This version of the dataset leaves the fillvalues in place
    """
    out_dir_path = Path(out_dir_path)
    fns = ["out_ncar_" + str(yr) + "_ndep.nc" for yr in range(first_yr, last_yr + 1)]
    ps = [out_dir_path.joinpath(fn) for fn in fns]

    # data_vars minimal
    # is critical to avoid
    # appending a time dimension to variables
    # that do not have it in a single file DataSet
    # iveg
    ds = xr.open_mfdataset(
        paths=ps,
        concat_dim="time",
        data_vars="minimal",
        # prevent the int32 iveg part to be automatically
        # translated to float32
        mask_and_scale=False,
        # use dask delayed
        parallel=True,
    )
    dd = ds.dims
    chunk_dict = {"time": dd["time"], "patch": dd["patch"], "land": 1}
    # We can rechunk the whole dataset (or could have given a chunk argument to  xr.open_mfdataset)
    ds = ds.chunk(chunk_dict)
    return ds


def reconstruct__first_day_A(fn: str) -> np.ndarray:
    ds = xr.open_dataset(
        fn,
        mask_and_scale=False,  # avoids the unintended conversion of iveg from int32 to float32
    )
    dd = ds.dims
    npool = sum(
        [dd[k] for k in ["plant_casa_pools", "soil_casa_pool", "litter_casa_pools"]]
    )
    npatch = dd["patch"]
    nland = dd["land"]
    mat_shape = (npool, npool, npatch, nland)
    A = np.zeros(mat_shape)
    for i in range(npool):
        A[i, i, :, :] = -1

    A[3:6, 0, :, :] = ds.fromLeaftoL[1, :, :, :]
    A[3:6, 1, :, :] = ds.fromRoottoL[1, :, :, :]
    A[3:6, 2, :, :] = ds.fromWoodtoL[1, :, :, :]
    A[6:9, 3, :, :] = ds.fromMettoS[1, :, :, :]
    A[6:9, 4, :, :] = ds.fromStrtoS[1, :, :, :]
    A[6:9, 5, :, :] = ds.fromCWDtoS[1, :, :, :]
    A[7, 6, :, :] = ds.fromSOMtoSOM[1, 0, :, :]
    A[8, 6, :, :] = ds.fromSOMtoSOM[1, 1, :, :]
    A[8, 7, :, :] = ds.fromSOMtoSOM[1, 2, :, :]
    return A


def single_year_ds(fn: str) -> xr.core.dataset.Dataset:
    return xr.open_dataset(
        fn,
        # prevent the int32 iveg part to be automatically
        mask_and_scale=False,
        # translated to float32
    )


def init_I_chunk(I, NPP, fracCalloc):
    # Input fluxes
    # Note that there is an implicit ordering of state variables
    # (leaf,fine_root,wood,metabolic_lit,structural_lit,cwd,fast_soil,slow_soil,passive_soil)

    I[:, 0, :, :] = NPP[:, :, :] * fracCalloc[:, 0, :, :]
    I[:, 1, :, :] = (
        NPP[:, :, :] * fracCalloc[:, 2, :, :]
    )  # note the different orderin in cable 1, 2
    I[:, 2, :, :] = (
        NPP[:, :, :] * fracCalloc[:, 1, :, :]
    )  # note the different orderin in cable 2, 1
    return I


def init_B_chunk(
    bc,
    # iveg,
    kplant,
    fromLeaftoL,
    fromRoottoL,
    fromWoodtoL,
    fromMettoS,
    fromStrtoS,
    fromCWDtoS,
    fromSOMtoSOM,
    C,
    # A,
    xktemp,
    xkwater,
    xkNlimiting,
):
    # Note that there is an implicit ordering of state variables
    # (leaf,fine_root,wood,metabolic_lit,structural_lit,cwd,fast_soil,slow_soil,passive_soil)
    #
    # numpy arrays elements can be assigned values, but dask arrays are immutable
    # and consequently can not be assigned values
    # Therefore we create the numpy variants of the chunks
    print(bc.shape)
    bn = np.array(bc)
    # valid

    A = np.array(bc)
    npool = 9
    for i in range(npool):
        A[:, i, i, :, :] = -1

    A[:, 3:6, 0, :] = fromLeaftoL[:, :, :]
    A[:, 3:6, 1, :] = fromRoottoL[:, :, :]
    A[:, 3:6, 2, :] = fromWoodtoL[:, :, :]
    A[:, 6:9, 3, :] = fromMettoS[:, :, :]
    A[:, 6:9, 4, :] = fromStrtoS[:, :, :]
    A[:, 6:9, 5, :] = fromCWDtoS[:, :, :]
    A[:, 7, 6, :] = fromSOMtoSOM[:, 0, :]
    A[:, 8, 6, :] = fromSOMtoSOM[:, 1, :]
    A[:, 8, 7, :] = fromSOMtoSOM[:, 2, :]
    # )# numpy arrays elements can be assigned values, dask arrays can not be assigned
    # a. Leaf turnover
    bn[:, 0, 0, :, :] = -kplant[:, 0, :, :]
    #
    # b. Root turnover
    bn[:, 1, 1, :, :] = -kplant[:, 2, :, :]  # note exchange of  1 and 2
    #
    # c. Wood turnover
    bn[:, 2, 2, :, :] = -kplant[:, 1, :, :]  # note exchange of  1 and 2
    #
    # d. Leaf to Metoblic litter
    bn[:, 3, 0, :, :] = fromLeaftoL[:, 0, :, :] * kplant[:, 0, :, :]
    #
    # e. Root to Metoblic litter
    bn[:, 3, 1, :, :] = fromRoottoL[0, 0, :, :] * kplant[:, 2, :, :]  # mm
    #
    # f. Metabolic turnover
    bn[:, 3, 3, :, :] = (
        -C[3, :, :] * xktemp[:, :, :] * xkwater[:, :, :] * xkNlimiting[:, :, :]
    )
    #
    # g. Leaf to Structural litter
    bn[:, 4, 0, :, :] = fromLeaftoL[:, 1, :, :] * kplant[:, 0, :, :]
    #
    # h. Root to Structural litter
    bn[:, 4, 1, :, :] = fromRoottoL[0, 1, :, :] * kplant[:, 2, :, :]
    #
    # i. Structural turnover
    bn[:, 4, 4, :, :] = (
        -C[4, :, :] * xktemp[:, :, :] * xkwater[:, :, :] * xkNlimiting[:, :, :]
    )
    #
    # j. Wood to CWD
    bn[:, 5, 2, :, :] = fromWoodtoL[:, 2, :, :] * kplant[:, 1, :, :]
    #
    # k. CWD turnover
    bn[:, 5, 5, :, :] = (
        -C[5, :, :] * xktemp[:, :, :] * xkwater[:, :, :] * xkNlimiting[:, :, :]
    )
    # l. Metabolic litter to Fast soil
    bn[:, 6, 3, :, :] = (
        A[:, 6, 3, :, :]
        * C[3, :, :]
        * xktemp[:, :, :]
        * xkwater[:, :, :]
        * xkNlimiting[:, :, :]
    )
    #
    # m. Structural litter to Fast soil
    bn[:, 6, 4, :, :] = (
        A[:, 6, 4, :, :]
        * C[4, :, :]
        * xktemp[:, :, :]
        * xkwater[:, :, :]
        * xkNlimiting[:, :, :]
    )
    #
    # n. CWD to Fast soil
    bn[:, 6, 5, :, :] = (
        A[:, 6, 5, :, :]
        * C[5, :, :]
        * xktemp[:, :, :]
        * xkwater[:, :, :]
        * xkNlimiting[:, :, :]
    )
    #
    # o. Fast soil turnover
    bn[:, 6, 6, :, :] = -C[6, :, :] * xktemp[:, :, :] * xkwater[:, :, :]
    #
    # p. Structural litter to Slow soil
    bn[:, 7, 4, :, :] = (
        A[:, 7, 4, :, :]
        * C[4, :, :]
        * xktemp[:, :, :]
        * xkwater[:, :, :]
        * xkNlimiting[:, :, :]
    )
    #
    # q. CWD to Slow soil
    bn[:, 7, 5, :, :] = (
        A[:, 7, 5, :, :]
        * C[5, :, :]
        * xktemp[:, :, :]
        * xkwater[:, :, :]
        * xkNlimiting[:, :, :]
    )
    #
    # r. Fast soil to Slow soil
    bn[:, 7, 6, :, :] = (
        A[:, 7, 6, :, :] * C[6, :, :] * xktemp[:, :, :] * xkwater[:, :, :]
    )
    #
    # s. Slow soil turnover
    bn[:, 7, 7, :, :] = -C[7, :, :] * xktemp[:, :, :] * xkwater[:, :, :]
    #
    # t. Slow soil to Passive soil
    bn[:, 8, 7, :, :] = (
        A[:, 8, 7, :, :] * C[7, :, :] * xktemp[:, :, :] * xkwater[:, :, :]
    )
    #
    # u. Passive soil turnover
    bn[:, 8, 8, :, :] = -C[8, :, :] * xktemp[:, :, :] * xkwater[:, :, :]
    # return dask.array.from_array(bn)
    return bn

def reconstruct_B(
    ifv,
    ffv,
    iveg,
    kplant,
    fromLeaftoL,
    fromRoottoL,
    fromWoodtoL,
    fromMettoS,
    fromStrtoS,
    fromCWDtoS,
    fromSOMtoSOM,
    xktemp,
    xkwater,
    xkNlimiting
) -> dask.array.core.Array:
    npool = 9
    ntime,_,npatch,nland=kplant.shape
    C_d=dask.array.from_array(
        c_diag_from_iveg(np.array(iveg),ifv),
        chunks=(npool,npatch,1)
    )
    B_chunk = (ntime,npool,npool,npatch)
    B_shape = B_chunk+(nland,)
    B_temp=dask.array.where(
        dask.array.isnan(iveg==ifv),
        dask.array.full(B_shape,ffv,chunks=B_chunk+(1,)),
        dask.array.zeros(B_shape,chunks=B_chunk+(1,)),
    )
    return B_temp.map_blocks(
        init_B_chunk,
        kplant,
        fromLeaftoL,
        fromRoottoL,
        fromWoodtoL,
        fromMettoS,
        fromStrtoS,
        fromCWDtoS,
        fromSOMtoSOM,
        C_d,
        xktemp,
        xkwater,
        xkNlimiting,
        dtype=np.float64
    )


def reconstruct_C(fn: str) -> np.ndarray:
    ds = single_year_ds(fn)
    ifv = ds.iveg.attrs["_FillValue"]
    C_diag = reconstruct_C_diag(fn)
    arr_shape = C_diag.shape
    npool = arr_shape[0]
    mat_shape = (npool,) + (arr_shape)
    C = np.where(ds.iveg.data == ifv, np.nan, np.zeros(mat_shape, dtype="float32"))
    for ipool in range(npool):
        C[ipool, ipool, :, :] = C_diag[ipool, :, :]
    return C


def c_diag_from_iveg(iveg_data: np.ndarray, ifv: np.int32) -> np.ndarray:
    """Reconstruct the output of the original ncl script
    `bgc_md2/src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_org/mkinitialmatrix.ncl`
    The original algorithm supposedly unintentionally produces `inf` values
    where actually `NaN` are appropriate because no `iveg` value is available.
    In the original script this happens by multiplying the (huge) fillValue and thereby 
    causing an overflow.  
    This function implements the expected behaviour, where lookups with a non
    available index just return a nan for the looked up value.
    The function can also be called on a iveg_data_chunk
    """

    npatch, nland = iveg_data.shape
    npool = 9
    arr_shape = (npool, npatch, nland)
    nans = [np.nan for i in range(npool)]
    tau = np.array(
        [
            [3.72837, 10, 20, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [1.65467, 10, 30, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [0.52343, 10, 20, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [0.50679, 10, 10, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [1.44000, 2, 4, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [0.2910, 0.28918, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [0.21420, 0.21404, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [0.54065, 0.54030, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [0.28935, 0.28935, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [0.37, 0.37000, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [1, 1, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [1, 1, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [1, 1, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [0.43293, 2, 5, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [1, 1, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [1, 1, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            [1, 1, 1, 0.04, 0.23, 0.824, 0.137, 5, 222],
            nans
            # original:
            # [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999]
        ]
    )

    xkoptlitter = np.array(
        [
            [1, 1, 1, 0.4, 0.4, 0.4, 1, 1, 1],
            [1, 1, 1, 0.4, 0.4, 0.4, 1, 1, 1],
            [1, 1, 1, 0.3, 0.3, 0.3, 1, 1, 1],
            [1, 1, 1, 0.6, 0.6, 0.6, 1, 1, 1],
            [1, 1, 1, 0.3, 0.3, 0.3, 1, 1, 1],
            [1, 1, 1, 0.3, 0.3, 0.3, 1, 1, 1],
            [1, 1, 1, 0.3, 0.3, 0.3, 1, 1, 1],
            [1, 1, 1, 0.2, 0.2, 0.2, 1, 1, 1],
            [1, 1, 1, 0.2, 0.2, 0.2, 1, 1, 1],
            [1, 1, 1, 0.4, 0.4, 0.4, 1, 1, 1],
            [1, 1, 1, 0.4, 0.4, 0.4, 1, 1, 1],
            [1, 1, 1, 0.4, 0.4, 0.4, 1, 1, 1],
            [1, 1, 1, 0.4, 0.4, 0.4, 1, 1, 1],
            [1, 1, 1, 2.0, 2.0, 2.0, 1, 1, 1],
            [1, 1, 1, 0.4, 0.4, 0.4, 1, 1, 1],
            [1, 1, 1, 0.4, 0.4, 0.4, 1, 1, 1],
            [1, 1, 1, 0.4, 0.4, 0.4, 1, 1, 1],
            nans
            # original:
            # [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999]
        ]
    )

    xkoptsoil = np.array(
        [
            [1, 1, 1, 1, 1, 1, 0.40, 0.40, 0.40],
            [1, 1, 1, 1, 1, 1, 0.40, 0.40, 0.40],
            [1, 1, 1, 1, 1, 1, 0.30, 0.30, 0.30],
            [1, 1, 1, 1, 1, 1, 0.60, 0.60, 0.60],
            [1, 1, 1, 1, 1, 1, 0.30, 0.30, 0.30],
            [1, 1, 1, 1, 1, 1, 0.3, 0.3, 0.3],
            [1, 1, 1, 1, 1, 1, 0.3, 0.3, 0.3],
            [1, 1, 1, 1, 1, 1, 0.2, 0.2, 0.2],
            [1, 1, 1, 1, 1, 1, 0.25, 0.3, 0.3],  # crop *1.25;1.5;1.5 of original number
            [1, 1, 1, 1, 1, 1, 0.25, 0.25, 0.25],
            [1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 0.65, 0.65, 0.65],
            [1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5],
            [1, 1, 1, 1, 1, 1, 2, 2, 2],
            [1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5],
            [1, 1, 1, 1, 1, 1, 1.0, 1.0, 1.0],
            [1, 1, 1, 1, 1, 1, 1.0, 1.0, 1.0],
            nans
            # original:
            # [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999]
        ]
    )

    fracLigninplant = np.array(
        [
            [0, 0, 0, 0, 0.25, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.2, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.2, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.2, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.2, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.15, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.15, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.15, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.15, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.15, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.25, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.1, 0, 0, 0, 0],
            nans
            # original: [-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999]
            # although the results computed from this line are masked in the
            # computation that used the array the masking happens AFTER they
            # have been processed which leads to errors and warnings in the
            # computation of values that will be thrown away afterwards.  In
            # other words np.where(cond(A),exp1(A,..),exp2(A,..) always
            # computes exp2(A) even if it is later discarded because cond(A)
            # was true for this part of the array.
            #
            # We could actually leave out this line of the array completely and
            # map the ivge _FillValue to any one  of the regular values
            # (0,1,...16) So we would always find a tau ,.. fracLigninpant and
            # could compute the result mapping the _FillValue to the 17 th line
            # of the arrays and putting a value of -9999 there can be seen as
            # some safety strategy since the results will look wired if we
            # forget to mask them.  unfortunately this (original) approch leads
            # to the expression  np.exp(-9999) to be evaluated which leads to a
            # value too big for float32 and thus to inf and warnings We use
            # np.nan instead
        ]
    )

    silt = np.array([0, 0, 0, 0, 0, 0, 0.33, 0, 0])
    clay = np.array([0, 0, 0, 0, 0, 0, 0.3, 0, 0])

    # since we want to use iveg as an index array with values between 0 and 17
    # to look up value in tau and the other arrays we have to convert the non
    # available parts attrs[_FillValue] to an integer (in this case 17 where 17
    # refers to the line with values that will never be used) this is
    # necessarry since numpy tries to lookup the values BEFORE masking them...
    iveg_m1 = np.where(iveg_data == ifv, 18, iveg_data) - 1
    flat_m1 = iveg_m1.flatten()

    C = np.zeros(arr_shape, dtype="float32")
    for ipool in range(npool):
        C[ipool, :, :] = (
            np.where(
                iveg_data.flatten == ifv,
                np.nan,
                1.0
                / tau[flat_m1, ipool]
                / 365.0
                * xkoptlitter[flat_m1, ipool]
                * xkoptsoil[flat_m1, ipool]
                * np.exp(-3.0 * fracLigninplant[flat_m1, ipool])
                * (1 - 0.75 * (silt[ipool] + clay[ipool])),
            )
        ).reshape((npatch, nland))
    return C


def reconstruct_C_diag(fn: str) -> np.ndarray:
    ds = single_year_ds(fn)
    return c_diag_from_iveg(ds.iveg.data, ds.iveg.attrs["_FillValue"])

def reconstruct_u(
    ifv,
    ffv,
    iveg,
    NPP,
    fracCalloc
)->dask.array.core.Array:
    npool = 9
    ntime,npatch,nland=NPP.shape
    I_chunk = (ntime,npool,npatch)
    I_shape = I_chunk+(nland,)
    I_temp = dask.array.where(
        dask.array.isnan(iveg==ifv),# this works since the last two dimensions of iveg and I match
        dask.array.full(I_shape,ffv,chunks=I_chunk+(1,)),
        dask.array.zeros(I_shape,chunks=I_chunk+(1,))
    )
    return  I_temp.map_blocks(
        init_I_chunk,
        NPP,
        fracCalloc,
        dtype=np.float64
    )

def write_vars_as_zarr(ds: xr.Dataset, dir_name: str):
    dir_p = Path(dir_name)

    dir_p.mkdir(parents=True, exist_ok=True)

    for var_name, v in ds.variables.items():
        zarr_dir_path = dir_p.joinpath(var_name)
        zarr_dir_name = str(zarr_dir_path)
        if not (zarr_dir_path.exists()):
            arr = dask.array.asarray(v.data)
            batchwise_to_zarr(arr, zarr_dir_name)


def batchwise_to_zarr(arr: dask.array.core.Array, zarr_dir_name: str, rm=False):
    dir_p = Path(zarr_dir_name)
    if rm & dir_p.exists():
        print("##########################################3")
        shutil.rmtree(dir_p)

    if arr.nbytes < 8 * 1024 ** 3:
        # if the array fits into memory
        # the direct call of the to_zarr method
        # is possible (allthough it seems to imply a compute()
        # for the whole array or at least a part that is too big
        # to handle for bigger arrays
        arr.to_zarr(zarr_dir_name)
    else:
        # if the array is bigger than memory we compute explicitly
        # a part of it and write it to the zarr array.
        # This takes longer but gives us control over the
        # memory usage
        z = zr.open(zarr_dir_name, mode="w", shape=arr.shape, chunks=arr.chunksize)
        ncores = 8
        slices = batchSlices(arr.shape[-1], ncores)
        for s in slices:
            print(s)
            z[..., s] = arr[..., s].compute()


def cable_dask_array_dict(dirpath_str):
    dir_p = Path(dirpath_str)
    paths = [p for p in dir_p.iterdir() if p.is_dir()]
    var_dict = {p.name: dask.array.from_zarr(str(p)) for p in paths}
    return var_dict


def reform(combi):
    return combi.map_blocks(
        lambda timeLine: np.expand_dims(timeLine, axis=-1), new_axis=len(combi.shape)
    )


def valid_combies_parallel(nz, mat):
    s = mat.shape
    print(len(s))
    ps, lps = nz
    ps.compute_chunk_sizes()
    lps.compute_chunk_sizes()
    rps = ps.rechunk(1,)
    rlps = lps.rechunk(1,)

    # we know the number of valid entries
    l_new = len(ps)
    # an so the shape of the new array.
    s_f = s[:-1][:-1]
    s_new = s_f + (l_new,)
    c_new = s_f + (1,)
    # now we create a template
    mt = dask.array.zeros(s_new, chunks=c_new)
    print(mt.shape)
    # we now construct the squezed dask array (without computing it)
    def f(m_c, rps_c, rlps_c):
        # notes
        # - rps_c and rlps are the chunks an have actually length 1
        # - we reference variable mat which is not part of the signature
        #   but part of the closure.
        return mat[..., rps_c[0], rlps_c[0]].reshape(m_c.shape)

    return mt.map_blocks(f, rps, rlps, dtype=np.float)


def valid_combies(nz, mat):
    """Although this function does not compute any array values
    it is very slow on large arrays because the perpetual recreation of the
    new dask array entails a lot of communication"""

    s = mat.shape
    print(len(s))
    ps, lps = nz
    ps.compute_chunk_sizes()
    lps.compute_chunk_sizes()

    i = 0
    mat_new = reform(mat[..., ps[i], lps[i]])
    print(mat_new.shape)
    # we now construct the squezed dask array (without computing it)
    while i < len(ps) - 1:
        print(i)
        i += 1
        mat_new = dask.array.concatenate(
            [mat_new, reform(mat[..., ps[i], lps[i]])], axis=len(s) - 2
        )
    return mat_new


# def write_valid_combies_as_zarr(
#    ds: xr.Dataset,
#    dir_name: str
# ):
#    dir_p = Path(dir_name)
#    if dir_p.exists():
#        shutil.rmtree(dir_p)
#
#    dir_p.mkdir(parents=True, exist_ok=True)
#
#    for var_name, v in ds.variables.items():
#        zarr_dir_name = str(dir_p.joinpath(var_name))
#        dask.array.asarray(v.data).to_zarr(zarr_dir_name)
#
