import dask
import shutil
import inspect
import xarray as xr
import zarr as zr
import numpy as np
from tqdm import tqdm  # Holger
from pathlib import Path
from sympy import Symbol, var
import CompartmentalSystems.helpers_reservoir as hr
from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
from bgc_md2.helper import batchSlices
from typing import Callable, Tuple, List, Union, Dict
from copy import copy, deepcopy
from functools import _lru_cache_wrapper

# fixme mm 10-14-2020
# This module contains hardcoded paths which will ultimately have to go.


def whoami(): 
    frame = inspect.currentframe()
    return inspect.getframeinfo(frame).function


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
    xkNlimiting,
) -> dask.array.core.Array:
    npool = 9
    ntime, _, npatch, nland = kplant.shape
    C_d = dask.array.from_array(
        c_diag_from_iveg(np.array(iveg), ifv), chunks=(npool, npatch, 1)
    )
    B_chunk = (ntime, npool, npool, npatch)
    B_shape = B_chunk + (nland,)
    B_temp = dask.array.where(
        dask.array.isnan(iveg == ifv),
        dask.array.full(B_shape, ffv, chunks=B_chunk + (1,)),
        dask.array.zeros(B_shape, chunks=B_chunk + (1,)),
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
        dtype=np.float64,
    )


def cache(
    zarr_dir_path: Path,
    name: str,
    arr: dask.array,
    rm: bool = False,
    batch_size: int = 1
):
    ''' 
    A wrapper to laod the desired result from a cache dir or
    compute it by calling the function with the given arguments
    '''
    zarr_dir_path.mkdir(exist_ok=True,parents=True)
    sub_dir_path = zarr_dir_path.joinpath(name)
    batchwise_to_zarr(
        arr,
        str(sub_dir_path),
        rm=rm,
        batch_size=batch_size
    )  # will do nothing if rm is False
    return dask.array.from_zarr(str(sub_dir_path))



# should be obsolete
#def reconstruct_B_u_x0_from_zarr(out: str):
#    outpath = Path(out)
#    # Allthough we will later use zarr arrays
#    # we open one of the original netcdf output files
#    # to get the meta information about fill values
#    # Since the file is small this is  fast
#    syds = single_year_ds(outpath.joinpath("out_ncar_1901_ndep.nc"))
#    tk = "_FillValue"
#    ifv = syds.iveg.attrs[tk]
#    ffv = syds.Cplant.attrs[tk]
#
#    zarr_dir_path = outpath.joinpath("zarr")
#    dad = cable_dask_array_dict(zarr_dir_path)
#
#    npool = sum([dad[name].shape[1] for name in ["Cplant", "Csoil", "Clitter"]])
#    s = dad["Cplant"].shape
#    npatch = s[2]
#    nland = s[3]
#    ntime = s[0]
#
#    for s in (
#        "leaf",
#        "wood",
#        "fine_root",
#        "metabolic_lit",
#        "structural_lit",
#        "cwd",
#        "fast_soil",
#        "slow_soil",
#        "passive_soil",
#    ):
#        var(s)
#
#    stateVariableTuple = (
#        leaf,
#        fine_root,
#        wood,
#        metabolic_lit,
#        structural_lit,
#        cwd,
#        fast_soil,
#        slow_soil,
#        passive_soil,
#    )
#    npool = len(stateVariableTuple)
#    var_arr_dict = {
#        leaf: dad["Cplant"][:, 0, :, :],
#        fine_root: dad["Cplant"][:, 2, :, :],
#        wood: dad["Cplant"][:, 1, :, :],
#        metabolic_lit: dad["Clitter"][:, 0, :, :],
#        structural_lit: dad["Clitter"][:, 1, :, :],
#        cwd: dad["Clitter"][:, 2, :, :],
#        fast_soil: dad["Csoil"][:, 0, :, :],
#        slow_soil: dad["Csoil"][:, 1, :, :],
#        passive_soil: dad["Csoil"][:, 2, :, :],
#    }
#    X = dask.array.stack([var_arr_dict[var] for var in stateVariableTuple], 1).rechunk(
#        (ntime, npool, npatch, 1)
#    )
#
#    X0 = X[0, ...]
#
#    B_res = reconstruct_B(
#        ifv,
#        ffv,
#        dad["iveg"],
#        dad["kplant"],
#        dad["fromLeaftoL"],
#        dad["fromRoottoL"],
#        dad["fromWoodtoL"],
#        dad["fromMettoS"],
#        dad["fromStrtoS"],
#        dad["fromCWDtoS"],
#        dad["fromSOMtoSOM"],
#        dad["xktemp"],
#        dad["xkwater"],
#        dad["xkNlimiting"],
#    )
#    u_res = reconstruct_u(ifv, ffv, dad["iveg"], dad["NPP"], dad["fracCalloc"])
#    return (B_res, u_res, X0)

def reconstruct_x(
    Cplant,
    Clitter,
    Csoil
) -> dask.array:
    # fixme mm:
    # obviosly the state variable tuple is fixed here
    # the reason that it is created here (and causes duplication)
    # is that the sister functions reconstruct_B and reconstruct_u
    # are not able to rearrange the variables yet.
    # and rely on the order of the state variables fixed.
    # actually the model description should get its variables from here
    npool = sum([var.shape[1] for var in [Cplant, Csoil, Clitter]])
    s = Cplant.shape
    npatch = s[2]
    nland = s[3]
    ntime = s[0]
    for s in (
        "leaf",
        "wood",
        "fine_root",
        "metabolic_lit",
        "structural_lit",
        "cwd",
        "fast_soil",
        "slow_soil",
        "passive_soil",
    ):
        var(s)

    stateVariableTuple = (
        leaf,
        fine_root,
        wood,
        metabolic_lit,
        structural_lit,
        cwd,
        fast_soil,
        slow_soil,
        passive_soil,
    )
    assert(npool == len(stateVariableTuple))
    var_arr_dict = {
        leaf: Cplant[:, 0, :, :],
        fine_root: Cplant[:, 2, :, :],
        wood: Cplant[:, 1, :, :],
        metabolic_lit: Clitter[:, 0, :, :],
        structural_lit: Clitter[:, 1, :, :],
        cwd: Clitter[:, 2, :, :],
        fast_soil: Csoil[:, 0, :, :],
        slow_soil: Csoil[:, 1, :, :],
        passive_soil: Csoil[:, 2, :, :],
    }
    X = dask.array.stack([var_arr_dict[var] for var in stateVariableTuple], 1).rechunk(
        (ntime, npool, npatch, 1)
    )
    return X


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


def reconstruct_u(ifv, ffv, iveg, NPP, fracCalloc) -> dask.array.core.Array:
    npool = 9
    ntime, npatch, nland = NPP.shape
    I_chunk = (ntime, npool, npatch)
    I_shape = I_chunk + (nland,)
    I_temp = dask.array.where(
        dask.array.isnan(
            iveg == ifv
        ),  # this works since the last two dimensions of iveg and I match
        dask.array.full(I_shape, ffv, chunks=I_chunk + (1,)),
        dask.array.zeros(I_shape, chunks=I_chunk + (1,)),
    )
    return I_temp.map_blocks(init_I_chunk, NPP, fracCalloc, dtype=np.float64)


#def write_vars_as_zarr(ds: xr.Dataset, dir_name: str, batch_size: int = 16):
#    dir_p = Path(dir_name)
#
#    dir_p.mkdir(parents=True, exist_ok=True)
#
#    for var_name, v in ds.variables.items():
#        zarr_dir_path = dir_p.joinpath(var_name)
#        zarr_dir_name = str(zarr_dir_path)
#        if not (zarr_dir_path.exists()):
#            arr = dask.array.asarray(v.data)
#            batchwise_to_zarr(arr, zarr_dir_name, batch_size)


def batchwise_to_zarr(
    arr: dask.array.core.Array,
    zarr_dir_name: str,
    batch_size: int = 1
):
    in_p = Path(zarr_dir_name + '.in_progress')
    dir_p = Path(zarr_dir_name)

    if dir_p.exists():
        raise Exception("zarr cache directory: " + str(dir_p) + "already exists.") 

    if in_p.exists():  
        #remove leftovers of previous attempts
        print("removing unfinished directory" + str(in_p))
        shutil.rmtree(in_p)

    # We compute explicitly a part of the array  and write it to the zarr
    # array.  This takes longer but gives us control over the memory usage
    z = zr.open(str(in_p), mode="w", shape=arr.shape, chunks=arr.chunksize)
    slices = batchSlices(arr.shape[-1], batch_size)
    print("result shape:", arr.shape)
    
    for s in tqdm(slices):  # Holger
             z[..., s] = arr[..., s].compute()
    # after we have written everything we rename the directory to its p
    # proper name
    shutil.move(in_p,dir_p)

def cable_dask_array_dict(dirpath_str):
    dir_p = Path(dirpath_str)
    paths = [p for p in dir_p.iterdir() if p.is_dir()]
    var_dict = {p.name: dask.array.from_zarr(str(p)) for p in paths}
    return var_dict


#def load_or_make_cable_dask_array_dict(
#    cable_run_output_dir: Union[str, Path],
#    zarr_dir_path: Union[str, Path],
#    names: List[str] = [
#		"iveg",
#		"kplant",
#		"fromLeaftoL",
#		"fromRoottoL",
#		"fromWoodtoL",
#		"fromMettoS",
#		"fromStrtoS",
#		"fromCWDtoS",
#		"fromSOMtoSOM",
#		"xktemp",
#		"xkwater",
#		"xkNlimiting",
#		"iveg",
#		"NPP",
#		"fracCalloc",
#		"Cplant",
#		"Clitter",
#		"Csoil"
#
#    rm: bool = False,
#    batch_size: int = 16,
#):
#    ds = cable_ds(cable_run_output_dir)
#
#    if not zarr_dir_path.exists():
#        zarr_dir_path.mkdir()
#
#    for var_name, v in ds.variables.items():
#        sub_dir_path = zarr_dir_path.joinpath(var_name)
#        zarr_dir_name = str(sub_dir_path)
#        arr = dask.array.asarray(v.data)
#        batchwise_to_zarr(arr, zarr_dir_name, rm=rm, batch_size=batch_size)
#
#    return cable_dask_array_dict(zarr_dir_path)
#

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
    if l_new == 0:
        raise Exception("no valid combies found.")
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


def pool_age_density_val(
    x0_c: np.ndarray,
    times: np.ndarray,
    age_bin_indices: np.ndarray,  # integer
    B_c: np.ndarray,
    U_c: np.ndarray,
    start_age_densities_of_x0_and_a: Callable[[np.ndarray, float], np.ndarray],
):
    # Note that the chunksize is supposed to be one
    # We have only one startvector
    npools, _ = x0_c.shape
    it_max = len(times)
    print(it_max)
    x0 = x0_c[:, 0]  # .reshape(9)

    Bs = [B_c[i, :, :, 0] + np.eye(9) for i in range(it_max)]
    Us = [U_c[i, :, 0] for i in range(it_max)]
    dmr = DMR.from_Bs_and_net_Us(x0, times, Bs, Us)

    def start_age_densities(a):
        return start_age_densities_of_x0_and_a(x0, a)

    start_age_densities_of_bin_index = hr.pool_wise_bin_densities_from_smooth_densities_and_index(
        start_age_densities, npools, dmr.dt
    )

    #
    res = np.zeros((len(age_bin_indices), len(times), npools, 1))

    p_dmr = dmr.pool_age_densities_func(start_age_densities_of_bin_index)
    pad = p_dmr(age_bin_indices)
    res[:, 1:, :, 0] = pad
    # sol=dmr.solve()[:-1,:].reshape(U_c.shape) #remove the last value
    # return sol
    return res


def valid_trajectory(x0_c, times, B_c, U_c):
    # Note that the chunksize is supposed to be one
    # We have only one startvector
    print("###########")
    print(x0_c.shape, type(x0_c))
    print(times.shape)
    # print(ages.shape)
    print(B_c.shape)
    print(U_c.shape)
    it_max = len(times)
    print(it_max)
    x0 = x0_c[:, 0]  # .reshape(9)
    Bs = [B_c[i, :, :, 0] + np.eye(9) for i in range(it_max)]
    Us = [U_c[i, :, 0] for i in range(it_max)]
    #    dmr = DMR.from_Bs_and_net_Us(x0,times,Bs,Us)
    #    sol=dmr.solve()[:-1,:].reshape(U_c.shape) #remove the last value
    dmr = DMR.from_Bs_and_net_Us(x0, times, Bs[:-1], Us[:-1])  # Holger
    sol = dmr.solve().reshape(U_c.shape)  # Holger
    return sol


def aggregate_xs(xs, nr_days):
    # Note that the chunksize is supposed to be one
    # We have only one startvector
    print("###########")
    print(x0_c.shape, type(x0_c))
    print(times.shape)
    # print(ages.shape)
    print(B_c.shape)
    print(U_c.shape)
    it_max = len(times)
    print(it_max)
    x0 = x0_c[:, 0]  # .reshape(9)
    Bs = [B_c[i, :, :, 0] + np.eye(9) for i in range(it_max)]
    Us = [U_c[i, :, 0] for i in range(it_max)]
    #    dmr = DMR.from_Bs_and_net_Us(x0,times,Bs,Us)
    #    sol=dmr.solve()[:-1,:].reshape(U_c.shape) #remove the last value
    dmr = DMR.from_Bs_and_net_Us(x0, times, Bs[:-1], Us[:-1])  # Holger
    sol = dmr.solve().reshape(U_c.shape)  # Holger
    return sol


# fixme mm 12-18-2020
# this function should not be necessarry any more
# fixme 6-30-2021
# deprecated dont use
#def load_or_make_B_u_x0_from_zarr(
#    out_path, zarr_sub_dir_name, rm=False, batch_size: int = 16
#):
#    zarr_dir_path = out_path.joinpath(zarr_sub_dir_name)
#    print(798, zarr_dir_path)
#    names = ("B", "u", "x0")
#    sub_dir_paths = [
#        #        zarr_dir_path.joinpath(name) for name in val_names
#        zarr_dir_path.joinpath(name)
#        for name in names  # Holger
#    ]
#    print(804, sub_dir_paths)
#    if all((p.exists() for p in sub_dir_paths)):
#        B, u, x0 = (dask.array.from_zarr(str(p)) for p in sub_dir_paths)
#        print(810)
#        return B, u, x0
#
#    else:
#        print(813)
#        #        B,u,x0 = reconstruct_B_u_x0_from_zarr(out_dir)
#        B, u, x0 = reconstruct_B_u_x0_from_zarr(out_path)  # Holger
#
#        #        for tup in zip((B,u_,x0),names):
#        for tup in zip((B, u, x0), names):  # Holger
#            arr, name = tup
#            sl = slice(0, None)  # Holger
#            #            suffix = "_Holger" # Holger
#            #            suffix = '_slice_' + str(sl.start) + '_' + str(sl.stop) # Holger
#            print(824)
#            suffix = ""  # Holger
#            batchwise_to_zarr(
#                arr[..., sl],
#                str(zarr_dir_path.joinpath(name + suffix)),
#                rm=rm,
#                batch_size=batch_size,
#            )
#        # redefine as based on the just written zarr array
#        # which is faster (due to better chunking)
#        # than the delayed object
#        # from which the zarr arrays are created
#        return load_or_make_B_u_x0_from_zarr(out_path, zarr_sub_dir_name)
#
## fixme mm
## this function is replaced by a
## deprecated
#def load_or_make_valid_B_u_x0(
#    out_path: Path,
#    zarr_sub_dir_name: str,
#    names: List[str] = ["B_val", "u_val", "x0_val"],
#    rm: bool = False,
#    batch_size: int = 16,
#) -> Tuple[dask.array.core.Array]:
#
#    zarr_dir_path = out_path.joinpath(zarr_sub_dir_name)
#    sub_dir_paths = [zarr_dir_path.joinpath(name) for name in names]
#
#    if all((p.exists() for p in sub_dir_paths)):
#        B_val, u_val, x0_val = (dask.array.from_zarr(str(p)) for p in sub_dir_paths)
#        print("loaded")
#        return B_val, u_val, x0_val
#
#    else:
#        print(867)
#        syds = single_year_ds(out_path.joinpath("out_ncar_1901_ndep.nc"))
#        tk = "_FillValue"
#        ifv = syds.iveg.attrs[tk]
#        ffv = syds.Cplant.attrs[tk]
#        dad = cable_dask_array_dict(zarr_dir_path)
#
#        iveg = dad["iveg"]
#        B, u, x0 = load_or_make_B_u_x0_from_zarr(out_path, zarr_sub_dir_name)
#        cond_1 = (iveg != ifv).compute()
#        cond_2 = (dad["Csoil"][0, 0, :, :] != 0).compute()
#        nz = dask.array.nonzero(cond_1 * cond_2)
#
#        B_val, u_val, x0_val = (valid_combies_parallel(nz, arr) for arr in (B, u, x0))
#        print(887)
#        for tup in zip([B_val, u_val, x0_val], sub_dir_paths):
#            arr, p = tup
#            batchwise_to_zarr(arr, str(p), rm=rm)
#        print(895)
#        return load_or_make_valid_B_u_x0(
#            out_path, zarr_sub_dir_name, names=names, rm=rm, batch_size=batch_size
#        )
#
# deprecated
#def load_or_make_valid_B_u_x0_slice(
#    out_path: Path, zarr_sub_dir_name: str, sl: slice, rm=False, batch_size: int = 16
#) -> Tuple[dask.array.core.Array]:
#
#    zarr_dir_path = out_path.joinpath(zarr_sub_dir_name)
#    suffix = "_slice_" + str(sl.start) + "_" + str(sl.stop)
#    names = ["B_val", "u_val", "x0_val"]
#    sub_dir_paths = [zarr_dir_path.joinpath(name + suffix) for name in names]
#
#    def load():
#        B_val, u_val, x0_val = (dask.array.from_zarr(str(p)) for p in sub_dir_paths)
#        print(922, "all exist")
#        return B_val, u_val, x0_val
#
#    if all((p.exists() for p in sub_dir_paths)):
#        return load()
#
#    else:
#        syds = single_year_ds(out_path.joinpath("out_ncar_1901_ndep.nc"))
#        tk = "_FillValue"
#        ifv = syds.iveg.attrs[tk]
#        ffv = syds.Cplant.attrs[tk]
#        dad = cable_dask_array_dict(zarr_dir_path)
#
#        iveg = dad["iveg"]
#        print(934)
#        B, u, x0 = load_or_make_B_u_x0_from_zarr(out_path, zarr_sub_dir_name)
#        cond_1 = (iveg != ifv).compute()
#        cond_2 = (dad["Csoil"][0, 0, :, :] != 0).compute()
#        nz = dask.array.nonzero(cond_1 * cond_2)
#
#        print(943)
#        B_val, u_val, x0_val = (valid_combies_parallel(nz, arr) for arr in (B, u, x0))
#
#        print(949)
#        for tup in zip([B_val, u_val, x0_val], sub_dir_paths):
#            arr, p = tup
#            print(952, sl, p)
#            batchwise_to_zarr(arr[..., sl], str(p), rm=rm)
#
#        print(959)
#        return load()
#
#

# fixme mm 6-30-2021
# never copy this function rather use parts
def load_or_make_B_u_x(
        dad: Dict,
        ifv,
        ffv,
        zarr_dir_path: Path,
        names: Tuple[str], 
        rm : bool = False,
        batch_size: int =1
):
    B = cache(
        zarr_dir_path=zarr_dir_path,
        name=names[0],
        arr=reconstruct_B(
            ifv,
            ffv,
            dad["iveg"],
            dad["kplant"],
            dad["fromLeaftoL"],
            dad["fromRoottoL"],
            dad["fromWoodtoL"],
            dad["fromMettoS"],
            dad["fromStrtoS"],
            dad["fromCWDtoS"],
            dad["fromSOMtoSOM"],
            dad["xktemp"],
            dad["xkwater"],
            dad["xkNlimiting"],
        ),
        rm=rm,
        batch_size=batch_size
    )
    u = cache(
        zarr_dir_path=zarr_dir_path,
        name=names[1],
        arr=reconstruct_u( ifv, ffv, dad["iveg"], dad["NPP"], dad["fracCalloc"]),
        rm=rm,
        batch_size=batch_size
    )
    
    x = cache(
        zarr_dir_path=zarr_dir_path,
        name=names[2],
        arr=reconstruct_x(
            dad["Cplant"],
            dad["Clitter"],
            dad["Csoil"]
        ),
        rm=rm,
        batch_size=batch_size
    )
    return (B,u,x)

def load_or_make_valid_B_u_x(
    B,
    u,
    x,
    nz,
    vcsp: Path,
    names: Tuple[str], 
    sl: slice = slice(None,None,None),
    rm: bool = False,
    batch_size: int =1
):    
    B_val = cache(
        zarr_dir_path=vcsp,
        name=names[0],
        arr=valid_combies_parallel(nz, B)[...,sl], 
        rm=rm,
        batch_size=batch_size
    )
    u_val = cache(
        zarr_dir_path=vcsp,
        name=names[1],
        arr=valid_combies_parallel(nz, u)[...,sl],
        rm=rm,
        batch_size=batch_size
    )
    x_val = cache(
        zarr_dir_path=vcsp,
        name=names[2],
        arr=valid_combies_parallel(nz, x)[...,sl],
        rm=rm,
        batch_size=batch_size
    )
    return (B_val, u_val, x_val)


def val_or_default(
    d: dict,
    key,
    default
):
    return d[key] if key in d.keys() else default


def cacheWrapper(
    cachable_func,
    **kwargs
):
    zarr_cache_path = kwargs['zarr_cache_path']
    name = cachable_func.__name__
    sub_dir_path = zarr_cache_path.joinpath(name)

    rm = val_or_default(kwargs,'rm',False) 
    rec_rm = val_or_default(kwargs,'rec_rm',False) 
    rm_pass = rm or rec_rm
    rec_kw_args={k:v for k,v in kwargs.items() if k != 'rm'}
    
    def create_cache():
        batch_size_pass = val_or_default(kwargs,'batch_size',1) 
        res = cachable_func(**rec_kw_args)
        # remove the rm flag since from the passed on kwargs
        # since we only want it to take effect on the highest level
        # which is here
        if isinstance(res,dask.array.core.Array):
            batchwise_to_zarr(
                res,
                str(sub_dir_path),
                batch_size=batch_size_pass
            )
        return res

    if sub_dir_path.exists():
        if rm_pass:
            print("##########################################")
            print("removing " + str(sub_dir_path))
            shutil.rmtree(sub_dir_path)
            return create_cache()
        else:
            print("##########################################")
            print("using existing zarr array "+ str(sub_dir_path))
            return dask.array.from_zarr(str(sub_dir_path))
        
    else:
        return create_cache()


def get_integer_fill_value(
        cable_data_set: xr.core.dataset.Dataset
) -> int:
    return cable_data_set['iveg'].attrs['_FillValue']


def get_float_fill_value(
        cable_data_set: xr.core.dataset.Dataset
) -> float:
    return cable_data_set['Cplant'].attrs['_FillValue']

