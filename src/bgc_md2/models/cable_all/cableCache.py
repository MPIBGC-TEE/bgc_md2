import dask.array
import numpy as np
from functools import reduce
from . import cableHelpers as cH
from testinfrastructure.helpers import pp


# iveg is not time dependent and so the function
# will not accept the time_slice argument
def iveg(cable_data_set, landpoint_slice=slice(None, None, None), **kwargs):
    arr = dask.array.asarray(cable_data_set["iveg"].data)
    return arr[..., landpoint_slice]


def time(cable_data_set, time_slice=slice(None, None, None), **kwargs):
    arr = dask.array.asarray(cable_data_set["time"].data)
    print(arr.shape)
    sl = arr[time_slice]
    print(sl.shape)
    return sl


# most of the other cable outputs are very similar
# to each other. So we create the functions that
# return them automatically from a template
def FuncMaker(name):
    def func(
        cable_data_set,
        time_slice=slice(None, None, None),
        landpoint_slice=slice(None, None, None),
        **kwargs
    ):
        arr = dask.array.asarray(cable_data_set[name].data)
        sl = arr[time_slice, ..., landpoint_slice]
        return sl

    func.__name__ = name
    return func


Clitter = FuncMaker("Clitter")
Cplant = FuncMaker("Cplant")
Csoil = FuncMaker("Csoil")
NPP = FuncMaker("NPP")
fracCalloc = FuncMaker("fracCalloc")
fromCWDtoS = FuncMaker("fromCWDtoS")
fromLeaftoL = FuncMaker("fromLeaftoL")
fromMettoS = FuncMaker("fromMettoS")
fromRoottoL = FuncMaker("fromRoottoL")
fromSOMtoSOM = FuncMaker("fromSOMtoSOM")
fromStrtoS = FuncMaker("fromStrtoS")
fromWoodtoL = FuncMaker("fromWoodtoL")
kplant = FuncMaker("kplant")
xkNlimiting = FuncMaker("xkNlimiting")
xktemp = FuncMaker("xktemp")
xkwater = FuncMaker("xkwater")


def FuncMaker0(name):
    def func(
        cable_data_set,
        time_slice=slice(0, None, None),
        landpoint_slice=slice(None, None, None),
        **kwargs
    ):
        start_time = time_slice.start
        pp("name", locals())
        arr = dask.array.asarray(cable_data_set[name].data)
        sl = arr[start_time, ..., landpoint_slice]
        return sl

    # important since cache_wrapper uses the function name
    func.__name__ = name + "0"
    return func


Clitter0 = FuncMaker0("Clitter")
Cplant0 = FuncMaker0("Cplant")
Csoil0 = FuncMaker0("Csoil")


def B_org(**kwargs) -> dask.array.core.Array:
    """The B matrix in the original cable shape"""
    cable_data_set = kwargs["cable_data_set"]

    res_ifv = cH.get_integer_fill_value(cable_data_set)
    res_ffv = cH.get_float_fill_value(cable_data_set)

    return cH.reconstruct_B(
        res_ifv,
        res_ffv,
        cH.cacheWrapper(iveg, **kwargs),
        cH.cacheWrapper(kplant, **kwargs),
        cH.cacheWrapper(fromLeaftoL, **kwargs),
        cH.cacheWrapper(fromRoottoL, **kwargs),
        cH.cacheWrapper(fromWoodtoL, **kwargs),
        cH.cacheWrapper(fromMettoS, **kwargs),
        cH.cacheWrapper(fromStrtoS, **kwargs),
        cH.cacheWrapper(fromCWDtoS, **kwargs),
        cH.cacheWrapper(fromSOMtoSOM, **kwargs),
        cH.cacheWrapper(xktemp, **kwargs),
        cH.cacheWrapper(xkwater, **kwargs),
        cH.cacheWrapper(xkNlimiting, **kwargs),
    )


def B_org_iveg(**kwargs):
    """The B matrix in the original cable shape but masked
    where the iveg array has fill values
    """
    res = cH.cacheWrapper(B_org, **kwargs)
    return createIvegMaskedArray(res, **kwargs)


def x_org(**kwargs) -> dask.array.core.Array:
    """The state vector in the original cable shape"""
    return cH.reconstruct_x(
        cH.cacheWrapper(Cplant, **kwargs),
        cH.cacheWrapper(Clitter, **kwargs),
        cH.cacheWrapper(Csoil, **kwargs),
    )


def iveg_mask(**kwargs):
    """Creates a mask where the iveg array has fillvalues"""
    cable_data_set = kwargs["cable_data_set"]
    ifv = cH.get_integer_fill_value(cable_data_set)
    mask = cH.cacheWrapper(iveg, **kwargs) == ifv
    return mask


def createIvegMaskedArray(res, **kwargs):
    mask= cH.cacheWrapper(iveg_mask, **kwargs) 
    full_mask = dask.array.broadcast_to(mask, res.shape)
    return dask.array.ma.masked_array(res, mask=full_mask)


def x_org_iveg(**kwargs):
    """The solution in the original cable shape
    but masked where the iveg array has fill values
    """
    res = cH.cacheWrapper(x_org, **kwargs)
    return createIvegMaskedArray(res, **kwargs)


def x0_org(**kwargs) -> dask.array.core.Array:
    """The state vector at time 0 in the original cable shape"""
    return cH.cacheWrapper(x_org, **kwargs)[0, ...]


def x0_org_iveg(**kwargs) -> dask.array.core.Array:
    """The state vector at time 0 in the original cable shape
    but masked where the iveg array has fill values
    """
    res = cH.cacheWrapper(x0_org, **kwargs)
    return createIvegMaskedArray(res, **kwargs)


def u_org(**kwargs) -> dask.array.core.Array:
    """The input vector in the original cable shape"""
    cable_data_set = kwargs["cable_data_set"]

    res_ifv = cH.get_integer_fill_value(cable_data_set)
    res_ffv = cH.get_float_fill_value(cable_data_set)

    return cH.reconstruct_u(
        res_ifv,
        res_ffv,
        cH.cacheWrapper(iveg, **kwargs),
        cH.cacheWrapper(NPP, **kwargs),
        cH.cacheWrapper(fracCalloc, **kwargs),
    )


def u_org_iveg(**kwargs) -> dask.array.core.Array:
    """The input vector at time 0 in the original cable shape
    but masked where the iveg array has fill values
    """
    res = cH.cacheWrapper(u_org, **kwargs)
    return createIvegMaskedArray(res, **kwargs)


def trajectory_amplitude_org_iveg(**kwargs) -> dask.array.core.Array:
    """The difference between the largest and smallest value"""
    return cH.cacheWrapper(trajectory_max_org_iveg, **kwargs) - cH.cacheWrapper(
        trajectory_min_org_iveg, **kwargs
    )


def trajectory_min_org_iveg(**kwargs) -> dask.array.core.Array:
    """The smallest value for all trajectories"""
    return dask.array.min(cH.cacheWrapper(x_org_iveg, **kwargs), axis=0)


def trajectory_max_org_iveg(**kwargs):
    """The largest value for all trajectories"""
    return dask.array.max(cH.cacheWrapper(x_org_iveg, **kwargs), axis=0)


def all_pools_vary_cond(**kwargs):
    """we compute the condition for patch,landpoint combies where
    all the pools change over time

    """
    amplitude = cH.cacheWrapper(trajectory_amplitude_org_iveg, **kwargs)
    n_pools = amplitude.shape[0]
    conds = (amplitude[i, ...] > 0 for i in range(n_pools))
    # Note that in python Boolean arrays are integer
    # arrays and therefore multiplication implements
    # the logical 'and'.
    return reduce(lambda acc, el: acc * el, conds)


def one_pool_varies_cond(**kwargs):
    """we compute the condition for patch,landpoint combies where
    at least one of the pools changes over time
    """
    amplitude = cH.cacheWrapper(trajectory_amplitude_org_iveg, **kwargs)
    n_pools = amplitude.shape[0]
    conds = (amplitude[i, ...] > 0 for i in range(n_pools))
    # Note that in python Boolean arrays are integer
    # arrays and therefore addition implements
    # the logical 'or'.
    return reduce(lambda acc, el: acc + el, conds)


def cond_patches(cond, **kwargs):
    # Note:
    # although the parallel version 
    #   res = dask.array.nonzero(cond)
    #   ps, lps = res
    #   ps.compute_chunk_sizes()
    # works
    # it suffers from the tiny chunksize (1,1) which in this case implies a lot
    # of communication because every chunk has its own thread to report if its
    # patch,landpoint combi fullfills the condition. For this application it is
    # much faster to:
    # 1.) convert the whole array to numpy (implying a complete compute)
    # 2.) call np.nonzero
    # 3.) create a properly chunked dask.array afterwards
    # This is only possible because the condition arrays are small enough to be
    # handled by a single worker even for the whole grid. 
    npc = np.array(cond)
    patches, _ = np.nonzero(npc)
    return dask.array.from_array(patches) 


def cond_landpoints(cond, **kwargs):
    # Note:
    # although the parallel version 
    #    res = dask.array.nonzero(cond)
    #    ps, lps = res
    #    lps.compute_chunk_sizes()
    #    return res
    # works
    # it suffers from the tiny chunksize (1,1) which in this case implies a lot
    # of communication because every chunk has its own thread to report if its
    # patch,landpoint combi fullfills the condition. For this application it is
    # much faster to:
    # 1.) convert the whole array to numpy (implying a complete compute)
    # 2.) call np.nonzero
    # 3.) create a properly chunked dask.array afterwards
    # This is only possible because the condition arrays are small enough to be
    # handled by a single worker even for the whole grid. 
    # we now know the number of valid entries
    npc = np.array(cond)
    _ , landpoints = np.nonzero(npc)
    return dask.array.from_array(landpoints) 


def all_pools_vary_cond_patches(**kwargs):
    cond = cH.cacheWrapper(all_pools_vary_cond, **kwargs)
    return cond_patches(cond, **kwargs)


def all_pools_vary_cond_landpoints(**kwargs):
    cond = cH.cacheWrapper(all_pools_vary_cond, **kwargs)
    return cond_landpoints(cond, **kwargs)


def all_pools_vary_cond_nz(**kwargs):
    # Note that this function can not be used within
    # the cacheWrapper since it returns a tuple
    return (
        cH.cacheWrapper(all_pools_vary_cond_patches, **kwargs),
        cH.cacheWrapper(all_pools_vary_cond_landpoints, **kwargs)
    )


def one_pool_varies_cond_patches(**kwargs):
    cond = cH.cacheWrapper(one_pool_varies_cond, **kwargs)
    return cond_patches(cond, **kwargs)


def one_pool_varies_cond_landpoints(**kwargs):
    cond = cH.cacheWrapper(one_pool_varies_cond, **kwargs)
    return cond_landpoints(cond, **kwargs)


def one_pool_varies_cond_nz(**kwargs):
    # Note that this function can not be used within
    # the cacheWrapper since it returns a tuple
    return (
        cH.cacheWrapper(one_pool_varies_cond_patches, **kwargs),
        cH.cacheWrapper(one_pool_varies_cond_landpoints, **kwargs)
    )




def B_val(**kwargs) -> dask.array.core.Array:
    return cH.valid_combies_parallel(
        cH.cacheWrapper(
			one_pool_varies_cond_nz,
			**kwargs
        ),
        cH.cacheWrapper(B_org, **kwargs)
    )


def u_val(**kwargs) -> dask.array.core.Array:
    return cH.valid_combies_parallel(
        cH.cacheWrapper(
			one_pool_varies_cond_nz,
			**kwargs
        ),
        cH.cacheWrapper(u_org, **kwargs)
    )


def x_val(**kwargs) -> dask.array.core.Array:
    return cH.valid_combies_parallel(
        cH.cacheWrapper(
			one_pool_varies_cond_nz,
			**kwargs
        ),
        cH.cacheWrapper(x_org, **kwargs)
    )


def x0_val(**kwargs) -> dask.array.core.Array:
    return cH.valid_combies_parallel(
        cH.cacheWrapper(
			one_pool_varies_cond_nz,
			**kwargs
        ),
        cH.cacheWrapper(x0_org, **kwargs)
    )


def sol_val(**kwargs) -> dask.array.core.Array:
    return dask.array.blockwise(
        cH.valid_trajectory,
        "ijk",
        cH.cacheWrapper(x0_val, **kwargs),
        "jk",
        cH.cacheWrapper(time, **kwargs),
        "i",
        cH.cacheWrapper(B_val, **kwargs),
        "ijjk",
        cH.cacheWrapper(u_val, **kwargs),
        "ijk",
        dtype="f8",
    )


def sol_org(**kwargs) -> dask.array.core.Array:
    return dask.array.blockwise(
        cH.trajectory_org,
        "ijkl",
        cH.cacheWrapper(x0_org, **kwargs),
        "jkl",
        cH.cacheWrapper(time, **kwargs),
        "i",
        cH.cacheWrapper(B_org, **kwargs),
        "ijjkl",
        cH.cacheWrapper(u_org, **kwargs),
        "ijkl",
        dtype="f8",
    )


def sol_org_iveg(**kwargs) -> dask.array.core.Array:
    return dask.array.blockwise(
        cH.trajectory_org,
        "ijkl",
        cH.cacheWrapper(x0_org_iveg, **kwargs),
        "jkl",
        cH.cacheWrapper(time, **kwargs),
        "i",
        cH.cacheWrapper(B_org_iveg, **kwargs),
        "ijjkl",
        cH.cacheWrapper(u_org_iveg, **kwargs),
        "ijkl",
        dtype="f8",
    )
