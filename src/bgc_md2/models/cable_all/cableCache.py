import dask.array
from . import cableHelpers as cH

# iveg is not time dependent and so the function
# will not accept the time_slice argument
def iveg(
    cable_data_set,
    landpoint_slice=slice(None,None,None),
    **kwargs
):
    arr = dask.array.asarray(cable_data_set['iveg'].data)
    return arr[...,landpoint_slice]

def time(
    cable_data_set,
    time_slice=slice(None,None,None),
    **kwargs
):
    arr = dask.array.asarray(cable_data_set['time'].data)
    print(arr.shape)
    sl = arr[time_slice]
    print(sl.shape)
    return sl

# most of the other cable outputs are very similar
# to each other. So we create the functions that
# return them automatically from a template
def get_singe_cable_var(
    name,
    cable_data_set,
    time_slice=slice(None,None,None),
    landpoint_slice=slice(None,None,None)
):
    arr = dask.array.asarray(cable_data_set[name].data)
    print(name)
    print(arr.shape)
    sl = arr[time_slice,...,landpoint_slice]
    print(sl.shape)
    return sl

def FuncMaker(name):
    def func(
        cable_data_set,
        time_slice=slice(None,None,None),
        landpoint_slice=slice(None,None,None),
        **kwargs
    ):
        return get_singe_cable_var(
            name=name, 
            cable_data_set=cable_data_set,
            time_slice=time_slice,
            landpoint_slice=landpoint_slice
        )
    
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
iveg = FuncMaker("iveg")
kplant = FuncMaker("kplant")
xkNlimiting = FuncMaker("xkNlimiting")
xktemp = FuncMaker("xktemp")
xkwater = FuncMaker("xkwater")

def FuncMaker0(name):
    def func(
        cable_data_set,
        time_slice=slice(0,None,None),
        landpoint_slice=slice(None,None,None),
        **kwargs
    ):
        start_time = time_slice.start
        return get_singe_cable_var(
            name=name, 
            cable_data_set=cable_data_set,
            time_slice=slice(start_time,start_time,1),
            landpoint_slice=landpoint_slice
        )
    
    func.__name__ = name
    return func

Clitter0 = FuncMaker0("Clitter")
Cplant0 = FuncMaker0("Cplant")
Csoil0 = FuncMaker0("Csoil")

def B_org(**kwargs) -> dask.array.core.Array:
    """The B matrix in the original cable shape 
    """
    cable_data_set = kwargs['cable_data_set']

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
        cH.cacheWrapper(xkNlimiting, **kwargs)
    )


def x_org(**kwargs) -> dask.array.core.Array:
    """The state vector in the original cable shape
    """
    return cH.reconstruct_x(
        cH.cacheWrapper(Cplant, **kwargs),
        cH.cacheWrapper(Clitter, **kwargs),
        cH.cacheWrapper(Csoil, **kwargs)
    )


def x0_org(**kwargs) -> dask.array.core.Array:
    """The state vector in the original cable shape
    """
    res_Cplant = cH.cacheWrapper(Cplant0, **kwargs)
    res_Clitter = cH.cacheWrapper(Clitter0, **kwargs)
    res_Csoil = cH.cacheWrapper(Csoil0, **kwargs)

    return cH.reconstruct_x(
            res_Cplant,
            res_Clitter,
            res_Csoil
    )


def u_org(**kwargs) -> dask.array.core.Array:
    """The input vector in the original cable shape
        """
    cable_data_set = kwargs['cable_data_set']

    res_ifv = cH.get_integer_fill_value(cable_data_set)
    res_ffv = cH.get_float_fill_value(cable_data_set)


    return reconstruct_u(
            res_ifv,
            res_ffv, 
    		cH.cacheWrapper(iveg, **kwargs),
    		cH.cacheWrapper(NPP, **kwargs),
    		cH.cacheWrapper(fracCalloc, **kwargs)
    )


def nz(**kwargs) -> dask.array.core.Array:
    cable_data_set = kwargs['cable_data_set']

    res_ifv = cH.get_integer_fill_value(cable_data_set)

    res_iveg = cH.cacheWrapper(iveg, **kwargs)
    res_Csoil = cH.cacheWrapper(Csoil, **kwargs)

    cond_1 = (res_iveg != res_ifv).compute()
    cond_2 = (res_Csoil[0, 0, :, :] != 0).compute()

    return dask.array.nonzero(cond_1 * cond_2)


def B_val(**kwargs) -> dask.array.core.Array:
    return cH.valid_combies_parallel(
    	cH.cacheWrapper(nz, **kwargs),
        cH.cacheWrapper(B_org, **kwargs)
    )


def u_val(**kwargs) -> dask.array.core.Array:
    return cH.valid_combies_parallel(
        cH.cacheWrapper(nz, **kwargs),
        cH.cacheWrapper(u_org, **kwargs)
    )


def x_val(**kwargs) -> dask.array.core.Array:
    return cH.valid_combies_parallel(
        cH.cacheWrapper(nz, **kwargs),
        cH.cacheWrapper(x_org, **kwargs)
    )


def x0_val(**kwargs) -> dask.array.core.Array:
    return cH.valid_combies_parallel(
        cH.cacheWrapper(nz, **kwargs),
        cH.cacheWrapper(x0_org, **kwargs)
    )


def sol_val(**kwargs) -> dask.array.core.Array:
    res_x0_val = cH.cacheWrapper(x0_val, **kwargs)
    res_times = cH.cacheWrapper(times_val, **kwargs)
    res_B_val = cH.cacheWrapper(B_val, **kwargs)
    res_u_val = cH.cacheWrapper(u_val, **kwargs)

    return dask.array.blockwise(
             cH.valid_trajectory, 'ijk',
             res_x0_val, 'jk',
             res_times, 'i',
             res_B_val, 'ijjk',
             res_u_val, 'ijk',
             dtype='f8'
    )
