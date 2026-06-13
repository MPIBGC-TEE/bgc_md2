import dask
import shutil
import inspect
from functools import _lru_cache_wrapper
from tqdm import tqdm  # Holger
from pathlib import Path
import zarr as zr

def batchSlices(nland, nproc):
    return [
        slice(i*nproc,min((i+1)*nproc,nland))
        for i in range(int(nland/nproc)+1)
    ]

def batchwise_to_zarr(
    arr: dask.array.core.Array,
    zarr_dir_name: str,
    rm: bool = False,
    batch_size: int = 1
):
    dir_p = Path(zarr_dir_name)
    if dir_p.exists():
        if rm:
            print("##########################################")
            print("removing " + str(dir_p))
            shutil.rmtree(dir_p)
        else:
            print("##########################################")
            print(str(dir_p) + " already exists")
            return

    if False:  # arr.nbytes < 8 * 1024 ** 3:
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
        #ncores = 32
        slices = batchSlices(arr.shape[-1], batch_size)
        print("result shape:", arr.shape)
        #        print(629, slices)
        #        for s in slices:

        for s in tqdm(slices):  # Holger
            z[..., s] = arr[..., s].compute()


def module_functions(cmod):
    def pred(a):
        return inspect.isfunction(a) or isinstance(a,_lru_cache_wrapper)
    return frozenset(
        [
            getattr(cmod, c)
            for c in cmod.__dir__()
            if pred(getattr(cmod, c))
        ]
    )

def val_or_default(
        d: dict
        ,key
        ,default 
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
    batch_size_pass = val_or_default(kwargs,'batch_size',1) 
    # remove the rm flag since from the passed on kwargs
    # since we only want it to take effect on the highest level
    # which is here
    rec_kw_args={k:v for k,v in kwargs.items() if k != 'rm'}
    res = cachable_func(**rec_kw_args)
    if isinstance(res,dask.array.core.Array):
        batchwise_to_zarr(
            res,
            str(sub_dir_path),
            rm=rm_pass,
            batch_size=batch_size_pass
        ) 
    return res 

def arg_names(org_func):
    s = inspect.signature(org_func)
    return [name for name in s.parameters.keys()]

def mod_func_names(cmod):
    return [f.__name__ for f in module_functions(cmod)]

def non_rec_arg_names(org_func, cmod):
    return frozenset.difference(
        frozenset(arg_names(org_func)),
        frozenset(mod_func_names(cmod))
    )

def rec_arg_funcs(org_func, cmod):
    return [
            func
            for func in module_functions(cmod)
            if func.__name__ in arg_names(org_func)
    ]

def fromCachedArgsMaker2(org_func, cmod):
    """This function takes an original function org_func from module cmod. 
    Function org_func does not look up the cache neither for its own value nor
    for its arguments.
    The function  returned by this function does:
    - look up its arguments in the cache
    - and is compatible to be called in cacheWrapper to be cached itself
    Note:
    The functions in cmod have to follow the convention that
    The names of the parameters in the signature are equal to the
    function in cmod that computes them.
    So:
    def C(A,B):
        ...
    in cmod implies that there are 
    functions wiht name A and B in cmod
    """

    def cachable(**kwargs):
        # we evaluate all the function arguments
        org_funcs = rec_arg_funcs(org_func,cmod)
        funcs = [fromCachedArgsMaker2(f,cmod) for f in org_funcs]
        def get_res(func):
            return cacheWrapper(
                func,
                **kwargs
            )


        results = map(get_res, funcs)
        # and call the original funtions with it and those 
        # kwargs that it excepts
        pass_on_kwargs = {
                k: v
                for k, v in kwargs.items()
                if k in non_rec_arg_names(org_func, cmod)
        }
        return org_func(*results,**pass_on_kwargs)

    cachable.__name__ = org_func.__name__
    return cachable

def cachedFromCachedArgsMaker2(org_func, cmod):
    """This function takes an original function org_func from module cmod.
    Function org_func does not look up the cache neither for its own value nor
    for its arguments.
    The function  returned by this function does:
    - loos up its own value in the cache
    - look up its arguments in the cache

    It looks for its dependencies in cmod. 
    
    Note:
    The functions in cmod have to follow the convention that
    The names of the parameters in the signature are equal to the
    function in cmod that computes them.
    So:
    def C(A,B):
        ...
    in cmod implies that there are 
    functions wiht name A and B in cmod
    """

    cachable = fromCachedArgsMaker2(org_func, cmod)

    def cached(**kwargs):
        return cacheWrapper(
                cachable,
                **kwargs
        )

    return cached
