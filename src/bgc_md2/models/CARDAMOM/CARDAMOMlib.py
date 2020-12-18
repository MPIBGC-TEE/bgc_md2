import dask
import time
import zarr 

import dask.array as da
import numpy as np

from dask.distributed import LocalCluster
from getpass import getuser

from bgc_md2.ModelStructure import ModelStructure
from bgc_md2.ModelDataObject_dict import ModelDataObject_dict as ModelDataObject
from bgc_md2.Variable import Variable
from bgc_md2.notebook_helpers import (
    write_to_logfile,
    write_header_to_logfile,
    custom_timeout
)


def prepare_cluster(n_workers, alternative_port=None):
    port_dict = {
        "cs": 8888        # change at will
        "mm": 8889,
        "hmetzler": 8890, # change at will
    }
    my_user_name = getuser()
    print("username:", my_user_name)

    if alternative_port is None:
        my_port = port_dict[my_user_name]
    else:
        myport = alternative_port

    print("notebook port:", my_port)
    
    # allow starting subprocesses --> allow controlled termination of the process at timeout
    dask.config.set({"distributed.worker.daemon": False})
    
    # dashboard needs a different port for accessing it remotely
    my_dashboard_port = my_port - 100

    # seem to be ignored, use dask.distributed.get_worker as the worker process is running
    # and set the values manually:
    # worker = get_worker()
    # worker.memory_limit ="2G"
    # ...
    worker_kwargs = {
        "memory_limit": "2G",
        "memory_target_fraction": 0.95,
        "memory_spill_fraction": 0.95,
        "memory_pause_fraction": 0.95,
        #"memory_terminate_fraction": False, # leads to errors if commented in
    }

    my_cluster = LocalCluster(
        dashboard_address="localhost:"+str(my_dashboard_port),
        n_workers=n_workers,
        threads_per_worker=1,
#        **worker_kwargs # seem to be ignored
    )
    print("default dashboard port (might change if busy):", my_dashboard_port)

    return my_cluster


def load_variables_from_zarr(zarr_data_path, non_data_variables):
    variable_paths = [p for p in zarr_data_path.iterdir() if p.is_dir()]

    variable_names = []
    variables_total = []
    for variable_path in variable_paths:
        name = variable_path.name
        if name not in non_data_variables:
            variable_names.append(name)
            variables_total.append(da.from_zarr(str(variable_path)))

    # manual reshaping for non_data_variables
    # necessary for map_blocks later on
    for k, name in enumerate(non_data_variables):
        shape = [1] * len(non_data_variables)
        shape[k] = -1
        variables_total.append(da.from_zarr(str(zarr_data_path.joinpath(name)), chunks=(1,)).reshape(shape))
        variable_names.append(name)

    return variable_names, variables_total


def select_slices_of_interest(
    variable_names,
    variables_total,
    slices,
    non_data_variables
):
    variables = []
    for nr, name, var in zip(range(len(variables_total)), variable_names, variables_total):
        if name not in non_data_variables:
            variables.append(variables_total[nr][slices["lat"], slices["lon"], slices["prob"], ...])

    for k, name in enumerate(non_data_variables):
        index = variable_names.index(name)
        v = variables_total[index]
        s = [slice(0, None, 1),] * v.ndim
        s[k] = slices[name]
        variables.append(v[tuple(s)])

    return variables


def batchwise_to_zarr(
    arr, # dask array
    z, # target zarr archive
    z_slices, # slices of interest
    batch_nrs,
    noisy=False
):  
    dims = {"lat": 0, "lon": 1, "prob": 2, "time": 3}

    def split_arr_to_slices(ns):
        res = dict()
        for name in ["lat", "lon", "prob"]:
            a = np.arange(0, arr.shape[dims[name]])
            l = np.array_split(a, ns[name])
            res[name] = [slice(x[0], x[-1]+1, 1) for x in l]
            
        return res
    
    arr_slices_dict = split_arr_to_slices(batch_nrs)
       
    arr_slices_tuples = [
        (slice_lat, slice_lon, slice_prob)
        for slice_lat in arr_slices_dict["lat"] 
        for slice_lon in arr_slices_dict["lon"]
        for slice_prob in arr_slices_dict["prob"]
    ]

    def split_slices(slices, ns):
        res = dict()
        for name, s in slices.items():
            if name != "time":
                a = np.arange(0, z.shape[dims[name]], 1)[s]
                l = np.array_split(a, ns[name])
                res[name] = [slice(x[0], x[-1]+1, s.step) for x in l]
        
        return res
          
    z_slices_dict = split_slices(z_slices, batch_nrs)
    
    z_slices_tuples = [
        (slice_lat, slice_lon, slice_prob)
        for slice_lat in z_slices_dict["lat"] 
        for slice_lon in z_slices_dict["lon"]
        for slice_prob in z_slices_dict["prob"]
    ]
   
    for s_arr, s_z in zip(arr_slices_tuples, z_slices_tuples):
        if noisy:
            print("CARDAMOMlib.batchwise_to_zarr")
            print("s_arr")
            print(s_arr)
            print("s_z")
            print(s_z)
            print()
        z[s_z[0], s_z[1], s_z[2], ...] = arr[s_arr[0], s_arr[1], s_arr[2], ...].compute()


# args[0] is the object on which map_blocks is called
# args[:-1] are the variables that get automatically chunked
# args[-1] is supposed to be a dictionary with additional parameters
def func_for_map_blocks(*args):
    additional_params = args[-1]

    func = additional_params["func"] # the actual function to be called

    v_names = additional_params["variable_names"]
    time_limit_in_min = additional_params["time_limit_in_min"]
    return_shape = additional_params["return_shape"]
    logfile_name = additional_params["logfile_name"]

    # create a dctionary that acts as a dataset replacement
    d = {v_names[i]: args[i].reshape((-1,)) for i in range(len(args[:-1]))}

    res = -np.inf * np.ones(return_shape)
    start_time = time.time()
    error_msg = ""
    try:
        res = custom_timeout(
            time_limit_in_min*60,
            func,
            d
        )
    except TimeoutError:
        duration = (time.time() - start_time) / 60
        error_msg = "Timeout after %2.2f min" % duration
        print(error_msg, flush=True)
    except Exception as e:
        res = np.nan * np.ones(return_shape)
        error_msg = str(e)
        print(error_msg, flush=True)

    write_to_logfile(
        logfile_name,
        "finished single,",
        "lat:", d["lat"],
        "lon:", d["lon"],
        "prob:", d["prob"],
        error_msg
    )

    return res.reshape(return_shape)


def compute_start_values(single_site_dict):
    mdo = _load_mdo(single_site_dict)
    xs = mdo.load_stocks()
    del mdo

    return xs.data.filled()[0, ...]


def compute_us(single_site_dict):
    mdo = _load_mdo(single_site_dict)
    us = mdo.load_us()
    del mdo

    return us


def compute_Bs(
    single_site_dict,
    integration_method='solve_ivp',
    nr_nodes=None
):
    mdo = _load_mdo(single_site_dict)
    Bs = mdo.load_Bs(
        integration_method,
        nr_nodes
    )
    del mdo

    return Bs


def get_timeout_sites(zarr_path, slices):
    timeout_zarr = zarr.open(str(zarr_path))
    fill_tup = (slice(0, 1, 1), ) * (timeout_zarr.ndim - 3)
    tup = (slices['lat'], slices['lon'], slices['prob']) + fill_tup
    timeout_da = da.from_zarr(timeout_zarr)[tup]
    timeout_coords = np.where(timeout_da.compute() == -np.inf)[:3]
    nr_timeout_sites = len(timeout_coords[0])

    return nr_timeout_sites, timeout_coords


def remove_timeouts(
    time_limit_in_min,
    zarr_path,
    variable_names,
    variables,
    non_data_variables,
    slices,
    task,
    additional_params
):
    logfile_name = additional_params["logfile_name"]

    nr_timeout_sites, timeout_coords = get_timeout_sites(zarr_path, slices)

    if nr_timeout_sites > 0:
        print('number of timeout sites:', nr_timeout_sites)
    else:
        print('timeout sites successfully removed')
        return

    # select timeout sites from variables
    timeout_variables = []
    for v, name in zip(variables, variable_names):
        if name not in non_data_variables:
            v_stack_list = []
            for lat, lon, prob in zip(*[timeout_coords[k] for k in range(3)]):
                v_stack_list.append(v[lat, lon, prob])

            timeout_variables.append(da.stack(v_stack_list))

    # add lat, lon, prob, time
    timeout_variables.append(da.from_array(timeout_coords[0].reshape(-1, 1), chunks=(1, 1))) # lat
    timeout_variables.append(da.from_array(timeout_coords[1].reshape(-1, 1), chunks=(1, 1))) # lon
    timeout_variables.append(da.from_array(timeout_coords[2].reshape(-1, 1), chunks=(1, 1))) # prob
    time_da = variables[variable_names.index('time')].reshape(1, -1).rechunk((1, 1152))
    timeout_variables.append(time_da)
#    x = timeout_variables[-1].rechunk((1, 1152))
#    print('x', x, flush=True)
#    timeout_variables[-1] = 'LECK MICH!'
#    timeout_variables[-1] = x
    # prepare the delayed computation
    additional_params['time_limit_in_min'] = time_limit_in_min

    _task = task.copy()
    _task['return_shape'] = task['return_shape'][2:]
    _task['meta_shape'] = (15, 1152, 6, 6)
    _task['drop_axis'] = []
    _task['new_axis'] = [2, 3]

    _additional_params = additional_params.copy()
    _additional_params['return_shape'] = _task["return_shape"]

    res_da = timeout_variables[0].map_blocks(
        func_for_map_blocks,
        *timeout_variables[1:], # variables[0] comes automatically as first argument
        _additional_params,
        drop_axis=_task["drop_axis"],
        new_axis=_task["new_axis"],
        chunks=_task["return_shape"],
        dtype=np.float64,
        meta=np.ndarray(task["meta_shape"][2:], dtype=np.float64)
    )

    # write header to logfile
    print(write_header_to_logfile(logfile_name, res_da, time_limit_in_min))
    print('starting, timeout (min) =', time_limit_in_min, flush=True)

    # do the computation
    res_computed = res_da.compute() # hopefully fits in memory

    # compute the z coords from the sliced Bs coords
    f_lat = lambda x: slices['lat'].start + x * slices['lat'].step
    f_lon = lambda x: slices['lon'].start + x * slices['lon'].step
    f_prob = lambda x: slices['prob'].start + x * slices['prob'].step
    timeout_coords_z_lat = np.array([f_lat(x) for x in timeout_coords[0]])
    timeout_coords_z_lon = np.array([f_lon(x) for x in timeout_coords[1]])
    timeout_coords_z_prob = np.array([f_prob(x) for x in timeout_coords[2]])

    timeout_coords_z = (timeout_coords_z_lat, timeout_coords_z_lon, timeout_coords_z_prob)

    # update zarr file
    z = zarr.open(str(zarr_path))
    for nr, lat, lon, prob in zip(range(nr_timeout_sites), *[timeout_coords_z[k] for k in range(3)]):
        z[lat, lon, prob, ...] = res_computed[nr, ...]

    write_to_logfile(logfile_name, 'done, timeout (min) =', time_limit_in_min)
    print('done, timeout (min) = ', time_limit_in_min, flush=True)


###############################################################################
#
# internal methods
#
###############################################################################


def _load_model_structure():
    # labile, leaf, root, wood, litter, and soil

    pool_structure = [
        {"pool_name": "Labile", "stock_var": "c_labile"},
        {"pool_name": "Leaf", "stock_var": "c_foliar"},
        {"pool_name": "Root", "stock_var": "c_root"},
        {"pool_name": "Wood", "stock_var": "c_wood"},
        {"pool_name": "Litter", "stock_var": "c_finelitter"},
        {"pool_name": "Soil", "stock_var": "c_som"},
    ]

    external_input_structure = {
        "Labile": ["gpp_to_labile"],
        "Leaf": ["gpp_to_leaf"],
        "Root": ["gpp_to_root"],
        "Wood": ["gpp_to_wood"],
    }

    horizontal_structure = {
        ("Labile", "Leaf"): ["labile_to_foliar"],
        ("Labile", "Litter"): ["fire_labile_to_litter"],
        ("Leaf", "Litter"): ["leaf_to_litter", "fire_foliar_to_litter"],
        ("Wood", "Soil"): ["wood_to_soilc", "fire_wood_to_som"],
        ("Root", "Litter"): ["root_to_litter", "fire_root_to_litter"],
        ("Litter", "Soil"): ["litter_to_som", "fire_litter_to_som"],
    }

    vertical_structure = {}

    external_output_structure = {
        "Labile": ["fire_em_labile"],
        "Leaf": ["fire_em_foliar"],
        "Root": ["fire_em_root"],
        "Wood": ["fire_em_wood"],
        "Litter": ["hetresp_litter", "fire_em_litter"],
        "Soil": ["hetresp_som", "fire_em_som"]
    }

    model_structure = ModelStructure(
        pool_structure=pool_structure,
        external_input_structure=external_input_structure,
        horizontal_structure=horizontal_structure,
        vertical_structure=vertical_structure,
        external_output_structure=external_output_structure,
    )

    return model_structure


def _load_mdo(ds_dict):
    ms = _load_model_structure()

    days_per_month = 31.0

    # no unit support for dictionary version
    time = Variable(
        name="time",
        data=np.arange(len(ds_dict['time'])) * days_per_month,
#        unit="d"
        unit="1"
    )

    mdo = ModelDataObject(
        model_structure=ms,
        dataset=ds_dict, 
#        stock_unit="gC/m2", 
        stock_unit="1", 
        time=time
    )

    return mdo



