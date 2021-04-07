import dask
import time
import traceback
import zarr 

import dask.array as da
import numpy as np

from pathlib import Path
from scipy.interpolate import interp1d
from sympy import symbols, Matrix
from tqdm import tqdm

from dask.distributed import LocalCluster
from getpass import getuser

from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
from CompartmentalSystems.discrete_model_run_14C import DiscreteModelRun_14C as DMR_14C
from CompartmentalSystems.helpers_reservoir import DECAY_RATE_14C_DAILY, ALPHA_14C
from CompartmentalSystems.pwc_model_run_fd import PWCModelRunFD
from bgc_md2.ModelStructure import ModelStructure
from bgc_md2.ModelDataObject import ModelDataObject
from bgc_md2.Variable import Variable
from bgc_md2.notebook_helpers import (
    write_to_logfile,
    write_header_to_logfile,
    custom_timeout,
    load_zarr_archive
)


def prepare_cluster(n_workers, alternative_dashboard_port=None):
    port_dict = {
        "cs": 8888,        # change at will
        "mm": 8889,
        "hmetzler": 8890 # change at will
    }
    my_user_name = getuser()
    print("username:", my_user_name)

    my_port = port_dict[my_user_name]
    print("notebook port:", my_port)
    
    # allow starting subprocesses --> allow controlled termination of the process at timeout
    dask.config.set({"distributed.worker.daemon": False})
    
    # dashboard needs a different port for accessing it remotely
    my_dashboard_port = my_port - 100
    if alternative_dashboard_port is not None:
        my_dashboard_port = alternative_dashboard_port

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

def load_params(time_resolution, delay_in_months):
    if time_resolution == "daily":
        return {
            "output_folder": "daily_output",
            "time_step_in_days": 1,
            "nr_time_steps": 10*12*31
        }
    elif time_resolution == "monthly":
        return {
            "output_folder": "monthly_output",
            "time_step_in_days": 31,
            "nr_time_steps": 10*12
        }
    elif time_resolution == "yearly":
        return {
            "output_folder": "yearly_%02d_output" % delay_in_months,
            "time_step_in_days": 12*31,
            "nr_time_steps": 10
        }
    else:
        raise(ValueError("Unknown time resolution"))

    
def check_data_consistency(ds, time_step_in_days):
    ms = load_model_structure()

    # data_test.py says tht for labile 31.0 is constantly right
    time = Variable(
        name="time",
        data=np.arange(len(ds.time)) * time_step_in_days,
        unit="d"
    )

    mdo = ModelDataObject(
        model_structure=ms,
        dataset=ds, 
        stock_unit="gC/m2", 
        time=time
    )

    abs_err, rel_err = mdo.check_data_consistency()
    return abs_err, rel_err


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


def linear_batchwise_to_zarr(
    res_da, # dask array
    z, # target zarr archive
    slices, # slices of interest,
    coords_tuples, # global coordinates (not sliced)
    batch_size,
):  
    nr_total = len(coords_tuples)
    nr_splits = nr_total // batch_size + (nr_total % batch_size > 0)
    batches = np.array_split(np.arange(nr_total), nr_splits)

    for batch in tqdm(batches):
        s = slice(batch[0], batch[-1]+1, 1)
        res_batch = res_da[s, ...].compute() # in parallel
       
        # update zarr file
        for nr, coords in enumerate(coords_tuples[s]):
            z[coords] = res_batch[nr, ...]


# args[0] is the object on which map_blocks is called
# args[:-1] are the variables that get automatically chunked
# args[-1] is supposed to be a dictionary with additional parameters
def func_for_map_blocks(*args):
    additional_params = args[-1]

    func = additional_params["func"] # the actual function to be called
    func_args = additional_params["func_args"]

    v_names = additional_params["variable_names"]
    time_limit_in_min = additional_params["time_limit_in_min"]
    return_shape = additional_params["return_shape"]
    logfile_name = additional_params["logfile_name"]

    # create a dctionary that acts as a dataset replacement
    d = {v_names[i]: args[i].reshape((-1,)) for i in range(len(args[:-1]))}

    res = -np.inf * np.ones(return_shape)
    start_time = time.time()
    error_msg = ""
    info = tuple()
    try:
        res, info = custom_timeout(
            time_limit_in_min*60,
            func,
            d,
            **func_args
        )
    except TimeoutError:
        duration = (time.time() - start_time) / 60
        error_msg = "Timeout after %2.2f min" % duration
        print(error_msg, flush=True)
    except Exception as e:
        tb = traceback.format_exc()
#        print(tb, flush=True)
        res = np.nan * np.ones(return_shape)
        error_msg = "Error: " + str(e)
        print(error_msg, flush=True)

    if error_msg == "":
        error_msg = "done"

    write_to_logfile(
        logfile_name,
        d["lat"],
        d["lon"],
        d["prob"],
        error_msg,
        *info
    )

    return res.reshape(return_shape)


# args[0] is the object on which map_blocks is called (Bs_da)
# args[:-1] are the variables that get automatically chunked
# args[-1] is supposed to be a dictionary with additional parameters
def func_for_map_blocks_with_mr(*args):
    additional_params = args[-1]

    model_type = additional_params["model_type"]
    computation = additional_params["computation"]
    nr_pools = additional_params["nr_pools"]
    time_step_in_days = additional_params["time_step_in_days"]
    return_shape = additional_params["return_shape"]
    func = additional_params["func"]
    func_args = additional_params["func_args"]
    time_limit_in_min = additional_params["time_limit_in_min"]
    logfile_name = additional_params["logfile_name"]
                                            
    Bs = args[0].reshape((-1, nr_pools, nr_pools))
    start_values = args[1].reshape(nr_pools)
    us = args[2].reshape((-1, nr_pools))
                                                            
    lat = args[3].reshape(-1)
    lon = args[4].reshape(-1)
    prob = args[5].reshape(-1)
    times = args[6].reshape(-1)

    nr_times = len(times)
    data_times = np.arange(0, nr_times, 1) * time_step_in_days
                                                                                        
    log_msg = "done"
    res = -np.inf * np.ones(return_shape)
    start_time = time.time()
    try:
        time_symbol = symbols('t')

        # using the custom_timeout function (even with timeout switched off) 
        # makes the worker report to the scheduler
        # and prevents frequent timeouts

        if model_type == "continuous":
            mr = custom_timeout(
                np.inf,
                PWCModelRunFD.from_Bs_and_us,
                time_symbol,
                data_times,
                start_values,
                Bs[:-1],
                us[:-1]
            )
        elif model_type == "discrete":       
            mr = custom_timeout(
                np.inf,
                DMR.from_Bs_and_net_Us,
                start_values,
                data_times,
                Bs[:-1],
                us[:-1]
            )
        else:
            raise(ValueError("model_type not recognized"))

        print("computing", computation, lat, lon, prob, flush=True)
        res = custom_timeout(
            time_limit_in_min*60,
            func,
            mr,
            **func_args
        )           
        print("done", lat, lon, prob, flush=True)
    except TimeoutError:
        duration = (time.time() - start_time) / 60
        log_msg = "Timeout after %2.2f min" % duration
        print(log_msg, flush=True)
    except Exception as e:
        tb = traceback.format_exc()
        res = np.nan * np.ones_like(res)
        print(str(e), flush=True)
        print(tb, flush=True)
        log_msg = "Error: " + str(e) + str(tb)

    write_to_logfile(
        logfile_name,
        lat,
        lon,
        prob,
        log_msg
    )
    
    return res.reshape(return_shape)


def compute_xs(single_site_dict, time_step_in_days):
    mdo = _load_mdo(single_site_dict, time_step_in_days, check_units=False)
    xs = mdo.load_stocks()
    del mdo

    info = tuple()
    return xs.data.filled(), info


def compute_start_values(single_site_dict, time_step_in_days):
    mdo = _load_mdo(single_site_dict, time_step_in_days, check_units=False)
    xs = mdo.load_stocks()
    del mdo

    info = tuple()
    return xs.data.filled()[0, ...], info


def compute_start_values_14C(
    mr,
    nr_time_steps,
    CARDAMOM_path=Path("/home/data/CARDAMOM/")
):
#    CARDAMOM_path = Path("/home/hmetzler/Desktop/CARDAMOM/")
    intcal20_path = CARDAMOM_path.joinpath("IntCal20_Year_Delta14C.csv")

    intcal20 = np.loadtxt(
        intcal20_path,
        delimiter=" ",
        skiprows=1,
        usecols=(1, 2)
    ).transpose()

    left_val = intcal20[1][np.argmin(intcal20[0])]
    F_atm_raw = interp1d(
        intcal20[0],
        intcal20[1],
#        kind="cubic",
        bounds_error=False,
        fill_value=(left_val, -1000.0) # no 14C oustide of dataset
    )

    t0 = 1920 # years AD
    lim = t0 + np.max(-intcal20[0]) # max. age rechable by intcal in years

    # input comes in days before t0 (age)
    # F_atm_raw input needs to be years AD 
    def F_atm(a):
        y = a / 31 / 12 # years before t0
        y = t0 - y # years AD

        return F_atm_raw(y)

    start_values_14C = mr.fake_eq_14C(
        nr_time_steps,
        F_atm,
        DECAY_RATE_14C_DAILY,
        lim * 31 * 12 * 2 # integration limit in days
    )

    return start_values_14C


def compute_solution_14C(dmr, nr_time_steps, Delta14C_atm_path):
    if not isinstance(dmr, DMR):
        raise(TypeError("wrong type of model run"))

    t0 = 1920
    Delta_14C = np.loadtxt(Delta14C_atm_path, delimiter=",", skiprows=1).transpose()

    F_atm = interp1d(
        Delta_14C[0],
        Delta_14C[1],
#        kind="cubic",
        bounds_error=False,
#       fill_value=(left_val, -1000.0) # no 14C oustide of dataset
        fill_value="extrapolate"
    )

    F_frac = lambda t: (F_atm(t)/1000+1)*ALPHA_14C
    F_frac_model = lambda t: F_frac(t0 + t/(31*12))

    # construct a 14C model run from the 12C model run

    net_Us_14C = np.array(
        [
            F_frac_model(dmr.times[ti]) * dmr.net_Us[ti] * np.exp(-DECAY_RATE_14C_DAILY * dmr.dt)
            for ti in range(len(dmr.times[:-1]))
        ]
    )

    start_values_14C = compute_start_values_14C(
        dmr,
        nr_time_steps
    )

    dmr_14C = DMR_14C(
        dmr,
        start_values_14C,
        net_Us_14C,
        DECAY_RATE_14C_DAILY
    )

    soln_dmr_14C = dmr_14C.solve()
    return soln_dmr_14C


def compute_us(single_site_dict, time_step_in_days):
    mdo = _load_mdo(single_site_dict, time_step_in_days, check_units=False)
    us = mdo.load_us()
    del mdo

    info = tuple()
    return us, info


def compute_Bs(
    single_site_dict,
    time_step_in_days,
    integration_method='solve_ivp',
    nr_nodes=None,
    check_success=True
):
    mdo = _load_mdo(single_site_dict, time_step_in_days, check_units=False)
    Bs, max_abs_err, max_rel_err = mdo.load_Bs(
        integration_method,
        nr_nodes,
        check_success
    )
    del mdo

    return Bs, (max_abs_err, max_rel_err)


def compute_Us(single_site_dict, time_step_in_days):
    mdo = _load_mdo(single_site_dict, time_step_in_days, check_units=False)
    _, Us, _, _ = mdo.load_xs_Us_Fs_Rs()

    nr_times = len(mdo.time_agg.data)
    nr_pools = mdo.model_structure.nr_pools
    data = np.nan * np.ones((nr_times, nr_pools))

    data[:-1, ...] = Us.data.filled()
    info = tuple()
    return data, info


def compute_Bs_discrete(single_site_dict, time_step_in_days):
    mdo = _load_mdo(single_site_dict, time_step_in_days, check_units=False)

    out = mdo.load_xs_Us_Fs_Rs()
    xs, Us, Fs, Rs = out
    times = mdo.time_agg.data.filled()
    nr_pools = mdo.model_structure.nr_pools

    data = np.nan * np.ones((len(times), nr_pools, nr_pools))

#        try:
    Bs = DMR.reconstruct_Bs(
        xs.data.filled(),
        Fs.data.filled(),
        Rs.data.filled()
    )
        
    data[:-1, ...] = Bs
#        except (DMRError, ValueError, OverflowError) as e:
#            error = str(e)
#            print(error, flush=True)

    info = tuple()
    return data, info


def compute_age_moment_vector_up_to(mr, nr_time_steps, up_to_order):
    start_age_moments = mr.fake_start_age_moments(nr_time_steps, up_to_order)
    res = mr.age_moment_vector_up_to(up_to_order, start_age_moments)
    return res


def compute_pool_age_quantile(mr, nr_time_steps, q, maxsize=None):
    if isinstance(mr, PWCModelRunFD):
        F0 = mr.fake_cumulative_start_age_distribution(nr_time_steps)

        mr.initialize_state_transition_operator_cache(
            maxsize,
            lru_stats=True,
            size=len(mr.times)
        )

        pool_age_quantiles = mr.pool_age_distributions_quantiles(
            quantile=q,
            F0=F0
        )

        cache_stats = mr._state_transition_operator_cache._cached_phi_tmax.cache_stats()
        print([tp[-1] for tp in cache_stats])

        return pool_age_quantiles
    elif isinstance(mr, DMR):
        P0 = mr.fake_cumulative_start_age_masses(nr_time_steps)

        mr.initialize_state_transition_operator_matrix_cache(
            maxsize
        )

        pool_age_quantiles = mr.pool_age_quantiles(
            q,
            P0
        )
        return pool_age_quantiles
    else:
        raise(TypeError("wrong type of model run"))


def compute_system_age_quantile(mr, nr_time_steps, q, maxsize=None):
    if isinstance(mr, PWCModelRunFD):
        F0 = mr.fake_cumulative_start_age_distribution(nr_time_steps)

        mr.initialize_state_transition_operator_cache(
            maxsize,
            size=len(mr.times) // 2
        )

        system_age_quantiles = mr.system_age_distribution_quantiles(
            quantile=q,
            F0=F0
        )
        return system_age_quantiles
    elif isinstance(mr, DMR):
        P0 = mr.fake_cumulative_start_age_masses(nr_time_steps)

        mr.initialize_state_transition_operator_matrix_cache(
            maxsize
        )

        system_age_quantiles = mr.system_age_quantiles(
            q,
            P0
        )
        return system_age_quantiles
    else:
        raise(TypeError("wrong type of model run"))


# needed because propert object PWCModelRunFD.external_output_vector cannot be pickled
def compute_external_output_vector(pwc_mr_fd):
    if not isinstance(pwc_mr_fd, PWCModelRunFD):
        raise(TypeError("wrong type of model run"))

    return pwc_mr_fd.external_output_vector


def compute_acc_net_external_output_vector(dmr):
    if not isinstance(dmr, DMR):
        raise(TypeError("wrong type of model run"))

    data = np.nan * np.ones((len(dmr.times), dmr.nr_pools))
    data[:-1] = dmr.acc_net_external_output_vector()
    return data


def compute_backward_transit_time_moment(pwc_mr_fd, nr_time_steps, order):
    if not isinstance(pwc_mr_fd, PWCModelRunFD):
        raise(TypeError("wrong type of model run"))

    start_age_moments = compute_start_age_moments(pwc_mr_fd, nr_time_steps, order)
    return pwc_mr_fd.backward_transit_time_moment(order, start_age_moments)


def compute_backward_transit_time_quantile(mr, nr_time_steps, q, maxsize=None):
    if isinstance(mr, PWCModelRunFD):
        F0 = mr.fake_cumulative_start_age_distribution(nr_time_steps)

        mr.initialize_state_transition_operator_cache(
            maxsize,
            size=len(mr.times)
        )

        btt_quantiles = mr.backward_transit_time_quantiles(
            q,
            F0
        )
        return btt_quantiles
    elif isinstance(mr, DMR):
        P0 = mr.fake_cumulative_start_age_masses(nr_time_steps)

        mr.initialize_state_transition_operator_matrix_cache(
            maxsize
        )

        data = np.nan * np.ones(len(mr.times))
        btt_quantiles = mr.backward_transit_time_quantiles(
            q,
            P0
        )
        del mr
        data[:-1] = btt_quantiles
        return data
    else:
        raise(TypeError("wrong type of model run"))


def get_complete_sites(z, slices):
    fill_tup = (slice(0, 1, 1), ) * (z.ndim - 3)
    tup = (slices['lat'], slices['lon'], slices['prob']) + fill_tup
    
    if isinstance(z, da.core.Array):
        sliced_da = z[tup]
    else:
        sliced_da = da.from_zarr(z)[tup]

    complete_coords = np.where(sliced_da.compute() != -np.inf)[:3]
    nr_complete_sites = len(complete_coords[0])

    return nr_complete_sites, complete_coords


def get_complete_non_nan_sites(z, slices):
    fill_tup = (slice(0, 1, 1), ) * (z.ndim - 3)
    tup = (slices['lat'], slices['lon'], slices['prob']) + fill_tup
    
    if isinstance(z, da.core.Array):
        sliced_da = z[tup]
    else:
        sliced_da = da.from_zarr(z)[tup]

    arr = sliced_da.compute()
    complete_sliced_non_nan_coords_linear = np.where(
        (arr != -np.inf) & (~np.isnan(arr))
    )[:3]
    
    complete_sliced_non_nan_coords_tuples = _convert_sliced_linear_coords_to_sliced_coords_tuples(
        *complete_sliced_non_nan_coords_linear
    )
    
    complete_non_nan_coords_tuples = _convert_sliced_linear_coords_to_global_coords_tuples(
        *complete_sliced_non_nan_coords_linear, # *(lat, lon, prob)
        slices
    )
    nr_complete_sites = len(complete_non_nan_coords_tuples)

    return nr_complete_sites, complete_non_nan_coords_tuples, complete_sliced_non_nan_coords_tuples


def get_nan_sites(z, slices):
    fill_tup = (slice(0, 1, 1), ) * (z.ndim - 3)
    tup = (slices['lat'], slices['lon'], slices['prob']) + fill_tup
    
    if isinstance(z, da.core.Array):
        sliced_da = z[tup]
    else:
        sliced_da = da.from_zarr(z)[tup]

    arr = sliced_da.compute()
    nan_sliced_coords_linear = np.where(np.isnan(arr))[:3]

    nan_coords_tuples = _convert_sliced_linear_coords_to_global_coords_tuples(
        *nan_sliced_coords_linear, # *(lat, lon, prob)
        slices
    )
    nr_nan_sites = len(nan_coords_tuples)

    return nr_nan_sites, nan_coords_tuples


def get_incomplete_sites(z, slices):
    fill_tup = (slice(0, 1, 1), ) * (z.ndim - 3)
    tup = (slices['lat'], slices['lon'], slices['prob']) + fill_tup
    
    if isinstance(z, da.core.Array):
        sliced_da = z[tup]
    else:
        sliced_da = da.from_zarr(z)[tup]

    incomplete_sliced_coords_linear = np.where(sliced_da.compute() == -np.inf)[:3]
    incomplete_sliced_coords_tuples = _convert_sliced_linear_coords_to_sliced_coords_tuples(
        *incomplete_sliced_coords_linear, # *(lat, lon, prob)
    )

    incomplete_coords_tuples = _convert_sliced_linear_coords_to_global_coords_tuples(
        *incomplete_sliced_coords_linear, # *(lat, lon, prob)
        slices
    )

    nr_incomplete_sites = len(incomplete_coords_tuples)
    return nr_incomplete_sites, incomplete_coords_tuples, incomplete_sliced_coords_tuples


def get_incomplete_site_tuples_for_mr_computation(
    start_values_zarr,
    us_zarr,
    Bs_zarr,
    z,
    slices
):
    _, complete_coords_svs, complete_sliced_coords_svs = get_complete_non_nan_sites(start_values_zarr, slices)
    _, complete_coords_us, complete_sliced_coords_us = get_complete_non_nan_sites(us_zarr, slices)
    _, complete_coords_Bs, complete_sliced_coords_Bs = get_complete_non_nan_sites(Bs_zarr, slices)
    _, incomplete_coords_z, incomplete_sliced_coords_z = get_incomplete_sites(z, slices)
    
    coords_tuples_list = [
        set(complete_coords_svs),
        set(complete_coords_us),
        set(complete_coords_Bs),
        set(incomplete_coords_z)
    ]
                                            
    coords_tuples = list(set.intersection(*coords_tuples_list))
    coords_tuples.sort()

    sliced_coords_tuples_list = [
        set(complete_sliced_coords_svs),
        set(complete_sliced_coords_us),
        set(complete_sliced_coords_Bs),
        set(incomplete_sliced_coords_z)
    ]
                                            
    sliced_coords_tuples = list(set.intersection(*sliced_coords_tuples_list))
    sliced_coords_tuples.sort()

    return len(coords_tuples), coords_tuples, sliced_coords_tuples


def get_nan_site_tuples_for_mr_computation(
    start_values_zarr,
    us_zarr,
    Bs_zarr,
    slices
):
    _, nan_coords_tuples_svs = get_nan_sites(start_values_zarr, slices)
    _, nan_coords_tuples_us = get_nan_sites(us_zarr, slices)
    _, nan_coords_tuples_Bs = get_nan_sites(Bs_zarr, slices)
                        
    nan_coords_tuples_list = [
        set(nan_coords_tuples_svs),
        set(nan_coords_tuples_us),
        set(nan_coords_tuples_Bs)
    ]
    nan_coords_tuples = list(set.union(*nan_coords_tuples_list))
    nan_coords_tuples.sort()

    return len(nan_coords_tuples), nan_coords_tuples


def compute_incomplete_sites(
    time_limit_in_min,
    z,
    nr_times,
    variable_names,
    variables,
    non_data_variables,
    slices,
    task,
    logfile_name,
):
    nr_incomplete_sites, incomplete_coords_tuples, incomplete_sliced_coords_tuples = get_incomplete_sites(z, slices)

    if nr_incomplete_sites > 0:
        print('number of incomplete sites:', nr_incomplete_sites)
    else:
        print('no incomplete sites remaining')
        return

    # select incomplete sites from variables
    incomplete_variables = []
    for v, name in zip(variables, variable_names):
        if name not in non_data_variables:
            v_stack_list = []
            for coords in incomplete_sliced_coords_tuples:
                v_stack_list.append(v[coords])

            incomplete_variables.append(da.stack(v_stack_list))

    # add lat, lon, prob, time
#    incomplete_variables.append(da.from_array(incomplete_coords[0].reshape(-1, 1), chunks=(1, 1))) # lat
    incomplete_variables.append(
        da.from_array(
            np.array([c[0] for c in incomplete_coords_tuples]).reshape(-1, 1),
            chunks=(1, 1)
        )
    )
#    incomplete_variables.append(da.from_array(incomplete_coords[1].reshape(-1, 1), chunks=(1, 1))) # lon
    incomplete_variables.append(
        da.from_array(
            np.array([c[1] for c in incomplete_coords_tuples]).reshape(-1, 1),
            chunks=(1, 1)
        )
    )
#    incomplete_variables.append(da.from_array(incomplete_coords[2].reshape(-1, 1), chunks=(1, 1))) # prob
    incomplete_variables.append(
        da.from_array(
            np.array([c[2] for c in incomplete_coords_tuples]).reshape(-1, 1),
            chunks=(1, 1)
        )
    )
    time_da = variables[variable_names.index('time')].reshape(1, -1).rechunk((1, nr_times))
    incomplete_variables.append(time_da)

    # prepare the delayed computation
    additional_params = {
        "func": task["func"],
        "func_args": task["func_args"],
        "variable_names": variable_names,
        "time_limit_in_min": time_limit_in_min,
        "return_shape": task["return_shape"],
        "logfile_name": logfile_name
    }
    meta_shape = list(task["meta_shape"])
    meta_shape[0] = nr_incomplete_sites
    meta_shape = tuple(meta_shape)

    res_da = incomplete_variables[0].map_blocks(
        func_for_map_blocks,
        *incomplete_variables[1:], # variables[0] comes automatically as first argument
        additional_params,
        drop_axis=task["drop_axis"],
        new_axis=task["new_axis"],
        chunks=task["return_shape"],
        dtype=np.float64,
        meta=np.ndarray(meta_shape, dtype=np.float64)
    )

    # write header to logfile
    print(write_header_to_logfile(logfile_name, res_da, time_limit_in_min))
    print('starting, timeout (min) = ', time_limit_in_min, flush=True)

    # do the computation
    linear_batchwise_to_zarr(
        res_da, # dask array
        z, # target zarr archive
        slices, # slices of interest,
        incomplete_coords_tuples,
        task["batch_size"]
    )

    write_to_logfile(logfile_name, 'done, timeout (min) = '+str(time_limit_in_min))
    print('done, timeout (min) = ', time_limit_in_min, flush=True)


def compute_incomplete_sites_with_mr(
    time_limit_in_min,
    z,
    nr_pools,
    time_step_in_days,
    times_da,
    start_values_zarr,
    us_zarr,
    Bs_zarr,
    slices,
    task,
    logfile_name,
):
    start_values_da = da.from_zarr(start_values_zarr)
    us_da = da.from_zarr(us_zarr)
    Bs_da = da.from_zarr(Bs_zarr)
#    B = Bs_da[7, 60, 10].compute()
#    print(B[:120].mean(axis=0))
    
    # copy nans from start_values, us, or Bs tu z
    nr_nan_sites, nan_coords_tuples = get_nan_site_tuples_for_mr_computation(
        start_values_zarr,
        us_zarr,
        Bs_zarr,
        slices
    )

    for nan_coords in nan_coords_tuples:
        z[nan_coords] = np.nan

    # identify non-nan computed sites in start_values, us, and Bs
    # combine with not yet computed sites in z
    nr_incomplete_sites, incomplete_coords_tuples, _ = get_incomplete_site_tuples_for_mr_computation(
        start_values_zarr,
        us_zarr,
        Bs_zarr,
        z,
        slices
    )

    if nr_incomplete_sites > 0:
        print('number of incomplete sites:', nr_incomplete_sites)
    else:
        print('no incomplete sites remaining')
        return True

    # select incomplete sites from variables
    incomplete_variables = []

    nr_times = len(times_da)
    shapes = [
        (nr_times, nr_pools, nr_pools),
        (1, nr_pools, 1),
        (nr_times, nr_pools, 1)
    ]
    for v, shape in zip([Bs_da, start_values_da, us_da], shapes):
        v_stack_list = []
        for ic in incomplete_coords_tuples:
            v_stack_list.append(v[ic].reshape(shape))

        incomplete_variables.append(da.stack(v_stack_list))

    # add lat, lon, prob
    for k, name in enumerate(["lat", "lon", "prob"]):
        incomplete_variables.append(
            da.from_array(
                np.array(
                    [ic[k] for ic in incomplete_coords_tuples]
                ).reshape(-1, 1, 1, 1),
                chunks=(1, 1, 1, 1)
            )
        )

    # add time
    incomplete_variables.append(times_da.reshape((1, -1, 1, 1)).rechunk((1, nr_times, 1, 1)))

    # prepare the delayed computation
    additional_params = {
        "model_type": task["model_type"],
        "computation": task["computation"],
        "nr_pools": nr_pools,
        "time_step_in_days": time_step_in_days,
        "return_shape": task["return_shape"],
        "func": task["func"],
        "func_args": task["func_args"],
        "time_limit_in_min": time_limit_in_min,
        "logfile_name": logfile_name
    }

    meta_shape = list(task["meta_shape"])
    meta_shape[0] = nr_incomplete_sites
    meta_shape = tuple(meta_shape)

    res_da = incomplete_variables[0].map_blocks(
        func_for_map_blocks_with_mr,
        *incomplete_variables[1:], # variables[0] comes automatically as first argument
        additional_params,
        drop_axis=task["drop_axis"],
        new_axis=task["new_axis"],
        chunks=task["return_shape"],
        dtype=np.float64,
        meta=np.ndarray(meta_shape, dtype=np.float64)
    )

    # write header to logfile
    print(write_header_to_logfile(logfile_name, res_da, time_limit_in_min))
    print('starting, timeout (min) = ', time_limit_in_min, flush=True)

    # do the computation
    linear_batchwise_to_zarr(
        res_da, # dask array
        z, # target zarr archive
        slices, # slices of interest,
        incomplete_coords_tuples,
        task["batch_size"]
    )

    write_to_logfile(logfile_name, 'done, timeout (min) = ' + str(time_limit_in_min))
    print('done, timeout (min) =', time_limit_in_min, flush=True)
    return False


def run_task_with_mr(
    project_path,
    task,
    nr_pools,
    time_step_in_days,
    times_da,
    start_values_zarr,
    us_zarr,
    Bs_zarr,
    slices
):
    print("task: computing", task["computation"])
    print()
            
    zarr_path = Path(project_path.joinpath(task["computation"]))
    print("zarr archive:", str(zarr_path))
    z = load_zarr_archive(
        zarr_path,
        task["result_shape"],
        task["result_chunks"],
        overwrite=task["overwrite"]
    )

#    nr_incomplete_sites, _ = get_incomplete_site_tuples_for_mr_computation(
#        start_values_zarr,
#        us_zarr,
#        Bs_zarr,
#        z,
#        slices
#    )
#    print("Number of incomplete sites:", nr_incomplete_sites)

    logfile_name = str(project_path.joinpath(task["computation"] + ".log"))
    print("Logfile:", logfile_name)

    for timeout in task["timeouts"]:
        done = False
        done = compute_incomplete_sites_with_mr(
            timeout,
            z,
            nr_pools,
            time_step_in_days,
            times_da,
            start_values_zarr,
            us_zarr,
            Bs_zarr,
            slices,
            task,
            logfile_name
        )
        if done:
            break

    if done:
        nr_incomplete_sites = 0
    else:
        nr_incomplete_sites, _, _ = get_incomplete_site_tuples_for_mr_computation(
            start_values_zarr,
            us_zarr,
            Bs_zarr,
            z,
            slices
        )

    write_to_logfile(logfile_name, nr_incomplete_sites, "incomplete sites remaining")
    print(nr_incomplete_sites, "incomplete sites remaining")
    print()


def load_model_structure():
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


###############################################################################
#
# internal methods
#
###############################################################################


def _load_mdo(ds_dict, time_step_in_days, check_units=True): # time step in days
    ms = load_model_structure()

    # no unit support for dictionary version
    time = Variable(
        name="time",
        data=np.arange(len(ds_dict['time'])) * time_step_in_days,
        unit="d"
#        unit="1"
    )

    mdo = ModelDataObject(
        model_structure=ms,
        dataset=ds_dict, 
        stock_unit="gC/m2", 
#        stock_unit="1", 
        time=time,
        check_units=check_units
    )

    return mdo


def _convert_sliced_linear_coords_to_sliced_coords_tuples(lats, lons, probs):
    sliced_coord_tuples = [
        (c[0], c[1], c[2]) for c in zip(lats, lons, probs)
    ]
    return sliced_coord_tuples


def _convert_sliced_linear_coords_to_global_coords_tuples(lats, lons, probs, slices):
    f_lat = lambda x: slices['lat'].start + x * slices['lat'].step
    f_lon = lambda x: slices['lon'].start + x * slices['lon'].step
    f_prob = lambda x: slices['prob'].start + x * slices['prob'].step

    coords_z_lat = [f_lat(x) for x in lats]
    coords_z_lon = [f_lon(x) for x in lons]
    coords_z_prob = [f_prob(x) for x in probs]

    coord_tuples = [
        (c[0], c[1], c[2]) for c in zip(coords_z_lat, coords_z_lon, coords_z_prob)
    ]

    return coord_tuples



