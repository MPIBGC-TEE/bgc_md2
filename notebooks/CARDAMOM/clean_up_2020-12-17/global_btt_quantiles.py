# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
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

# # Compute global backward transit time quantiles

# +
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from functools import lru_cache
from pathlib import Path
from tqdm import tqdm

from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
from CompartmentalSystems.helpers_reservoir import generalized_inverse_CDF

from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask.distributed import Client
# -

my_cluster = CARDAMOMlib.prepare_cluster(
    n_workers=10,
    alternative_dashboard_port=8792
)
Client(my_cluster)

# +
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
netCDF_filestem = "sol_acc_age_btt"

dc = ("monthly", None, "discrete")
time_resolution, delay_in_months, model_type = dc
params = CARDAMOMlib.load_params(time_resolution, delay_in_months)
output_path = data_path.joinpath(data_path.joinpath(params["output_folder"]))
project_path = output_path.joinpath(model_type)

ds_path = project_path.joinpath(netCDF_filestem)
print(dc, ds_path)
ds = xr.open_mfdataset(str(ds_path) + "*.nc")
ds

# +
data_vars = {
    "start_values": ds.start_values,
    "Us": ds.Us,
    "Bs": ds.Bs
#    "solution": ds.solution,
#    "acc_net_external_output_vector": ds.acc_net_external_output_vector
}

small_ds = xr.Dataset(
    data_vars=data_vars,
    attrs=ds.attrs
)
small_ds


# -

# sub_ds contains 
# start_values (lat, lon, prob, pool)
# Us (lat, lon, prob, time, pool)
# Bs (lat, lon, prob, time, pool_to, pool_from)
def compute_global_btt_quantile(ds, name, q, nr_time_steps, time_step_in_days, nr_sites=None, sto_cache_size=None):
    # load and adapt area and landfrac datases
    
    ds_area_lf = xr.open_dataset(data_path.joinpath("LSFRAC2.nc"))
    ds_area_lf_adapted = ds_area_lf.sel(lat=ds.lat, lon=ds.lon)
    ds_area_lf_adapted.area_sphere.attrs["units"] = "m^2"
    
    ds_area_lf_adapted.coords["lat"] = np.float64(ds_area_lf_adapted.coords["lat"])
    ds_area_lf_adapted.coords["lon"] = np.float64(ds_area_lf_adapted.coords["lon"])
    
    for var_name, var in ds_area_lf_adapted.data_vars.items():
        ds_area_lf_adapted[var_name].coords["lat"] = np.float64(var.coords["lat"])
        ds_area_lf_adapted[var_name].coords["lon"] = np.float64(var.coords["lon"])

    # collect all sites with useful data
    coords_linear = np.where(
        (~np.isnan(ds.start_values[:, :, 0])\
        & (ds.start_values[:, :, 0] != -np.inf))\
        & (~np.isnan(ds.Bs[:, :, 0, 0, 0]))\
    )

    coords_tuples = [(x, y) for x, y in zip(*coords_linear)]
#    print(len(coords_tuples), "sites available")

    if nr_sites is None:
        nr_sites = len(coords_tuples)
        
    nr_times = len(ds.time)

    F_btt_svs = []

    weights = np.nan * np.ones((nr_sites, nr_times-1))
    total_outflux_C = np.zeros(nr_times-1)
    times = np.arange(len(ds.time)) * time_step_in_days

    if nr_sites is None:
        nr_sites = len(coords_tuples)

    def func_maker(dmr):
        F0 = dmr.fake_cumulative_start_age_masses(nr_time_steps)
        F_sv = dmr.cumulative_pool_age_masses_single_value(F0)
        rho = 1 - dmr.Bs.sum(1)
        F_btt_sv = lambda ai, ti: (rho[ti] * F_sv(ai, ti)).sum() 

        return F_btt_sv

    for nr_coord, coords in enumerate(tqdm(coords_tuples[:nr_sites])):
        lat_idx, lon_idx = coords
        sub_ds = ds.isel(lat=lat_idx, lon=lon_idx)

        # load DMR
        start_values = sub_ds["start_values"].values
        Bs = sub_ds["Bs"].values
        Us = sub_ds["Us"].values

        dmr = DMR.from_Bs_and_net_Us(
            start_values,
            times,
            Bs[:nr_times-1],
            Us[:nr_times-1]
        )

        # create cumulative BTT distribution single value
        F_btt_sv = func_maker(dmr)

        land_sub_ds = ds_area_lf_adapted.sel(
            lat=ds.lat[lat_idx],
            lon=ds.lon[lon_idx]
        )
        weight = np.float(land_sub_ds.area_sphere * land_sub_ds.landfrac) / 1e12

        F_btt_svs.append(F_btt_sv)

        weights[nr_coord, :] = weight
        R = dmr.acc_net_external_output_vector()
        total_outflux_C += weight * R.sum(axis=1)

        # initialize a state transition operator cache for this dmr
        if sto_cache_size:
            dmr.initialize_state_transition_operator_matrix_cache(
                sto_cache_size
            )

    # the global quantile function is a weighted average of the local quantile functions
    # weights are the respective masses of outflux (are*landfrac*R.sum) in each time step
    @lru_cache()
    def F_btt_sv_global(ai, ti):
        res = np.sum(
            np.array(
                [
                    weight * F_btt_sv(int(ai), ti) 
                        for weight, F_btt_sv in zip(weights[:, ti], F_btt_svs)
                ]
            )
        )
        return res

    global_btt_quantiles = np.nan * np.ones(nr_times, dtype=np.float64)

    quantile_ai = 0
    for ti in tqdm(range(len(times)-1)):
        quantile_ai = generalized_inverse_CDF(
            lambda ai: F_btt_sv_global(ai, ti),
            q * total_outflux_C[ti],
            x1=quantile_ai
        )

        if F_btt_sv_global(quantile_ai, ti) > q * total_outflux_C[ti]:
            if quantile_ai > 0:
                quantile_ai = quantile_ai - 1

        # save result for timestep ti
        global_btt_quantiles[ti] = quantile_ai * time_step_in_days
        
    # prepare return dataset
    data_vars = dict()
    data_vars[name] = xr.DataArray(
        data=global_btt_quantiles,
        dims=["time"]
    )
    
    res_ds = xr.Dataset(
        data_vars=data_vars,
        coords = {"time": ds.time.data}    
    )
    sub_ds.close()
    
    return res_ds


# +
chunk_dict = {"prob": 5}

sub_ds = small_ds.isel(prob=slice(0, None, 1), time=slice(0, None, 1)).chunk(chunk_dict)
#sub_ds = ds.chunk(chunk_dict)

q = 0.5
nr_time_steps = 10 * 12

name = "global_btt_median"
fake_data = np.zeros((len(sub_ds.prob), len(sub_ds.time)))

fake_array = xr.DataArray(
    data=fake_data,
    dims=['prob', "time"]
)

fake_coords = {
    "prob": sub_ds.prob.data,
    "time": sub_ds.time.data
}

fake_ds = xr.Dataset(
    data_vars={name: fake_array},
    coords=fake_coords
).chunk(chunk_dict)


# -

def func_chunk(ds, func, **kwargs):
    res = ds.groupby("prob").apply(func, **kwargs)
    return res


def compute_global_btt_quantile_complete_ensemble(sub_ds, name, q, nr_time_steps, sto_cache_size):
    chunk_dict = {"prob": 5}

    sub_ds = small_ds.isel(prob=slice(0, None, 1), time=slice(0, None, 1)).chunk(chunk_dict)

    fake_data = np.zeros((len(sub_ds.prob), len(sub_ds.time)))

    fake_array = xr.DataArray(
        data=fake_data,
        dims=['prob', "time"]
    )

    fake_coords = {
        "prob": sub_ds.prob.data,
        "time": sub_ds.time.data
    }

    fake_ds = xr.Dataset(
        data_vars={name: fake_array},
        coords=fake_coords
    ).chunk(chunk_dict)
    
    ds_res = xr.map_blocks(
        func_chunk,
        sub_ds,
        args=(compute_global_btt_quantile, ),
        kwargs = {
            "name": name,
            "q": q,
            "nr_time_steps": nr_time_steps,
            "time_step_in_days": params["time_step_in_days"],
            "nr_sites": None,
            "sto_cache_size": sto_cache_size
        },
        template=fake_ds
    )

    return ds_res


task_list = [
    {
        "name": "global_btt_median",
        "q": 0.5,
        "sto_cache_size": 7_500
    },
    {
        "name": "global_btt_quantile_05",
        "q": 0.05,
        "sto_cache_size": 10_500
    },
    {
        "name": "global_btt_quantile_95",
        "q": 0.95,
        "sto_cache_size": 7_500
    },
]

# +
# %%time

res_list = []
for task in task_list[2:3]:
    print("starting", task)
    res = compute_global_btt_quantile_complete_ensemble(
        small_ds,
        task["name"],
        task["q"],
        nr_time_steps,
        task["sto_cache_size"]
    )
    current_ds = res.compute()
    for name, var in current_ds.data_vars.items():
        current_ds[name] = var / 31 / 12 # convert age from days to years
        
    current_ds.to_netcdf(project_path.joinpath(task["name"]+".nc"))
    res_list.append(current_ds)
    print("finished", task)
# -

res_ds = xr.merge(res_list)

res_ds.to_netcdf(project_path.joinpath("global_btt_quantiles.nc"))

# ## Load saved single datasets, merge them, and store the merged dataset

res_list = [xr.open_dataset(project_path.joinpath(task["name"]+".nc")) for task in task_list]
xr.merge(res_list).to_netcdf(project_path.joinpath("global_btt_quantiles.nc"))

# # Load merged dataset

ds_btt = xr.open_dataset(project_path.joinpath("global_btt_quantiles.nc"))

# +
fig, axes = plt.subplots(nrows=len(ds_btt.data_vars), figsize=(18, 8*len(ds_btt.data_vars)))

for (name, var), ax in zip(ds_btt.data_vars.items(), axes):
    var.rolling(time=12).mean().plot.line(ax=ax, x="time", add_legend=False, c="blue", alpha=0.2, lw=4)
    var.mean(dim="prob").rolling(time=12).mean().plot(ax=ax, c="red", label="mean")

    ax.set_ylabel("yr")
    ax.set_title(var.name + " (rolling mean over 12 months)")
    ax.set_xlim([ds_btt.time[0], ds_btt.time[-1]])
    ax.legend()

fig.tight_layout()
plt.draw()
# -

# ## Load area and landfrac data

ds_area_lf = xr.open_dataset(data_path.joinpath("LSFRAC2.nc"))

# +
ds_area_lf_adapted = ds_area_lf.sel(lat=ds.lat, lon=ds.lon)
ds_area_lf_adapted.area_sphere.attrs["units"] = "m^2"

ds_area_lf_adapted.coords["lat"] = np.float64(ds_area_lf_adapted.coords["lat"])
ds_area_lf_adapted.coords["lon"] = np.float64(ds_area_lf_adapted.coords["lon"])

for var_name, var in ds_area_lf_adapted.data_vars.items():
    ds_area_lf_adapted[var_name].coords["lat"] = np.float64(var.coords["lat"])
    ds_area_lf_adapted[var_name].coords["lon"] = np.float64(var.coords["lon"])

# +
fig, ax = plt.subplots(figsize=(18, 8))
global_mean_btt_wrong_weights = ds["mean_btt"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])

weights = (ds_area_lf.area_sphere * ds_area_lf.landfrac * ds.acc_net_external_output_vector.sum(dim="pool")).fillna(0)
global_mean_btt = ds["mean_btt"].weighted(weights).mean(dim=["lat", "lon"])

global_mean_btt_wrong_weights.rolling(time=12).mean().plot.line(ax=ax, x="time", c="orange", add_legend=False, alpha=0.2)
global_mean_btt.rolling(time=12).mean().plot.line(ax=ax, x="time", c="green", add_legend=False, alpha=0.2)
ds_btt.global_btt_median.rolling(time=12).mean().plot.line(ax=ax, x="time", c="red", add_legend=False, alpha=0.2)

from matplotlib.lines import Line2D
colors = ["orange", "green", "red"]
lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors]
labels = ["mean BTT (weights missing outflux)", "mean BTT (corrected weights)", "global BTT median"]
_ = ax.legend(lines, labels)
# -

fig, ax = plt.subplots(figsize=(18, 8))
var = np.log(2)*global_mean_btt - ds_btt.global_btt_median
var.rolling(time=12).mean().plot.line(ax=ax, x="time", add_legend=False, c="blue", alpha=0.2)
_ = ax.set_title("[log(2) * (Global mean BTT)] - global btt median")


