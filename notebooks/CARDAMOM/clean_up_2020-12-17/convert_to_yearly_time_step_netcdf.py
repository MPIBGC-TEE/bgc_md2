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

# # Convert CARDAMOM netcdf data to yearly time steps
#
# This notebook loads the CARDAMOM netcdf data files and saves them as a single netcdf file in yearly time steps.

# +
import xarray as xr

from pathlib import Path
from tqdm import tqdm

from bgc_md2.models.CARDAMOM import CARDAMOMlib
from dask.distributed import Client
# -

my_cluster = CARDAMOMlib.prepare_cluster(n_workers=48)
Client(my_cluster)

data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
target_path = data_path.joinpath("yearly_ds.nc")
ds = xr.open_mfdataset(str(data_path) + "/SUM*.nc")

ds_sub = ds.isel(lat=9, lon=26, prob=0)
ds_sub.compute()

nr_months = 12

yearly_coords = {k: v for k,v in ds.coords.items()}
yearly_coords["time"] = ds["time"].isel(time=slice(0, None, nr_months))

# +
stock_variable_names = ["c_finelitter", "c_foliar", "c_labile", "c_root", "c_som", "c_wood"]

data_vars = dict()
for variable_name, variable in ds.data_vars.items():
    if variable_name in stock_variable_names:
        yearly_stock_variable = variable.isel(time=slice(0, None, nr_months))
        data_vars[variable_name] = yearly_stock_variable
    else:
        # flux variable
        yearly_flux_variable = variable.shift(time=-1).coarsen( # first u is useless in CARDAMOM and ELM
            time=nr_months,
            boundary="pad",
            coord_func="min",
            keep_attrs=True
        ).mean().shift(time=1) # unit is gC/m2/d, create useless first u again
        data_vars[variable_name] = yearly_flux_variable

# +
# %%time

yearly_ds.to_netcdf(target_path, compute=True)
# -











check_ds = xr.open_dataset(target_path)
check_ds

sub_ds = check_ds.isel(lat=9, lon=26, prob=0)
sub_ds.compute()

CARDAMOMlib.check_data_consistency(sub_ds.compute(), time_step_in_days=31*12)



# +

variable_name = stock_variable_names[3]
x_da = ds_sub.data_vars[variable_name]
x = x_da.compute().data

fin_da = ds_sub.data_vars["gpp_to_root"]
fin = fin_da.compute().data
fout_da = ds_sub.data_vars["root_to_litter"] + ds_sub.data_vars["fire_root_to_litter"] + ds_sub.data_vars["fire_em_root"]
fout = fout_da.compute().data

# +
# x_da.shift?
# -

n = 8
x[..., n:(n+1)] - (x[..., (n-1):n] + (fin[..., n:(n+1)]*31.-fout[..., n:(n+1)]*31.))
2

yearly_x_da = x_da.isel(time=slice(0, None, 12))
yearly_x = yearly_x_da.compute().data
yearly_x_da

yearly_fin_da = fin_da[..., 1:].coarsen(
    time=12,
    boundary="pad",
#    coord_func="mean",
    keep_attrs=True
).mean().shift(time=1)
yearly_fout_da = fout_da[..., 1:].coarsen(
    time=12,
    boundary="pad",
#    coord_func="mean",
    keep_attrs=True
).mean().shift(time=1)
yearly_fin = yearly_fin_da.compute().data
yearly_fout = yearly_fout_da.compute().data
yearly_fin_da



np.max(np.abs(yearly_x[..., 1:] - (yearly_x[..., :-1] + (yearly_fin-yearly_fout)[..., 1:]*31*12)))

yearly_x.shape


# convert variable to daily data (increasing time dimension by factor of 31)
def c(variable_name, variable):   
    source_zarr_path = source_path.joinpath(variable_name)
    z_source = zarr.open(str(source_zarr_path))

    target_zarr_path = target_path.joinpath(variable_name)
    print("starting", target_zarr_path, flush=True)
    
    target_shape = z_source.shape[:-1] + (z_source.shape[-1]*31,)
    target_chunks = (1, 1, 1, -1)
    z_target = load_zarr_archive(
        target_zarr_path,
        target_shape,
        target_chunks,
        overwrite=True
    )
    
    lat_subs = np.array_split(np.arange(z_source.shape[0]), 5)
    lon_subs = np.array_split(np.arange(z_source.shape[1]), 5)
    prob_subs = np.array_split(np.arange(z_source.shape[2]), 1)
    
    coord_tuples = [(lat, lon, prob) for lat in lat_subs for lon in lon_subs for prob in prob_subs]
    for coord_tuple in tqdm(coord_tuples):
        lat, lon, prob = coord_tuple
        s0 = slice(lat[0], lat[-1]+1, 1)
        s1 = slice(lon[0], lon[-1]+1, 1)
        s2 = slice(prob[0], prob[-1]+1, 1)
        z_target[s0, s1, s2, :] = np.repeat(z_source[s0, s1, s2, :], 31, -1)
                
    print("done", target_zarr_path, flush=True)
    return 1


results = []
for variable_name, variable in ds.data_vars.items():
    y = delayed(c)(variable_name, variable)
    results.append(y)                      
total = delayed(sum)(results)
total.visualize()
# +
# %%time

print(total.compute(), "data variables converted")
# -

for variable_name in ["lat", "lon", "prob"]:
    variable = ds[variable_name]
    zarr_path = target_path.joinpath(variable_name)
    if zarr_path.exists():
        shutil.rmtree(zarr_path)
    da.asarray(variable.data).to_zarr(str(zarr_path))

variable_name = "time"
time = ds["time"]
time_in_days = np.array(time[0], dtype="datetime64[D]") + np.arange(len(time)*31)
time_in_days

zarr_path = target_path.joinpath("time")
if zarr_path.exists():
    shutil.rmtree(zarr_path)
da.from_array(time_in_days).to_zarr(str(zarr_path))


