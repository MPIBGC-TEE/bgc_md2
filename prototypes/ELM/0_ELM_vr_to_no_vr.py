# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# # Convert ELM vertical resolved data to data accumulated over vertical dimension

# +
import xarray as xr
import numpy as np

from pathlib import Path

#from pint import Quantity as Q_ # has no gC unit
from ACGCA import Q_


# +
data_path = Path("../../../../SOIL-R/Manuscripts/Berkeley/2019/Data/Holger_veg_soil/")
ds = xr.open_mfdataset(str(data_path.joinpath("*.nc")))

target_path = data_path.joinpath("no_vr")
target_path.mkdir(exist_ok=True)

ds
# -

ds_depth = xr.open_dataset("../../../../SOIL-R/Manuscripts/Berkeley/2019/Data/DZSOI.nc")
ds_depth

q_dz = Q_(ds_depth.DZSOI[:, 0].data, "m")
q_dz

# +
# list of necessary variables

var_names = [
    # vegetation component

    "LEAFC_TOT",
    "LIVESTEMC_TOT",
    "DEADSTEMC_TOT",
    "FROOTC_TOT",
    "LIVECROOTC_TOT",
    "DEADCROOTC_TOT",
    
    # soil component
    "CWDC_vr",
    "LITR1C_vr",
    "LITR2C_vr",
    "LITR3C_vr",
    "SOIL1C_vr",
    "SOIL2C_vr",
    "SOIL3C_vr",
    
    # fluxes from outside the model into vegetation
    "PSN_TO_LEAFC_TOT",
    "PSN_TO_LIVESTEMC_TOT",
    "PSN_TO_DEADSTEMC_TOT",
    "PSN_TO_FROOTC_TOT",
    "PSN_TO_LIVECROOTC_TOT",
    "PSN_TO_DEADCROOTC_TOT",
    
    # internal horizontal fluxes
    "LIVESTEMC_TO_DEADSTEMC",
    "LIVECROOTC_TO_DEADCROOTC",

    # vegetation component to soil component
    "LEAFC_TO_LITTER_MET_vr",
    "LEAFC_TO_LITTER_CEL_vr",
    "LEAFC_TO_LITTER_LIG_vr",

    "LIVESTEMC_TO_LITTER_MET_vr",
    "LIVESTEMC_TO_CWDC_vr",

    "DEADSTEMC_TO_LITTER_MET_vr",
    "DEADSTEMC_TO_CWDC_vr",

    "FROOTC_TO_LITTER_MET_vr",
    "FROOTC_TO_LITTER_CEL_vr",
    "FROOTC_TO_LITTER_LIG_vr",

    "LIVECROOTC_TO_LITTER_MET_vr",
    "LIVECROOTC_TO_CWDC_vr",

    "DEADCROOTC_TO_LITTER_MET_vr",
    "DEADCROOTC_TO_CWDC_vr",

    # soil component
    "CWDC_TO_LITR2C_vr",
    "CWDC_TO_LITR3C_vr",
    "LITR1C_TO_SOIL1C_vr",
    "LITR2C_TO_SOIL1C_vr",
    "LITR3C_TO_SOIL2C_vr",
    "SOIL1C_TO_SOIL2C_vr",
    "SOIL1C_TO_SOIL3C_vr",
    "SOIL2C_TO_SOIL1C_vr",
    "SOIL2C_TO_SOIL3C_vr",
    "SOIL3C_TO_SOIL1C_vr",

    # fluxes leaving the model
    "LEAF_MR", "LEAF_GR",
    "LIVESTEM_MR", "LIVESTEM_GR",
    "DEADSTEM_GR",
    "FROOT_MR", "FROOT_GR",
    "LIVECROOT_MR", "LIVECROOT_GR",
    "DEADCROOT_GR",

    # soil component
    "CWD_HR_L2_vr", "CWD_HR_L3_vr", "M_CWDC_TO_FIRE_vr",
    "LITR1_HR_vr", "M_LITR1C_TO_FIRE_vr",
    "LITR2_HR_vr", "M_LITR2C_TO_FIRE_vr",
    "LITR3_HR_vr", "M_LITR3C_TO_FIRE_vr",
    "SOIL1_HR_S2_vr", "SOIL1_HR_S3_vr",
    "SOIL2_HR_S1_vr", "SOIL2_HR_S3_vr",
    "SOIL3_HR_vr"
]

print(len(var_names))

# +
# remove "_vr" by accumulating over depth dimension if necessary

new_data_vars = dict()
for var_name in var_names:
    var = ds[var_name]
    var_name_no_vr = var_name.replace("_vr", "")
    if  var_name_no_vr == var_name:
        print(f"Not vertically resolved: {var_name}")
        new_data_vars[var_name_no_vr] = var
    else:
        if var_name_no_vr in ds.data_vars.keys():
            new_data_vars[var_name_no_vr] = ds[var_name_no_vr]
            print(f"Vertically resolved, already there: {var_name_no_vr}")
        else:
            print(f"{var_name} --> {var_name_no_vr}")
            q_flux_vr = Q_(var.data, var.attrs["units"].replace("m3", "m^3"))
            q_flux = np.einsum("ijkl,j->ikl", q_flux_vr, q_dz).compute()
            
            print(q_flux.units)
            unit_str = str(q_flux.units)
            unit_str = unit_str.replace("g_carbon", "gC")
            unit_str = unit_str.replace("meter ** 2", "m^2")
            unit_str = unit_str.replace("second", "s")
            print(unit_str)
            
            new_data_vars[var_name_no_vr] = xr.DataArray(
                data=q_flux.data,
                dims=["time", "lat", "lon"],
                attrs={"units": unit_str}
            )
        
# -

ds_no_vr = xr.Dataset(new_data_vars)
target_filename = str(target_path.joinpath("global_spinup_yearly.nc"))

ds_no_vr.to_netcdf(target_filename)
Path(target_filename).absolute()

# ## Check success

ds_no_vr = xr.open_dataset(target_filename)
ds_no_vr

for var_name in var_names:
    var_name_no_vr = var_name.replace("_vr", "")
    if not (var_name_no_vr in ds_no_vr.data_vars.keys()):
        print(f"{var_name} --> {var_name_no_vr}")

ds.close()
ds_no_vr.close()


