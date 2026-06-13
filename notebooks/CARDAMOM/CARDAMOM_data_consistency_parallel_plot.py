# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import xarray as xr
import matplotlib.pyplot as plt

# +
data_folder = "/home/data/CARDAMOM/"  # matagorda, antakya, raka..
output_folder = "output/"
dataset_name = "Greg_2020_10_26/"

ds = xr.open_dataset(data_folder + dataset_name + output_folder + "data_consistency.nc")

# +
# find nonnull coordinates (land coordinates, basically)
# (just an exercise here)

da_stacked = ds.abs_err.stack(notnull=['lat', 'lon', 'prob'])
da_stacked[da_stacked.notnull()]

# +
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12,6))

ds.abs_err.mean(dim='prob').plot(
    ax=ax1,
    cbar_kwargs={"label": '$gC/m^2$'},
    robust=True
)
ax1.set_title('mean absolute error')

ds.rel_err.mean(dim='prob').plot(
    ax=ax2,
    cbar_kwargs={"label": '%'},
    robust=True
)
ax2.set_title('mean relative error')

plt.suptitle('CARDAMOM - data consistency (robust version)')

plt.tight_layout()
plt.draw()
# -

ds.close()
del ds
