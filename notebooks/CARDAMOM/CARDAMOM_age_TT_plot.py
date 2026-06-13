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

# +
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

xr.set_options(display_style='html')
# -


age_output_dir = '/home/data/CARDAMOM/output/age_output/'
#filestem = "cardamom_for_holger_10_ensembles"

ds_age_global_clean = xr.open_mfdataset(age_output_dir + '*.nc')
ds_age_global_clean

# +
fig, ax = plt.subplots()
ds_age_global_clean.mean_age_vector.mean(dim=['ens', 'lat', 'lon']).plot.line(ax=ax, hue='pool')
    
plt.show()
# -

ds_age_global_clean.mean_system_age.mean(dim=['lat', 'lon']).plot()

ds_age_global_clean.mean_system_age.mean(dim=['ens', 'lat', 'lon']).plot()

diff = (ds_age_global_clean.mean_age_vector/ds_age_global_clean.age_standard_deviation_vector)
diff.mean(dim=['ens', 'lat', 'lon']).plot.line(hue='pool')

diff = ds_age_global_clean.mean_system_age/ds_age_global_clean.system_age_standard_deviation
diff.mean(dim=['ens', 'time']).plot()


