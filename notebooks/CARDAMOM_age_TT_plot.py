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

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd


age_output_dir = 'age_output/'
filestem = "cardamom_for_holger_10_ensembles"

ds_age_global_clean = xr.open_mfdataset('age_output/*.nc')
ds_age_global_clean
