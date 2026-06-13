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

# The original ncl script `mm_mkinitialmatrices_A_C_1901_to_2004.ncl`
# only uses data of the first day of each year to create only one A per year.
# We produce a completely time dependent A that changes daily but do not 
# store it but use it directly in the computation of B.

from getpass import getuser
import xarray as xr
import numpy as np
import dask
import zarr
from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
import shutil
import matplotlib.pyplot as plt
import bgc_md2.models.cable_all.cableHelpers as cH
from bgc_md2.helper import batchSlices
from time import time
#%load_ext autoreload 
#%autoreload 2


from bgc_md2.sitespecificHelpers import get_client, getCluster



if __name__ == '__main__':
    #cluster = getCluster()
    #client = Client(cluster)
    out_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4'
    arrs = cH.reconstruct_B_u_x0_from_zarr(out_dir) 

    for t in zip(arrs,['B','u','x0']):
        cH.batchwise_to_zarr(t[0], str(zarr_dir_path.joinpath(t[1])))
