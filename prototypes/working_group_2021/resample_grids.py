#!/usr/bin/env python
# %load_ext autoreload
# %autoreload 2
import numpy as np
import general_helpers as gh
import matplotlib.pyplot as plt
from pathlib import Path
import netCDF4 as nc


# +
model_names={
    "ab_classic":"CLASSIC",  
    "clm5":"CLM5.0",
    "kv_ft_dlem": "DLEM", 
    "bian_ibis2":"IBIS",  # exclude IBIS to get results with deserts  
    "cj_isam": "ISAM",    
    "isba-ctrip":"ISBA-CTRIP",    
    "jsbach":"JSBACH",
    "yz_jules": "JULES-ES-1p0",    
    "lpj-guess":"LPJ-GUESS",
    "lpjwsl":"LPJ",
    "lpx-bern":"LPX-Bern",
    "ORCHIDEE-V2":"OCN",    
    "ORCHIDEE":"ORCHIDEE",
    "ORCHIDEE-CNP":"ORCHIDEE-CNP",    
    "ORCHIDEEv3":"ORCHIDEEv3",
    "Aneesh_SDGVM":"SDGVM",
    "kv_visit2": "VISIT",
    "jon_yib": "YIBs"    
}

global_mask=gh.globalMask(file_name="common_mask_all_models.nc")    
#global_mask=gh.globalMask(file_name="common_mask_all_models_w_deserts.nc")   
# -

# Compute temporal average + difference(last-first) and resample all NetCDF4 data streams for each model to the global mask
gh.average_and_resample_nc(
    model_names=model_names,
    experiment_names=['S2','S3'],
    target_mask=global_mask,
    method="nearest",
    radius_of_influence=500000, 
    )

# +
# Resample all NetCDF4 data streams for each model to the global mask (takes time and creates large .nc files)
# gh.resample_nc(
#     model_names=model_names,
#     experiment_names=['S2'],#,'S3'],
#     target_mask=global_mask,
#     method="nearest",
#     radius_of_influence=500000, 
#     )
# -

gh.add_gridded_vars (
        model_names=model_names,
        experiment_names=['S2','S3'],
        global_mask=global_mask,       
        )

gh.uncertainty_grids(
    model_names=model_names,
    experiment_names=['S2','S3'],
    global_mask=global_mask,
    output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

gh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2','S3'],
    global_mask=global_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

a=[1,2,3,4,5]
start=len(a)-3
a[-1]-a[start]

# +
#change
# -



