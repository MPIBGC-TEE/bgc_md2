#!/usr/bin/env python
# %load_ext autoreload
# %autoreload 2
import numpy as np
import general_helpers as gh
import MIP_output_helpers as moh
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

from pyproj import Geod, Proj, transform

from pyproj.aoi import AreaOfUse

# Compute temporal average + difference(last-first) and resample all NetCDF4 data streams for each model to the global mask
moh.average_and_resample_nc(
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

moh.add_gridded_vars (
        model_names=model_names,
        experiment_names=['S2','S3'],
        global_mask=global_mask,       
        )

moh.uncertainty_grids(
    model_names=model_names,
    experiment_names=['S2','S3'],
    global_mask=global_mask,
    output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

plt.rcParams.update({'font.size': 15})
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=global_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

gh.nc_classes_2_masks(
    FilePath = "C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\global_biomes_geotiff\\biomes_05deg_rasterized.nc", 
    var_name = 'biomes', 
    classes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14], 
    global_mask = global_mask)

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_1.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_2.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_3.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_4.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_5.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_6.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_7.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_8.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_9.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_10.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_11.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_12.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_13.nc")
moh.grid_attribution(
    #model_names=model_names,
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )


