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
    #"bian_ibis2":"IBIS",  # exclude IBIS to get results with deserts  
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

#global_mask=gh.globalMask(file_name="common_mask_all_models.nc")    
global_mask=gh.globalMask(file_name="common_mask_all_models_w_deserts.nc")   
# -

# ### Resampling and cleaning model output data

moh.subset_and_resample_nc(
    model_names=model_names,    
    experiment_names=['S2','S3'],
    target_mask=global_mask,
    method="nearest",
    radius_of_influence=500000, 
    )

new_global_mask=gh.globalMask(file_name="common_mask_expanded.nc")  

moh.find_global_outliers (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=new_global_mask,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
        )

final_mask=gh.globalMask(file_name="Final_mask.nc")  

# ### Computing uncertainty

moh.ensamble_uncertainty (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=final_mask,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
        )

moh.ensamble_uncertainty (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=new_global_mask,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
        )

outlier_veg_mask=gh.globalMask(file_name="Outlier_veg_mask.nc") 

moh.ensamble_uncertainty (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=outlier_veg_mask,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_veg",
        )

outlier_soil_mask=gh.globalMask(file_name="Outlier_soil_mask.nc") 

moh.ensamble_uncertainty (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=outlier_soil_mask,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_soil",
        )

# ### Uncertainty attribtution

moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=final_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Global"
    )

moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=final_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Global"
    )

moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=new_global_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="Global"
    )

moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=new_global_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="Global"
    )

moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=outlier_veg_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_veg",
    biome="Veg_Outliers"
    )

moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=outlier_veg_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Veg_Outliers"
    )

moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=outlier_soil_mask,
    data_path="C:\\Users\\KV248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_soil",
    biome="Soil_Outliers"
    )

moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=outlier_soil_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Soil_Outliers"
    )

# ### Biomes

gh.nc_classes_2_masks(
    FilePath = "C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\global_biomes_geotiff\\biomes_05deg_rasterized.nc", 
    var_name = 'biomes', 
    classes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14], 
    global_mask = new_global_mask)

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_1.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TUNDRA Permafrost"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TUNDRA Permafrost"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_2.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TUNDRA_Interfrost"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TUNDRA_Interfrost"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_3.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="BOREAL_Semi-arid"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="BOREAL_Semi-arid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_4.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="BOREAL_Humid"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="BOREAL_Humid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_5.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TEMPERATE_Semi-arid"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TEMPERATE_Semi-arid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_6.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TEMPERATE_Humid"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TEMPERATE_Humid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_7.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="MEDITERRANEAN_Warm"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="MEDITERRANEAN_Warm"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_8.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="MEDITERRANEAN_Cold"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="MEDITERRANEAN_Cold"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_9.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="DESERT_Tropical"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="DESERT_Tropical"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_10.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="DESERT_Temperate"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="DESERT_Temperate"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_11.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="DESERT_Cold"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="DESERT_Cold"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_12.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TROPICAL_Semi-arid"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TROPICAL_Semi-arid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_13.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TROPICAL_Humid"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
    biome="TROPICAL_Humid"
    )

text_displacement=[[2,0.05],[3,-0.05]]
text_displacement

pct_displacement=[[2,0.05],[3,-0.08]]
pct_displacement

for i in text_displacement:
    print(i[1])

cont_var_names = [
            "$NPP_{0}$", "Δ $NPP$",
            "$τ_{veg_0}$", "Δ $τ_{veg}$",
            "$τ_{soil_0}$", "Δ $τ_{soil}$"
            ]   
cont_var_names

uncert_var_name = "$Δ X$" 
uncert_var_name

[uncert_var_name]+cont_var_names

3.7e9

subset="last_5_years"

if subset in ["last_5_years", "diff_total"]:
    print('yo')
else: 
    print('ne')


