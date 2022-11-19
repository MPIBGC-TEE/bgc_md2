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

new_mask=moh.subset_and_resample_nc(
    model_names=model_names,    
    experiment_names=['S2','S3'],
    target_mask=global_mask,
    method="nearest",
    radius_of_influence=500000, 
    )

new_global_mask=gh.globalMask(file_name="common_mask_expanded.nc")    
moh.ensamble_uncertainty (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=new_global_mask,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
        )

plt.rcParams.update(plt.rcParamsDefault)

moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=new_global_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Global"
    )

gh.nc_classes_2_masks(
    FilePath = "C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\global_biomes_geotiff\\biomes_05deg_rasterized.nc", 
    var_name = 'biomes', 
    classes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14], 
    global_mask = new_global_mask)

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_1.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="TUNDRA Permafrost"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_2.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="TUNDRA_Interfrost"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_3.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="BOREAL_Semi-arid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_4.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="BOREAL_Humid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_5.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="TEMPERATE_Semi-arid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_6.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="TEMPERATE_Humid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_7.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2'],#'S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="MEDITERRANEAN_Warm"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_8.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="MEDITERRANEAN_Cold"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_9.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="DESERT_Tropical"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_10.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="DESERT_Temperate"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_11.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="DESERT_Cold"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_12.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="TROPICAL_Semi-arid"
    )

biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_13.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome=" TROPICAL_Humid"
    )

np.mean(a[mid:mid+5])

last5 = a.shape[0]-5
last5

np.mean(a[last5:last5+5])

a[last5+4]

arr = np.array([[11, -1, 15, -7, 12, 14],
               [1,   2,  3,  4,  5,  7]])
arr.flatten()

np.where(arr < 0)

arr[np.where(arr < 0)]=-50
arr

np.mean(arr)
np.std(arr)
arr[np.where(abs(arr) > np.mean(arr)+2*np.std(arr))]=-50
arr

arr1=np.ma.array(data=[(1, 1.0), (2, 2.0)],
             mask=[(0, 0), (1, 0)],)

arr1 

mask=[(0, 1), (1, 0)]
arr1=np.ma.array(arr1,
             mask=mask)
arr1

arr_m=arr1.data[np.where(arr1.mask==0)]
arr_m

# +
#np.random.seed(42)
data = np.random.normal(size=259200)
mask=np.random.randint(2, size=259200)
data_masked = np.ma.array(data,mask=mask)
n, bins, patches = plt.hist(data_masked, density=True, bins=[-3. , -2.5, -2. , -1.5, -1. , -0.5,  0. ,  0.5,  1. ,  1.5,  2. ,
        2.5], facecolor='green', alpha=0.75)
import matplotlib.mlab as mlab
import scipy.stats as stat 
# best fit of data
(mu, sigma) = stat.norm.fit(data_masked)
y = stat.norm.pdf ( bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=2)



# plt.hist(x, density=False, bins=20)  # density=False would make counts
# plt.ylabel('Probability')
# plt.xlabel('Data');
# -

n

bins

patches

# +
# Empirical average and variance are computed
avg = np.mean(data)
var = np.var(data)
# From that, we know the shape of the fitted Gaussian.
pdf_x = np.linspace(np.min(data),np.max(data),100)
pdf_y = 1.0/np.sqrt(2*np.pi*var)*np.exp(-0.5*(pdf_x-avg)**2/var)

# Then we plot :
plt.figure()
plt.hist(data,30,normed=True)
plt.plot(pdf_x,pdf_y,'k--')
#plt.legend(("Fit","Data"),"best")
plt.show()
# -



H, bins =np.histogram(data, bins=10, range=None, normed=None, weights=None, density=None)
H

plt.bar(bins[:-1],H,width=1,color="green")

a=np.array([1,2,3,4,5,6,7,8,9])
a

b=np.array([2,3,4])
b

a[np.where(np.isin(a,b))]=50
a

all_data = [np.random.normal(0, std, 100) for std in range(6, 10)]

all_data.append(np.random.normal(0, 3, 100))

np.array(all_data).shape

all_data=list()

all_data

a=["1","2"]
a

b=a+["3"]
b


