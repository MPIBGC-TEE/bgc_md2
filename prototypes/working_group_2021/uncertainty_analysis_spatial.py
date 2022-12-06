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

#global_mask=gh.globalMask(file_name="common_mask_all_models.nc")    
global_mask=gh.globalMask(file_name="common_mask_all_models_w_deserts.nc")   
# -

# ### Resampling and cleaning model output data

moh.subset_and_resample_nc_2(
    model_names=model_names,    
    experiment_names=['S2','S3'],
    target_mask=global_mask,
    method="nearest",
    radius_of_influence=500000, 
    )

moh.mask_inverse_sum(model_names=model_names,
                     experiment_names=['S2','S3'],
                     global_mask=global_mask,
                     output_path=".",
                     threshold=14
                    )

combined_global_mask=gh.globalMask(file_name="combined_global_mask.nc") 
inverse_mask_sum=gh.globalMask(file_name="inv_mask_sum.nc") 

high_outlier_npp, high_outlier_veg, high_outlier_soil=moh.ensamble_uncertainty_2 (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=combined_global_mask,
        inverse_mask_sum=inverse_mask_sum,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
        )

combined_global_mask=gh.globalMask(file_name="mask_corrected.nc") 

print(high_outlier_npp)
print(high_outlier_veg)
print(high_outlier_soil)

# ### Resampling and cleaning model output data

# +
# moh.subset_and_resample_nc(
#     model_names=model_names,    
#     experiment_names=['S2','S3'],
#     target_mask=global_mask,
#     method="nearest",
#     radius_of_influence=500000, 
#     )

# +
#new_global_mask=gh.globalMask(file_name="combined_global_mask.nc")  
# -

moh.find_global_outliers (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=combined_global_mask,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
        high_outlier_npp=high_outlier_npp,
        high_outlier_veg=high_outlier_veg,
        high_outlier_soil=high_outlier_soil
        )

no_outliers_mask=gh.globalMask(file_name="No_outliers_mask.nc")  

moh.ensamble_uncertainty_2 (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=no_outliers_mask,
        inverse_mask_sum=inverse_mask_sum,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_excluded",
        )

outlier_mask=gh.globalMask(file_name="Outlier_mask.nc") 

moh.ensamble_uncertainty (
        model_names=model_names,
        experiment_names=["S2","S3"],
        global_mask=outlier_mask,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_only",
        )

# ### Uncertainty attribution

moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=combined_global_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Global"
    )

moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=combined_global_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Global"
    )

moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=no_outliers_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_excluded",
    biome="Global"
    )

moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=no_outliers_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_excluded",
    biome="Global"
    )

moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=outlier_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_only",
    biome="Outliers"
    )

moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=outlier_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_only",
    biome="Outliers"
    )

# +
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=outlier_soil_mask,
#     data_path="C:\\Users\\KV248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_soil",
#     biome="Soil_Outliers"
#     )

# +
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=outlier_soil_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_soil",
#     biome="Soil_Outliers"
#     )
# -

# ### Biomes

moh.biome_masks_aggregated(
    mask_file = "C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\global_biomes_geotiff\\biomes_05deg_rasterized.nc",
    output_path = ".",
    var_name = 'biomes', 
    classes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14], 
    global_mask = combined_global_mask)

# +
# gh.nc_classes_2_masks(
#     FilePath = "C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\global_biomes_geotiff\\biomes_05deg_rasterized.nc", 
#     var_name = 'biomes', 
#     classes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14], 
#     global_mask = new_global_mask)
# -

biome_mask=gh.globalMask(file_name="mask_tundra.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Tundra"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Tundra"
    )

biome_mask=gh.globalMask(file_name="mask_boreal.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Boreal"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Boreal"
    )

biome_mask=gh.globalMask(file_name="mask_temperate.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Temperate"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Temperate"
    )

biome_mask=gh.globalMask(file_name="mask_mediterranean.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Mediterranean"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Mediterranean"
    )

biome_mask=gh.globalMask(file_name="mask_tropical.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Tropical"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Tropical"
    )

biome_mask=gh.globalMask(file_name="mask_desert.nc")
moh.C_sink_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Desert"
    )
moh.C_storage_uncertainty_attribution(
    experiment_names=['S2','S3'],
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Desert"
    )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_1.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TUNDRA Permafrost"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TUNDRA Permafrost"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_2.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TUNDRA_Interfrost"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TUNDRA_Interfrost"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_3.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="BOREAL_Semi-arid"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="BOREAL_Semi-arid"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_4.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="BOREAL_Humid"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="BOREAL_Humid"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_5.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TEMPERATE_Semi-arid"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TEMPERATE_Semi-arid"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_6.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TEMPERATE_Humid"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TEMPERATE_Humid"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_7.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2'],#'S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="MEDITERRANEAN_Warm"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="MEDITERRANEAN_Warm"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_8.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="MEDITERRANEAN_Cold"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="MEDITERRANEAN_Cold"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_9.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="DESERT_Tropical"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="DESERT_Tropical"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_10.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="DESERT_Temperate"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="DESERT_Temperate"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_11.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="DESERT_Cold"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="DESERT_Cold"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_12.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TROPICAL_Semi-arid"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TROPICAL_Semi-arid"
#     )

# +
# biome_mask=gh.globalMask(file_name="biomes_05deg_rasterized.nc_mask_13.nc")
# moh.C_sink_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TROPICAL_Humid"
#     )
# moh.C_storage_uncertainty_attribution(
#     experiment_names=['S2','S3'],
#     global_mask=biome_mask,
#     data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble\\Outliers_included",
#     biome="TROPICAL_Humid"
#     )
# -

# ### Difference between scenarios / time periods

moh.difference_computation(
    model_names=model_names,
    global_mask=global_mask
)

moh.difference_uncertainty (
        model_names=model_names,
        global_mask=combined_global_mask,
        inverse_mask_sum=inverse_mask_sum,
        output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
        )

moh.C_sink_difference_uncertainty_attribution(
    global_mask=combined_global_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Global"
    )

biome_mask=gh.globalMask(file_name="mask_tundra.nc")
moh.C_sink_difference_uncertainty_attribution(
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Tundra"
    )

biome_mask=gh.globalMask(file_name="mask_boreal.nc")
moh.C_sink_difference_uncertainty_attribution(
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Boreal"
    )

biome_mask=gh.globalMask(file_name="mask_temperate.nc")
moh.C_sink_difference_uncertainty_attribution(
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Temperate"
    )

biome_mask=gh.globalMask(file_name="mask_mediterranean.nc")
moh.C_sink_difference_uncertainty_attribution(
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Mediterranean"
    )

biome_mask=gh.globalMask(file_name="mask_tropical.nc")
moh.C_sink_difference_uncertainty_attribution(
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Tropical"
    )

biome_mask=gh.globalMask(file_name="mask_desert.nc")
moh.C_sink_difference_uncertainty_attribution(
    global_mask=biome_mask,
    data_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    biome="Desert"
    )



moh.get_uncertainty_mask(
    path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble",
    output_path='.',
    global_mask=global_mask
)

# +
import matplotlib.pyplot as plt
import numpy as np
import math

def abc_to_rgb(A=0.0,B=0.0,C=0.0):
    ''' Map values A, B, C (all in domain [0,1]) to
    suitable red, green, blue values.'''
    #return min(A*3,1), min(B*3,1), min(C*3,1)
    return min(A/0.65,1), min(B/0.65,1), min(C/0.65,1)
    #return 1 / max(A,B,C) * A , 1 / max(A,B,C) * B, 1 / max(A,B,C) * C    
    #return (A,B,C)
    #return (min(B+C,1.0),min(A+C,1.0),min(A+B,1.0))
    ''' Plots a legend for the colour scheme
    given by abc_to_rgb. Includes some code adapted
    from http://stackoverflow.com/a/6076050/637562'''

    # Basis vectors for triangle
basis = np.array([[0.0, 1.0], [-1.5/np.sqrt(3), -0.5],[1.5/np.sqrt(3), -0.5]])


    # Plot points
a, b, c = np.mgrid[0.0:1.0:50j, 0.0:1.0:50j, 0.0:1.0:50j]
a, b, c = a.flatten(), b.flatten(), c.flatten()
abc = np.dstack((a,b,c))[0]
#abc = list(filter(lambda x: x[0]+x[1]+x[2]==1, abc)) # remove points outside triangle
abc = list(map(lambda x: x/sum(x), abc)) # or just make sure points lie inside triangle ...
data = np.dot(abc, basis)
colours = [abc_to_rgb(A=point[0],B=point[1],C=point[2]) for point in abc]
fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')

ax.scatter(data[:,0], data[:,1],marker=',',edgecolors='none',facecolors=colours)

    # Plot triangle
ax.plot([basis[_,0] for _ in range(3)],[basis[_,1] for _ in range(3)],**{'color':'black','linewidth':3})
ax.plot([1.5/np.sqrt(3),0],[-0.5,1],**{'color':'black','linewidth':3})
    # Plot labels at vertices
offset = 0.25
fontsize = 32
ax.text(basis[0,0]*(1+offset), basis[0,1]*(1+offset), '$A$', horizontalalignment='center',
        verticalalignment='center', fontsize=fontsize)
ax.text(basis[1,0]*(1+offset), basis[1,1]*(1+offset), '$B$', horizontalalignment='center',
        verticalalignment='center', fontsize=fontsize)
ax.text(basis[2,0]*(1+offset), basis[2,1]*(1+offset), '$C$', horizontalalignment='center',
        verticalalignment='center', fontsize=fontsize)    

ax.set_frame_on(False)
ax.set_xticks(())
ax.set_yticks(())
plt.show()
# -



a=np.log(0.5)

type(a)

a=np.array([1,3,5,7,9])
b=np.array([8,6,4,2,0])
a,b

c=np.maximum(a,b)
c

0.8/0.65


0.8*1.25

np.max((1,2,3))

A=0.8
np.min((A,1))


