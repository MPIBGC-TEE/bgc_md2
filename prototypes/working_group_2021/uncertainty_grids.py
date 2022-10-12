#!/usr/bin/env python
# ### Computing multi-model ensemble uncertainty from gridded output and writing uncertainty maps to NetCDF

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import netCDF4 as nc
import general_helpers as gh # function repository for the project


# Selecting models and variables
model_names={ # dictionary - model_folder : model_name
    "ab_classic":"CLASSIC",  
    "clm5":"CLM5.0",
    "kv_ft_dlem": "DLEM", 
    "bian_ibis2":"IBIS", 
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
var_names=['cVeg', 'cSoil_total', 'gpp','ra','rh','f_v2s','nep','RT_veg', 
           'RT_soil', 'dist', 'X', 'RT', 'X_c', 'X_p', 'X_c_veg', 'X_p_veg', 
           'X_c_soil', 'X_p_soil']  


# function for computing and writing uncertainty 
def uncertainty_grids (
        model_names,
        experiment_names,
        var_names,
        global_mask,
        output_path,
        ):
    model_folders=[(m) for m in model_names]
    m_names=list(model_names.values())
    g_mask=global_mask.index_mask
    for experiment in experiment_names:
        print('\033[1m'+'. . . Computing uncertainty for '+experiment+' experiment . . .')
        print('\033[0m')
        for vn in var_names:
            # initializing variables
            var_sum_zero=np.zeros_like(g_mask)
            var_sum=np.ma.array(var_sum_zero,mask = g_mask)
            var_diff_sqr_zero=np.zeros_like(g_mask)
            var_diff_sqr=np.ma.array(var_diff_sqr_zero,mask = g_mask)
            var_abs_diff_zero=np.zeros_like(g_mask)
            var_abs_diff=np.ma.array(var_abs_diff_zero,mask = g_mask)                       
            
            ### LOADING DATA AND COMPUTING ENSEMBLE MEAN ###      
            k=0 # model counter        
            for mf in model_folders:
                # reading gata from NetCDF 
                experiment_name=m_names[k]+"_"+experiment+"_"
                conf_dict = gh.confDict(mf)
                dataPath=Path(conf_dict["dataPath"])       
                file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc") 
                ds = nc.Dataset(str(file_path))
                var=ds.variables[vn][:, :].data
                # accumulating sum
                var_sum=var_sum+var
                k+=1 # model counter
                ds.close()                
            # computing mean
            mean=var_sum/len(model_folders) 
            
            ### COMPUTING UNCERTAINTY MEASURES: STANDARD DEVIATION AND AVERAGE DEVIATION ###                  
            k=0 # model counter        
            for mf in model_folders:
                # reading data again
                experiment_name=m_names[k]+"_"+experiment+"_"
                conf_dict = gh.confDict(mf)
                dataPath=Path(conf_dict["dataPath"])                
                file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc")    
                ds = nc.Dataset(str(file_path))
                var=ds.variables[vn][:, :].data           
                # computing differences from the mean
                var_diff_sqr = var_diff_sqr + (var-mean)**2
                var_abs_diff = var_abs_diff + np.abs(var-mean)
                k+=1 # model counter 
                ds.close()             
            variance = var_diff_sqr / (len(model_folders)-1)
            st_dev=np.sqrt(variance) # standard deviation
            ave_dev=var_abs_diff / len(model_folders) # average deviation
            
            ### WRITING UNCERTAINTY MAP IN A NEW NetCDF FILE ###
            # final masking
            var_mean_final = np.ma.array(mean,mask = g_mask)
            var_sd_final = np.ma.array(st_dev,mask = g_mask)
            var_avd_final = np.ma.array(ave_dev,mask = g_mask)
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=Path(output_path).joinpath(experiment+"_"+vn+"_uncertainty.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var_mean = ds_new.createVariable(vn+'_mean', "float32", ["lat", "lon"])
            var_mean[:, :] = var_mean_final
            # standard deviation
            var_sd = ds_new.createVariable(vn+'_sd', "float32", ["lat", "lon"])
            var_sd[:, :] = var_sd_final            
            # average deviation 
            var_avd = ds_new.createVariable(vn+'_avd', "float32", ["lat", "lon"])
            var_avd[:, :] = var_avd_final
            # relative average deviation in %
            var_avd_relative = ds_new.createVariable(vn+'_avd_relative', "float32", ["lat", "lon"])
            relative_uncertainty = var_avd_final / abs(var_mean_final) *100 
            var_avd_relative[:, :] = relative_uncertainty              
            # latitudes and longitudes
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()
            # checking the mean of each grid
            print('\033[1m'+vn+' - '+'\033[0m'+ 
                  'mean: '+str(round(np.ma.mean(var_mean_final),3))+
                  '; std: '+str(round(np.ma.mean(var_sd_final),3))+
                  '; avd: '+str(round(np.ma.mean(var_avd_final),3)))
        
    print('\033[1m'+'Done!')    


# applying the function to selected models, variables and experiments
uncertainty_grids(
    model_names=model_names,
    experiment_names=['S2','S3'],
    var_names=var_names,    
    global_mask=gh.globalMask(file_name="common_mask_all_models.nc")     
    output_path="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\Matrix MIP data\\TRENDY\\Ensemble"
    )


