#!/usr/bin/env python
# %load_ext autoreload
# %autoreload 2
import numpy as np
import general_helpers as gh
import matplotlib.pyplot as plt
from pathlib import Path
import netCDF4 as nc


model_names={
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
model_folders=[(k) for k in model_names]
m_names=list(model_names.values())
experiment_names_S2=list()
for name in m_names:
    experiment_names_S2.append(name + "_S2_") 
experiment_names_S3=list()
for name in m_names:
    experiment_names_S3.append(name + "_S3_") 

global_mask=gh.globalMask(file_name="common_mask.nc")
# gm = nc.Dataset("common_mask_all.nc")
# mask=gm.variables["mask"][:, :].data
# global_mask = gh.CoordMask(
#     index_mask=mask,
#     tr=gh.globalMaskTransformers(mask)
# )
# f = plt.figure(figsize=(20,10))
# ax = f.add_subplot()
# gm.plot_dots(ax)

##################### making temporal averages of nc files
k=0
for mf in model_folders:
    print(mf)
    conf_dict = gh.confDict(mf)
    dataPath=Path(conf_dict["dataPath"])
    model_mask=gh.msh(mf).spatial_mask(dataPath=Path(conf_dict["dataPath"]))    
    experiment_name=experiment_names_S2[k]    
    for vn in gh.msh(mf).data_str._fields:      
        if vn=="npp_nlim": file_path = dataPath.joinpath(gh.msh(mf).nc_file_name("npp", experiment_name=experiment_name))
        else: file_path = dataPath.joinpath(gh.msh(mf).nc_file_name(vn, experiment_name=experiment_name))
        print(file_path)
        ds = nc.Dataset(str(file_path))
        var=ds.variables[vn][:, :, :].data
        if len(var[:,0,0])>500: years_10=10*12 # guessing if data is monthly or yearly
        else: years_10=10
        var_avg=np.ma.mean(var, axis=0)    
        var_first_10=np.ma.mean(var[0:years_10], axis=0) # average of first 10 years
        var_last_10=np.ma.mean(var[-years_10:], axis=0) # average of last 10 years
                
        var_final_avg = np.ma.array(var_avg, mask = model_mask.index_mask)
        var_final_first_10 = np.ma.array(var_first_10, mask = model_mask.index_mask)    
        var_final_last_10 = np.ma.array(var_last_10, mask = model_mask.index_mask)          
              
        s = model_mask.index_mask.shape
        n_lats, n_lons = s
        new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_avg.nc")
        ds_new = nc.Dataset(str(new_path), "w", persist=True)
        
        lat = ds_new.createDimension("lat", size=n_lats)
        lon = ds_new.createDimension("lon", size=n_lons)
                      
        avg = ds_new.createVariable(vn, "float32", ["lat", "lon"])
        first_10 = ds_new.createVariable(str(vn)+"_first_10", "float32", ["lat", "lon"])
        last_10 = ds_new.createVariable(str(vn)+"_last_10", "float32", ["lat", "lon"])        
        avg[:, :] = var_final_avg
        first_10[:, :] = var_final_first_10
        last_10[:, :] = var_final_last_10        

        lats = ds_new.createVariable("lat", "float32", ["lat"])
        lats[:] = list(map(model_mask.tr.i2lat, range(n_lats)))
        lons = ds_new.createVariable("lon", "float32", ["lon"])
        lons[:] = list(map(model_mask.tr.i2lon, range(n_lons)))        
        
        ds.close()        
        ds_new.close()
        #################### same for S3
    experiment_name=experiment_names_S3[k]    
    for vn in gh.msh(mf).data_str._fields:      
        if vn=="npp_nlim": file_path = dataPath.joinpath(gh.msh(mf).nc_file_name("npp", experiment_name=experiment_name))
        else: file_path = dataPath.joinpath(gh.msh(mf).nc_file_name(vn, experiment_name=experiment_name))
        print(file_path)
        ds = nc.Dataset(str(file_path))
        var=ds.variables[vn][:, :, :].data
        if len(var[:,0,0])>500: years_10=10*12 # guessing if data is monthly or yearly
        else: years_10=10
        var_avg=np.ma.mean(var, axis=0)    
        var_first_10=np.ma.mean(var[0:years_10], axis=0) # average of first 10 years
        var_last_10=np.ma.mean(var[-years_10:], axis=0) # average of last 10 years
                
        var_final_avg = np.ma.array(var_avg, mask = model_mask.index_mask)
        var_final_first_10 = np.ma.array(var_first_10, mask = model_mask.index_mask)    
        var_final_last_10 = np.ma.array(var_last_10, mask = model_mask.index_mask)          
              
        s = model_mask.index_mask.shape
        n_lats, n_lons = s
        new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_avg.nc")
        ds_new = nc.Dataset(str(new_path), "w", persist=True)
        
        lat = ds_new.createDimension("lat", size=n_lats)
        lon = ds_new.createDimension("lon", size=n_lons)
                      
        avg = ds_new.createVariable(vn, "float32", ["lat", "lon"])
        first_10 = ds_new.createVariable(str(vn)+"_first_10", "float32", ["lat", "lon"])
        last_10 = ds_new.createVariable(str(vn)+"_last_10", "float32", ["lat", "lon"])        
        avg[:, :] = var_final_avg
        first_10[:, :] = var_final_first_10
        last_10[:, :] = var_final_last_10        

        lats = ds_new.createVariable("lat", "float32", ["lat"])
        lats[:] = list(map(model_mask.tr.i2lat, range(n_lats)))
        lons = ds_new.createVariable("lon", "float32", ["lon"])
        lons[:] = list(map(model_mask.tr.i2lon, range(n_lons)))        
        
        ds.close()        
        ds_new.close()         
    k+=1
print("Done!")

################### resampling (to be completed) ########################
k=0
for mf in model_folders:
    print(mf)
    experiment_name=experiment_names_S2[k]
    conf_dict = gh.confDict(mf)
    dataPath=Path(conf_dict["dataPath"])
    model_mask=gh.msh(mf).spatial_mask(dataPath=Path(conf_dict["dataPath"]))    
    for vn in ["cVeg"]:#gh.msh(mf).data_str._fields:
        print("Resampling "+vn)        
        file_path = dataPath.joinpath(gh.msh(mf).nc_file_name(vn, experiment_name=experiment_name))
        ds = nc.Dataset(str(file_path))
        var=ds.variables[vn][:, :, :].data
#         mask=model_mask.index_mask
#         mask[mask==1]=float("nan")        
        mask = global_mask.index_mask
        zero_array=np.zeros((var.shape[0],mask.shape[0],mask.shape[1]))
        gm=zero_array.copy()
        for i in range(gm.shape[0]):
            gm[i,:,:]=mask 
        var_final=np.ma.array(zero_array,mask = gm)
        print(var_final.shape)
        for i in range(var.shape[0]):
            var_current=var[i,:,:]
        #var_avg=np.mean(var, axis=0)     
        #var_first_10=np.mean(var_arr[0:10], axis=0)
        #var_last_10=np.mean(var_arr[-10:], axis=0)
 
            #var_masked = var_current+mask
            var_masked = np.ma.array(var_current, mask = model_mask.index_mask)    
            var_resampled=gh.resample_grid (
                source_coord_mask=model_mask, 
                target_coord_mask=global_mask, 
                var=var_masked, 
                method="nearest"
            )
            #var_final[i,:,:] = np.ma.array(var_resampled.index_mask, mask = global_mask.index_mask)
            var_final[i,:,:] = var_resampled.index_mask

            if i//100 == i/100:
                print(str(i+1)+" out of "+str(var.shape[0])+" time steps completed")
        print(type(var_final))        
        #print(var_final)
        #print(var_final.shape)
        #print(var_final.shape[0,:,:])
        #s = var_final.shape[0,:,:]
        s = global_mask.index_mask.shape
        n_lats, n_lons = s
        new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_r.nc")
        ds_new = nc.Dataset(str(new_path), "w", persist=True)
        # mask = ds.createDimension('mask',size=s)
        lat = ds_new.createDimension("lat", size=n_lats)
        lon = ds_new.createDimension("lon", size=n_lons)
        
        source_times=ds.variables["time"][:].data        
        time = ds_new.createDimension("time", size=len(source_times))
                
        test = ds_new.createVariable(vn, "float32", ["time", "lat", "lon"])
        test[:, :, :] = var_final

        lats = ds_new.createVariable("lat", "float32", ["lat"])
        lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
        lons = ds_new.createVariable("lon", "float32", ["lon"])
        lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))        
        
        times = ds_new.createVariable ("time", "float32", ["time"])
        times[:] = source_times
        
        ds.close()        
        ds_new.close() 
    k+=1


