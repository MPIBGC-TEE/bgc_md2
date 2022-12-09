#!/usr/bin/env python
# ### Computing multi-model ensemble uncertainty from gridded output and writing uncertainty maps to NetCDF

from osgeo import gdal  # run this notebook in a separate environment which includes gdal package
from pathlib import Path
import netCDF4 as nc
import numpy as np

def TIFF_2_NetCDF (FilePath, NewPath, VarName='var'): 
    ds = gdal.Open(FilePath)
    a = ds.ReadAsArray()
    nlat, nlon = np.shape(a)
    b = ds.GetGeoTransform()  # bbox, interval
    # creating and writing a new NetCDF file 
    ds_new = nc.Dataset(str(NewPath), "w", persist=True)
    # creating dimensions
    lat = ds_new.createDimension("lat", size=nlat)
    lon = ds_new.createDimension("lon", size=nlon)         
    # creating variables                         
    var = ds_new.createVariable(VarName, "float32", ["lat", "lon"])
    var[:, :] = a               
    lats = ds_new.createVariable("lat", "float32", ["lat"])
    lats[:] = np.arange(nlat)*b[5]+b[3]
    lons = ds_new.createVariable("lon", "float32", ["lon"])
    lons[:] = np.arange(nlon)*b[1]+b[0]   
    # closing NetCDF file      
    ds_new.close() 
    print('File written as '+NewPath)


# +
# TIFF_2_NetCDF(
#     FilePath="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\global_biomes_geotiff\\biomes_05deg_rasterized.tif",
#     NewPath="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\global_biomes_geotiff\\biomes_05deg_rasterized.nc",
#     VarName='biomes'
# )
# -

TIFF_2_NetCDF(
    FilePath="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\global_biomes_geotiff\\wwf_biomes\\wwf_biomes_rasterized2.tif",
    NewPath="C:\\Users\\kv248\\OneDrive - Cornell University\\Data\\global_biomes_geotiff\\wwf_biomes\\wwf_biomes_rasterized2.nc",
    VarName='biomes'
)
