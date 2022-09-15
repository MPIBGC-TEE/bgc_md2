import sys
import json 
from pathlib import Path
from collections import namedtuple
import netCDF4 as nc
import numpy as np
from sympy import Symbol
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
        TimeStepIterator2,
)
from copy import copy
from typing import Callable
from functools import reduce

sys.path.insert(0,'..') # necessary to import general_helpers
import general_helpers as gh

def spatial_mask(dataPath)->'CoorMask':
    mask=nc.Dataset(dataPath.joinpath("JSBACH_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
    sym_tr= gh.SymTransformers(
        itr=make_model_index_transforms(),
        ctr=make_model_coord_transforms()
    )
    return gh.CoordMask(
        mask,
        sym_tr
    )

def make_model_coord_transforms():
    """ This function can is used to achieve a target grid LAT,LON with
    - LAT ==   0 at the equator with 
    - LAT == -90 at the south pole,
    - LAT== +90 at the north pole,
    - LON ==   0 at Greenich and 
    - LON is counted positive eastwards from -180 to 180
    """
    return gh.CoordTransformers(
            lat2LAT=lambda lat: lat,
            LAT2lat=lambda LAT: LAT,
            lon2LON=lambda lon: lon,
            LON2lon=lambda LON: LON,
    )
    

# def make_model_index_transforms():
    # return gh.transform_maker(
    # lat_0 = -89,
    # lon_0 = -180,
    # step_lat = 1.875,
    # step_lon = 1.875,
 # )
 
def make_model_index_transforms():
    # returns a tuple of functions to describe the grid 
    # index_to_latitude
    # latitude_to_index
    # index_to_longitude
    # longitude_to_index
    # 
    #
    # These function represent the indexing
    # scheme used in the  netcdf4 files.
    # So they are specific to the dataset and model
    # (in this case trendy and yibs )
    # In this case the indexing is somewhat different  
    # from other models (which can use a general_helpers function)
    # since -180 and +180 are among the latitudes
    # These half-pixels are different since there lat value
    # is not in the center but on the boundary

    # The netcdf variables contain lats and lons as
    # arrays (possibly under a different name)
    # if we choose an index i we want:
    # lats[i]==index_to_latitude(i) and
    # latitude_to_index(lats[i])==i
    #
    # lons[i]=index_to_longitude(i) and
    # longitude_to_index[lons[i])==i

    n_lat = 96
    n_lon = 192
    lat_0 = -88.57216851
    lon_0 = -180
    step_lat=1.864677231789474 #special case n.e. 180/nlatssince 
    step_lon=360.0/n_lon
    def i2lat(i_lat):
        if i_lat > (n_lat-1):
            raise IndexError("i_lat > n_lat; with i_lat={}, n_lat={}".format(i_lat,n_lat))
        return lat_0+(step_lat*i_lat)
    
    #def i2lat_min_max(i):
    #    #compute the lat boundaries of pixel i
    #    center=i2lat(i)
    #    lat_min = center if center==-90 else center - step_lat/2 
    #    lat_max= center if center==90 else center + step_lat/2 
    #    return lat_min,lat_max
    #
    #def lat2i(lat):
    #    # the inverse finds the indices of the pixel containing
    #    # the point with the given coordinates
    #    # we cant use round since we want ir=3.5 to be already in pixel 4
    #    ir=(lat-lat_0)/step_lat
    #    ii=int(ir)
    #    d=ir-ii
    #    return ii if d<0.5 else ii+1
    #
    #def i2lon_min_max(i):
    #    #compute the lon boundaries of pixel i
    #    center=i2lon(i)
    #    lon_min = center - step_lon/2 
    #    lon_max=  center + step_lon/2 
    #    return lon_min,lon_max


    def i2lon(i_lon):
        if i_lon > (n_lon-1):
            raise IndexError("i_lon > n_lon; with i_lon={0}, n_lon={1}".format(i_lon,n_lon))
        return lon_0+(step_lon*i_lon)
    
        
    #def lon2i(lon):
    #    # we cant use round since we want ir=3.5 to be already in pixel 4
    #    ir=(lon-lon_0)/step_lon
    #    ii=int(ir)
    #    d=ir-ii
    #    return ii if d<0.5 else ii+1
    return gh.Transformers(
            i2lat=i2lat,
            #i2lat_min_max=i2lat_min_max,
            #lat2i=lat2i,
            i2lon=i2lon,
            #i2lon_min_max=i2lon_min_max,
            #lon2i=lon2i,
        ) 
 

def start_date():
    ## this function is important to syncronise our results
    ## because our data streams start at different times the first day of 
    ## a simulation day_ind=0 refers to different dates for different models
    ## we have to check some assumptions on which this calculation is based
    ## Here is how to get these values
    #ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("VISIT_S2_gpp.nc")))
    #times = ds.variables["time"]
    #tm = times[0] #time of first observation in Months_since_1860-01 # print(times.units)
    #td = int(tm *30)  #in days since_1860-01-01 
    #import datetime as dt
    #ad = dt.date(1, 1, 1) # first of January of year 1 
    #sd = dt.date(1860, 1, 1)
    #td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    ## from td_aD (days since 1-1-1) we can compute the year month and day
    return gh.date(
        year=1700, 
        month=1,
        day=16
    )
