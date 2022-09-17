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
    mask=nc.Dataset(dataPath.joinpath("CLASSIC_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
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
    # lat_0 = -87.86380005,
    # lon_0 = -180,
    # step_lat = 2.8125,
    # step_lon = 2.8125,
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

    n_lat = 64
    n_lon = 128
    lat_0 = -87.86380005
    lon_0 = -180
    step_lat=2.789326985714286 #special case n.e. 180/nlatssince
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
        year=1701, 
        month=1,
        day=31
    )

################ function for computing global mean for custom data streams ###################
    
# def get_global_mean_vars_all(experiment_name="CLASSIC_S2_"):
    
    # def nc_file_name(nc_var_name, experiment_name="CLASSIC_S2_"):
        # return experiment_name+"{}.nc".format(nc_var_name)

    # def nc_global_mean_file_name(nc_var_name, experiment_name="CLASSIC_S2_"):
        # return experiment_name+"{}_gm_all.nc".format(nc_var_name)

    # data_str = namedtuple( # data streams available in the model
        # 'data_str',
        # ["cVeg", "cLitter", "cSoil", "gpp", "npp", "ra", "rh"]
        # )
        
    # names = data_str._fields
    # conf_dict = gh.confDict("ab_classic")
    # # with Path('config.json').open(mode='r') as f:
        # # conf_dict = frozendict(json.load(f))
    # dataPath=Path(conf_dict["dataPath"])    
    
    # if all([dataPath.joinpath(nc_global_mean_file_name(vn, experiment_name=experiment_name)).exists() for vn in names]):
        # print(""" Found cached global mean files. If you want to recompute the global means
            # remove the following files: """
        # )
        # for vn in names:
            # print( dataPath.joinpath(nc_global_mean_file_name(vn,experiment_name=experiment_name)))

        # def get_cached_global_mean(vn):
            # gm = gh.get_cached_global_mean(dataPath.joinpath(nc_global_mean_file_name(vn,experiment_name=experiment_name)),vn)
            # return gm * 86400 if vn in ["gpp", "npp", "rh", "ra"] else gm

        # #map variables to data
        # output=gh.data_streams(*map(get_cached_global_mean, data_str._fields))
        # return (
            # output
        # )

    # else:
        # gm=gh.globalMask()
        # # load an example file with mask
        # template = nc.Dataset(dataPath.joinpath("CLASSIC_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
        # gcm=gh.project_2(
                # source=gm,
                # target=gh.CoordMask(
                    # index_mask=np.zeros_like(template),
                    # tr=gh.SymTransformers(
                        # ctr=make_model_coord_transforms(),
                        # itr=make_model_index_transforms()
                    # )
                # )
        # )

        # print("computing means, this may take some minutes...")

        # def compute_and_cache_global_mean(vn):
            # path = dataPath.joinpath(nc_file_name(vn, experiment_name=experiment_name))
            # ds = nc.Dataset(str(path))
            # vs=ds.variables
            # lats= vs["latitude"].__array__()
            # lons= vs["longitude"].__array__()
            # print(vn)
            # var=ds.variables[vn]
            # # check if we have a cached version (which is much faster)
            # gm_path = dataPath.joinpath(nc_global_mean_file_name(vn, experiment_name=experiment_name))

            # gm=gh.global_mean_var(
                    # lats,
                    # lons,
                    # gcm.index_mask,
                    # var
            # )
            # gh.write_global_mean_cache(
                    # gm_path,
                    # gm,
                    # vn
            # )
            # return gm * 86400 if vn in ["gpp", "npp", "rh", "ra"] else gm
        
        # #map variables to data
        # output=data_str(*map(compute_and_cache_global_mean, data_str._fields))
        # return (
            # gh.data_streams( # required data streams
                # cVeg=output.cVeg,
                # cSoil=output.cLitter+output.cSoil,
                # gpp=output.gpp,
                # npp=output.npp,
                # ra=output.ra,
                # rh=output.rh,
            # )
        # )
        

################ function for computing global mean for custom data streams ###################  
    
def get_global_mean_vars_all(experiment_name="CLASSIC_S2_"):
    def nc_file_name(nc_var_name, experiment_name="CLASSIC_S2_"):
        return experiment_name+"{}.nc".format(nc_var_name)

    def nc_global_mean_file_name(experiment_name="CLASSIC_S2_"):
        return experiment_name+"gm_all_vars.nc"

    data_str = namedtuple(
        'data_str',
        ["cVeg", "cLitter", "cSoil", "gpp", "npp", "ra", "rh"]
        )      
    names = data_str._fields
    conf_dict = gh.confDict("ab_classic")
    dataPath=Path(conf_dict["dataPath"])    
    
    if dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name)).exists():
        print(""" Found cached global mean files. If you want to recompute the global means
            remove the following files: """
        )
        print( dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name)))

        def get_cached_global_mean(vn):
            gm = gh.get_cached_global_mean(dataPath.joinpath(
                nc_global_mean_file_name(experiment_name=experiment_name)),vn)
            return gm

        #map variables to data
        output=gh.data_streams(*map(get_cached_global_mean, gh.data_streams._fields))      
        return (
            output
        )

    else:
        gm=gh.globalMask()
        # load an example file with mask
        template = nc.Dataset(dataPath.joinpath("CLASSIC_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
        gcm=gh.project_2(
                source=gm,
                target=gh.CoordMask(
                    index_mask=np.zeros_like(template),
                    tr=gh.SymTransformers(
                        ctr=make_model_coord_transforms(),
                        itr=make_model_index_transforms()
                    )
                )
        )

        print("computing means, this may take some minutes...")

        def compute_and_cache_global_mean(vn):
            path = dataPath.joinpath(nc_file_name(vn, experiment_name=experiment_name))
            ds = nc.Dataset(str(path))
            vs=ds.variables
            lats= vs["latitude"].__array__()
            lons= vs["longitude"].__array__()
            print(vn)
            var=ds.variables[vn]
            # check if we have a cached version (which is much faster)
            gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))

            gm=gh.global_mean_var(
                    lats,
                    lons,
                    gcm.index_mask,
                    var
            )
            return gm * 86400 if vn in ["gpp", "npp", "rh", "ra"] else gm
        
        #map variables to data
        output=data_str(*map(compute_and_cache_global_mean, data_str._fields)) 
        cVeg=output.cVeg if output.cVeg.shape[0]<500 else gh.avg_timeline(output.cVeg, 12)
        cLitter=output.cLitter if output.cLitter.shape[0]<500 else gh.avg_timeline(output.cLitter, 12)
        cSoil=output.cSoil if output.cSoil.shape[0]<500 else gh.avg_timeline(output.cSoil, 12)        
        gpp=output.gpp if output.gpp.shape[0]<500 else gh.avg_timeline(output.gpp, 12)
        npp=output.npp if output.npp.shape[0]<500 else gh.avg_timeline(output.npp, 12)
        ra=output.ra if output.ra.shape[0]<500 else gh.avg_timeline(output.ra, 12)
        rh=output.rh if output.rh.shape[0]<500 else gh.avg_timeline(output.rh, 12)
        output_final=gh.data_streams(
            cVeg=cVeg,
            cSoil=cLitter+cSoil,
            gpp=gpp, 
            npp=npp,
            ra=ra,
            rh=rh,
            )      
        gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))
        gh.write_data_streams_cache(gm_path, output_final)        
        return (
            output_final
        )