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
    mask=nc.Dataset(dataPath.joinpath("LPJ-GUESS_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
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
    

def make_model_index_transforms():
    return gh.transform_maker(
    lat_0 = -89.75,
    lon_0 = -179.75,
    step_lat = 0.5,
    step_lon = 0.5,
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

data_str = namedtuple( # data streams available in the model
    'data_str',
    ["cVeg", "cLitter", "cSoil", "gpp", "npp", "ra", "rh"]
    )

def get_global_mean_vars_all(experiment_name):
        return(
            gh.get_global_mean_vars_all(model_folder="lpj_guess", 
                            experiment_name=experiment_name,
                            lat_var="latitude",
                            lon_var="longitude",
                            ) 
        )  

################ function for computing global mean for custom data streams ###################
    
# def get_global_mean_vars_all(experiment_name="LPJ-GUESS_S2_"):
    
    # def nc_file_name(nc_var_name, experiment_name="LPJ-GUESS_S2_"):
        # if nc_var_name in ["rh"]: return experiment_name+"{}_annual.nc".format(nc_var_name)
        # else: return experiment_name+"{}.nc".format(nc_var_name)

    # def nc_global_mean_file_name(nc_var_name, experiment_name="LPJ-GUESS_S2_"):
        # return experiment_name+"{}_gm_all.nc".format(nc_var_name)

    # data_str = namedtuple( # data streams available in the model
        # 'data_str',
        # ["cVeg", "cLitter", "cSoil", "gpp", "npp", "ra", "rh"]
        # )
        
    # names = data_str._fields
    # conf_dict = gh.confDict("lpj-guess")
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
        # template = nc.Dataset(dataPath.joinpath("LPJ-GUESS_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
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