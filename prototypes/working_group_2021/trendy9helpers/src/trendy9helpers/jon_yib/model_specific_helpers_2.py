import sys
import json 
from pathlib import Path
from collections import namedtuple 
from importlib import import_module
import netCDF4 as nc
import numpy as np
from sympy import Symbol
from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.ArrayDictResult import ArrayDictResult
from copy import copy
from typing import Callable, Tuple, Dict
from functools import partial, reduce

model_mod='bgc_md2.models.jon_yib'
cp_mod=import_module(f"{model_mod}.CachedParameterization")
make_func_dict=import_module(f"{model_mod}.CachedParameterization").make_func_dict
Drivers=cp_mod.Drivers
CachedParameterization=cp_mod.CachedParameterization

from .. import general_helpers as gh

def spatial_mask(dataPath)->'CoorMask':
    # the yibs data set has 
    # 1.) ONE file with a mask YIBs_S2_Monthly_tas.nc
    # 2.) missing masks for all other files
    # 3.) trajectories that brake of after some time and are filled wiht NAN
    #     without beeing masked
    # We therefore create a mask by checking for the NANs
    # we now check if any of the arrays has a time lime containing nan values 
    # APART FROM values that are already masked by the fillvalue
    
    # 1.)
    f_mask=nc.Dataset(dataPath.joinpath("YIBs_S2_Monthly_tas.nc")).variables['tas'][0,:,:].mask
    
    # 2.) 
    print("computing masks to exclude pixels with nan entries, this may take some minutes...")
    
    def f(vn):
        path = dataPath.joinpath(nc_file_name(vn))
        ds = nc.Dataset(str(path))
        var =ds.variables[vn]
        ##return after assessing NaN data values
        return gh.get_nan_pixel_mask(var)

    o_names=Observables._fields
    d_names=Drivers._fields
    names = o_names + d_names 

    masks=[ f(name)    for name in names ]
    # We compute the common mask so that it yields valid pixels for ALL variables 
    combined_mask= reduce(lambda acc,m: np.logical_or(acc,m),masks,f_mask)
    
    sym_tr= gh.SymTransformers(
        itr=make_model_index_transforms(),
        ctr=make_model_coord_transforms()
    )
    return gh.CoordMask(
        combined_mask,
        sym_tr
    )

def make_model_coord_transforms():
    return gh.identicalTransformers()

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

    n_lat = 181
    n_lon = 360
    lat_0 = -90
    lon_0 = -179.5
    step_lat=1.0 #special case n.e. 180/nlatssince yibs includes -180 and 180
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

        

Observables_annual = namedtuple(
    'Observables_annual',
    ["cVeg", "cSoil"]
)
Observables_monthly = namedtuple(
    'Observables_monthly',
    ["rh"] #,"ra"]
)
Observables = namedtuple(
    "Observables",
    Observables_annual._fields+Observables_monthly._fields
)



# when downloading data make sure model names match TRENDY server names:
# "CABLE-POP","CLASSIC","CLM5","DLEM","IBIS","ISAM","ISBA_CTRIP",
# "JSBACH","JULES-ES-1.0","LPJ-GUESS","LPJwsl","LPX-Bern","OCN",
# "ORCHIDEE","ORCHIDEE-CNP","ORCHIDEEv3","ORCHIDEEv3_0.5deg"
# "SDGVM","VISIT","YIBs"

#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['YIBs'],
        variables =Observables._fields + Drivers._fields
    )

def get_example_site_vars(dataPath):
    # pick up 1 site   
    s = slice(None, None, None)  # this is the same as :
    t = s, 74, 118  # [t] = [:,49,325]
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them 
        ds = nc.Dataset(str(path))
        #check for npp/gpp/rh/ra to convert from kg/m2/s to kg/m2/day
        if vn in ["npp","gpp","rh","ra"]:
            #for name, variable in ds.variables.items():            
            #    for attrname in variable.ncattrs():
            #        print("{} -- {}".format(attrname, getattr(variable, attrname)))
            return (ds.variables[vn][t]*24*60*60)
        else:
            #for name, variable in ds.variables.items():            
            #    for attrname in variable.ncattrs():
            #        print("{} -- {}".format(attrname, getattr(variable, attrname)))
            return ds.variables[vn][t]

    # Link symbols and data:
    # YIBS has annual vs monthly file names so they are linked separately
    # If all your data is similarly named you can do this in one step

    # Create annual file names (single step if files similarly named)
    o_names=[(f,"YIBs_S2_Annual_{}.nc".format(f)) for f in Observables_annual._fields]

    # Create monthly file names (can remove if done in one step above)
    monthly_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Observables_monthly._fields]
    # Extend name list with monthly names
    o_names.extend(monthly_names)

    # create file names for Drivers
    d_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Drivers._fields]

    # Link symbols and data for Observables/Drivers
    return (Observables(*map(f, o_names)),Drivers(*map(f,d_names)))


# +
# deprecated because it uses an already deprecated global_mean_JULES function

def nc_file_name(nc_var_name, experiment_name="YIBs_S2_"):
    if nc_var_name in ["cVeg", "cSoil"]:
        return experiment_name+"Annual_{}.nc".format(nc_var_name)
    else:
        return experiment_name+"Monthly_{}.nc".format(nc_var_name)


# +
def nc_global_mean_file_name(nc_var_name, experiment_name="YIBs_S2_"):
    return experiment_name+"{}_gm.nc".format(nc_var_name)

def nc_clip_file_name(nc_var_name):
    return experiment_name+"{}_clipped.nc".format(nc_var_name)


def compute_global_mean_arr_var_dict(dataPath):
    ds = nc.Dataset(dataPath.joinpath("YIBs_S2_Monthly_tas.nc"))
    vs=ds.variables
    lats= vs["latitude"].__array__()
    lons= vs["longitude"].__array__()
    template =vs['tas'][0,:,:].mask
    def var(vn):
        return nc.Dataset(
                str(
                    dataPath.joinpath(
                        nc_file_name(vn)
                    )
                 )
             ).variables[vn]

    
    def gm_func (var,time_slice=slice(None,None,None)):
        return gh.global_mean_var_with_resampled_mask(
            template=template,
            ctr=make_model_coord_transforms(), 
            itr=make_model_index_transforms(),
            lats=lats,
            lons=lons,
            var=var,
            time_slice=time_slice
        )

    
    scaled=["npp","gpp","rh","ra"]
    arr_dict = {
        **{vn: gm_func(var(vn) )
            for vn in set(Observables._fields + Drivers._fields).difference(scaled)
        }, 
        **{vn: gm_func(var(vn))*86400 # kg/m2/s kg/m2/day;
            for vn in scaled 
        } 
    }
    #from IPython import embed;embed()
    return arr_dict

def get_global_mean_vars(dataPath, targetPath=None, flash_cache=False):
    if targetPath is None:
        targetPath = dataPath

    arr_dict= gh.cached_var_dict(
        dataPath,
        targetPath,
        nc_global_mean_file_name,
        compute_global_mean_arr_var_dict,
        names=Observables._fields + Drivers._fields,
        flash_cache=flash_cache
    )
    obs = Observables(*(arr_dict[k] for k in Observables._fields))
    dvs = Drivers(*(arr_dict[k] for k in Drivers._fields))
    
    return (obs,dvs)

def get_global_mean_vars_2(conf_dict,targetPath=None):
    o_names = Observables._fields
    d_names = Drivers._fields
    names = o_names + d_names
    dataPath = Path(conf_dict["dataPath"])
    if not(all([dataPath.joinpath(nc_file_name(vn)).exists() for vn in names])):
        download_my_TRENDY_output(conf_dict)
    return get_global_mean_vars(dataPath,targetPath)    

def make_da_iterator(
        mvs,
        X_0, #: StartVector,
        par_dict,
        func_dict,
        delta_t_val=1 # defaults to 1day timestep
    ):
    # this function has to be modelspecific but because some models 
    # don't fall in the general case that we can use here
    return gh.rh_iterator(
            mvs,
            X_0,
            par_dict,
            func_dict,
            delta_t_val
    )
    return mit




def make_weighted_cost_func(
        obs: Observables
    ) -> Callable[[Observables],np.float64]:
    
    # Define cost function 
    def costfunction(out_simu: Observables) -> np.float64:
        stream_cost = partial(gh.single_stream_cost,obs,out_simu)
        
        J_new = (
            stream_cost("cVeg")
            + stream_cost("cSoil")
            + stream_cost("rh")
        ) # + J_obj4)

        return J_new
   
    return costfunction


def start_date():
    ## this function is important to syncronise our results
    ## because our data streams start at different times the first day of 
    ## a simulation day_ind=0 refers to different dates for different models
    ## we have to check some assumptions on which this calculation is based
    ## Here is how to get these values
    #ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("YIBs_S2_Monthly_npp.nc")))
    #times = ds.variables["time"]
    #tm = times[0] #time of first observation in Months_since_1860-01 # print(times.units)
    #td = int(tm *30)  #in days since_1700-01-01 
    #import datetime as dt
    #ad = dt.date(1, 1, 1) # first of January of year 1 
    #sd = dt.date(1700, 1, 1)
    #td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    ## from td_aD (days since 1-1-1) we can compute the year month and day
    return gh.date(
        year=1700, 
        month=1,
        day=1
    )

def make_sim_day_2_day_since_a_D(conf_dict):
    # this function is extremely important to syncronise our results
    # because our data streams start at different times the first day of 
    # a simulation day_ind=0 refers to different dates for different models
    # we have to check some assumptions on which this calculation is based
    # for jules the data points are actually spaced monthly with different numbers of days
    ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("YIBs_S2_Monthly_gpp.nc")))
    times = ds.variables["time"]

    # we have to check some assumptions on which this calculation is based
    tm = times[0] #time of first observation in Months_since_1860-01 # print(times.units)
    td = int(tm *31)  #in days since_1700-01-01 
    #NOT assuming a 30 day month...
    import datetime as dt
    ad = dt.date(1, 1, 1) # first of January of year 1 
    sd = dt.date(1700, 1, 1)
    td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    
    def f(day_ind: int)->int:
        return day_ind+td_aD

    return f



data_str = namedtuple( # data streams available in the model
    'data_str',
    ["cVeg", "cSoil", "gpp", "npp", "ra", "rh"]
    )
    

# fixme mm 7-4 2023:
# model_folder="..." is self referential and would
# break if we rename the model folder. This is a design flaw (anty pattern)
# If there is any model specific information it should be computed
# by model specific functions here and included transmitted in the function call
def get_global_mean_vars_all(experiment_name):
        return(
            gh.get_global_mean_vars_all(model_folder="jon_yib", 
                            experiment_name=experiment_name,
                            lat_var="latitude",
                            lon_var="longitude",
                            ) 
        )       


