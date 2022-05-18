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
Drivers=namedtuple(
    "Drivers",
    ["npp","tas","gpp"]
)

Constants = namedtuple(
    "Constants",
    [
        'npp_0',       #Initial input/pools
        'rh_0',
        #'ra_0',
        'c_veg_0',
        'c_soil_0',
        'clay',        #Constants like clay
        'silt',
        'nyears',       #Run time (years for my model)        
    ]
)
EstimatedParameters = namedtuple(
    'EstimatedParameters', 
    [
        'beta_leaf',
        'beta_root',
        #'r_c_leaf_rh',
        #'r_c_root_rh',
        #'r_c_wood_rh',
        'r_c_leaf_2_c_lit_met',
        'r_c_leaf_2_c_lit_str',
        'r_c_root_2_c_soil_met',
        'r_c_root_2_c_soil_str',
        'r_c_wood_2_c_lit_cwd',
        'c_leaf_0',               
        'c_root_0',
        'r_c_lit_cwd_rh',
        'r_c_lit_met_rh',
        'r_c_lit_str_rh',
        'r_c_lit_mic_rh',
        'r_c_soil_met_rh',
        'r_c_soil_str_rh',
        'r_c_soil_mic_rh',
        'r_c_soil_slow_rh',
        'r_c_soil_passive_rh',
        'r_c_lit_cwd_2_c_lit_mic',
        'r_c_lit_cwd_2_c_soil_slow',
        'r_c_lit_met_2_c_lit_mic',
        'r_c_lit_str_2_c_lit_mic',
        'r_c_lit_str_2_c_soil_slow',
        'r_c_lit_mic_2_c_soil_slow',
        'r_c_soil_met_2_c_soil_mic',
        'r_c_soil_str_2_c_soil_mic',
        'r_c_soil_str_2_c_soil_slow',
        'r_c_soil_mic_2_c_soil_slow',
        'r_c_soil_mic_2_c_soil_passive',
        'r_c_soil_slow_2_c_soil_mic',
        'r_c_soil_slow_2_c_soil_passive',
        'r_c_soil_passive_2_c_soil_mic',
        'c_lit_cwd_0',
        'c_lit_met_0',
        'c_lit_str_0',
        'c_lit_mic_0',
        'c_soil_met_0',
        'c_soil_str_0',
        'c_soil_mic_0',
        'c_soil_slow_0',
    ]
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
# def get_global_mean_vars(dataPath):
    
#     #define function to average variables
#     def f(tup):
#         #define parts of function from nc file
#         vn, fn = tup
#         path = dataPath.joinpath(fn)
#         ds = nc.Dataset(str(path))
#         lats = ds.variables["latitude"]
#         lons = ds.variables["longitude"]
        
#         #check for npp/gpp/rh/ra to convert from kg/m2/s to kg/m2/day
#         if vn in ["npp","gpp","rh","ra"]:
#             #for name, variable in ds.variables.items():            
#             #    for attrname in variable.ncattrs():
#             #        print("{} -- {}".format(attrname, getattr(variable, attrname)))
#             return (gh.global_mean_JULES(lats, lons, ds.variables[vn])*24*60*60)
#         else:
#             #for name, variable in ds.variables.items():            
#             #    for attrname in variable.ncattrs():
#             #        print("{} -- {}".format(attrname, getattr(variable, attrname)))
#             return (gh.global_mean_JULES(lats, lons, ds.variables[vn]))

#     # Link symbols and data:
#     # YIBS has annual vs monthly file names so they are linked separately
#     # If all your data is similarly named you can do this in one step

#     # Create annual file names (single step if files similarly named)
#     o_names=[(f,"YIBs_S2_Annual_{}.nc".format(f)) for f in Observables_annual._fields]

#     # Create monthly file names (can remove if done in one step above)
#     monthly_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Observables_monthly._fields]
#     # Extend name list with monthly names
#     o_names.extend(monthly_names)

#     # create file names for Drivers
#     d_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Drivers._fields]

#     # Link symbols and data for Observables/Drivers
#     return (Observables(*map(f, o_names)),Drivers(*map(f,d_names)))
# -

experiment_name="YIBs_S2_"
def nc_file_name(nc_var_name):
    if nc_var_name in ["cVeg", "cSoil"]:
        return experiment_name+"Annual_{}.nc".format(nc_var_name)
    else:
        return experiment_name+"Monthly_{}.nc".format(nc_var_name)


# +
def nc_global_mean_file_name(nc_var_name):
    return experiment_name+"{}_gm.nc".format(nc_var_name)

def nc_clip_file_name(nc_var_name):
    return experiment_name+"{}_clipped.nc".format(nc_var_name)


# -

def get_global_mean_vars(dataPath):
    # Define function to select geospatial cell and scale data
    #def f(tup):
    #    vn, fn = tup
    #    path = dataPath.joinpath(fn)
    #    # Read NetCDF data but only at the point where we want them
    #    ds = nc.Dataset(str(path))
    #    lats = ds.variables["latitude"].__array__()
    #    lons = ds.variables["longitude"].__array__()
    #    
    #    #check for npp/gpp/rh/ra to convert from kg/m2/s to kg/m2/day
    #    if vn in ["npp","gpp","rh","ra"]:
    #        return (gh.global_mean(lats, lons, ds.variables[vn].__array__())*24*60*60)
    #    else:
    #        return (gh.global_mean(lats, lons, ds.variables[vn].__array__()))
    #
    ## Link symbols and data:
    #
    ## Create annual file names (single step if files similarly named)
    #o_names=[(f,"YIBs_S2_Annual_{}.nc".format(f)) for f in Observables_annual._fields]
    #
    ## Create monthly file names (can remove if done in one step above)
    #monthly_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Observables_monthly._fields]
    ## Extend name list with monthly names
    #o_names.extend(monthly_names)
    #
    ## create file names for Drivers
    #d_names=[(f,"YIBs_S2_Monthly_{}.nc".format(f)) for f in Drivers._fields]
    #
    ## Link symbols and data for observables/drivers
    ## print(o_tuples)
    #return (
    #    Observables(*map(f, o_names)),
    #    Drivers(*map(f, d_names))
    #)
    o_names=Observables._fields
    d_names=Drivers._fields
    names = o_names + d_names 

    if all([dataPath.joinpath(nc_global_mean_file_name(vn)).exists() for vn in names]):
        print(""" Found cached global mean files. If you want to recompute the global means
            remove the following files: """
        )
        for vn in names:
            print( dataPath.joinpath(nc_global_mean_file_name(vn)))

        def get_cached_global_mean(vn):
            return gh.get_cached_global_mean(dataPath.joinpath(nc_global_mean_file_name(vn)),vn)
    
        return (
            Observables(*map(get_cached_global_mean, o_names)),
            Drivers(*map(get_cached_global_mean,d_names))
        )

    else:
        # we now check if any of the arrays has a time lime containing nan values 
        # APART FROM values that are already masked by the fillvalue
        print("computing masks to exclude pixels with nan entries, this may take some minutes...")
        def f(vn):
            path = dataPath.joinpath(nc_file_name(vn))
            #ds = nc.Dataset(str(path))
            #scale fluxes vs pools
            #var =ds.variables[vn]
            
            ########### apply shape file mask of the world to each data file
            #import needed librarys - these have to be added to conda env
            #use "conda install --channel conda-forge geopandas rioxarray xarray shapely gdal=3.2.1"
            #for some reason the newest gdal created an unknown conflict and couldn't be imported
            #also this will break general_helpers "from forzendict import frozendict" so run_my_tests.py fails
            import geopandas
            import rioxarray
            import xarray
            from shapely.geometry import mapping
            #use xarray to open dataset
            nc_x = xarray.open_dataset(path, decode_times=False)
            #set spatial dimensions by name
            nc_x.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)
            #set projection
            nc_x.rio.write_crs("epsg:4326", inplace=True)
            #use geopandas to read in map of world to use as clip
            wrld_shp = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'), crs="epsg:4326")
            #clip the imported array by the projected world map
            nc_clip = nc_x.rio.clip(wrld_shp.geometry.apply(mapping), wrld_shp.crs, drop=False)

            # we can clip/mask again by Greenland,Antarctica, or other shapefile if we want to here

            #set save directory and filename
            dataDIR = dataPath.joinpath(nc_clip_file_name(vn))
            #save to file - inefficient but only option, xarray can't create netcdf4 python object, can only save
            #This is fine though b/c it allows inspection of clipped netcdf files
            #I do this to keep all further code intact that uses netcdf4
            nc_clip.to_netcdf(dataDIR)
            #read netcdf file back in
            var = nc.Dataset(str(dataDIR)).variables[vn]
            ###########
            #return after assessing NaN data values
            return gh.get_nan_pixel_mask(var)

        masks=[ f(name)    for name in names ]
        # We compute the common mask so that it yields valid pixels for ALL variables 
        combined_mask= reduce(lambda acc,m: np.logical_or(acc,m),masks)
        print("computing means, this may also take some minutes...")

        def compute_and_cache_global_mean(vn):
            path = dataPath.joinpath(nc_file_name(vn))
            ds = nc.Dataset(str(path))
            vs=ds.variables
            lats= vs["latitude"].__array__()
            lons= vs["longitude"].__array__()
            print(vn)
            var=ds.variables[vn]
            # check if we have a cached version (which is much faster)
            gm_path = dataPath.joinpath(nc_global_mean_file_name(vn))

            gm=gh.global_mean_var(
                    lats,
                    lons,
                    combined_mask,
                    var
            )
            
            #scale per second to per day before caching
            if vn in ["npp","gpp","rh","ra"]:
                gm = gm * 86400
                
            gh.write_global_mean_cache(
                    gm_path,
                    gm,
                    vn
            )
            return gm
    
        #map variables to data
        return (
            Observables(*map(compute_and_cache_global_mean, o_names)),
            Drivers(*map(compute_and_cache_global_mean, d_names))
        )


def make_iterator_sym(
        mvs,
        V_init, 
        par_dict,
        func_dict,
        delta_t_val=1 # defaults to 1day timestep
    ):
    B_func, u_func = gh.make_B_u_funcs_2(mvs,par_dict,func_dict,delta_t_val)  

    sv=mvs.get_StateVariableTuple()
    n=len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct 
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr=np.array(
        [V_init.__getattribute__(str(v)) for v in sv]+
        [V_init.rh] #, V_init.ra]
    ).reshape(n+1,1) #reshaping is neccessary for matmul (the @ in B @ X)

    numOutFluxesBySymbol={
        k: gh.numfunc(expr_cont, mvs, delta_t_val=delta_t_val, par_dict=par_dict, func_dict=func_dict) 
        for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
    } 
    def f(it,V):
        X = V[0:n]
        b = u_func(it,X)
        B = B_func(it,X)
        X_new = X + b + B @ X
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep
        
        rh_flux=[
            numOutFluxesBySymbol[Symbol(k)](it,*X.reshape(n,))
            for k in ["c_lit_cwd","c_lit_met","c_lit_str","c_lit_mic","c_soil_met","c_soil_str","c_soil_mic","c_soil_slow","c_soil_passive"] 
            if Symbol(k) in numOutFluxesBySymbol.keys()
        ]
        
        #ra_flux=[
        #    numOutFluxesBySymbol[Symbol(k)](it,*X.reshape(n,))
        #    for k in ["c_leaf","c_root","c_wood"] 
        #    if Symbol(k) in numOutFluxesBySymbol.keys()
        #]
        
        rh = np.array(rh_flux).sum()
        #ra = np.array(ra_flux).sum()
        
        V_new = np.concatenate(
            (
                X_new.reshape(n,1),
                np.array([rh]).reshape(1,1) #,
                #np.array([ra]).reshape(1,1)
            )
            , axis=0
        )
        return V_new

    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )

def make_StartVector(mvs):
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()]+
        ["rh"] #,"ra"]
    ) 

def make_func_dict(mvs,dvs,cpa,epa):
    
    def make_temp_func(dvs):
        def temp_func(day):
            month=gh.day_2_month_index_vm(day)
            # kg/m2/s kg/m2/day;
            return (dvs.tas[month])
        return temp_func

    def make_npp_func(dvs):
        def npp_func(day):
            month=gh.day_2_month_index_vm(day)
            # kg/m2/s kg/m2/day;
            return (dvs.npp[month])
        return npp_func

    def make_gpp_func(dvs):
        def gpp_func(day):
            month=gh.day_2_month_index_vm(day)
            # kg/m2/s kg/m2/day;
            return (dvs.gpp[month])
        return gpp_func

    def make_xi_func_leaf(dvs):
        t_ref = 273.15 + 24
        t_half = 273.15 + 33
        t_exp = 1.8
        tf_frac = 0.2
        def xi_func_leaf(day):
            month = gh.day_2_month_index_vm(day)
            s_t = t_exp ** ((dvs.tas[month] - t_ref)/10)
            s_f = (1 + np.exp(tf_frac * (dvs.tas[month]-t_half)))
            return s_t / s_f 
        return xi_func_leaf
    
    def make_xi_func_soil(dvs):
        t_ref = 273.15 + 28
        t_half = 273.15 + 0
        t_exp = 1.9
        def xi_func_soil(day):
            month = gh.day_2_month_index_vm(day)
            s_t = t_exp ** ((dvs.tas[month] - t_ref)/10)
            s_f = 1 / (1 + np.exp(t_half - dvs.tas[month]))
            return s_t * s_f 
        return xi_func_soil
    
    return {
        "GPP": make_gpp_func(dvs),
        "NPP": make_npp_func(dvs),
        "xi_leaf": make_xi_func_leaf(dvs),
        "xi_soil": make_xi_func_soil(dvs),
        "temp": make_temp_func(dvs)
    }


def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers,
    ) -> Callable[[np.ndarray], np.ndarray]:
    
    # Build dictionary of model parameters
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
     
    # Create namedtuple for initial values
    StartVector=make_StartVector(mvs)
    
    # Define actual forward simulation function
    def param2res(pa):
        
        # Parameter vector
        epa=EstimatedParameters(*pa)
        
        # Build input and environmental scaler functions
        func_dict = make_func_dict(mvs,dvs,cpa,epa)
        
        # Parameter dictionary for the iterator
        apa = {**cpa._asdict(),**epa._asdict()}
        model_par_dict = {
            Symbol(k):v for k,v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }
        
        # Create a startvector for the iterator 
        V_init = StartVector(
            c_leaf = apa['c_leaf_0'],
            c_root = apa['c_root_0'],
            c_wood = apa['c_veg_0'] - (
                apa['c_leaf_0'] + 
                apa['c_root_0']
            ),
            c_lit_cwd = apa['c_lit_cwd_0'],
            c_lit_met = apa['c_lit_met_0'],
            c_lit_str = apa['c_lit_str_0'],
            c_lit_mic = apa['c_lit_mic_0'],
            c_soil_met = apa['c_soil_met_0'],
            c_soil_str = apa['c_soil_str_0'],
            c_soil_mic = apa['c_soil_mic_0'],
            c_soil_slow = apa['c_soil_slow_0'],
            c_soil_passive = apa['c_soil_0'] - (
                apa['c_lit_cwd_0'] +
                apa['c_lit_met_0'] +
                apa['c_lit_str_0'] +
                apa['c_lit_mic_0'] +
                apa['c_soil_met_0'] +
                apa['c_soil_str_0'] +
                apa['c_soil_mic_0'] +
                apa['c_soil_slow_0']
            ),
            rh = apa['rh_0'] #,
            #ra = apa['ra_0']
        )
        
        # define time step and iterator
        delta_t_val=15 
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
        )
    
        # Build functions to sum pools
        def cVegF(V):
            return float(V.c_leaf+V.c_wood+V.c_root)
        
        def cSoilF(V): 
            return float(
                V.c_lit_cwd+
                V.c_lit_met+
                V.c_lit_str+
                V.c_lit_mic+ 
                V.c_soil_met+
                V.c_soil_str+
                V.c_soil_mic+
                V.c_soil_slow+
                V.c_soil_passive
            )
        
        #empty arrays for saving data
        cVeg_arr=np.zeros(cpa.nyears)
        cSoil_arr=np.zeros(cpa.nyears)
        rh_arr=np.zeros(cpa.nyears*12)
        #ra_arr=np.zeros(cpa.nyears*12)
        
        #set days per month, month counter, and step counts
        dpm = 30                      
        im = 0
        steps_per_month=int(dpm/delta_t_val)
        steps_per_year=int(dpm/delta_t_val)*12
        
        # forward simulation by year
        for y in range(cpa.nyears):
            cVeg_avg= 0    
            cSoil_avg = 0
            for m in range(12):
                rh_avg=0
                #ra_avg=0
                for d in range(steps_per_month):    
                    V = StartVector(*it_sym.__next__())                  
                    rh_avg += V.rh
                    #ra_avg += V.ra
                    cVeg_avg += cVegF(V)
                    cSoil_avg += cSoilF(V)
                rh_arr[im] = rh_avg/steps_per_month
                #ra_arr[im] = ra_avg/steps_per_month
                im += 1
            cVeg_arr[y] = cVeg_avg/steps_per_year
            cSoil_arr[y] = cSoil_avg/steps_per_year
            #if y == 100:
            #    print(V)
        return Observables(
            cVeg = cVeg_arr,
            cSoil = cSoil_arr,
            rh = rh_arr) #,
            #ra = ra_arr)
    return param2res

def make_weighted_cost_func(
        obs: Observables
    ) -> Callable[[Observables],np.float64]:
    
    # Define cost function 
    def costfunction(out_simu: Observables) -> np.float64:
   
        #calculate costs for each data stream
        J_obj1 = (100/obs.cVeg.shape[0]) * np.sum((out_simu.cVeg - obs.cVeg)**2, axis=0) / (obs.cVeg.mean(axis=0)**2)
        J_obj2 = (100/obs.cSoil.shape[0]) * np.sum((out_simu.cSoil -  obs.cSoil)**2, axis=0) / (obs.cSoil.mean(axis=0)**2)
        J_obj3 = (100/obs.rh.shape[0]) * np.sum((out_simu.rh - obs.rh)**2, axis=0) / (obs.rh.mean(axis=0)**2)
        #J_obj4 = (100/obs.ra.shape[0]) * np.sum((out_simu.ra - obs.ra)**2, axis=0) / (obs.ra.mean(axis=0)**2)
        
        # sum costs
        J_new = 100 * (J_obj1 + J_obj2 + J_obj3) # + J_obj4)
        return J_new
    return costfunction

def make_param_filter_func(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters 
        ) -> Callable[[np.ndarray], bool]:

    # find position of beta_leaf and beta_wood
    beta_leaf_ind=EstimatedParameters._fields.index("beta_leaf")
    beta_root_ind=EstimatedParameters._fields.index("beta_root")

    def isQualified(c):
        cond1 =  (c >= c_min).all() 
        cond2 =  (c <= c_max).all() 
        cond3 =  c[beta_leaf_ind]+c[beta_root_ind] <= 0.99  
        return (cond1 and cond2) # and cond3)
        
    return isQualified


def make_traceability_iterator(mvs,dvs,cpa,epa):
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "c_leaf": apa['c_leaf_0'],     
        "c_root": apa['c_root_0'],     
        "c_wood": apa['c_veg_0'] - (apa['c_leaf_0'] +  apa['c_root_0']),  
        "c_lit_cwd": apa['c_lit_cwd_0'],
        "c_lit_met": apa['c_lit_met_0'],
        "c_lit_str": apa['c_lit_str_0'],
        "c_lit_mic": apa['c_lit_mic_0'],
        "c_soil_met": apa['c_soil_met_0'],
        "c_soil_str": apa['c_soil_str_0'],
        "c_soil_mic": apa['c_soil_mic_0'],
        "c_soil_slow": apa['c_soil_slow_0'],
        "c_soil_passive": apa['c_soil_0'] - (
                              apa['c_lit_cwd_0'] 
                            + apa['c_lit_met_0'] 
                            + apa['c_lit_str_0'] 
                            + apa['c_lit_mic_0'] 
                            + apa['c_soil_met_0'] 
                            + apa['c_soil_str_0'] 
                            + apa['c_soil_mic_0'] 
                            + apa['c_soil_slow_0']
                        )
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    fd=make_func_dict(mvs,dvs,cpa,epa)
    V_init = gh.make_InitialStartVectorTrace(
            X_0,mvs,
            par_dict=par_dict,
            func_dict=fd
    )
    it_sym_trace = gh.make_daily_iterator_sym_trace(
        mvs,
        V_init=V_init,
        par_dict=par_dict,
        func_dict=fd
    )
    return it_sym_trace


# Define start and end dates of the simulation
import datetime as dt
start_date=dt.date(1700, 1, 1)
end_date = dt.date(2019, 11, 30)

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


def numeric_X_0(mvs,dvs,cpa,epa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "c_leaf": apa['c_leaf_0'],     
        "c_root": apa['c_root_0'],     
        "c_wood": apa['c_veg_0'] - (apa['c_leaf_0'] +  apa['c_root_0']),  
        "c_lit_cwd": apa['c_lit_cwd_0'],
        "c_lit_met": apa['c_lit_met_0'],
        "c_lit_str": apa['c_lit_str_0'],
        "c_lit_mic": apa['c_lit_mic_0'],
        "c_soil_met": apa['c_soil_met_0'],
        "c_soil_str": apa['c_soil_str_0'],
        "c_soil_mic": apa['c_soil_mic_0'],
        "c_soil_slow": apa['c_soil_slow_0'],
        "c_soil_passive": apa['c_soil_0'] - (
                              apa['c_lit_cwd_0'] 
                            + apa['c_lit_met_0'] 
                            + apa['c_lit_str_0'] 
                            + apa['c_lit_mic_0'] 
                            + apa['c_soil_met_0'] 
                            + apa['c_soil_str_0'] 
                            + apa['c_soil_mic_0'] 
                            + apa['c_soil_slow_0']
                        )
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    return X_0
