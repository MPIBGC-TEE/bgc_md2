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
            return gh.get_nc_array(dataPath.joinpath(nc_global_mean_file_name(vn)),vn)
    
        return (
            Observables(*map(get_cached_global_mean, o_names)),
            Drivers(*map(get_cached_global_mean,d_names))
        )

    else:
        print("computing means, this may also take some minutes...")

        gm=gh.globalMask()
        # load an example file with mask
        template = nc.Dataset(
                        dataPath.joinpath("YIBs_S2_Monthly_tas.nc")
                    ).variables['tas'][0,:,:].mask
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
                    gcm.index_mask,
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

# deprecated
def make_func_dict_old(mvs,dvs,cpa,epa):
    
    def make_temp_func(dvs):
        def temp_func(day):
            month=gh.day_2_month_index(day)
            # kg/m2/s kg/m2/day;
            return (dvs.tas[month])
        return temp_func

    def make_npp_func(dvs):
        def npp_func(day):
            month=gh.day_2_month_index(day)
            # kg/m2/s kg/m2/day;
            return (dvs.npp[month])
        return npp_func

    def make_gpp_func(dvs):
        def gpp_func(day):
            month=gh.day_2_month_index(day)
            # kg/m2/s kg/m2/day;
            return (dvs.gpp[month])
        return gpp_func

    def make_xi_func_leaf(dvs):
        t_ref = 273.15 + 24
        t_half = 273.15 + 33
        t_exp = 1.8
        tf_frac = 0.2
        def xi_func_leaf(day):
            month = gh.day_2_month_index(day)
            s_t = t_exp ** ((dvs.tas[month] - t_ref)/10)
            s_f = (1 + np.exp(tf_frac * (dvs.tas[month]-t_half)))
            return s_t / s_f 
        return xi_func_leaf
    
    def make_xi_func_soil(dvs):
        t_ref = 273.15 + 28
        t_half = 273.15 + 0
        t_exp = 1.9
        def xi_func_soil(day):
            month = gh.day_2_month_index(day)
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


def make_func_dict(dvs, **kwargs):

    def xi_leaf(tas):
        t_ref = 273.15 + 24
        t_half = 273.15 + 33
        t_exp = 1.8
        tf_frac = 0.2
        s_t = t_exp ** ((tas - t_ref)/10)
        s_f = (1 + np.exp(tf_frac * (tas-t_half)))
        return s_t / s_f 

    def xi_soil(tas):
        t_ref = 273.15 + 28
        t_half = 273.15 + 0
        t_exp = 1.9
        s_t = t_exp ** ((tas - t_ref)/10)
        s_f = 1 / (1 + np.exp(t_half - tas))
        return s_t * s_f 

    gpp_func, npp_func, tas_func = map(
        gh.make_interpol_of_t_in_days,
        (dvs.gpp, dvs.npp, dvs.tas)
    )

    return {
        "temp": tas_func,
        "GPP": gpp_func,
        "NPP": npp_func,
        "xi_leaf": lambda t: xi_leaf(tas_func(t)),
        "xi_soil": lambda t: xi_soil(tas_func(t))
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
        func_dict = make_func_dict(dvs)
        
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
        delta_t_val=30 
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

# this function is deprecated - see general helpers traceability_iterator
# def make_traceability_iterator(mvs,dvs,cpa,epa):

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

data_str = namedtuple( # data streams available in the model
    'data_str',
    ["cVeg", "cSoil", "gpp", "npp", "ra", "rh"]
    )
    
def get_global_mean_vars_all(experiment_name):
        return(
            gh.get_global_mean_vars_all(model_folder="jon_yib", 
                            experiment_name=experiment_name,
                            lat_var="latitude",
                            lon_var="longitude",
                            ) 
        )       
################ function for computing global mean for custom data streams ###################
    
# def get_global_mean_vars_all(experiment_name="YIBs_S2_"):

    # def nc_file_name(nc_var_name, experiment_name="YIBs_S2_"):
        # if nc_var_name in ["cVeg", "cSoil"]: return experiment_name+"Annual_{}.nc".format(nc_var_name)
        # else: return experiment_name+"Monthly_{}.nc".format(nc_var_name)

    # def nc_global_mean_file_name(nc_var_name, experiment_name="YIBs_S2_"):
        # return experiment_name+"{}_gm_all.nc".format(nc_var_name)

    # data_str = namedtuple( # data streams available in the model
        # 'data_str',
        # ["cVeg", "cSoil", "gpp", "npp", "ra", "rh"]
        # )

    # names = data_str._fields
    # conf_dict = gh.confDict("jon_yib")
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
            # gm = gh.get_nc_array(dataPath.joinpath(nc_global_mean_file_name(vn,experiment_name=experiment_name)),vn)
            # return gm * 86400 if vn in ["gpp", "npp", "rh", "ra"] else gm

        # #map variables to data
        # output=gh.data_streams(*map(get_cached_global_mean, data_str._fields))
        # return (
            # output
        # )

    # else:
        # gm=gh.globalMask()
        # # load an example file with mask
        # template = nc.Dataset(dataPath.joinpath("YIBs_S2_Monthly_tas.nc")).variables['tas'][0,:,:].mask
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
                # cSoil=output.cSoil,
                # gpp=output.gpp,
                # npp=output.npp,
                # ra=output.ra,
                # rh=output.rh,
            # )
        # )

# ################## remove after
    # if all([dataPath.joinpath(nc_global_mean_file_name(vn)).exists() for vn in names]):
        # print(""" Found cached global mean files. If you want to recompute the global means
            # remove the following files: """
        # )
        # for vn in names:
            # print( dataPath.joinpath(nc_global_mean_file_name(vn)))

        # def get_cached_global_mean(vn):
            # return gh.get_nc_array(dataPath.joinpath(nc_global_mean_file_name(vn)),vn)
    
        # return (
            # Observables(*map(get_cached_global_mean, o_names)),
            # Drivers(*map(get_cached_global_mean,d_names))
        # )

    # else:
        # print("computing means, this may also take some minutes...")

        # gm=gh.globalMask()
        # # load an example file with mask
        # template = nc.Dataset(
                        # dataPath.joinpath("YIBs_S2_Monthly_tas.nc")
                    # ).variables['tas'][0,:,:].mask
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
        # def compute_and_cache_global_mean(vn):
            # path = dataPath.joinpath(nc_file_name(vn))
            # ds = nc.Dataset(str(path))
            # vs=ds.variables
            # lats= vs["latitude"].__array__()
            # lons= vs["longitude"].__array__()
            # print(vn)
            # var=ds.variables[vn]
            # # check if we have a cached version (which is much faster)
            # gm_path = dataPath.joinpath(nc_global_mean_file_name(vn))

            # gm=gh.global_mean_var(
                    # lats,
                    # lons,
                    # gcm.index_mask,
                    # var
            # )
            
            # #scale per second to per day before caching
            # if vn in ["npp","gpp","rh","ra"]:
                # gm = gm * 86400
                
            # gh.write_global_mean_cache(
                    # gm_path,
                    # gm,
                    # vn
            # )
            # return gm
    
        # #map variables to data
        # return (
            # Observables(*map(compute_and_cache_global_mean, o_names)),
            # Drivers(*map(compute_and_cache_global_mean, d_names))
        # )
