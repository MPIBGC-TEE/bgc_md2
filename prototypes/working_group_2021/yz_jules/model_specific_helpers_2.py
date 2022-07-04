import sys
# import json
from pathlib import Path
from collections import namedtuple
import netCDF4 as nc
import numpy as np
from sympy import Symbol
# from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (TimeStepIterator2)
# from copy import copy
from typing import Callable

# from general_helpers import month_2_day_index, monthly_to_yearly
from functools import reduce

sys.path.insert(0, '..')  # necessary to import general_helpers
import general_helpers as gh

# dictionaries to link file names to var names
file_name_from_var_name = {
   "npp_nlim": "JULES-ES-1p0_S2_npp.nc",
   **{
        vn: "JULES-ES-1p0_S2_{}.nc".format(vn) 
        for vn in [ "mrsos",
	 "tsl",
	 "landCoverFrac",
	 "cVeg",
	 "cSoil",
	 "rh",
	    "fVegSoil" ]
   }
}
# global dictionary to link var names to var names in files
d_name2varname_in_file = {
    "npp": 'npp_nlim',
    **{
        vn: vn
        for vn in ["fVegSoil", "mrsos", "tsl", "landCoverFrac", "cVeg", "cSoil", "rh"]
    }

}

def spatial_mask(dataPath)->'CoorMask':
    cache_path=dataPath.joinpath('mask.nc')
    

    # We read the mask of a file and also create a masks by checking for the NANs
    # we now check if any of the arrays has a time lime containing nan values 
    # APART FROM values that are already masked by the fillvalue
    
    print("computing masks to exclude pixels with nan entries, this may take some minutes...")
    

    # We compute the common mask so that it yields valid pixels for ALL variables 
    def f(d_name):
        vn_in_file = d_name2varname_in_file[d_name]
        file_name = file_name_from_var_name[vn_in_file]
        path = dataPath.joinpath(file_name)
        ds = nc.Dataset(str(path))
        var =ds.variables[vn_in_file]
        #return after assessing NaN data values
        return gh.get_nan_pixel_mask(var)

    o_names=Observables._fields
    d_names=Drivers._fields
    names = o_names + d_names 
    
    masks=[ f(name)    for name in names ]
    combined_mask = reduce(lambda acc,m: np.logical_or(acc,m),masks)

    sym_tr= gh.SymTransformers(
        itr=make_model_index_transforms(),
        ctr=make_model_coord_transforms()
    )
    return gh.CoordMask(
        combined_mask,
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
            lon2LON=lambda lon: -180+ lon-180 if lon > 180 else lon,
            LON2lon=lambda LON: 360+LON if LON < 0 else LON
    )
    
def make_model_index_transforms():
    return gh.transform_maker(
    lat_0 = -89.375,
    lon_0 = 0.9375,
    step_lat = 1.25,
    step_lon = 1.875,
 )

Observables = namedtuple(
    'Observables',
    ["cVeg", "cSoil", "rh", "fVegSoil"]
)
# Driver data streams on TRENDY server
Drivers = namedtuple(
    "Drivers",
    ["npp", "mrsos", "tsl", "landCoverFrac"]
)

Constants = namedtuple(
    "Constants",
    [
        'npp_0',  # Initial input/pools
        'rh_0',
        'c_veg_0',
        'c_soil_0',
        'fVegSoil_0',  # Total carbon mass from vegetation directly into the soil
        'nyears'  # Run time (years for my model)
    ]
)

EstimatedParameters = namedtuple(
    'EstimatedParameters',
    [
        'c_leaf_0',  # Names: c_poolname_0
        'c_wood_0',  # Only initial pools that are estimated
        'c_DPM_0',
        'c_RPM_0',
        'c_BIO_0',

        'beta_leaf',
        'beta_wood',
        'Mw',
        'Ms',
        'Topt',
        'Tcons',

        'r_c_DPM_rh',
        'r_c_RPM_rh',
        'r_c_BIO_rh',
        'r_c_HUM_rh',

        'r_c_leaf_2_c_DPM',
        'r_c_leaf_2_c_RPM',
        'r_c_wood_2_c_DPM',
        'r_c_wood_2_c_RPM',
        'r_c_root_2_c_DPM',
        'r_c_root_2_c_RPM',
        'r_c_DPM_2_c_BIO',
        'r_c_DPM_2_c_HUM',
        'r_c_RPM_2_c_BIO',
        'r_c_RPM_2_c_HUM',
        'r_c_BIO_2_c_HUM',
        'r_c_HUM_2_c_BIO',
    ]
)
def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),  # platform independent path desc. (Windows vs. linux)
        models=["JULES-ES"],
        variables=Observables._fields + Drivers._fields
    )

def get_example_site_vars(dataPath):
    # Define single geospatial cell from (3840, 144, 192)
    slayer = slice(None, None, None)  # this is the same as :
    slat = 120
    slon = 50
    t = slayer, slat, slon  # a site in South America

    # Define function to select geospatial cell and scale data
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them
        ds = nc.Dataset(str(path))

        if vn in ["npp_nlim", "gpp", "rh", "ra", "fVegSoil"]:  # (3840, 144, 192), kg m-2 s-1
            # f_veg2soil: Total carbon mass from vegetation directly into the soil
            print("reading ", vn, ", size is ", ds.variables[vn].shape)
            return ds.variables[vn][t] * 86400  # convert from kg/m2/s to kg/m2/day
        elif vn in ["tsl"]:  # Temperature of Soil - top layer, 192x144x4 (depth) x3840, 'K'
            print("reading ", vn, ", size is ", ds.variables[vn].shape)  ## (3840, 4, 144, 192)
            tmp = ds.variables[vn][:, 0, slat, slon]
            # tmp = np.mean(ds.variables[vn][:, :, slat, slon], axis=1)
            return tmp  # average soil temperature at different depth
            print("converted size is ", tmp.shape)
        elif vn in ["landCoverFrac"]:  # Plant Functional Type Grid Fraction, 192x144x17 (vegtype) x3840
            print("reading ", vn, ", size is ", ds.variables[vn].shape)  ## (3840, 17, 144, 192)
            # from IPython import embed;embed()
            var = ds.variables[vn]
            sh = var.shape
            tmp = np.sum(var[:, 0:13, slat, slon], axis=1)
            print("converted size is ", tmp.shape)
            #'0.BdlDcd; 1.BdlEvgTrop; 2.BdlEvgTemp; 3.NdlDcd; 4.NdlEvg; 5.c3grass; 6.c3crop; 7.c3pasture; 8.c4grass; 9.c4crop; 10.c4pasture; 11.shrubDcd; 12.shrubEvg; 13.urban; 14.lake; 15.soil; 16.ice'
            print("Forest cover (t=0) is ", np.sum(var[1, 0:5, slat, slon]))
            print("Grass + crop + shrub (t=0) cover is ", np.sum(var[1, 5:13, slat, slon]))
            print("Non-vegetation (t=0) cover is ", np.sum(var[1, 13:17, slat, slon]))
            return tmp  # sum the vegetation coverages
        else:
            print("reading ", vn, ", size is ", ds.variables[vn].shape)
            return ds.variables[vn][t]

    ## Link symbols and data:

    o_tuples = [
        (
            d_name2varname_in_file[f],
            file_name_from_var_name[d_name2varname_in_file[f]]
        )
        for f in Observables._fields
    ]

    d_tuples = [
        (
            d_name2varname_in_file[f],
            file_name_from_var_name[d_name2varname_in_file[f]]
        )
        for f in Drivers._fields
    ]

    # Link symbols and data for observables/drivers
    # print(o_tuples)
    return (
        Observables(*map(f, o_tuples)),
        Drivers(*map(f, d_tuples))
    )


def get_global_mean_vars(dataPath):
    # Define function to select geospatial cell and scale data
    gm=gh.globalMask()
    # load an example file with mask
    template = nc.Dataset(dataPath.joinpath("JULES-ES-1p0_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
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
    def f(tup):
        vn, fn = tup
        path = dataPath.joinpath(fn)
        # Read NetCDF data but only at the point where we want them
        ds = nc.Dataset(str(path))
        lats = ds.variables["latitude"].__array__()
        lons = ds.variables["longitude"].__array__()

        if vn in ["npp_nlim", "gpp", "rh", "ra", "fVegSoil"]:  # (3840, 144, 192), kg m-2 s-1
            # f_veg2soil: Total carbon mass from vegetation directly into the soil
            #print("reading ", vn, ", size is ", ds.variables[vn].shape)
            ma=ds.variables[vn][:,:,:].data * 86400# convert from kg/m2/s to kg/m2/day
        
        elif vn in ["tsl"]:  # Temperature of Soil - top layer, 192x144x4 (depth) x3840, 'K'
            print("reading ", vn, ", size is ", ds.variables[vn].shape)  ## (3840, 4, 144, 192)
            ma=ds.variables[vn][:, 0, :, :].data
            #print("converted size is ", ma.shape)
            # ma= np.mean(ds.variables[vn][:, :, slat, slon].data, axis=1)
            # average soil temperature at different depth
    
        elif vn in ["landCoverFrac"]:  # Plant Functional Type Grid Fraction, 192x144x17 (vegtype) x3840
            #print("reading ", vn, ", size is ", ds.variables[vn].shape)  ## (3840, 17, 144, 192)
            # from IPython import embed;embed()
            var = ds.variables[vn]
            sh = var.shape
            ma= np.zeros(shape = var[:, 0, :, :].shape)
            #print("creating a zero arry, shape is ", var[:, 0, :, :].shape)
            for i in range(13):
                ma = ma + var[:, i, :, :].data
                # print("converted size is ", ma.shape)
            #'0.BdlDcd; 1.BdlEvgTrop; 2.BdlEvgTemp; 3.NdlDcd; 4.NdlEvg; 5.c3grass; 6.c3crop; 7.c3pasture; 8.c4grass; 9.c4crop; 10.c4pasture; 11.shrubDcd; 12.shrubEvg; 13.urban; 14.lake; 15.soil; 16.ice'
            # print("Forest cover (t=0) is ", np.sum(var[1, 0:5, :, :]))
            # print("Grass + crop + shrub (t=0) cover is ", np.sum(var[1, 5:13, slat, slon]))
            # print("Non-vegetation (t=0) cover is ", np.sum(var[1, 13:17, slat, slon]))
            # sum the vegetation coverages
            
        
        else:
            #print("reading ", vn, ", size is ", ds.variables[vn].shape)
            ma=ds.variables[vn][:,:,:].data
       

        return gh.global_mean(
            lats, 
            lons, 
            np.ma.masked_array(
                data=ma,
                # stack the 2d gcm 
                mask=np.stack(
                    [
                        gcm.index_mask 
                        for i in range(ma.shape[0])
                    ],
                    axis=0
                )
            )                
         )  

    # Link symbols and data:

    # Create file names (single step if files similarly named)


    o_tuples = [
        (
            d_name2varname_in_file[f],
            file_name_from_var_name[d_name2varname_in_file[f]]
        )
        for f in Observables._fields
    ]

    d_tuples = [
        (
            d_name2varname_in_file[f],
            file_name_from_var_name[d_name2varname_in_file[f]]
        )
        for f in Drivers._fields
    ]

    # Link symbols and data for observables/drivers
    # print(o_tuples)
    return (
        Observables(*map(f, o_tuples)),
        Drivers(*map(f, d_tuples))
    )


def make_StartVector(mvs):
    # fixme mm:
    # add the fraction!
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()] +
        ["rh","fVegSoil"]
    )


def make_xi_func(tsl, Mw, Ms, Topt, Tcons, mrsos, landCoverFrac):
    def xi_func(day):
        mi = gh.day_2_month_index(day)
        # alternative FT
            # Q10 function (this is not what Clark et al 2011 Fig. 2 presented, the equation must be wrong)
        #FT = 2.0 ** ((tsl[mi] - 298.15) / 10)  # temperature rate modifier
            # RothC temperature function (Jenkinson 1990)
        FT = Tcons / (1 + np.exp(106/(tsl[mi] - 273.1 + Topt)))
        FV = 0.6 + 0.4 * (1 - landCoverFrac[mi] / 100)  # effect of vegetation cover
        # Mw is soil moisture at wilting point as a fraction of saturation
        # Ms is soil moisture content at saturation
        S0 = 0.5 * (1 + Mw)  # optimum soil moisture
        Smin = 1.7 * Mw  # lower threshold soil moisture for soil respiration
        if S0 < mrsos[mi]/Ms:
            FS = 1 - 0.8 * (mrsos[mi]/Ms - S0)  # effect of soil moisture
        if (Smin < mrsos[mi]/Ms) and (mrsos[mi]/Ms <= S0):
            FS = 0.2 + 0.8 * (mrsos[mi]/Ms - Smin) / (S0 - Smin)
        if mrsos[mi]/Ms <= Smin:
            FS = 0.2
        # print("FT,FV,FS", FT, FV, FS)
        rh_factor = FT * FV * FS
        return rh_factor # 1.0     # Set to 1 if no scaler implemented
        # return 1.0

    return xi_func


def make_npp_func(dvs):
    def func(day):
        month = gh.day_2_month_index(day)
        # kg/m2/s kg/m2/day;
        return (dvs.npp[month])

    return func


def make_func_dict(mvs, dvs, cpa, epa):
    return {
        "NPP": make_npp_func(dvs),
        "xi": make_xi_func(dvs.tsl, epa.Mw, epa.Ms, epa.Topt, epa.Tcons, dvs.mrsos, dvs.landCoverFrac)
    }


def make_iterator_sym(
        mvs,
        V_init: "StartVector",
        par_dict,
        func_dict,
        delta_t_val=1  # defaults to 1day timestep
):
    B_func, u_func = gh.make_B_u_funcs_2(mvs, par_dict, func_dict, delta_t_val)
    sv = mvs.get_StateVariableTuple()
    n = len(sv)
    # build an array in the correct order of the StateVariables which in our case is already correct
    # since V_init was built automatically but in case somebody handcrafted it and changed
    # the order later in the symbolic formulation....
    V_arr = np.array(
        [V_init.__getattribute__(str(v)) for v in sv] +
        [V_init.rh, V_init.fVegSoil]
    ).reshape(n + 2, 1)  # reshaping is neccessary for matmul (the @ in B @ X)

    # numOutFluxesBySymbol={
    #    k:numfunc(expr_cont,delta_t_val=delta_t_val)
    #    for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
    # }
    numOutFluxesBySymbol = {
        k: gh.numfunc(
            expr_cont, 
            mvs,
            delta_t_val=delta_t_val,
            par_dict=par_dict,
            func_dict=func_dict
        )
        for k, expr_cont in mvs.get_OutFluxesBySymbol().items()
    }
    v2sfl = {
        k:v for k,v in mvs.get_InternalFluxesBySymbol().items()
        # included DPM and RPM as donor pools 
        # becaue the previous list didn't have mass balance between increase of cVeg (~1kg) and (NPP - fVegSoil) (0.09kg)
        if k[0] in map(Symbol,["c_leaf","c_wood","c_root"])  # ,"c_DPM","c_RPM"
    }
    numv2sfl = {
        k: gh.numfunc(
            expr_cont,
            mvs,delta_t_val=delta_t_val,
            par_dict=par_dict,
            func_dict=func_dict
        )
        for k, expr_cont in v2sfl.items()
    }

    def f(it, V):
        X = V[0:n]
        b = u_func(it, X)
        B = B_func(it, X)
        X_new = X + b + B @ X
        # we also compute the autotrophic and heterotrophic respiration in every (daily) timestep

        rh=np.array([
            f(it, *X.reshape(n, ))
            for k,f in numOutFluxesBySymbol.items()
            if str(k) in ["c_DPM", "c_RPM", "c_BIO", "c_HUM"]
        ]).sum()
        
        fVegSoil=np.array(
            [
                f(it, *(X.reshape(n, )))
                for f in numv2sfl.values()
            ]
        ).sum()


        V_new = np.concatenate(
            (
                X_new.reshape(n,1),
                np.array([rh]).reshape(1,1),
                np.array([fVegSoil]).reshape(1,1)
            )
            , axis=0
        )
        return V_new
    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )

def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers,
) -> Callable[[np.ndarray], np.ndarray]:
    # Build iterator
    # Need dictionary of numeric values for all parameters that are not state variables/time
    srm = mvs.get_SmoothReservoirModel()
    model_par_dict_keys = srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))] +
        list(mvs.get_StateVariableTuple())
    )
    StartVector = make_StartVector(mvs)

    ########### Ask Markus here
    # Build environmental scaler function  ############### day or monthly, 
    # monthly inputs here, Mw and Ms are the to-be-estimated parameters

    # Define actual forward simulation function
    def param2res(pa):
        # Define function dictionary
        # Parameter vector
        epa = EstimatedParameters(*pa)

        func_dict = make_func_dict(mvs, dvs, cpa, epa)
        # Create a startvector for the iterator
        V_init = StartVector(
            c_leaf=epa.c_leaf_0,
            c_wood=epa.c_wood_0,
            c_root=cpa.c_veg_0 - (epa.c_leaf_0 + epa.c_wood_0),

            c_DPM=epa.c_DPM_0,
            c_RPM=epa.c_RPM_0,
            c_BIO=epa.c_BIO_0,
            c_HUM=cpa.c_soil_0 - (epa.c_DPM_0 + epa.c_RPM_0 + epa.c_BIO_0),

            rh=cpa.rh_0,
            fVegSoil=cpa.fVegSoil_0
        )

        # Parameter dictionary for the iterator
        apa = {**cpa._asdict(), **epa._asdict()}
        model_par_dict = {
            Symbol(k): v for k, v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }

        # size of the timestep in days
        # We could set it to 30 o
        # it makes sense to have a integral divisor of 30 (15,10,6,5,3,2)
        delta_t_val = 15
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val=delta_t_val
        )
        # Now that we have the iterator we can start to compute.
        # the first thing to notice is that we don't need to store all values (daili yi delta_t_val=1)
        # since the observations are recorded monthly while our iterator possibly has a smaller timestep.
        # - for the pool contents we only need to record the last day of the month
        # - for the respiration(s) ra and rh we want an accumulated value (unit mass)
        #   have to sum up the products of the values*delta_t over a month
        #
        # Note: check if TRENDY months are like this...
        # days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        
        # fixme mm:
        # - We want monthly averages instead of the presently implemented last value of the month
        # - We want to compute the flux from Veg 2 Soil.
        sols = []
        dpm = 30  #
        n = len(V_init)
        for m in range(cpa.nyears*12):
            # dpm = days_per_month[ m % 12]
            mrh = 0
            for d in range(int(dpm / delta_t_val)):
                v = it_sym.__next__().reshape(n, )
                # actually the original datastream seems to be a flux per area (kg m^-2 ^-s )
                # at the moment the iterator also computes a flux but in kg^-2 ^day
            V = StartVector(*v)
            # from IPython import embed;embed()
            o = Observables(
                cVeg=float(V.c_leaf + V.c_wood + V.c_root),
                cSoil=float(V.c_DPM + V.c_RPM + V.c_BIO + V.c_HUM),
                rh=V.rh,
                fVegSoil=V.fVegSoil #Total carbon mass from vegetation directly into the soil
            )
            sols.append(o)

        # quick conversion of the list of observables into an  Observables tuple of arrays 
        def arr(field):
            n=len(sols)
            a=np.zeros(n)
            for i in range(n):
                a[i]=sols[i].__getattribute__(field)
            return a

        return Observables(*( arr(field) for field in Observables._fields))
        # sol = np.stack(sols)
        # convert to yearly output if necessary (the monthly pool data looks very "yearly")
        # sol_yr=np.zeros(int(cpa.number_of_months/12)*sol.shape[1]).reshape([int(cpa.number_of_months/12),sol.shape[1]])
        # for i in range(sol.shape[1]):
        #   sol_yr[:,i]=monthly_to_yearly(sol[:,i])
        # sol=sol_yr
        # return sol

    return param2res

def make_weighted_cost_func(
        obs: Observables
) -> Callable[[Observables], np.float64]:
    # first unpack the observation array into its parts
    # cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
    def costfunction(out_simu: Observables) -> np.float64:
        # fixme
        # loop over the fields and move to general_helpers

        J_obj1 = np.mean((out_simu.cVeg - obs.cVeg) ** 2) / (2 * np.var(obs.cVeg))
        J_obj2 = np.mean((out_simu.cSoil - obs.cSoil) ** 2) / (2 * np.var(obs.cSoil))

        J_obj3 = np.mean((out_simu.rh - obs.rh) ** 2) / (2 * np.var(obs.rh))
        J_obj4 = np.mean((out_simu.fVegSoil - obs.fVegSoil) ** 2) / (2 * np.var(obs.fVegSoil))

        J_new = (J_obj1 + J_obj2 + J_obj3*2) #  + J_obj4 removed because of the issues about fVegSoil definition and mass balance 
        print("cVeg J1 = ", J_obj1, "c Soil J2 = ", J_obj2, "rh J3 = ", J_obj3*2)
        return J_new

    return costfunction


def make_weighted_cost_func_2(
        obs: Observables
) -> Callable[[Observables], np.float64]:
    # first unpack the observation array into its parts
    # cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
    
    denominator_cVeg = np.sum((obs.cVeg - np.mean(obs.cVeg)) ** 2)
    denominator_cSoil = np.sum((obs.cSoil - np.mean(obs.cSoil)) ** 2)
    denominator_rh = np.sum((obs.rh - np.mean(obs.rh)) ** 2)
    # now we compute a scaling factor per observable stream
    # fixme mm 10-28-2021
    #   The denominators in this case are actually the TEMPORAL variances of the data streams

    #   The desired effect of automatically adjusting weight could be achieved
    #   by the mean itself.
    # dominators = means
    def costfunction(out_simu: Observables) -> np.float64:
        J1 = np.sum((obs.cVeg - out_simu.cVeg) ** 2, axis=0) / denominator_cVeg
        J2 = np.sum((obs.cSoil - out_simu.cSoil) ** 2, axis=0) / denominator_cSoil
        J3 = np.sum((obs.rh - out_simu.rh) ** 2, axis=0) / denominator_rh
        cost = J1 + J2 + J3
        print("cVeg J1 = ", J1, "c Soil J2 = ", J2, "rh J3 = ", J3)
        #cost = np.mean(
        #    np.sum((obs - mod) ** 2, axis=0) / denominators * 100
        #)
        return cost
    return costfunction


def make_traceability_iterator(mvs,dvs,cpa,epa):
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "c_leaf": apa['c_leaf_0'],     
        "c_wood": apa['c_wood_0'],     
        "c_root": apa['c_veg_0'] - (apa['c_leaf_0'] +  apa['c_wood_0']),  
        "c_DPM": apa['c_DPM_0'],
        "c_RPM": apa['c_RPM_0'],
        "c_BIO": apa['c_BIO_0'],
        "c_HUM": apa['c_soil_0'] - (apa['c_DPM_0'] + apa['c_RPM_0'] + apa['c_BIO_0'])
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    fd = make_func_dict(mvs, dvs, cpa, epa)
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

def numeric_X_0(mvs,dvs,cpa,epa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "c_leaf": apa['c_leaf_0'],     
        "c_wood": apa['c_wood_0'],     
        "c_root": apa['c_veg_0'] - (apa['c_leaf_0'] +  apa['c_wood_0']),  
        "c_DPM": apa['c_DPM_0'],
        "c_RPM": apa['c_RPM_0'],
        "c_BIO": apa['c_BIO_0'],
        "c_HUM": apa['c_soil_0'] - (apa['c_DPM_0'] + apa['c_RPM_0'] + apa['c_BIO_0'])
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    return X_0

def make_tuple_traceability_iterator(mvs,dvs,cpa,epa):

    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0 = numeric_X0(mvs, dvs, cpa, epa)
    fd = make_func_dict(mvs, dvs, cpa, epa)
    V_init = gh.make_InitialTraceTuple(
            X_0,
            mvs,
            par_dict=par_dict,
            func_dict=fd
    )
    it= gh.make_trace_tuple_iterator(
        mvs,
        V_init=V_init,
        par_dict=par_dict,
        func_dict=fd
    )
    return it


def start_date():
    ## this function is important to syncronise our results
    ## because our data streams start at different times the first day of 
    ## a simulation day_ind=0 refers to different dates for different models
    ## we have to check some assumptions on which this calculation is based
    ## for jules the data points are actually spaced monthly with different numbers of days
    ## 
    ## Here is how to get these values
    #ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("JULES-ES-1p0_S2_cVeg.nc")))
    #times = ds.variables["time"]
    ## we have to check some assumptions on which this calculation is based
    ## for jules the data points are actually spaced with different numbers of days between monthly
    ## data point
    ## we can see this by looking at the first 24 months
    ## for i in range(24):
    ##     print((times[i + 1] - times[i])/(3600 * 24))

    #ts = times[0] #time of first observation in seconds_since_2010_01_01_00_00_00
    #td = int(ts / (3600 * 24)) #in days since_2010_01_01_00_00_00
    #import datetime as dt
    #ad = dt.date(1, 1, 1) # first of January of year 1 
    #sd = dt.date(2010, 1, 1)
    #td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    ## from td_aD (days since 1-1-1) we can compute the year month and day
    return gh.date(
        year=1700, 
        month=1,
        day=16
    )
