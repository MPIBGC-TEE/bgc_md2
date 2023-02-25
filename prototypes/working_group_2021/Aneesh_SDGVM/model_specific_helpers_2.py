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
import general_helpers as gh
from functools import reduce

sys.path.insert(0,'..') # necessary to import general_helpers
import general_helpers as gh

def spatial_mask(dataPath)->'CoorMask':
    mask=nc.Dataset(dataPath.joinpath("SDGVM_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
    sym_tr= gh.SymTransformers(
        itr=make_model_index_transforms(),
        ctr=make_model_coord_transforms()
    )
    return gh.CoordMask(
        mask,
        sym_tr
    )

def make_model_coord_transforms():
    return gh.identicalTransformers()

def make_model_index_transforms():
    return gh.transform_maker(
    lat_0 = -89.5,
    lon_0 = -179.5,
    step_lat = 1,
    step_lon = 1,
 )

Observables = namedtuple(
    'Observables',
    ["cVeg", "cRoot", "cLitter", "cSoil", "rh"]
)
Drivers=namedtuple(
    "Drivers",
    ["npp"]
)    

Constants = namedtuple(
    "Constants",
    [
        "cLitter_0",
        "cSoil_0",
        "cRoot_0",
        "cVeg_0",
        "npp_0",
        "rh_0",
        "number_of_months" # necessary to prepare the output in the correct lenght
    ]
)

EstimatedParameters = namedtuple(
    "EstimatedParameters",[
         'beta_leaf',
         'beta_wood',
         #beta_root,
         'r_C_leaf2abvstrlit',
         'r_C_abvmetlit2surface_microbe',
         'r_C_abvstrlit2slowsom',
         'r_C_abvstrlit2surface_microbe',
         'r_C_belowmetlit2soil_microbe',
         'r_C_belowstrlit2slowsom',
         'r_C_belowstrlit2soil_microbe',
         #'r_C_leached',
         'r_C_leaf2abvmetlit',
         'r_C_passsom2soil_microbe',
         'r_C_root2belowmetlit',
         'r_C_root2belowstrlit',
         'r_C_slowsom2passsom',
         'r_C_slowsom2soil_microbe',
         'r_C_soil_microbe2passsom',
         'r_C_soil_microbe2slowsom',
         'r_C_surface_microbe2slowsom',
         'r_C_wood2abvmetlit',
         'r_C_wood2abvstrlit',
         'r_C_abvstrlit_rh',
         'r_C_abvmetlit_rh',
         'r_C_belowstrlit_rh',
         'r_C_belowmetlit_rh',
         'r_C_surface_microbe_rh',
         'r_C_slowsom_rh',
         'r_C_passsom_rh',
         'r_C_soil_microbe_rh',
         'C_leaf_0',
        #'C_root_0',
         'C_abvstrlit_0',
         'C_abvmetlit_0',
         'C_blwstrlit_0',
         'C_surfacemic_0',
         'C_soilmic_0',
         'C_slow_0'
    ]
)
# note that the observables are from the point of view of the mcmc also considered to be constant (or unestimated)
# parameters. In this case we may use only the first entry e.g. to derive startvalues.
# The destinction is only made for the data assimilation to isolate those parameters
# that the mcmc will try to optimise

#create a small model specific function that will later be stored in the file model_specific_helpers.py
def download_my_TRENDY_output(conf_dict):
    gh.download_TRENDY_output(
        username=conf_dict["username"],
        password=conf_dict["password"],
        dataPath=Path(conf_dict["dataPath"]),#platform independent path desc. (Windows vs. linux)
        models=['SDGVM'],
        variables = Observables._fields + Drivers._fields
    )

def get_example_site_vars(dataPath):
    # According to the netcdf metadata the datasets are not uniform
    # - npp and rh start at 360h (15 days) after 01-01-1900 and are recorded every 30 days
    #   these are interpreted as mid-monthly
    # - C_litter, C_soil, C_veg, C_root start at 4320 h = 180 days = 6 months after 01-01-1700
    #   These can be seen at midyearly values since there are 6 (midmonth) times of the npp and rh measurements after the last (midyear)
    #   measurement of C_litter, C_soil, C_veg, C_root

    # To combine these streams into a consistent array of observations we will:
    # 1. Make C_litter, C_soil, C_veg, C_root refer to hours after 01/01/1900 (as npp and rh)
    #
    # 2. cut away everything before 1900 from them (cutting of the first 200y)
    #
    # Note:
    #    We will have to adapt the costfunction and param2res later to accommodate the different
    #    resolution of the C_pool and rh observations.



    # 1.)
    # pick one of the 1700 yearly example ds to get at the times
    # convert time to refer to the same starting point (from h after 1700 to h after 1900)
    hs_from_1900=nc.Dataset(dataPath.joinpath('SDGVM_S2_cLitter.nc')).variables['time'][:]-200*12*30*24

    #2.)
    # find the index after which we are after 01/01 1900
    ind_start = 200
    c_h_from_1900_after_1900=hs_from_1900[ind_start:]
    c_h_from_1900_after_1900.shape

    # pick up 1 site   wombat state forest for the spacial selection
    s_rh  = slice(None, None, None)  # this is the same as :
    s_c  = slice(ind_start, None, None)  # this is the same as ind_start:
    #t = s, 50, 33  # [t] = [:,49,325]
    loc=(-25,16)
    t_rh = s_rh,*loc
    t_c = s_c, *loc
    print(t_c)

    # Read NetCDF data and slice out our site
    arr_dict={
        **{
            vn:nc.Dataset(str(dataPath.joinpath(fn))).variables[vn][t_c]
            for vn,fn in  {
                'cLitter': 'SDGVM_S2_cLitter.nc',
                'cSoil': 'SDGVM_S2_cSoil.nc',
                'cVeg': 'SDGVM_S2_cVeg.nc',
                'cRoot': 'SDGVM_S2_cRoot.nc',
            }.items()
        },
        **{
            vn:nc.Dataset(str(dataPath.joinpath(fn))).variables[vn][t_rh]*86400   # kg/m2/s kg/m2/day;
            for vn,fn in {
                'npp': 'SDGVM_S2_npp.nc',
                'rh': 'SDGVM_S2_rh.nc'
            }.items()
        }
    }

    return (
        Observables(*(arr_dict[k] for k in Observables._fields)),
        Drivers(*(arr_dict[k] for k in Drivers._fields))
    )

def get_global_mean_vars(dataPath):
    # According to the netcdf metadata the datasets are not uniform
    # - npp and rh start at 360h (15 days) after 01-01-1900 and are recorded every 30 days
    #   these are interpreted as mid-monthly
    # - C_litter, C_soil, C_veg, C_root start at 4320 h = 180 days = 6 months after 01-01-1700
    #   These can be seen at midyearly values since there are 6 (midmonth) times of the npp and rh measurements after the last (midyear)
    #   measurement of C_litter, C_soil, C_veg, C_root

    # To combine these streams into a consistent array of observations we will:
    # 1. Make C_litter, C_soil, C_veg, C_root refer to hours after 01/01/1900 (as npp and rh)
    #
    # 2. cut away everything before 1900 from them (cutting of the first 200y)
    #
    # Note:
    #    We will have to adapt the costfunction and param2res later to accommodate the different
    #    resolution of the C_pool and rh observations.



    # 1.)
    # pick one of the 1700 yearly example ds to get at the times
    # convert time to refer to the same starting point (from h after 1700 to h after 1900)
    hs_from_1900=nc.Dataset(dataPath.joinpath('SDGVM_S2_cLitter.nc')).variables['time'][:]-200*12*30*24

    #2.)
    # find the index after which we are after 01/01 1900
    ind_start = 200
    c_h_from_1900_after_1900=hs_from_1900[ind_start:]
    c_h_from_1900_after_1900.shape

    # pick up 1 site   wombat state forest for the spacial selection
    s_rh  = slice(None, None, None)  # this is the same as :
    s_c  = slice(ind_start, None, None)  # this is the same as ind_start:
    #t = s, 50, 33  # [t] = [:,49,325]
    #loc=(-25,16)
    t_rh = s_rh
    t_c = s_c
    print(t_c)
    ds = nc.Dataset(str(dataPath.joinpath('SDGVM_S2_cLitter.nc')))
    lats = ds.variables["latitude"].__array__()
    lons = ds.variables["longitude"].__array__()
    gm = gh.global_mean
    
    # Read NetCDF data and slice out our site
    arr_dict={
        **{
            vn:gm(lats,lons,nc.Dataset(str(dataPath.joinpath(fn))).variables[vn][t_c])
            for vn,fn in  {
                'cLitter': 'SDGVM_S2_cLitter.nc',
                'cSoil': 'SDGVM_S2_cSoil.nc',
                'cVeg': 'SDGVM_S2_cVeg.nc',
                'cRoot': 'SDGVM_S2_cRoot.nc',
            }.items()
        },
        **{
            vn:gm(lats,lons, nc.Dataset(str(dataPath.joinpath(fn))).variables[vn].__array__())*86400   # kg/m2/s kg/m2/day;
            for vn,fn in {
                'npp': 'SDGVM_S2_npp.nc',
                'rh': 'SDGVM_S2_rh.nc'
            }.items()
        }
    }

    return (
        Observables(*(arr_dict[k] for k in Observables._fields)),
        Drivers(*(arr_dict[k] for k in Drivers._fields))
    )

## deprecated
#def make_npp_func(dvs):
#    def func(day):
#        month=gh.day_2_month_index(day)
#        # kg/m2/s kg/m2/day;
#        return (dvs.npp[month])
#
#    return func
#
#
## deprecated
#def make_xi_func(dvs):
#    def func(day):
#        return 1.0 # preliminary fake for lack of better data... 
#    return func
#
#
## deprecated
#def make_func_dict_old(mvs,dvs,cpa,epa):
#    return {
#        "NPP": make_npp_func(dvs),
#        "xi": make_xi_func(dvs)
#    }
#

def make_func_dict(dvs, **kwargs):
    return {
        "NPP": gh.make_interpol_of_t_in_days(dvs.npp),
        "xi": lambda t: 1
    }

def make_StartVector(mvs):
    return namedtuple(
        "StartVector",
        [str(v) for v in mvs.get_StateVariableTuple()]+
        ["rh"]
    ) 

from typing import Callable
from functools import reduce


def make_param2res_sym(
        mvs,
        cpa: Constants,
        dvs: Drivers
    ) -> Callable[[np.ndarray], np.ndarray]: 
    # To compute numeric solutions we will need to build and iterator 
    # as we did before. As it will need numeric values for all the parameters 
    # we will have to create a complete dictionary for all the symbols
    # exept those for the statevariables and time.
    # This set of symbols does not change during the mcmc procedure, since it only
    # depends on the symbolic model.
    # Therefore we create it outside the mcmc loop and bake the result into the 
    # param2res function.
    # The iterator does not care if they are estimated or not it only 
    # wants a dictionary containing key: value pairs for all
    # parameters that are not state variables or the symbol for time
    StartVector=make_StartVector(mvs)
    srm=mvs.get_SmoothReservoirModel()
    model_par_dict_keys=srm.free_symbols.difference(
        [Symbol(str(mvs.get_TimeSymbol()))]+
        list(mvs.get_StateVariableTuple())
    )
    
    # the time dependent driver function for gpp does not change with the estimated parameters
    # so its enough to define it once as in our test
    
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # create a startvector for the iterator from both estimated and fixed parameters 
        # The order is important and corresponds to the order of the statevariables
        # Our namedtuple StartVector takes care of this       
        V_init = StartVector(
            C_leaf = epa.C_leaf_0,
            C_root= cpa.cRoot_0,
            C_wood = cpa.cVeg_0 - epa.C_leaf_0 - cpa.cRoot_0,
            C_abvstrlit = epa.C_abvstrlit_0,
            C_abvmetlit = epa.C_abvmetlit_0,
            C_belowstrlit = epa.C_blwstrlit_0,
            C_belowmetlit = cpa.cLitter_0 - epa.C_abvstrlit_0 - epa.C_abvmetlit_0 - epa.C_blwstrlit_0,
            C_surface_microbe = epa.C_surfacemic_0,
            C_soil_microbe = epa.C_soilmic_0,
            C_slowsom = epa.C_slow_0,
            C_passsom= cpa.cSoil_0 - epa.C_surfacemic_0 - epa.C_soilmic_0 - epa.C_slow_0,
            rh=cpa.rh_0
        )
        # next we create the parameter dict for the iterator
        # The iterator does not care if they are estimated or not so we look for them
        # in the combination
        apa = {**cpa._asdict(),**epa._asdict()}
        model_par_dict = {
            Symbol(k):v for k,v in apa.items()
            if Symbol(k) in model_par_dict_keys
        }
        # Beside the par_dict the iterator also needs the python functions to replace the symbolic ones with
        # our fake xi_func could of course be defined outside of param2res but in general it
        # could be build from estimated parameters and would have to live here...
        func_dict = make_func_dict(dvs)
        
        delta_t_val=1
        it_sym = make_iterator_sym(
            mvs,
            V_init=V_init,
            par_dict=model_par_dict,
            func_dict=func_dict,
            delta_t_val = delta_t_val
        )
        
        # Now that we have the iterator we can start to compute.
        # the first thing to notice is that we don't need to store all daily values,
        # since the observations are recorded monthly while our iterator has a
        # daily timestep.
        # - for the pool contents we only need to record the last day of the month
        # - for the respiration(s) ra and rh we have to sum up the daily values 
        #   over a month
        # 
        # Note: check if TRENDY months are like this...
        #days_per_month = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
        rhs=np.zeros(cpa.number_of_months)
        number_of_years=int(cpa.number_of_months/12)
        cVegs=np.zeros(number_of_years)
        cRoots=np.zeros(number_of_years)
        cLitters=np.zeros(number_of_years)
        cSoils=np.zeros(number_of_years)
        sols=[]
        m_id=0
        #y_id=0
        dpm=30
        #dpy=30*12
        steps_per_month = int(dpm / delta_t_val)
        steps_per_year = steps_per_month*12
        for y in range(number_of_years):
            cVeg=0
            cRoot = 0
            cLitter= 0
            cSoil= 0
            for m in range(12):
                #dpm = days_per_month[ m % 12]  
                #mra=0
                mrh=0
                for d in range(steps_per_month):
                    v = it_sym.__next__()
                    #mra +=v[9,0]
                    V=StartVector(*v)
                    print(V)
                    cVeg+=float(V.C_leaf+V.C_wood+V.C_root)
                    cRoot+= float(V.C_root)
                    cLitter+=float(V.C_abvstrlit + V.C_abvmetlit + V.C_belowstrlit + V.C_belowmetlit)
                    cSoil+=float(V.C_surface_microbe + V.C_soil_microbe + V.C_slowsom + V.C_passsom)
                    mrh +=V.rh
                rhs[m_id]=mrh/steps_per_month
                m_id+=1
            cVegs[y]=cVeg/steps_per_year
            cRoots[y] = cRoot/steps_per_year
            cLitters[y]= cLitter/steps_per_year 
            cSoils[y]= cSoil/steps_per_year
            #y_id+= 1
        return Observables(
            cVeg=cVegs,
            cRoot=cRoots,
            cLitter=cLitters,
            cSoil=cSoils,
            rh=rhs#/(60*60*24)
        )
            
        
    return param2res

from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
        TimeStepIterator2,
)
from copy import copy

# def make_daily_iterator_sym(
        # mvs,
        # V_init,
        # par_dict,
        # func_dict,
        # delta_t_val=1
    # ):
    # B_func, u_func = gh.make_B_u_funcs_2(mvs,par_dict,func_dict,delta_t_val)
    # def numfunc(expr_cont,delta_t_val):
        # # build the discrete expression (which depends on it,delta_t instead of
        # # the continius one that depends on t (TimeSymbol))
        # it=Symbol("it")           #arbitrary symbol for the step index )
        # t=mvs.get_TimeSymbol()
        # delta_t=Symbol('delta_t')
        # expr_disc = expr_cont.subs({t:delta_t*it})
        # return hr.numerical_function_from_expression(
            # expr=expr_disc.subs({delta_t:delta_t_val}),
            # tup=(it, *mvs.get_StateVariableTuple()),
            # parameter_dict=par_dict,
            # func_set=func_dict
        # )

    # sv=mvs.get_StateVariableTuple()
    # n=len(sv)
    # # build an array in the correct order of the StateVariables which in our case is already correct
    # # since V_init was built automatically but in case somebody handcrafted it and changed
    # # the order later in the symbolic formulation....
    # V_arr=np.array(
        # [V_init.__getattribute__(str(v)) for v in sv]+
        # #[V_init.ra,V_init.rh]
        # [V_init.rh]
    # ).reshape(n+1,1) #reshaping is neccessary for matmux
    # numOutFluxesBySymbol={
        # k:numfunc(expr_cont,delta_t_val=delta_t_val)
        # for k,expr_cont in mvs.get_OutFluxesBySymbol().items()
    # }

    # def f(it,V):
        # X = V[0:n]
        # b = u_func(it,X)
        # B = B_func(it,X)
        # outfluxes = B @ X
        # X_new = X + b + outfluxes
        # rh = np.sum(
            # [
                # numOutFluxesBySymbol[Symbol(k)](it,*X)
                # for k in [
                    # "C_abvstrlit",
                    # "C_abvmetlit",
                    # "C_belowstrlit",
                    # "C_belowmetlit",
                    # "C_surface_microbe",
                    # "C_soil_microbe",
                    # "C_slowsom",
                    # "C_passsom"
                # ]
                # if Symbol(k) in numOutFluxesBySymbol.keys()
            # ]
        # )

        # V_new = np.concatenate((X_new.reshape(n,1),np.array([rh]).reshape(1,1)), axis=0)
        # return V_new

    # return TimeStepIterator2(
        # initial_values=V_arr,
        # f=f,
    # )


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
        [V_init.rh]
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
        
        # rh_flux=[
            # numOutFluxesBySymbol[Symbol(k)](it,*X.reshape(n,))
            # for k in [
                # "C_abvstrlit",
                # "C_abvmetlit",
                # "C_belowstrlit",
                # "C_belowmetlit",
                # "C_surface_microbe",
                # "C_soil_microbe",
                # "C_slowsom",
                # "C_passsom"] 
            # if Symbol(k) in numOutFluxesBySymbol.keys()
        # ]
        
        # rh = np.array(rh_flux).sum()
        rh = np.sum(
            [
                numOutFluxesBySymbol[Symbol(k)](it,*X)
                for k in [
                    "C_abvstrlit",
                    "C_abvmetlit",
                    "C_belowstrlit",
                    "C_belowmetlit",
                    "C_surface_microbe",
                    "C_soil_microbe",
                    "C_slowsom",
                    "C_passsom"]
                if Symbol(k) in numOutFluxesBySymbol.keys()
            ]
        )
        
        V_new = np.concatenate(
            (
                X_new.reshape(n,1),
                np.array([rh]).reshape(1,1)
            )
            , axis=0
        )
        return V_new

    return TimeStepIterator2(
        initial_values=V_arr,
        f=f,
    )


def make_weighted_cost_func(
        obs: Observables
    ) -> Callable[[Observables],np.float64]:
    # first unpack the observation array into its parts
    #cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
    def costfunction(out_simu: np.ndarray) ->np.float64:
        # fixme
        #   as indicated by the fact that the function lives in this
        #   model-specific module it is not apropriate for (all) other models.
        #   There are model specific properties:
        #   1.) The weight for the different observation streams
        #
        number_of_ys=out_simu.cVeg.shape[0]
        number_of_ms=out_simu.rh.shape[0]

        J_obj1 = np.mean (( out_simu.cVeg - obs.cVeg )**2)/(2*np.var(obs.cVeg))
        J_obj2 = np.mean (( out_simu.cRoot - obs.cRoot )**2)/(2*np.var(obs.cRoot))
        J_obj3 = np.mean (( out_simu.cLitter - obs.cLitter )**2)/(2*np.var(obs.cLitter))
        J_obj4 = np.mean (( out_simu.cSoil -  obs.cSoil )**2)/(2*np.var(obs.cSoil))

        J_obj5 = np.mean (( out_simu.rh - obs.rh )**2)/(2*np.var(obs.rh))

        J_new = (J_obj1 + J_obj2 + J_obj3 + J_obj4)/200 + J_obj5/4
        # to make this special costfunction comparable (in its effect on the
        # acceptance rate) to the general costfunction proposed by Feng we
        # rescale it by a factor
        return J_new*400
    return costfunction


def make_param_filter_func(
        c_max: EstimatedParameters,
        c_min: EstimatedParameters 
        ) -> Callable[[np.ndarray], bool]:

    # find position of beta_leaf and beta_wood
    beta_leaf_ind=EstimatedParameters._fields.index("beta_leaf")
    beta_wood_ind=EstimatedParameters._fields.index("beta_wood")

    def isQualified(c):
        cond1 =  (c >= c_min).all() 
        cond2 =  (c <= c_max).all() 
        cond3 =  c[beta_leaf_ind]+c[beta_wood_ind] <= 0.99  
        return (cond1 and cond2 and cond3)
        
    return isQualified

# this function is deprecated - see general helpers traceability_iterator
# def make_traceability_iterator(mvs,dvs,cpa,epa):


def numeric_X_0(mvs,dvs,cpa,epa):
    # This function creates the startvector for the pools
    # It can be used inside param_2_res and for other iterators that
    # track all carbon stocks
    apa = {**cpa._asdict(), **epa._asdict()}
    par_dict=gh.make_param_dict(mvs,cpa,epa)
    X_0_dict={
        "C_leaf": apa['C_leaf_0'],     
        "C_root": apa['cRoot_0'],     
        "C_wood": apa['cVeg_0'] - (apa['C_leaf_0'] +  apa['cRoot_0']),  
        "C_abvstrlit": apa['C_abvstrlit_0'],
        "C_abvmetlit": apa['C_abvmetlit_0'],
        "C_belowstrlit": apa["C_blwstrlit_0"],
        "C_belowmetlit":apa["cLitter_0"]- apa["C_abvstrlit_0"] - apa["C_abvmetlit_0"] - apa["C_blwstrlit_0"],
        "C_surface_microbe":apa["C_surfacemic_0"],
        "C_soil_microbe": apa["C_soilmic_0"],
        "C_slowsom":apa["C_slow_0"],
        "C_passsom":apa["cSoil_0"] - apa["C_surfacemic_0"] - apa["C_soilmic_0"] - apa["C_slow_0"]
    }
    X_0= np.array(
        [
            X_0_dict[str(v)] for v in mvs.get_StateVariableTuple()
        ]
    ).reshape(len(X_0_dict),1)
    return X_0


def nc_file_name(nc_var_name,experiment_name="SDGVM_S2_"):
    return experiment_name+"{}.nc".format(nc_var_name)


def nc_global_mean_file_name(nc_var_name):
    return experiment_name+"{}_gm.nc".format(nc_var_name)

# +
# def get_global_mean_vars(dataPath):
#     o_names=Observables._fields
#     d_names=Drivers._fields
#     names = o_names + d_names 
    
    
#     if all([dataPath.joinpath(nc_global_mean_file_name(vn)).exists() for vn in names]):
#         print(""" Found cached global mean files. If you want to recompute the global means
#             remove the following files: """
#         )
#         for vn in names:
#             print( dataPath.joinpath(nc_global_mean_file_name(vn)))

#         def get_cached_global_mean(vn):
#             return gh.get_nc_array(dataPath.joinpath(nc_global_mean_file_name(vn)),vn)
    
#         return (
#             Observables(*map(get_cached_global_mean, o_names)),
#             Drivers(*map(get_cached_global_mean,d_names))
#         )

#     else:
#         # we now check if any of the arrays has a time lime containing nan values 
#         # APART FROM values that are already masked by the fillvalue
#         print("computing masks to exclude pixels with nan entries, this may take some minutes...")
#         def f(vn):
#             path = dataPath.joinpath(nc_file_name(vn))
#             ds = nc.Dataset(str(path))
#             #scale fluxes vs pools
#             var =ds.variables[vn]
#             return gh.get_nan_pixel_mask(var)

#         masks=[ f(name)    for name in names ]
#         # We compute the common mask so that it yields valid pixels for ALL variables 
#         combined_mask= reduce(lambda acc,m: np.logical_or(acc,m),masks)
#         print("computing means, this may also take some minutes...")

#         def compute_and_cache_global_mean(vn):
#             path = dataPath.joinpath(nc_file_name(vn))
#             ds = nc.Dataset(str(path))
#             vs=ds.variables
#             lats= vs["latitude"].__array__()
#             lons= vs["longitude"].__array__()
#             print(vn)
#             var=ds.variables[vn]
#             # check if we have a cached version (which is much faster)
#             gm_path = dataPath.joinpath(nc_global_mean_file_name(vn))

#             gm=gh.global_mean_var(
#                     lats,
#                     lons,
#                     combined_mask,
#                     var
#             )
#             gh.write_global_mean_cache(
#                     gm_path,
#                     gm,
#                     vn
#             )
#             return gm #* 86400 if vn in ["npp", "rh"] else gm
    
#         #map variables to data
#         return (
#             Observables(*map(compute_and_cache_global_mean, o_names)),
#             Drivers(*map(compute_and_cache_global_mean, d_names))
#         )

# +
# def get_global_mean_vars(dataPath):
#     o_names=Observables._fields
#     d_names=Drivers._fields
#     names = o_names + d_names 
    
#     if all([dataPath.joinpath(nc_global_mean_file_name(vn)).exists() for vn in names]):
#         print(""" Found cached global mean files. If you want to recompute the global means
#             remove the following files: """
#         )
#         for vn in names:
#             print( dataPath.joinpath(nc_global_mean_file_name(vn)))

#         def get_cached_global_mean(vn):
#             gm = gh.get_nc_array(dataPath.joinpath(nc_global_mean_file_name(vn)),vn)
#             return gm * 86400 if vn in ["npp", "rh"] else gm

#         #map variables to data
#         odvs=Drivers(*map(get_cached_global_mean, d_names))
#         obss=Observables(*map(get_cached_global_mean, o_names))
#         dvs=Drivers(
#             npp=odvs.npp
#             #rso=odvs.mrso,
#             #as=odvs.tas
#         )
        
#         return (
#             obss,
#             dvs
#             #Observables(*map(get_cached_global_mean, o_names)),
#             #OrgDrivers(*map(get_cached_global_mean,d_names))
#         )

#     else:
#         # we now check if any of the arrays has a time lime containing nan values 
#         # APART FROM values that are already masked by the fillvalue
#         print("computing masks to exclude pixels with nan entries, this may take some minutes...")
#         def f(vn):
#             path = dataPath.joinpath(nc_file_name(vn))
#             ds = nc.Dataset(str(path))
#             #scale fluxes vs pools
#             var =ds.variables[vn]
#             return gh.get_nan_pixel_mask(var)

#         masks=[ f(name)    for name in names ]
#         # We compute the common mask so that it yields valid pixels for ALL variables 
#         combined_mask= reduce(lambda acc,m: np.logical_or(acc,m),masks)
#         print("computing means, this may also take some minutes...")

#         def compute_and_cache_global_mean(vn):
#             path = dataPath.joinpath(nc_file_name(vn))
#             ds = nc.Dataset(str(path))
#             vs=ds.variables
#             lats= vs["latitude"].__array__()
#             lons= vs["longitude"].__array__()
#             print(vn)
#             var=ds.variables[vn]
#             # check if we have a cached version (which is much faster)
#             gm_path = dataPath.joinpath(nc_global_mean_file_name(vn))

#             gm=gh.global_mean_var(
#                     lats,
#                     lons,
#                     combined_mask,
#                     var
#             )
#             gh.write_global_mean_cache(
#                     gm_path,
#                     gm,
#                     vn
#             )
#             return gm * 86400 if vn in ["npp", "rh"] else gm
        
#         #map variables to data
#         odvs=Drivers(*map(compute_and_cache_global_mean, d_names))
#         obss=Observables(*map(compute_and_cache_global_mean, o_names))
#         dvs=Drivers(
#             npp=odvs.npp
#             #mrso=odvs.mrso,
#             #tas=odvs.tas
#         )
    
#         return (
#             obss,
#             dvs
#             #Observables(*map(compute_and_cache_global_mean, o_names)),
#             #Drivers(*map(compute_and_cache_global_mean, d_names))            
#         )

# +
# def make_sim_day_2_day_since_a_D(conf_dict):
#     # this function is extremely important to syncronise our results
#     # because our data streams start at different times the first day of 
#     # a simulation day_ind=0 refers to different dates for different models
#     # we have to check some assumptions on which this calculation is based
#     # for jules the data points are actually spaced monthly with different numbers of days
#     ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("SDGVM_S2_cVeg.nc")))
#     times = ds.variables["time"]
#     # we have to check some assumptions on which this calculation is based
#     # for jules the data points are actually spaced with different numbers of days between monthly
#     # data point
#     # we can see this by looking at the first 24 months
#     # for i in range(24):
#     #     print((times[i + 1] - times[i])/(3600 * 24))


#     ts = times[0] #time of first observation in hours since 1900-01-01 00:00:00
#     td = ts / 24 #in days since 1900-01-01 00:00:00
    
#     import datetime as dt
#     ad = dt.date(1, 1, 1) # first of January of year 1 
#     sd = dt.date(1900, 1, 16)
#     td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    

#     def f(day_ind: int)->int:
#         return day_ind+td_aD

#     return f
# -

def start_date():
    ## this function is important to syncronise our results
    ## because our data streams start at different times the first day of 
    ## a simulation day_ind=0 refers to different dates for different models
    ## we have to check some assumptions on which this calculation is based
    ## Here is how to get these values
    #ds=nc.Dataset(str(Path(conf_dict['dataPath']).joinpath("SDGVM_S2_npp.nc")))
    #times = ds.variables["time"]
    #tm = times[0] #time of first observation in Months_since_1900-01 # print(times.units)
    #td = ts/24  #in days since_1900-01-01 
    #import datetime as dt
    #ad = dt.date(1, 1, 1) # first of January of year 1 
    #sd = dt.date(1900, 1, 1)
    #td_aD = td+(sd - ad).days #first measurement in days_since_1_01_01_00_00_00
    ## from td_aD (days since 1-1-1) we can compute the year month and day
    return gh.date(
        year=1900, 
        month=1,
        day=16
    )

data_str = namedtuple( # data streams available in the model
        'data_str',
        ["cVeg", "cLitter", "cSoil", "gpp", "npp", "ra", "rh"]
        )
        
def get_global_mean_vars_all(experiment_name):
        return(
            gh.get_global_mean_vars_all(model_folder="Aneesh_SDGVM", 
                            experiment_name=experiment_name,
                            lat_var="latitude",
                            lon_var="longitude",
                            ) 
        )
        
################ function for computing global mean for custom data streams ###################
    
# def get_global_mean_vars_all(experiment_name="SDGVM_S2_"):
    
    # def nc_file_name(nc_var_name, experiment_name="SDGVM_S2_"):
        # return experiment_name+"{}.nc".format(nc_var_name)

    # def nc_global_mean_file_name(nc_var_name, experiment_name="SDGVM_S2_"):
        # return experiment_name+"{}_gm_all.nc".format(nc_var_name)

    # data_str = namedtuple( # data streams available in the model
        # 'data_str',
        # ["cVeg", "cLitter", "cSoil", "gpp", "npp", "ra", "rh"]
        # )
        
    # names = data_str._fields
    # conf_dict = gh.confDict("Aneesh_SDGVM")
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
        # template = nc.Dataset(dataPath.joinpath("SDGVM_S2_cSoil.nc")).variables['cSoil'][0,:,:].mask
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
