from typing import Callable
import netCDF4 as nc
import numpy as np
# to do
# 1.) make a namedtuple for the yycable data or use xarray to create a multifile dataset
def get_variables_from_files(dataPath):
    # Read NetCDF data  ******************************************************************************************************************************
    path = dataPath.joinpath("npp_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc")
    ds = nc.Dataset(str(path))
    var_npp = ds.variables['npp'][:,:,:]
    ds.close()
    
    path = dataPath.joinpath("rh_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc")
    ds = nc.Dataset(str(path))
    var_rh = ds.variables['rh'][:,:,:]
    ds.close()
    
    path = dataPath.joinpath("cLeaf_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc")
    ds = nc.Dataset(str(path))
    var_cleaf = ds.variables['cLeaf'][:,:,:]
    ds.close()
    
    path = dataPath.joinpath("cRoot_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc")
    ds = nc.Dataset(str(path))
    var_croot = ds.variables['cRoot'][:,:,:]
    ds.close()
    
    
    path = dataPath.joinpath("cVeg_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc")
    ds = nc.Dataset(str(path))
    var_cveg = ds.variables['cVeg'][:,:,:]
    ds.close()
    
    
    path = dataPath.joinpath("cSoil_Emon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc")
    ds = nc.Dataset(str(path))
    var_csoil = ds.variables['cSoil'][:,:,:]
    lat_data=ds.variables['lat'][:].data
    lon_data=ds.variables['lon'][:].data
    ds.close()
    
    path = dataPath.joinpath("cLitter_Lmon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc")
    ds = nc.Dataset(str(path))
    var_clitter = ds.variables['cLitter'][:,:,:]
    ds.close()

    return (var_npp, var_rh, var_cleaf, var_croot, var_cveg, var_csoil, var_clitter)       

def get_example_site_vars(dataPath):
    var_npp, var_rh, var_cleaf, var_croot, var_cveg, var_csoil, var_clitter = get_variables_from_files(dataPath)       
    # pick up 1 site   62.8125 W, 17.5S
    s = slice(None,None,None) # this is the same as : 
    t = s,58,159 # [t] = [:,58,159]
    npp= var_npp[t]* 86400   #   kg/m2/s kg/m2/day; 
    rh= var_rh[t]*86400;   # per s to per day  
    (
        clitter,
        csoil,
        cveg,
        cleaf,
        croot
    ) = map(
        lambda var: var[t],
        (
            var_clitter,
            var_csoil, 
            var_cveg,
            var_cleaf,
            var_croot
        )
    )
    cwood = cveg - cleaf - croot; 
    return (npp, rh, clitter, csoil, cveg, cleaf, croot, cwood)

def make_param_filter_func(
        c_max: np.ndarray,
        c_min: np.ndarray
        ) -> Callable[[np.ndarray], bool]:

    def isQualified(c):
        # fixme
        #   this function is model specific: It discards parameter proposals
        #   where beta1 and beta2 are >0.99
        paramNum = len(c)
        flag = True
        for i in range(paramNum):
           if(c[i] > c_max[i] or c[i] < c_min[i]):
              flag = False
              break
           if(c[0] + c[1] > 0.99):
              flag = False
              break
        return flag
    
    return isQualified

def make_weighted_cost_func(
        obs: np.ndarray
    ) -> Callable[[np.ndarray],np.float64]:
    # first unpack the observation array into its parts
    cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
    def costfunction(out_simu: np.ndarray) ->np.float64:
        # fixme 
        #   as indicated by the fact that the function lives in this  
        #   model-specific module it is not apropriate for other models 
        #   It is determined by the relationship of the 6 observabal data streams
        #   and the available computed results (e.g. the combination of the 
        #   3 computed litter pools into one or the 3 computed soil pools into
        #   one. 

        # as well as the 
        n_pools=9
        tot_len=out_simu.shape[0]
        # we assume the x part to be in the first n_pools columns
        # this convention has to be honored by the forwar_simulation as well
        x_simu = out_simu[:,0:n_pools] 
        
        # we assume the rh  part to be in the remaining columns again
        # this convention has to be honored by the forwar_simulation as well
        rh_simu = out_simu[:,n_pools:]
        #from IPython import embed; embed()

        J_obj1 = np.mean (( x_simu[:,0] - cleaf[0:tot_len] )**2)/(2*np.var(cleaf[0:tot_len]))
        J_obj2 = np.mean (( x_simu[:,1] - croot[0:tot_len] )**2)/(2*np.var(croot[0:tot_len]))
        J_obj3 = np.mean (( x_simu[:,2] - cwood[0:tot_len] )**2)/(2*np.var(cwood[0:tot_len]))
        J_obj4 = np.mean (( np.sum(x_simu[:,3:6],axis=1)- clitter[0:tot_len] )**2)/(2*np.var(clitter[0:tot_len]))
        J_obj5 = np.mean (( np.sum(x_simu[:,6:9],axis=1)- csoil[0:tot_len] )**2)/(2*np.var(csoil[0:tot_len]))
        
        J_obj6 = np.mean (( rh_simu[:,0] - rh[0:tot_len] )**2)/(2*np.var(rh[0:tot_len]))
        
        J_new = (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 )/200+ J_obj6/4
        return J_new
    return costfunction    

    
