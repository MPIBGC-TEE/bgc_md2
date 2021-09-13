from tqdm import tqdm
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
        #   model-specific module it is not apropriate for (all) other models.
        #   There are model specific properties:
        #   1.) The weight for the different observation streams
        #   
        tot_len=out_simu.shape[0]
        # we assume the model output to be in the same shape and order 
        # as the obeservation
        # this convention has to be honored by the forwar_simulation as well
        # which in this instance already compresses the 3 different litter pools
        # to c_litter and the 3 different soil pools to one
        c_simu = out_simu[:,0:5] 
        
        # we assume the rh  part to be in the remaining columns again
        # this convention has to be honored by the forwar_simulation as well
        rh_simu = out_simu[:,5:]
        #from IPython import embed; embed()

        J_obj1 = np.mean (( c_simu[:,0] - cleaf[0:tot_len] )**2)/(2*np.var(cleaf[0:tot_len]))
        J_obj2 = np.mean (( c_simu[:,1] - croot[0:tot_len] )**2)/(2*np.var(croot[0:tot_len]))
        J_obj3 = np.mean (( c_simu[:,2] - cwood[0:tot_len] )**2)/(2*np.var(cwood[0:tot_len]))
        J_obj4 = np.mean (( c_simu[:,3]-  clitter[0:tot_len] )**2)/(2*np.var(clitter[0:tot_len]))
        J_obj5 = np.mean (( c_simu[:,4]-  csoil[0:tot_len] )**2)/(2*np.var(csoil[0:tot_len]))
        
        J_obj6 = np.mean (( rh_simu[:,0] - rh[0:tot_len] )**2)/(2*np.var(rh[0:tot_len]))
        
        J_new = (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 )/200+ J_obj6/4
        return J_new
    return costfunction     


def make_matrix_simu(
            cleaf_0,
            croot_0,
            cwood_0,
            clitter_0,
            csoil_0,
            npp,
            tot_len,
            clay, 
            silt,
            lig_wood
    ) -> Callable[[np.ndarray], np.ndarray]: 

    def matrix_simu(pa):
        days=[31,28,31,30,31,30,31,31,30,31,30,31]
        # Construct B matrix 
        beta1=pa[0]; beta2=pa[1]; beta3= 1- beta1- beta2
        B = np.array([beta1, beta2, beta3, 0, 0, 0, 0,0,0]).reshape([9,1])   # allocation
        # Now construct A matrix
        lig_leaf = pa[2]
    
        f41 = pa[3]; f42 = pa[4]; f51 = 1-f41; f52 = 1-f42; f63 = 1;
        f74 = 0.45; f75 = 0.45*(1-lig_leaf); 
        f85 = 0.7*lig_leaf; f86 = 0.4*(1-lig_wood);
        f96 = 0.7*lig_wood;  
        f87=(0.85 - 0.68 * (clay+silt))* (0.997 - 0.032*clay)
        f97=(0.85 - 0.68 * (clay+silt))* (0.003 + 0.032*clay)
        f98=0.45 * (0.003 + 0.009 *clay)
    
        A = np.array([  -1,   0,   0,   0,   0,   0,   0,   0,   0,
                         0,  -1,   0,   0,   0,   0,   0,   0,   0,
                         0,   0,  -1,   0,   0,   0,   0,   0,   0,
                       f41, f42,   0,  -1,   0,   0,   0,   0,   0,
                       f51, f52,   0,   0,  -1,   0,   0,   0,   0,
                         0,   0, f63,   0,   0,  -1,   0,   0,   0,
                         0,   0,   0, f74, f75,   0,  -1,   0,   0,
                         0,   0,   0,   0, f85, f86, f87,  -1,   0,
                         0,   0,   0,   0,   0, f96, f97, f98,  -1 ]).reshape([9,9])   # tranfer
    
        #turnover rate per day of pools: 
        temp = [pa[5],pa[6],pa[7], pa[8],pa[8]/(5.75*np.exp(-3*pa[2])), pa[8]/20.6, pa[9],pa[10], pa[11]]
        K = np.zeros(81).reshape([9, 9])
        for i in range(0, 9):
            K[i][i] = temp[i]
          
        x_fin=np.zeros((tot_len,9))
        rh_fin=np.zeros((tot_len,1))
        # leaf, root , wood, metabolic, structural, CWD, microbial, slow, passive 
        x_init = np.array([cleaf_0,croot_0,cwood_0,pa[12],pa[13],clitter_0-pa[12]-pa[13],pa[14],csoil_0- pa[14] - pa[15], pa[15]]).reshape([9,1])   # Initial carbon pool size
        X=x_init   # initialize carbon pools 
        jj=0
        for m in tqdm(np.arange(0,tot_len)):
          npp_in = npp[jj] 
          co2_rh = 0  
          for d in np.arange(0,days[m%12]):
              X=X + B*npp_in + np.array(A@K@X).reshape([9,1])
              co2_rate = [0,0,0, (1-f74)*K[3,3],(1-f75-f85)*K[4,4],(1-f86-f96)*K[5,5], (1- f87 - f97)*K[6,6], (1-f98)*K[7,7], K[8,8] ]
              co2=np.sum(co2_rate*X.reshape(1,9))
              co2_rh = co2_rh + co2/days[m%12]   # monthly average rh 
          x_fin[jj,:]=X.reshape(1,9)
          rh_fin[jj,0]=co2_rh
          jj= jj+1
             
        # We create an output that has the same shape
        # as the obvervations to make the costfunctions 
        # easier. 
        # To this end we project our 10 output variables
        # onto the 6 data streams
        c_litter = np.sum(x_fin[:,3:6],axis=1).reshape(tot_len,1)
        c_soil = np.sum(x_fin[:,6:9],axis=1).reshape(tot_len,1)
        #from IPython import embed; embed()
        out_simu = np.concatenate(
            [
                x_fin[:,0:3], # the first 3 pools are used as they are
                c_litter,
                c_soil,
                rh_fin
            ]
            ,axis=1
        )
        return out_simu

    return matrix_simu
