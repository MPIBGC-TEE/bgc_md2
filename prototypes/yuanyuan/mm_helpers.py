from bgc_md2.resolve.mvars import NumericStartValueDict
from sympy import Symbol
from CompartmentalSystems.smooth_model_run import SmoothModelRun
import netCDF4 as nc
import numpy as np
from pathlib import Path

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

def one_step_matrix_simu(pa,X,npp_in):
    clay = 0.2028
    silt = 0.2808
    lig_wood = 0.4
    # Construct B matrix 
    beta1 = pa.beta_leaf
    beta2 = pa.beta_root
    beta3 = 1- beta1- beta2
    B = np.array([beta1, beta2, beta3, 0, 0, 0, 0,0,0]).reshape([9,1])   # allocation
    # Now construct A matrix
    lig_leaf = pa.lig_leaf

    f41 = pa.f_leaf2metlit
    f42 = pa.f_root2metlit
    f51 = 1 - f41
    f52 = 1 - f42
    f63 = 1
    f74 = 0.45
    f75 = 0.45 * (1 - lig_leaf)
    f85 = 0.7 * lig_leaf
    f86 = 0.4 * (1 - lig_wood)
    f96 = 0.7 * lig_wood
    f87=(0.85 - 0.68 * (clay+silt)) * (0.997 - 0.032 * clay)
    f97=(0.85 - 0.68 * (clay+silt)) * (0.003 + 0.032 * clay)
    f98=0.45 * (0.003 + 0.009 * clay)

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
    temp = [
        pa.k_leaf,
        pa.k_root,
        pa.k_wood,
        pa.k_metlit,
        pa.k_metlit/(5.75*np.exp(-3*pa.lig_leaf)),
        pa.k_metlit/20.6,
        pa.k_mic,
        pa.k_slowsom,
        pa.k_passsom
    ]
    K = np.zeros(81).reshape([9, 9])
    for i in range(0, 9):
        K[i][i] = temp[i]
      
    # 1 leaf
    # 2 root 
    # 3 wood
    # 4 metlit
    # 5 structural
    # 6 CWD
    # 7 mic
    # 8 slow
    # 9 passive 
    X=X + B*npp_in + np.array(A@K@X).reshape([9,1])
    co2_rate = [0,0,0, (1-f74)*K[3,3],(1-f75-f85)*K[4,4],(1-f86-f96)*K[5,5], (1- f87 - f97)*K[6,6], (1-f98)*K[7,7], K[8,8] ]
         
    co2=np.sum(co2_rate*X.reshape(1,9))
    return X, co2 

def tot_len(nyears):
    return nyears * 12


def mcmc(data,start_pa, nyears):

    npp, rh, clitter, csoil, cveg, cleaf, croot, cwood = data
    #===
    c_min=np.array([0.09, 0.09,0.09,0.01,0.01,  1/(2*365), 1/(365*10), 1/(60*365), 0.1/(0.1*365), 0.06/(0.137*365), 0.06/(5*365), 0.06/(222.22*365),clitter[0]/100,clitter[0]/100,csoil[0]/100,csoil[0]/2])
    c_max=np.array([1   ,    1,0.21,   1,   1,1/(0.3*365),1/(0.8*365),      1/365,   1/(365*0.1),  0.6/(365*0.137),  0.6/(365*5),  0.6/(222.22*365),    clitter[0],    clitter[0],  csoil[0]/3,csoil[0]])
    
    def GenerateParamValues(c_op):
       flag = True
       while (flag):
          c_new = c_op + (np.random.random((paramNum)) - 0.5)*(c_max - c_min)/10.0
          #c_new = c_op + (np.random.normal(0, 1, paramNum))*(c_max - c_min)/15.0
          if (isQualified(c_new)):
             flag = False
       return c_new
    
    def isQualified(c):
       flag = True
       for i in range(paramNum):
          if(c[i] > c_max[i] or c[i] < c_min[i]):
             flag = False
             break
          if(c[0] + c[1] > 0.99):
             flag = False
             break
       return flag

    def matrix_simu(pa):
        tl = tot_len(nyears) 
        days=[31,28,31,30,31,30,31,31,30,31,30,31]
        x_fin=np.zeros((tl,9))
        rh_fin=np.zeros((tl,1))
        jj=0
        x_init = np.array([cleaf[0],croot[0],cwood[0],pa[12],pa[13],clitter[0]-pa[12]-pa[13],pa[14],csoil[0]- pa[14] - pa[15], pa[15]]).reshape([9,1])   # Initial carbon pool size
        X=x_init   # initialize carbon pools 
        for y in np.arange(0,nyears):
           for m in np.arange(0,12):
             npp_in = npp[jj] 
             co2_rh = 0  
             for d in np.arange(0,days[m]):
                 X,co2 = one_step_matrix_simu(pa,X,npp_in)
                 co2_rh = co2_rh + co2/days[m]   # monthly average rh 
             x_fin[jj,:]=X.reshape(1,9)
             rh_fin[jj,0]=co2_rh
             jj= jj+1
             
        return x_fin, rh_fin     


    np.random.seed(seed=10)
    
    paramNum=len(start_pa)
    nsimu    = 20000
    
    upgraded=0;
    J_last = 400
    C_op = start_pa
    
    C_upgraded = np.zeros((paramNum, nsimu))
    J_upgraded = np.zeros((1, nsimu))
    
    for simu in range(nsimu):
        c_new = GenerateParamValues(C_op)
    
        x_simu,rh_simu = matrix_simu(c_new)
    
        tl=tot_len(nyears)
        J_obj1 = np.mean (( x_simu[:,0] - cleaf[0:tl] )**2)/(2*np.var(cleaf[0:tl]))
        J_obj2 = np.mean (( x_simu[:,1] - croot[0:tl] )**2)/(2*np.var(croot[0:tl]))
        J_obj3 = np.mean (( x_simu[:,2] - cwood[0:tl] )**2)/(2*np.var(cwood[0:tl]))
        J_obj4 = np.mean (( np.sum(x_simu[:,3:6],axis=1)- clitter[0:tl] )**2)/(2*np.var(clitter[0:tl]))
        J_obj5 = np.mean (( np.sum(x_simu[:,6:9],axis=1)- csoil[0:tl] )**2)/(2*np.var(csoil[0:tl]))
        
        J_obj6 = np.mean (( rh_simu[:,0] - rh[0:tl] )**2)/(2*np.var(rh[0:tl]))
        
        J_new= (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 )/200+ J_obj6/4
    
        delta_J =  J_last - J_new;
        
        randNum = np.random.uniform(0, 1)
        if (min(1.0, np.exp(delta_J)) > randNum):
                C_op=c_new;
                J_last=J_new;
                C_upgraded[:,upgraded]=C_op;
                J_upgraded[:,upgraded]=J_last; 
                upgraded=upgraded+1;
    
    df=pd.DataFrame(C_upgraded)
    df_j=pd.DataFrame(J_upgraded)
    #df.columns = ["std_tomcat","r_tomcat","std_lmdz","r_lmdz","std_jamstec","r_jamstec"]
    #df.index = ['all', 'ha', 'ar', 'sa','ds','hu'] 
    return df, df_j
     
def matrix_simul_from_symbolic(pa,X0,npp_in):
    from sympy import var
    
    from bgc_md2.models.cable_yuanyuan.source import mvs 
    symbol_names = mvs.get_BibInfo().sym_dict.keys()   
    for name in symbol_names:
        var(name)
   
    srm = mvs.get_SmoothReservoirModel()
    # we create a parameterdict for the fixed values
    # and extend it by the parameters provided 
    parDict = {
        clay: 0.2028,
        silt: 0.2808,
        lig_wood: 0.4,
        f_wood2CWD: 1,
        f_metlit2mic: 0.45,
        NPP: npp_in
    }
    #from IPython import embed; embed()
    model_params = {Symbol(k): v for k,v in pa._asdict().items()}

    parDict.update(model_params)

    nsv1 = {
        Symbol(k): v 
        for k,v in X0._asdict().items()
    }
    start_values=np.array(
        [
            nsv1[k] for k in mvs.get_StateVariableTuple()
        ]
    )
    #from IPython import embed; embed()
    smr = SmoothModelRun(
            srm,
            parameter_dict=parDict,
            start_values=start_values,
            times=np.array([0,1,2,3,4]),
            func_set={}
    )
    smr.solve()
