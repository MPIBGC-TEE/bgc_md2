from bgc_md2.resolve.mvars import NumericStartValueDict
from sympy import Symbol, var
from functools import lru_cache
from CompartmentalSystems.smooth_model_run import SmoothModelRun
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
import CompartmentalSystems.helpers_reservoir as hr
import netCDF4 as nc
import numpy as np
from pathlib import Path
from collections import namedtuple
from functools import reduce
from bgc_md2.models.cable_yuanyuan.source import mvs 

# The unestimated parameters  
UnEstimatedParameters = namedtuple(
    "UnEstimatedParameters",
    [
        'cleaf_0',
        'croot_0',
        'cwood_0',
        'clitter_0',
        'csoil_0',
        'clay',
        'silt',
        'lig_wood',
        'f_wood2CWD',
        'f_metlit2mic',
        'npp',
        'number_of_months'
    ]
)

# The estimated parameters include start values
# for some of the pools
EstimatedParameters = namedtuple(
    "EstiamatedParameters",
    [
        "beta_leaf",
        "beta_root",
        "lig_leaf",
        "f_leaf2metlit",
        "f_root2metlit",
        "k_leaf",
        "k_root",
        "k_wood",
        "k_metlit",
        "k_mic",
        "k_slowsom",
        "k_passsom",
        "cmet_init",
        "cstr_init",
        "cmic_init",
        "cpassive_init"
    ]
)

# the model parameters are parameters are just
# constants in the model
ModelParameters = namedtuple(
    "ModelParameters",
    [
        "beta_leaf",
        "beta_root",
        "lig_leaf",
        "f_leaf2metlit",
        "f_root2metlit",
        "k_leaf",
        "k_root",
        "k_wood",
        "k_metlit",
        "k_mic",
        "k_slowsom",
        "k_passsom",
    ]
)
# This namedtuple defines the order of the parameters
InitialValues = namedtuple(
    'IntitialValues',
    [
        'C_leaf',
        'C_root',
        'C_wood',
        'C_metlit',
        'C_stlit',
        'CWD',
        'C_mic',
        'C_slowsom',
        'C_passsom',
    ]
)

def construct_matrix_func_sym(pa):
    # we create a parameterdict for the fixed values
    # and extend it by the parameters provided 
    symbol_names = mvs.get_BibInfo().sym_dict.keys()   
    for name in symbol_names:
        var(name)
    parDict = {
        clay: 0.2028,
        silt: 0.2808,
        lig_wood: 0.4,
        f_wood2CWD: 1,
        f_metlit2mic: 0.45,
    #    NPP: npp_in
    }
    model_params = {Symbol(k): v for k,v in pa._asdict().items()}
    parDict.update(model_params)
    B_func = hr.numerical_array_func(
            state_vector = mvs.get_StateVariableTuple(),
            time_symbol=mvs.get_TimeSymbol(),
            expr=mvs.get_CompartmentalMatrix(),
            parameter_dict=parDict,
            func_dict={}
    )
    # in the general nonautonomous nonlinear B_func is a function of t,x
    # although for this example it does not depend on either t, nor x.
    return B_func 




#def n_day_matrix_simu(pa, Xnt, npp_in, n):
#    # Construct b vector for allocation
#    b_func  = construct_allocation_vector_func(pa)
#    # Now construct B matrix B=A*K
#    B_func  = construct_matrix_func(pa)
#    def f(it,V):
#        X = V[0:9]
#        co2 = V[9]
#        b = b_func(it,X)
#        B = B_func(it,X)
#        X_new = X + b*npp_in + B@X
#        co2_new = respiration_from_compartmental_matrix(B,X) #accumulate respiration
#        V_new = np.concatenate((X_new,np.array([co2_new]).reshape(1,1)), axis=0)
#        return V_new
#
#
#    X_0 = np.array(Xnt).reshape(9,1) #transform the tuple to a matrix
#   
#    co2_0 = np.array([0]).reshape(1,1)
#    tsi = TimeStepIterator2(
#            initial_values=np.concatenate((X_0, co2_0), axis=0),
#            f=f,
#            max_it=n
#    )
#    
#    RES = np.concatenate([ np.transpose(val) for val in tsi], axis=0)  
#    return RES

def monthly_matrix_simu(mpa, Xnt, npp, ns):
    # here we want to store only monthly values
    # although the computation takes place on a dayly basis
    # Construct b vector for allocation
    b_func  = construct_allocation_vector_func(mpa)
    # Now construct B matrix B=A*K
    B_func  = construct_matrix_func(mpa)
    def f(it,V):
        X = V[0:9]
        co2 = V[9]
        b = b_func(it,X)
        B = B_func(it,X)
        X_new = X + b*npp[it] + B@X
        # here we add the co2 up in every step which allows us not to save it
        #co2_new = co2 + respiration_from_compartmental_matrix(B,X) #accumulate respiration
        co2_new =  respiration_from_compartmental_matrix(B,X) #accumulate respiration
        V_new = np.concatenate((X_new,np.array([co2_new]).reshape(1,1)), axis=0)
        return V_new


    X_0 = np.array(Xnt).reshape(9,1) #transform the tuple to a matrix
   
    co2_0 = np.array([0]).reshape(1,1)
    
    day_indices = month_2_day_index(ns)
    tsi = TimeStepIterator2(
            initial_values=np.concatenate((X_0, co2_0), axis=0),
            f=f,
            max_it=max(day_indices)+1
    )
    
   
    #values = [ np.transpose(val) for val in tsi if tsi.i in day_indices]
    #values =  [v for i, v in enumerate(tsi) if i in day_indices]
    def g(acc, el):
        xs,co2s,acc_co2 = acc
        i, v = el
        d_pools = v[0:-1,:]
        d_co2=v[-1:,:]
        acc_co2 += d_co2
        if i in day_indices:
            xs += [d_pools] 
            co2s +=[acc_co2]
            acc_co2=np.array([0.0]).reshape(1,1)
            
        acc = (xs,co2s,acc_co2)
        return acc
    xs,co2s,_ =  reduce(g,enumerate(tsi),([],[],0))
            
    def h(tup):
        x, co2 = tup
        return np.transpose(np.concatenate([x,co2]))

    values_with_monthly_co2 = [v for v in  map(h,zip(xs,co2s))]
    RES = np.concatenate(values_with_monthly_co2 , axis=0)  
    return RES

            


def one_day_matrix_simu(pa,day,Xnt,npp_in):
    # Construct b vector for allocation
    b_func  = construct_allocation_vector_func(pa)
    # Now construct B matrix B=A*K
    B_func  = construct_matrix_func(pa)
    X=np.array(Xnt).reshape(9,1) #transform the tuple to a matrix
    B = B_func(day,X)
    b = b_func(day,X)
    X=X + b*npp_in + B@X
    co2 = respiration_from_compartmental_matrix(B,X)
    return X, co2 


def number_of_months(nyears):
    return nyears * 12

def matrix_simu_maker_org(cpa):
    # this function produces the function f of the estimated parameters 
    # which is called by the mcmc function to produce the forward solution

    def f(pa):
        days=[31,28,31,30,31,30,31,31,30,31,30,31]
        x_fin=np.zeros((cpa.number_of_months,9))
        rh_fin=np.zeros((cpa.number_of_months,1))
        x_init = np.array([
            cpa.cleaf_0,
            cpa.croot_0,
            cpa.cwood_0,
            pa[12],
            pa[13],
            cpa.clitter_0-pa[12]-pa[13],
            pa[14],
            cpa.csoil_0 - pa[14] - pa[15],
            pa[15]
        ]).reshape([9,1])   # Initial carbon pool size
        X=x_init   # initialize carbon pools 
        for m in np.arange(0,cpa.number_of_months):
          npp_in = cpa.npp[m] 
          co2_rh = 0  
          for d in np.arange(0,days[m]):
              X,co2 = one_day_matrix_simu(pa,d,X,npp_in)
              co2_rh = co2_rh + co2/days[m%12]   # monthly average rh 
          x_fin[m,:]=X.reshape(1,9)
          rh_fin[m,0]=co2_rh
             
        return x_fin, rh_fin     

    return f

def matrix_simu_maker(cpa):
    # this function produces the function f(epa) of the estimated parameters epa
    # which called by the mcmc function

    def matrix_simu(epa):
        x_init = np.array([
            cpa.cleaf_0,
            cpa.croot_0,
            cpa.cwood_0,
            epa[12],
            epa[13],
            cpa.clitter_0-epa[12]-epa[13],
            epa[14],
            cpa.csoil_0 - epa[14] - epa[15],
            epa[15]
        ]).reshape([9,1])   # Initial carbon pool size
        X=x_init   # initialize carbon pools 
        RES = monthly_matrix_simu()

    return matrix_simu

def mcmc(data,start_pa, nyears):
    #from IPython import embed; embed()

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
        tl = number_of_months(nyears) 
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
                 X,co2 = one_day_matrix_simu(pa,d,X,npp_in)
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
        # converte to a named tuple of type EstimatedParameters
        # to be able to use the parameter names instead of the 
        # positions inside matrix_simu
        pa_new = EstimatedParameters(*c_new)
        x_simu,rh_simu = matrix_simu(pa_new)
    
        tl=number_of_months(nyears)
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
     
def matrix_simul_from_symbolic(
    # for this model we also want the accumulated respiration
        pa,
        X0,
        npp_in,
        times
    ):
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
    smr = SmoothModelRun(
            srm,
            parameter_dict=parDict,
            start_values=start_values,
            times=times,
            func_set={}
    )
    sol = smr.solve()
    RESP_vec = smr.acc_gross_external_output_vector()
    RESP = np.sum(RESP_vec,axis=1)
    # in order to attach it to the solution we 
    # have to have equal length. (in timesteps)
    # since the accumulated respiration vector
    # does not return the zeros for the t=0
    # but the solution contains the start values
    # we add the zero at the first time step 
    RESP_w0 = np.concatenate([np.array([0]),RESP]).reshape(sol.shape[0],1)
    # have to add zeros at the start (
    result = np.concatenate([sol,RESP_w0], axis=1)
    return result 
