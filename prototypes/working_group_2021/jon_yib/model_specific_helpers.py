from typing import Callable
import netCDF4 as nc
import numpy as np
from collections import namedtuple
from functools import reduce
from general_helpers import day_2_month_index, month_2_day_index, months_by_day_arr, TimeStepIterator2, respiration_from_compartmental_matrix

pseudo_days_per_month = 30
# fixme:
# Your parameters will most likely differ but you can still use the
# destinctions between different sets of parameters. The aim is to make
# the different tasks in the code more obvious. In principal you could
# have a lot of overlapping sets and just have to keep them consistent. 
# 'namedtuples' can be used just as normal tuples by functions
# that are not aware of the names. They can still use the positions like 
# in the original code

# @Kostia and the 'R'tists: 
# It is not necessary to replicate the complete functionality of the #
# namedtuple classes. A simple 'R'proximation is a list with named entries 
# pa=list(C_leaf_0=2,...)


# This set is used by the functions that produce the 
# specific ingredients (functions) that will be run by
# mcmc alg.
UnEstimatedParameters = namedtuple(
    "UnEstimatedParameters",
    [
        'C_soil_0',
        'C_veg_0',
        'rh_0',
        'ra_0',
        'npp',
        'clay',
        'silt',
        'nyears'
    ]
)

# This set is used (as argument) by the functions that are called
# inside the mcmc
EstimatedParameters = namedtuple(
    "EstiamatedParameters",
    [
        "beta_leaf",        # 0         
        "beta_root",        # 1      
        "k_leaf",           # 2      
        "k_root",           # 3         
        "k_wood",           # 4
        "k_cwd",            # 5      
        "k_samet",          # 6      
        "k_sastr",          # 7      
        "k_samic",          # 8      
        "k_slmet",          # 9      
        "k_slstr",          # 10      
        "k_slmic",          # 11      
        "k_slow",           # 12      
        "k_arm",            # 13      
        "f_samet_leaf",     # 14      
        "f_slmet_root",     # 15      
        "f_samic_cwd",      # 16     
        "C_leaf_0",         # 17      
        "C_root_0",         # 18      
        "C_cwd_0",          # 19      
        "C_samet_0",        # 20      
        "C_sastr_0",        # 21      
        "C_samic_0",        # 22      
        "C_slmet_0",        # 23      
        "C_slstr_0",        # 24      
        "C_slmic_0",        # 25      
        "C_slow_0"          # 26      
    ]
)

# This is the set off all 
_Parameters=namedtuple(
        "Parameters",
        EstimatedParameters._fields + UnEstimatedParameters._fields
)
class Parameters(_Parameters):
    @classmethod
    def from_EstimatedParametersAndUnEstimatedParameters(
            cls,
            epa :EstimatedParameters,
            cpa :UnEstimatedParameters
        ):
        return cls(*(epa + cpa))


# This set defines the order of the c pools
# The order is crucial for the compatibility
# with the matrices (B and b) If you change ist
# the matrix changes
StateVariables = namedtuple(
    'StateVariables',
    [
        'leaf',
        'root',
        'wood',
        'cwd',
        'samet',
        'sastr',
        'samic',
        'slmet',
        'slstr',
        'slmic',
        'slow',
        'arm'
    ]
)
Observables = namedtuple(
    'Observables',
    [
        'c_veg',
        'c_soil',
       #'a_respiration',
        'h_respiration'
    ]
)

# We define another set of parameters which describes
# the parameters of the matrices A,K and the vector b
# and drivers like npp (in form of arrays)
# but does not include start values and hyperparameters like the 'number_of_months'
# This distinction is helpful for the forward simulation where the
# distinction between estimated and constant is irrelevant.
ModelParameters = namedtuple(
    "ModelParameters",
    [
        name for name in Parameters._fields 
        if name not in [
            'C_leaf_0',
            'C_root_0',
            'C_wood_0',
            'clitter_0',
            'csoil_0',
            "C_metlit_0",
            "CWD_0",
            "C_mic_0",
            "C_passom_0",
            'number_of_months'
        ] 
    ]
)
# to do
# 1.) make a namedtuple for the yycable data or use xarray to create a multifile dataset
def monthly_to_yearly(monthly):
    #TRENDY specific - months weighted like all months are 30 days
    if len(monthly.shape) > 1:
        sub_arrays=[monthly[i*12:(i+1)*12,:,:] for i in range(int(monthly.shape[0]/12))]
    else:
        sub_arrays=[monthly[i*12:(i+1)*12,] for i in range(int(monthly.shape[0]/12))]
    return np.stack(list(map(lambda sa:sa.mean(axis=0), sub_arrays)), axis=0)

def pseudo_daily_to_yearly(daily):
    # compute a yearly average from pseudo daily data
    # for one data point
    pseudo_days_per_year = pseudo_days_per_month*12 
    sub_arrays=[daily[i*pseudo_days_per_year:(i+1)*pseudo_days_per_year,:] for i in range(int(daily.shape[0]/pseudo_days_per_year))]
    return np.stack(list(map(lambda sa:sa.mean(axis=0), sub_arrays)), axis=0)

def get_variables_from_files(dataPath):
    # Read NetCDF data  ******************************************************************************************************************************
    path = dataPath.joinpath("YIBs_S2_Monthly_npp.nc")
    ds = nc.Dataset(str(path))
    var_npp = ds.variables['npp'][:,:,:]
    ds.close()
    
    path = dataPath.joinpath("YIBs_S2_Monthly_rh.nc")
    ds = nc.Dataset(str(path))
    var_rh = ds.variables['rh'][:,:,:]
    ds.close()
    
    path = dataPath.joinpath("YIBs_S2_Monthly_ra.nc")
    ds = nc.Dataset(str(path))
    var_ra = ds.variables['ra'][:,:,:]
    ds.close()

    path = dataPath.joinpath("YIBs_S2_Annual_cVeg.nc")
    ds = nc.Dataset(str(path))
    var_cveg = ds.variables['cVeg'][:,:,:]
    ds.close()
    
    path = dataPath.joinpath("YIBs_S2_Annual_cSoil.nc")
    ds = nc.Dataset(str(path))
    var_csoil = ds.variables['cSoil'][:,:,:]
    ds.close()
    
    return (var_npp, var_rh, var_ra, var_cveg, var_csoil)       

def get_example_site_vars(dataPath):
    var_npp, var_rh, var_ra, var_cveg, var_csoil = get_variables_from_files(dataPath)       
    # pick up 1 site   62.8125 W, 17.5S
    s = slice(None,None,None) # this is the same as : 
    t = s,74,118 # [t] = [:,58,159]
    npp= var_npp[t]*86400   #   kg/m2/s kg/m2/day; 
    rh= var_rh[t]*86400;   # per s to per day 
    ra= var_ra[t]*86400;
    (
        csoil,
        cveg,
    ) = map(
            lambda var: var[t],
        (
            var_csoil, 
            var_cveg,
        )
    ) 
    return (npp, rh, ra, csoil, cveg)

def make_param_filter_func(
        c_max: np.ndarray,
        c_min: np.ndarray,
        c_veg,
        c_soil 
        ) -> Callable[[np.ndarray], bool]:

    def isQualified(c):
        # fixme
        #   this function is model specific: It discards parameter proposals
        #   where beta1 and beta2 are > 0.99
        #   checks for soil/veg sums as well
        #print("Is Qualified was called")
        paramNum = len(c)
        flag = True
        for i in range(paramNum):
            if(c[i] > c_max[i] or c[i] < c_min[i]):
                flag = False
                break
        if(c[0] + c[1] > 0.99):
            flag = False
        if(np.sum(c[17:18]) > 0.99*c_veg):
            flag = False
        if(np.sum(c[19:26]) > 0.99*c_soil):
            flag = False
        if((c[14] > 0.99) or (c[15] > 0.99) or (c[16] > 0.99)):
            flag = False
        return flag
    
    return isQualified

# def make_weighted_cost_func(
#         obs: np.ndarray
#     ) -> Callable[[np.ndarray],np.float64]:
#     # first unpack the observation array into its parts
#     cleaf, croot, cwood, clitter, csoil, rh = np.split(obs,indices_or_sections=6,axis=1)
#     def costfunction(out_simu: np.ndarray) ->np.float64:
#         # fixme 
#         #   as indicated by the fact that the function lives in this  
#         #   model-specific module it is not apropriate for (all) other models.
#         #   There are model specific properties:
#         #   1.) The weight for the different observation streams
#         #   
#         tot_len=out_simu.shape[0]
#         # we assume the model output to be in the same shape and order 
#         # as the obeservation
#         # this convention has to be honored by the forwar_simulation as well
#         # which in this instance already compresses the 3 different litter pools
#         # to c_litter and the 3 different soil pools to one
#         c_simu = out_simu[:,0:5] 
#         
#         # we assume the rh  part to be in the remaining columns again
#         # this convention has to be honored by the forwar_simulation as well
#         rh_simu = out_simu[:,5:]
#         #from IPython import embed; embed()
# 
#         J_obj1 = np.mean (( c_simu[:,0] - cleaf[0:tot_len] )**2)/(2*np.var(cleaf[0:tot_len]))
#         J_obj2 = np.mean (( c_simu[:,1] - croot[0:tot_len] )**2)/(2*np.var(croot[0:tot_len]))
#         J_obj3 = np.mean (( c_simu[:,2] - cwood[0:tot_len] )**2)/(2*np.var(cwood[0:tot_len]))
#         J_obj4 = np.mean (( c_simu[:,3]-  clitter[0:tot_len] )**2)/(2*np.var(clitter[0:tot_len]))
#         J_obj5 = np.mean (( c_simu[:,4]-  csoil[0:tot_len] )**2)/(2*np.var(csoil[0:tot_len]))
#         
#         J_obj6 = np.mean (( rh_simu[:,0] - rh[0:tot_len] )**2)/(2*np.var(rh[0:tot_len]))
#         
#         J_new = (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 )/200+ J_obj6/4
#         return J_new
#     return costfunction     


def make_param2res(
        cpa
    ) -> Callable[[np.ndarray], np.ndarray]: 
    """The returned function 'param2res' is used by the metropolis hastings algorithm
    to produce a result from a proposed parameter(tuple).
    'param2res' will be called only with the parameters to be estimated 'pe'.
    The constant parameters are the parameters of the present function 'pc' and are automatically available in the returned 'param2res' (closure).
    
    In the case of this model 'param2res' has to perform the following 
    tasks.
    -   use both sets of parameters 'pc' and 'pe' to build a forward model
    -   run the forward model to produce output for the times 
        when obeservations are available.(Here monthly)
    -   project the output of the forward model to the obeserved variables.
        (here summation of all different soil-pools to compare to the single
        observed soil-pool and the same for the litter pools)
        
    In this version of the function all tasks are performed at once as in 
    the original script.
    See the alternative implementation to see how all the tasks can be 
    delegated to smaller functions.
    """
    # fixme
    # This function is model-specific in several ways:
    # 0. Which are the fixed and which are the estimated variables
    # 1. The matrices for the forward simulation 
    # 2. the driver (here monthly npp)
    # 3. the projection of the simulated pool values to observable variables
    #    (here summation of   
    def param2res(pa):
        # pa is a numpy array when pa comes from the predictor
        # so we transform it to be able to use names instead of positions 
        epa = EstimatedParameters(*pa)
        days = 30 #for trendy all months have 30days
        tot_len = cpa.nyears * pseudo_days_per_month * 12 #total length is now nyears b/c monthly -> annual

        #define matrix dimensions
        abh = 12  #height of A and B matrix
        aw = 12   #width of A matrix
        bw = 1    #width of B matrix
        
        #assign betas (allocation) from input parameters
        beta_leaf=epa.beta_leaf
        beta_root=epa.beta_root
        beta_wood=1-beta_leaf-beta_root
        
        #create B matrix
        b = np.zeros(abh).reshape([abh,bw])  #empty matrix of zeros
        b[0] = beta_leaf   #fill matrix with beta values
        b[1] = beta_root   #remember python n starts at zero
        b[2] = beta_wood

        # define transfer values for A matrix    
        f_samet_leaf = epa.f_samet_leaf * epa.k_leaf     #optimized: leaf -> surface metabolic
        f_sastr_leaf = (1-f_samet_leaf) * epa.k_leaf     #remainder to surface str 
        f_slmet_root = epa.f_slmet_root * epa.k_root     #optimized: root -> soil metabolic
        f_slstr_root = (1-f_slmet_root) * epa.k_root     #remainder to soil str 
        f_cwd_wood = 0.4 * epa.k_wood                    #fixed - made up values so far
        f_samic_cwd = epa.f_samic_cwd * epa.k_cwd        #optimized: cwd ->  surface microbial
        f_slow_cwd = (1-f_samic_cwd) * epa.k_cwd         #remainder to slow soil 
        f_samic_samet = 0.3 * epa.k_samet                #fixed - made up values so far
        f_samic_sastr = 0.1 * epa.k_sastr                #fixed - made up values so far
        f_slow_sastr = 0.1 * epa.k_sastr                 #fixed - made up values so far
        f_slow_samic = 0.1 * epa.k_samic                 #fixed - made up values so far
        f_slmic_slmet = 0.4 * epa.k_slmet                #fixed - made up values so far
        f_slmic_slstr = 0.3 * epa.k_slstr                #fixed - made up values so far
        f_slow_slstr = 0.2 * epa.k_slstr                 #fixed - made up values so far
        f_slow_slmic = 0.4 * epa.k_slmic                 #fixed - made up values so far
        f_arm_slmic = 0.4 * epa.k_slmic                  #fixed - made up values so far
        f_slmic_slow = 0.10 * epa.k_slow                 #fixed - made up values so far
        f_arm_slow = 0.45*(0.003+0.009*cpa.clay) * epa.k_slow	#copied CABLES slow to passive
        f_slmic_arm = 0.10 * epa.k_arm                   #fixed - made up values so far
        
        #create A matrix
        A = np.zeros(abh*aw).reshape([abh,aw])  #create empty A matrix
        np.fill_diagonal(A,-1)                  #replace diagonal with -1s
        A[4,0] = f_samet_leaf   	#f5_1 #add in transfers - A is a zero ordered matrix
        A[5,0] = f_sastr_leaf   	#f6_1 #minus one from matrix position for [row,col]
        A[7,1] = f_slmet_root   	#f8_2 #f8_2, to pool8 from pool2 matrix pos [7,1]
        A[8,1] = f_slstr_root   	#f9_2
        A[3,2] = f_cwd_wood             #f4_3
        A[6,3] = f_samic_cwd            #f7_4
        A[10,3] = f_slow_cwd    	#f11_4
        A[6,4] = f_samic_samet          #f7_5
        A[6,5] = f_samic_sastr          #f7_6
        A[10,5] = f_slow_sastr          #f11_6
        A[10,6] = f_slow_samic          #f11_7
        A[9,7] = f_slmic_slmet          #f10_8
        A[9,8] = f_slmic_slstr          #f10_9
        A[10,8] = f_slow_slstr          #f11_9
        A[10,9] = f_slow_slmic          #f11_10
        A[11,9] = f_arm_slmic   	#f12_10
        A[9,10] = f_slmic_slow          #f10_11
        A[11,10] = f_arm_slow  	        #f12_11
        A[9,11] =f_slmic_arm    	#f10_12
        
        #turnover rate per day of pools:
        #leaf(p1),root(p2),wood(p3),cwd(p4),samet(p5),sastr(p6),samic(p7),
        #slmet(p8),slstr(p9),slmic(p10),slow(p11),arm(p12) 
        k_val = np.array(
            [
                epa.k_leaf,
                epa.k_root,
                epa.k_wood, 
                epa.k_cwd,
                epa.k_samet, 
                epa.k_sastr, 
                epa.k_samic,
                epa.k_slmet, 
                epa.k_slstr, 
                epa.k_slmic, 
                epa.k_slow, 
                epa.k_arm
            ]
        )
        K = np.zeros(abh*aw).reshape([12,12])
        np.fill_diagonal(K, k_val)
          
        x_fin=np.zeros((cpa.nyears,12))
        rh_fin=np.zeros((cpa.nyears,1))
        #ra_fin=np.zeros((cpa.nyears,1))
        #leaf(p1),root(p2),wood(p3),cwd(p4),samet(p5),sastr(p6),samic(p7),
        #slmet(p8),slstr(p9),slmic(p10),slow(p11),arm(p12)
        #x_init = np.array([cleaf[0],croot[0],cwood[0],cwd[0],pa[17],pa[18],
        #clitter[0]-pa[17]-pa[18],pa[19],pa[20],pa[21],pa[22],
        #csoil[0]-pa[19]-pa[20]-pa[21]-pa[22]]).reshape([12,1])
        x_init = np.array(
            [
                epa.C_leaf_0,
                epa.C_root_0,
                cpa.C_veg_0 - epa.C_leaf_0 - epa.C_root_0,
                epa.C_cwd_0,
                epa.C_samet_0,
                epa.C_sastr_0,
                epa.C_samic_0,
                epa.C_slmet_0,
                epa.C_slstr_0,
                epa.C_slmic_0,
                epa.C_slow_0,
                cpa.C_soil_0 - epa.C_slmet_0 - epa.C_slstr_0 - epa.C_slmic_0 - epa.C_slow_0 - epa.C_cwd_0 - epa.C_samet_0 - epa.C_sastr_0 - epa.C_samic_0
            ]
        ).reshape([12,1])
        X=x_init   # initialize carbon pools 
        i_m = 0
        i_pd =0
        for y in np.arange(0,cpa.nyears):
            x_year_avg=0
            ra_year_avg=0
            rh_year_avg=0
            for m in np.arange(0,12):
                npp_in = cpa.npp[i_m]
                co2_rh = 0
                co2_ra = 0
                for pd in np.arange(0,pseudo_days_per_month):
                    co2_hrate = [
                        0,
                        0,
                        0, 
                        (K[3,3]-f_samic_cwd-f_slow_cwd), 
                        (K[4,4]-f_samic_samet), 
                        (K[5,5]-f_samic_sastr-f_slow_sastr), 
                        (K[6,6]-f_slow_samic), 
                        (K[7,7]-f_slmic_slmet), 
                        (K[8,8]-f_slmic_slstr-f_slow_slstr), 
                        (K[9,9]-f_slow_slmic-f_arm_slmic), 
                        (K[10,10]-f_slmic_slow-f_arm_slow), 
                        (K[11,11]-f_slmic_arm)
                    ]
                    #similarly calculate autotrophic respiration ra
                    #co2_arate = [
                    #    (K[0,0]-f_samet_leaf-f_sastr_leaf), 
                    #    (K[1,1]-f_slmet_root-f_slstr_root),
                    #    (K[2,2]-f_cwd_wood),
                    #    0,
                    #    0,
                    #    0,
                    #    0,
                    #    0,
                    #    0,
                    #    0,
                    #    0,
                    #    0
                    #]
                    X=X + b*npp_in + np.array(A@X).reshape([12,1])
                    x_year_avg += X.reshape(1,12)/(pseudo_days_per_month*12)
                    co2h=np.sum(co2_hrate*X.reshape(1,12))
                    rh_year_avg += co2h/(pseudo_days_per_month*12)
                    
                    #co2a=np.sum(co2_arate*X.reshape(1,12))
                    #ra_year_avg += co2a/(pseudo_days_per_month*12)
                    i_pd += 1
                i_m += 1
            x_fin[y,:] = x_year_avg
            #ra_fin[y,:] = ra_year_avg
            rh_fin[y,:] = rh_year_avg

        # end od part I (run the nodel dayly
        # part II projection to yearly values 
        #x_fin = pseudo_daily_to_yearly(x) 
        #ra_fin = pseudo_daily_to_yearly(ra) 
        #rh_fin = pseudo_daily_to_yearly(rh) 
        # part IIb projection to combinded varaible 
        # We create an output that has the same shape
        # as the obvervations to make the costfunctions 
        # easier. 
        # above ground and belowground   
        # 'leaf',
        # 'root',
        # 'wood',
        # 'cwd',
        # 'samet',
        # 'sastr',
        # 'samic',
        # 'slmet',
        # 'slstr',
        # 'slmic',
        # 'slow',
        # 'arm'
        # x_veg = Leaf + wood + root + ?
        x_veg = np.sum(x_fin[:,0:6],axis=1).reshape(cpa.nyears,1)
        # x_soil = C_slmet + C_slstr + C_slmic + C_slow +C_arm
        x_soil = np.sum(x_fin[:,7:11],axis=1).reshape(cpa.nyears,1)
            
        out_simu = np.concatenate(
            [
                x_veg,
                x_soil,
                #ra_fin,
                rh_fin
            ]
            ,axis=1
        )
        return out_simu

    return param2res

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#
## alternative implementation and helper functions. If you do not want to use it # comment it.
#def make_param2res_2(
#        cpa: UnEstimatedParameters
#    ) -> Callable[[np.ndarray], np.ndarray]: 
#    # This function is an alternative implementation with the same results as 
#    # make_param2res 
#    # Internally it uses a slightly higher level of abstraction and devides the task 
#    # into more numerous but smaller testable and reusable units.
#    # Although both properties are not necessary to replicate the original
#    # code they provide a different perspective which is easier to reuse the code
#    # in  a wider range of models e.g. models with daily driver data, time
#    # dependent or even nonlinear matrices and so on. 
#    # Thus it allows to test small parts of the code against the respective parts of the 
#    # symbolic framework 
#    # 
#    # It differs in the following respects:
#    # 0.)   We use named tuples for the parameters (just for readability)
#    # 1.)   We use a seperate function to costruct the matrix which will be
#    #       a function B(i,X)  of the iteration it  and the present poolvalues X 
#    #       (althouhg in this case unnecessary since B is constant but
#    #       the dependence on 'i' becomes useful if you want to include 
#    #       an environmental scalar and the dependence on X makes it possible
#    #       to include nonlinearities.
#    # 2.)   We build a daily advancing model that can provide output for an arbitrary 
#    #       selection of days.  To do so we provide all driver data as
#    #       functions of the smalles timestep (here day with index i), which in
#    #       the case of this model means that we provide the same gpp value for
#    #       all the days in a given month. 
#    # 3.)   We translate the index of a given month to the appropriate day index
#    #       and apply the dayly model of 2.). Again this seems cumbersome for this
#    #       example but allows us to reuse the daily model for all applications.
#    #       This is espacially usefull for testing since we only need some dayly timesteps.
#    #       It makes it also easier to keep an overview over the appropriate 
#    #       units: If the smallest timestep is a day, then all time related parameters
#    #       have to be provided in the corresponding  units, regardless of the
#    #       number of available datapoints 
#    #
#    def param2res(pa):
#        epa=EstimatedParameters(*pa)
#        # here we want to store only monthly values
#        # although the computation takes place on a daily basis
#        V_init = construct_V0(cpa,epa)
#        # compute the days when we need the results
#        # to be able to compare to the monthly output
#        day_indices = month_2_day_index(range(cpa.number_of_months)) 
#        print("day_indices=",day_indices)
#        apa = Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)
#
#        mpa = ModelParameters(
#            **{
#                k:v for k,v in apa._asdict().items() 
#                if k in ModelParameters._fields
#            }
#        )
#        full_res = run_forward_simulation(
#            V_init=V_init,
#            day_indices=day_indices,
#            mpa=mpa
#        )
#        
#        # project the litter and soil pools
#        tot_len=cpa.number_of_months
#        c_litter = np.sum(full_res[:,3:6],axis=1).reshape(tot_len,1)
#        c_soil = np.sum(full_res[:,6:9],axis=1).reshape(tot_len,1)
#        out_simu = np.concatenate(
#            [
#                full_res[:,0:3], # the first 3 pools are used as they are
#                c_litter,
#                c_soil,
#                full_res[:,9:10]
#            ]
#            ,axis=1
#        )
#        return out_simu
#        
#    return param2res
#
#
#def run_forward_simulation(
#        V_init,
#        day_indices,
#        mpa
#    ):
#        tsi = make_daily_iterator(
#            V_init,
#            mpa=mpa
#        )
#
#        def g(acc, i):
#            xs,co2s,acc_co2,acc_days = acc
#            v = tsi.__next__()
#            d_pools = v[0:-1,:]
#            d_co2=v[-1:,:]
#            acc_co2 += d_co2
#            acc_days += 1
#            if i in day_indices:
#                xs += [d_pools] 
#                co2s +=[acc_co2/acc_days]
#                acc_co2=np.array([0.0]).reshape(1,1)
#                acc_days = 0
#                
#            acc = (xs,co2s,acc_co2,acc_days)
#            return acc
#        xs,co2s,acc_days,_ =  reduce(g,range(max(day_indices)+1),([],[],0,0))
#                
#        def h(tup):
#            x, co2 = tup
#            return np.transpose(np.concatenate([x,co2]))
#    
#        values_with_accumulated_co2 = [v for v in  map(h,zip(xs,co2s))]
#        RES = np.concatenate(values_with_accumulated_co2 , axis=0)  
#        return RES
#
#def make_daily_iterator(
#        V_init,
#        mpa
#    ):
#         
#        # Construct gpp(day)
#        # in general this function can depend on the day i and the state_vector X
#        # e.g. typically the size fo X.leaf...
#        # In this case it only depends on the day i 
#        def gpp_func(day,X):
#            return mpa.gpp[day_2_month_index(day)] 
#
#        # b (b vector for partial allocation) 
#        beta_wood = 1- mpa.beta_leaf- mpa.beta_root
#        b = np.array(
#            [
#                mpa.beta_leaf,
#                mpa.beta_root,
#                beta_wood,
#                0,
#                0,
#                0,
#                0,
#                0,
#                0
#            ],
#        ).reshape(9,1)
#        # Now construct B matrix B=A*K 
#        B_func  = make_compartmental_matrix_func(
#            mpa=mpa,
#        )
#        # Build the iterator which is the representation of a dynamical system
#        # for looping forward over the results of a difference equation
#        # X_{i+1}=f(X_{i},i)
#        # This is the discrete counterpart to the initial value problem 
#        # of the continuous ODE  
#        # dX/dt = f(X,t) and initial values X(t=0)=X_0
#        
#        def f(it,V):
#            X = V[0:9]
#            co2 = V[9]
#            gpp  = gpp_func(it,X)
#            B = B_func(it,X)
#            X_new = X + gpp * b + B@X
#
#            # we also compute the respired co2 in every (daily) timestep
#            # and use this part of the solution later to sum up the monthly amount
#            co2_new =  respiration_from_compartmental_matrix(B,X) 
#            
#            V_new = np.concatenate((X_new,np.array([co2_new]).reshape(1,1)), axis=0)
#            return V_new
#    
#    
#        #X_0 = np.array(x_init).reshape(9,1) #transform the tuple to a matrix
#        #co2_0 = np.array([0]).reshape(1,1)
#        #V_0 = np.concatenate((X_0, co2_0), axis=0)
#
#        return TimeStepIterator2(
#                initial_values=V_init,
#                f=f,
#                #max_it=max(day_indices)+1
#        )
#        
#
##def make_compartmental_matrix_func(
##        mpa
##    ):
##    # Now construct A matrix
##    # diagonal 
##    # make sure to use 1.0 instead of 1 otherwise it will be an interger array
##    # and round your values....
##    A =np.diag([-1.0 for i in range(9)]) 
##    # because of python indices starting at 0 we have A[i-1,j-1]=fij
##    #
##    #A = np.array([  -1,   0,   0,   0,   0,   0,   0,   0,   0,
##    #                 0,  -1,   0,   0,   0,   0,   0,   0,   0,
##    #                 0,   0,  -1,   0,   0,   0,   0,   0,   0,
##    #               f41, f42,   0,  -1,   0,   0,   0,   0,   0,
##    #               f51, f52,   0,   0,  -1,   0,   0,   0,   0,
##    #                 0,   0, f63,   0,   0,  -1,   0,   0,   0,
##    #                 0,   0,   0, f74, f75,   0,  -1,   0,   0,
##    #                 0,   0,   0,   0, f85, f86, f87,  -1,   0,
##    #                 0,   0,   0,   0,   0, f96, f97, f98,  -1 ]).reshape([9,9])   # tranfer
##
##    # note that the matrix description of a model is implicitly related to the
##    # order of pools in the state_vector and not independent from it.
##
##    A[3,0] = mpa.f_leaf2metlit
##    A[3,1] = mpa.f_root2metlit
##    A[4,0] = 1.0 - mpa.f_leaf2metlit
##    A[4,1] = 1.0 - mpa.f_root2metlit
##    A[5,2] = mpa.f_wood2CWD 
##    A[6,3] = mpa.f_metlit2mic 
##    A[6,4] = mpa.f_metlit2mic * (1 - mpa.lig_leaf)
##    A[7,4] = 0.7 * mpa.lig_leaf
##    A[7,5] = 0.4 * (1 - mpa.lig_wood)
##    A[8,5] = 0.7 * mpa.lig_wood
##    A[7,6] = (0.85 - 0.68 * (mpa.clay+mpa.silt)) * (0.997 - 0.032 * mpa.clay)
##    A[8,6] = (0.85 - 0.68 * (mpa.clay+mpa.silt)) * (0.003 + 0.032 * mpa.clay)
##    A[8,7] = 0.45 * (0.003 + 0.009 * mpa.clay)
##
##    #turnover rate per day of pools: 
##    K = np.diag([
##        mpa.k_leaf,
##        mpa.k_root,
##        mpa.k_wood,
##        mpa.k_metlit,
##        mpa.k_metlit/(5.75*np.exp(-3*mpa.lig_leaf)),
##        mpa.k_metlit/20.6,
##        mpa.k_mic,
##        mpa.k_slowsom,
##        mpa.k_passsom
##    ])
##    # in the general nonautonomous nonlinear case
##    # B will depend on an it,X (althouh in this case it does not depend on either
##    def B_func(it,X): 
##        return A@K
##
##    return B_func
##
#def run_forward_simulation_sym(
#        V_init,
#        day_indices,
#        mpa
#    ):
#        # Construct gpp(day)
#        # in general this function can depend on the day i and the state_vector X
#        # e.g. typically the size fo X.leaf...
#        # In this case it only depends on the day i 
#        def gpp_func(day,X):
#            return mpa.gpp[day_2_month_index(day)] 
#
#        func_dict = {Symbol('gpp'):gpp_func}
#        tsi = make_daily_iterator_sym(
#            V_init,
#            mpa=mpa,
#            func_dict=func_dict
#        )
#
#        def g(acc, i):
#            xs,co2s,acc_co2,acc_days = acc
#            v = tsi.__next__()
#            d_pools = v[0:-1,:]
#            d_co2=v[-1:,:]
#            acc_co2 += d_co2
#            acc_days += 1
#            if i in day_indices:
#                xs += [d_pools] 
#                co2s +=[acc_co2/acc_days]
#                acc_co2=np.array([0.0]).reshape(1,1)
#                acc_days = 0
#                
#            acc = (xs,co2s,acc_co2,acc_days)
#            return acc
#        xs,co2s,acc_days,_ =  reduce(g,range(max(day_indices)+1),([],[],0,0))
#                
#        def h(tup):
#            x, co2 = tup
#            return np.transpose(np.concatenate([x,co2]))
#    
#        values_with_accumulated_co2 = [v for v in  map(h,zip(xs,co2s))]
#        RES = np.concatenate(values_with_accumulated_co2 , axis=0)  
#        return RES

def make_daily_iterator_sym(
        V_init,
        mpa,
        func_dict
    ):
         

        # b (b vector for partial allocation) 
        beta_wood = 1- mpa.beta_leaf- mpa.beta_root
        b = np.array(
            [
                mpa.beta_leaf,
                mpa.beta_root,
                beta_wood,
                0,
                0,
                0,
                0,
                0,
                0
            ],
        ).reshape(9,1)
        # Now construct B matrix B=A*K 
        B_func  = construct_matrix_func_sym(
            pa=mpa,
        )
        # Build the iterator which is the representation of a dynamical system
        # for looping forward over the results of a difference equation
        # X_{i+1}=f(X_{i},i)
        # This is the discrete counterpart to the initial value problem 
        # of the continuous ODE  
        # dX/dt = f(X,t) and initial values X(t=0)=X_0
        
        def f(it,V):
            X = V[0:9]
            co2 = V[9]
            gpp  = func_dict[Symbol('gpp')](it,X)
            B = B_func(it,X)
            X_new = X + gpp * b + B@X

            # we also compute the respired co2 in every (daily) timestep
            # and use this part of the solution later to sum up the monthly amount
            co2_new = -np.sum(B @ X) # fixme add computer for respirattion
            
            V_new = np.concatenate((X_new,np.array([co2_new]).reshape(1,1)), axis=0)
            return V_new
    
    
        #X_0 = np.array(x_init).reshape(9,1) #transform the tuple to a matrix
        #co2_0 = np.array([0]).reshape(1,1)
        #V_0 = np.concatenate((X_0, co2_0), axis=0)

        return TimeStepIterator2(
                initial_values=V_init,
                f=f#,
                #max_it=max(day_indices)+1
        )

def construct_V0(
        cpa :UnEstimatedParameters,
        epa :EstimatedParameters
    ) -> np.ndarray:
    """Construct the initial values for the forward simulation
    from constant and eveluated parameters

    param: cpa : constant parameeters
    param: epa : estimated parameters 
    """
    # to make sure that we are in the right order we use the 
    # StateVariables namedtuple 
    X_0 = StateVariables( 
        leaf=  epa.C_leaf_0,
        root=  epa.C_root_0,
        wood=  cpa.C_veg_0 - epa.C_leaf_0 - epa.C_root_0 ,
        cwd=   epa.C_cwd_0,
        samet= epa.C_samet_0,
        sastr= epa.C_sastr_0,
        samic= epa.C_samic_0,
        slmet= epa.C_slmet_0,
        slstr= epa.C_slstr_0,
        slmic= epa.C_slmic_0,
        slow=  epa.C_slow_0,
        arm=   cpa.C_soil_0 - epa.C_slmet_0 - epa.C_slstr_0 - epa.C_slmic_0 - epa.C_slow_0
    )
    # add the respiration start value to the tuple
    return np.array(X_0).reshape(12,1)   

def make_param2res_sym(
        cpa: UnEstimatedParameters
    ) -> Callable[[np.ndarray], np.ndarray]: 
    # This function is an alternative implementation with the same results as 
    # make_param2res but uses the symbolic model in bgc_md2
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # here we want to store only monthly values
        # although the computation takes place on a daily basis
        V_init = construct_V0(cpa,epa)
        # compute the days when we need the results
        # to be able to compare to the monthly output
        day_indices = month_2_day_index(range(cpa.number_of_months)) 
        print("day_indices=",day_indices)
        apa = Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)

        mpa = ModelParameters(
            **{
                k:v for k,v in apa._asdict().items() 
                if k in ModelParameters._fields
            }
        )
        full_res = run_forward_simulation_sym(
            V_init=V_init,
            day_indices=day_indices,
            mpa=mpa
        )
        
        # project the litter and soil pools
        tot_len=cpa.number_of_months
        c_litter = np.sum(full_res[:,3:6],axis=1).reshape(tot_len,1)
        c_soil = np.sum(full_res[:,6:9],axis=1).reshape(tot_len,1)
        out_simu = np.concatenate(
            [
                full_res[:,0:3], # the first 3 pools are used as they are
                c_litter,
                c_soil,
                full_res[:,9:10]
            ]
            ,axis=1
        )
        return out_simu
        
    return param2res
