from typing import Callable
from collections import namedtuple
import numpy as np
from functools import reduce
from general_helpers import day_2_month_index, month_2_day_index, months_by_day_arr, TimeStepIterator2, respiration_from_compartmental_matrix

# fixme:
# Your parameters will most likely differ but you can still use the
# destinctions between different sets of parameters. The aim is to make
# the different tasks in the code more obvious. In principal you could
# have a lot of overlapping sets and just have to keep them consistent. 
# 'namedtuples' can be used just as normal tuples by functions
# that are not aware of the names. They can still use the positions like 
# in the original code

# This set is used by the functions that produce the 
# specific ingredients (functions) that will be run by
# mcmc alg.
UnEstimatedParameters = namedtuple(
    "UnEstimatedParameters",
    [
        'C_leaf_0',
        'C_root_0',
        'C_wood_0',
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

# This set is used (as argument) by the functions that are called
# inside the mcmc
EstimatedParameters = namedtuple(
    "EstiamatedParameters",
    [
        "beta_leaf",    #  0 (indices uses in original code) 
        "beta_root",    #  1
        "lig_leaf",     #  2
        "f_leaf2metlit",#  3
        "f_root2metlit",#  4
        "k_leaf",       #  5
        "k_root",       #  6
        "k_wood",       #  7
        "k_metlit",	#  8
        "k_mic",	#  9
        "k_slowsom",	# 10
        "k_passsom",	# 11
        "C_metlit_0",	# 12
        "CWD_0",	# 13
        "C_mic_0",	# 14
        "C_passom_0"    # 15
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
        'C_leaf',
        'C_root',
        'C_wood',
        'C_metlit',
        'C_stlit',
        'CWD',
        'C_mic',
        'C_slowsom',
        'C_passsom'
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
def project(pa,target_cls):
        return target_cls(
                **{ k: v for k,v in pa._asdict().items() if k in target_cls._fields }
            )

def make_param2res(
        cpa: UnEstimatedParameters
    ) -> Callable[[np.ndarray], np.ndarray]: 
    #        cleaf_0,
    #        croot_0,
    #        cwood_0,
    #        clitter_0,
    #        csoil_0,
    #        npp,
    #        tot_len,
    #        clay, 
    #        silt,
    #        lig_wood,
    #        f_wood2CWD,
    #        f_metlit2mic
    # this function is an alternative implementation with the same results as 
    # make_param2res in yy_cable_specific_helpers.py
    # Internally it uses a slightly higher level of abstraction and devides the task 
    # into more numerous but smaller testable and reusable units.
    # Although both properties are not necessary to replicate the original
    # code they provide a different perspective which is easier to reuse the code
    # in  a wider range of models e.g. models with daily driver data, time
    # dependent or even nonlinear matrices and so on. 
    # Thus it allows to test small parts of the code against the respective parts of the 
    # symbolic framework 
    # 
    # It differs in the following respects:
    # 0.)   We use named tuples for the parameters (just for readability)
    # 1.)   We use a seperate function to costruct the matrix which will be
    #       a function B(i,X)  of the iteration it  and the present poolvalues X 
    #       (althouhg in this case unnecessary since B is constant but
    #       the dependence on 'i' becomes useful if you want to include 
    #       an environmental scalar and the dependence on X makes it possible
    #       to include nonlinearities.
    # 2.)   We build a daily advancing model that can provide output for an arbitrary 
    #       selection of days.  To do so we provide all driver data as
    #       functions of the smalles timestep (here day with index i), which in
    #       the case of this model means that we provide the same npp value for
    #       all the days in a given month. 
    # 3.)   We translate the index of a given month to the appropriate day index
    #       and apply the dayly model of 2.). Again this seems cumbersome for this
    #       example but allows us to reuse the daily model for all applications.
    #       This is espacially usefull for testing since we only need some dayly timesteps.
    #       It makes it also easier to keep an overview over the appropriate 
    #       units: If the smallest timestep is a day, then all time related parameters
    #       have to be provided in the corresponding  units, regardless of the
    #       number of available datapoints 
    #
    def param2res(pa):
        epa=EstimatedParameters(*pa)
        # here we want to store only monthly values
        # although the computation takes place on a daily basis
        #x_init = np.array([cleaf_0,croot_0,cwood_0,epa[12],epa[13],clitter_0-epa[12]-epa[13],epa[14],csoil_0- epa[14] - epa[15], epa[15]]).reshape([9,1])   
        # Initial carbon pool size
        x_init = construct_X0(cpa,epa)
        # compute the days when we need the results
        # to be able to compare to the monthly output
        day_indices = month_2_day_index(range(cpa.number_of_months)) 
        apa = Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)

        mpa = ModelParameters(
            **{
                k:v for k,v in apa._asdict().items() 
                if k in ModelParameters._fields
            }
        )
        full_res = run_forward_simulation(
            x_init=x_init,
            day_indices=day_indices,
            **mpa._asdict()
            #npp=npp,
            #clay=clay, 
            #silt=silt,
            #lig_wood=lig_wood,
            #lig_leaf=epa.lig_leaf,
            #beta_leaf=epa.beta_leaf,
            #beta_root=epa.beta_root,              
            #f_leaf2metlit=epa.f_leaf2metlit,
            #f_root2metlit=epa.f_root2metlit,
            #f_wood2CWD=f_wood2CWD,
            #f_metlit2mic=f_metlit2mic,
            #k_leaf=epa.k_leaf,
            #k_root=epa.k_root,
            #k_wood=epa.k_wood,
            #k_metlit=epa.k_metlit,
            #k_mic=epa.k_mic,
            #k_slowsom=epa.k_slowsom,
            #k_passsom=epa.k_passsom
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


def run_forward_simulation(
        x_init,
        day_indices,
        npp,
        clay, 
        silt,
        lig_wood,
        lig_leaf,
        beta_leaf,
        beta_root,
        f_leaf2metlit,
        f_root2metlit,
        f_wood2CWD,
        f_metlit2mic,
        k_leaf,
        k_root,
        k_wood,
        k_metlit,
        k_mic,
        k_slowsom,
        k_passsom      
    ):
        tsi = make_daily_iterator(
            x_init,
            npp,
            clay, 
            silt,
            lig_wood,
            lig_leaf,
            beta_leaf,
            beta_root,
            f_leaf2metlit,
            f_root2metlit,
            f_wood2CWD,
            f_metlit2mic,
            k_leaf,
            k_root,
            k_wood,
            k_metlit,
            k_mic,
            k_slowsom,
            k_passsom      
        )
        

        def g(acc, i):
            xs,co2s,acc_co2 = acc
            v = tsi.__next__()
            d_pools = v[0:-1,:]
            d_co2=v[-1:,:]
            acc_co2 += d_co2
            if i in day_indices:
                xs += [d_pools] 
                co2s +=[acc_co2]
                acc_co2=np.array([0.0]).reshape(1,1)
                
            acc = (xs,co2s,acc_co2)
            return acc
        xs,co2s,_ =  reduce(g,range(max(day_indices)+1),([],[],0))
                
        def h(tup):
            x, co2 = tup
            return np.transpose(np.concatenate([x,co2]))
    
        values_with_monthly_co2 = [v for v in  map(h,zip(xs,co2s))]
        RES = np.concatenate(values_with_monthly_co2 , axis=0)  
        return RES

def make_daily_iterator(
        x_init,
        npp,
        clay, 
        silt,
        lig_wood,
        lig_leaf,
        beta_leaf,
        beta_root,
        f_leaf2metlit,
        f_root2metlit,
        f_wood2CWD,
        f_metlit2mic,
        k_leaf,
        k_root,
        k_wood,
        k_metlit,
        k_mic,
        k_slowsom,
        k_passsom      
    ):
         
        # Construct npp(day)
        # in general this function can depend on the day i and the state_vector X
        # e.g. typically the size fo X.leaf...
        # In this case it only depends on the day i 
        def npp_func(day,X):
            return npp[day_2_month_index(day)] * b

        # b (b vector for partial allocation) 
        beta_wood = 1- beta_leaf- beta_root
        b = np.array(
            [
                beta_leaf,
                beta_root,
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
        B_func  = make_compartmental_matrix_func(
            clay, 
            silt,
            lig_wood,
            lig_leaf,
            f_leaf2metlit,
            f_root2metlit,
            f_wood2CWD,
            f_metlit2mic,
            k_leaf,
            k_root,
            k_wood,
            k_metlit,
            k_mic,
            k_slowsom,
            k_passsom
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
            npp  = npp_func(it,X)
            B = B_func(it,X)
            X_new = X + npp * b + B@X

            # we also compute the respired co2 in every (daily) timestep
            # and use this part of the solution later to sum up the monthly amount
            co2_new =  respiration_from_compartmental_matrix(B,X) 
            
            V_new = np.concatenate((X_new,np.array([co2_new]).reshape(1,1)), axis=0)
            return V_new
    
    
        X_0 = np.array(x_init).reshape(9,1) #transform the tuple to a matrix
        co2_0 = np.array([0]).reshape(1,1)
        V_0 = np.concatenate((X_0, co2_0), axis=0)

        return TimeStepIterator2(
                initial_values=V_0,
                f=f,
                #max_it=max(day_indices)+1
        )
        

def make_compartmental_matrix_func(
        clay, 
        silt,
        lig_wood,
        lig_leaf,
        f_leaf2metlit,
        f_root2metlit,
        f_wood2CWD, 
        f_metlit2mic,
        k_leaf,
        k_root,
        k_wood,
        k_metlit,
        k_mic,
        k_slowsom,
        k_passsom      
    ):
    # Now construct A matrix
    # diagonal 
    # make sure to use 1.0 instead of 1 otherwise it will be an interger array
    # and round your values....
    A =np.diag([-1.0 for i in range(9)]) 
    #A[i-1,j-1]=fij
    #A = np.array([  -1,   0,   0,   0,   0,   0,   0,   0,   0,
    #                 0,  -1,   0,   0,   0,   0,   0,   0,   0,
    #                 0,   0,  -1,   0,   0,   0,   0,   0,   0,
    #               f41, f42,   0,  -1,   0,   0,   0,   0,   0,
    #               f51, f52,   0,   0,  -1,   0,   0,   0,   0,
    #                 0,   0, f63,   0,   0,  -1,   0,   0,   0,
    #                 0,   0,   0, f74, f75,   0,  -1,   0,   0,
    #                 0,   0,   0,   0, f85, f86, f87,  -1,   0,
    #                 0,   0,   0,   0,   0, f96, f97, f98,  -1 ]).reshape([9,9])   # tranfer

    # note that the matrix description of a model is implicitly related to the
    # order of pools in the state_vector and not independent from it.

    A[3,0] = f_leaf2metlit
    A[3,1] = f_root2metlit
    A[4,0] = 1.0 - f_leaf2metlit
    A[4,1] = 1.0 - f_root2metlit
    A[5,2] = f_wood2CWD 
    A[6,3] = f_metlit2mic 
    A[6,4] = f_metlit2mic * (1 - lig_leaf)
    A[7,4] = 0.7 * lig_leaf
    A[7,5] = 0.4 * (1 - lig_wood)
    A[8,5] = 0.7 * lig_wood
    A[7,6] = (0.85 - 0.68 * (clay+silt)) * (0.997 - 0.032 * clay)
    A[8,6] = (0.85 - 0.68 * (clay+silt)) * (0.003 + 0.032 * clay)
    A[8,7] = 0.45 * (0.003 + 0.009 * clay)

    #turnover rate per day of pools: 
    K = np.diag([
        k_leaf,
        k_root,
        k_wood,
        k_metlit,
        k_metlit/(5.75*np.exp(-3*lig_leaf)),
        k_metlit/20.6,
        k_mic,
        k_slowsom,
        k_passsom
    ])
    # in the general nonautonomous nonlinear case
    # B will depend on an it,X (althouh in this case it does not depend on either
    def B_func(it,X): 
        return A@K

    return B_func

def construct_X0(
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
        C_leaf=cpa.C_leaf_0,
        C_root=cpa.C_root_0,
        C_wood=cpa.C_wood_0,
        C_metlit=epa.C_metlit_0,
        C_stlit=epa.CWD_0,
        CWD=cpa.clitter_0-epa.C_metlit_0-epa.CWD_0,
        C_mic=epa.C_mic_0,
        C_slowsom=cpa.csoil_0- epa.C_mic_0 - epa.C_passom_0, 
        C_passsom=epa.C_passom_0
    )
    return np.array(X_0).reshape([9,1])   
