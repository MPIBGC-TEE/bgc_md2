import numpy as np
from sympy import Symbol, var
from functools import lru_cache, reduce
from typing import Callable

from CompartmentalSystems import helpers_reservoir as hr
from CompartmentalSystems.TimeStepIterator import (
        TimeStepIterator2,
)
from general_helpers import (
        day_2_month_index, 
        #months_by_day_arr,
        #respiration_from_compartmental_matrix,
        month_2_day_index,
        plot_solutions,
        make_B_u_funcs,
)

from ParameterMappings import (
    Observables,
    StateVariables,
    UnEstimatedParameters,
    EstimatedParameters,
    ModelParameters,
    Parameters
)
from model_specific_helpers import construct_V0

def run_forward_simulation_sym(
        V_init,
        day_indices,
        mpa
    ):

        func_dict = {Symbol('NPP'):make_npp_func(mpa)}
        tsi = make_daily_iterator_sym(
            V_init,
            mpa=mpa,
            func_dict=func_dict
        )

        def g(acc, i):
            xs,co2s,acc_co2,acc_days = acc
            v = tsi.__next__()
            d_pools = v[0:-1,:]
            d_co2=v[-1:,:]
            acc_co2 += d_co2
            acc_days += 1
            if i in day_indices:
                xs += [d_pools] 
                co2s +=[acc_co2/acc_days]
                acc_co2=np.array([0.0]).reshape(1,1)
                acc_days = 0
                
            acc = (xs,co2s,acc_co2,acc_days)
            return acc
        xs,co2s,acc_days,_ =  reduce(g,range(max(day_indices)+1),([],[],0,0))
                
        def h(tup):
            x, co2 = tup
            return np.transpose(np.concatenate([x,co2]))
    
        values_with_accumulated_co2 = [v for v in  map(h,zip(xs,co2s))]
        RES = np.concatenate(values_with_accumulated_co2 , axis=0)  
        return RES



def make_daily_iterator_sym(
        V_init,
        mpa,
        func_dict
    ):
        from bgc_md2.models.cable_yuanyuan.source import mvs 
        B_func, u_func = make_B_u_funcs(mvs,mpa,func_dict)  
        
        def f(it,V):
            X = V[0:9]
            co2 = V[9]
            b = u_func(it,X)
            B = B_func(it,X)
            X_new = X + b + B@X

            # we also compute the respired co2 in every (daily) timestep
            # and use this part of the solution later to sum up the monthly amount
            co2_new = -np.sum(B @ X) # fixme add computer for respirattion
            
            V_new = np.concatenate((X_new,np.array([co2_new]).reshape(1,1)), axis=0)
            return V_new
    
        return TimeStepIterator2(
                initial_values=V_init,
                f=f#,
                #max_it=max(day_indices)+1
        )
        
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



def make_npp_func(mpa):
        # Construct npp(day)
        # in general this function can depend on the day i and the state_vector X
        # e.g. typically the size fo X.leaf...
        # In this case it only depends on the day i 
    def npp_func(day):
        return mpa.npp[day_2_month_index(day)] 

    return npp_func


#def construct_matrix_func_sym(pa):
#    # we create a parameterdict for the fixed values
#    # and extend it by the parameters provided 
#    from bgc_md2.models.cable_yuanyuan.source import mvs 
#    symbol_names = mvs.get_BibInfo().sym_dict.keys()   
#    for name in symbol_names:
#        var(name)
#    parDict = {
#        clay: 0.2028,
#        silt: 0.2808,
#        lig_wood: 0.4,
#        f_wood2CWD: 1,
#        f_metlit2mic: 0.45,
#    #    NPP: npp_in
#    }
#    model_params = {Symbol(k): v for k,v in pa._asdict().items()}
#    parDict.update(model_params)
#    B_func = hr.numerical_array_func(
#            state_vector = mvs.get_StateVariableTuple(),
#            time_symbol=mvs.get_TimeSymbol(),
#            expr=mvs.get_CompartmentalMatrix(),
#            parameter_dict=parDict,
#            func_dict={}
#    )
#    # in the general nonautonomous nonlinear B_func is a function of t,x
#    # although for this example it does not depend on either t, nor x.
#    return B_func 
