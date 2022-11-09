# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Exporting yearly traceable components to a .csv

# ### Loading required packages  and functions

# %load_ext autoreload
# %autoreload 2
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# from functools import lru_cache
import general_helpers as gh
import MIP_output_helpers as moh
# from bgc_md2.resolve.mvars import (
#     CompartmentalMatrix,
#     InputTuple,
#     StateVariableTuple
# )

# +
delta_t_val=30

model_names={
    "yz_jules": "JULES",
    "kv_visit2": "VISIT",
}
model_folders=[(k) for k in model_names]
test_arg_list=gh.get_test_arg_list(model_folders)

model_cols={
    "yz_jules": "blue",
    "kv_visit2": "orange",
    "jon_yib": "green",
    "kv_ft_dlem": "red",
    "Aneesh_SDGVM":"yellow",
    "cj_isam": "purple",
    "bian_ibis2":"magenta",
    "ORCHIDEE-V2":"teal",
}

var_names={
    "x": "X",
    "x_c":"X_c",
    "x_p": "X_p",
    "u": "C input",
    "rt": "Residence time",
}


# -

# this is to export yearly traceable components for use outside of framework (experimenting with Enqing's code) 
def write_yearly_components(
    model_names,  # dictionary (folder name : model name)
    test_arg_list,  # a list of test_args from all models involved
    var_names,  # dictionary (trace_tuple name : descriptive name)
    delta_t_val,  # model time step
    model_cols,  # dictionary (folder name :color)
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging=12,  # number of iterator steps over which to average results. 1 for no averaging
):
    if (part < 0) | (part > 1):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        )
    model_folders = [(k) for k in model_names]
    variables = [(k) for k in var_names]
    n = len(variables)

    Xs=np.array((0,))
    X_cs=np.array((0,))
    X_ps=np.array((0,))
    us=np.array((0,))
    RTs=np.array((0,))
 
    for i, name in enumerate(variables):
        k = 0
        for mf in model_folders:
            itr = moh.traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
            start_min, stop_max = gh.min_max_index(
                test_arg_list[k],
                delta_t_val,
                *gh.t_min_tmax_overlap(test_arg_list, delta_t_val)
            )
            # if we do not want the whole interval but look at a smaller part to observe the dynamics
            # start,stop = start_min, int(start_min+(stop_max-start_min)*part)
            start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
            times = (
                gh.times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
                / gh.days_per_year()
            )
            print("Gathering " + str(name) + " data for " + str(mf) + " model")
            vals = itr[start:stop]
            vals_output = gh.avg_timeline(vals.__getattribute__(name), averaging)
            
            if name == "x":  
                Xs=np.concatenate((Xs,vals_output))
            if name == "x_c":  
                X_cs=np.concatenate((X_cs,vals_output))  
            if name == "x_p":  
                X_ps=np.concatenate((X_ps,vals_output)) 
            if name == "u":  
                us=np.concatenate((us,vals_output))                  
            if name == "rt":  
                RTs=np.concatenate((RTs,vals_output)) 
            k += 1       
        
    # Saving yearly components to a csv for further analysis
    test_keys = ["Model", "Year", "X", "X_c", "X_p", "u", "RT"]
        
    model = np.repeat(model_folders[0], 160).reshape(160,1) 
    models=np.array((model)).reshape(160,1)
    for i in range (len(model_folders)-1):
        model=np.repeat(model_folders[i+1], 160).reshape(160,1) 
        models=np.concatenate((models,model))
        
    print("Model")
    print(models.shape)
    #print(models)  
        
    year = np.array(range(1860,2020)).reshape(160,1)
    years=np.array((year)).reshape(160,1)
    for i in range (len(model_folders)-1):
        years=np.concatenate((years,year))
        
    print("Year")
    print(years.shape)
    
    test_values = [models.flatten(), years.flatten(), Xs.flatten()[1:], X_cs.flatten()[1:], X_ps.flatten()[1:], us.flatten()[1:], RTs.flatten()[1:]]
    
    template = {test_keys[i]: test_values[i] for i in range(len(test_keys))}
    output = pd.DataFrame(template)    
    
    pd.DataFrame(output).to_csv('traceable_components_output.csv', sep=',')
    
    return output  


output=write_yearly_components(model_names=model_names,
                        test_arg_list=test_arg_list,   
                        var_names=var_names,
                        delta_t_val=delta_t_val,
                        model_cols=model_cols,
                        part=1,
                        averaging=12*30//delta_t_val # yearly averaging
                       )


output


