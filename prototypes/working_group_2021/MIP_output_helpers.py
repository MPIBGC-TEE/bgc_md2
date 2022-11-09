import numpy as np
from copy import copy
import matplotlib.pyplot as plt
from pathlib import Path
import json
import netCDF4 as nc    
import general_helpers as gh

def average_and_resample_nc(
            model_names, # dictionary e.g. "ab_classic":"CLASSIC"
            experiment_names, # e.g. ['S2', 'S3']
            target_mask,
            method="nearest",
            radius_of_influence=500000,
            ):
    for experiment in experiment_names:
        print('\033[1m'+'Resampling data for '+experiment+' experiment...')
        k=0 # model counter
        model_folders=[(m) for m in model_names] 
        m_names=list(model_names.values())  
        for mf in model_folders:
            print('\033[1m'+m_names[k])
            print('\033[0m')
            experiment_name=m_names[k]+"_"+experiment+"_"
            conf_dict = gh.confDict(mf)
            dataPath=Path(conf_dict["dataPath"])
            model_mask=gh.msh(mf).spatial_mask(dataPath=Path(conf_dict["dataPath"]))       
            for vn in gh.msh(mf).data_str._fields:      
                if vn=="npp_nlim": file_path = dataPath.joinpath(gh.msh(mf).nc_file_name("npp", experiment_name=experiment_name))
                else: file_path = dataPath.joinpath(gh.msh(mf).nc_file_name(vn, experiment_name=experiment_name))
                print(file_path)
                ds = nc.Dataset(str(file_path))
                var=ds.variables[vn][:, :, :].data
                
                # finding start value (differs for different models)
                if var.shape[0]>1000: # monthly data
                    start=var.shape[0]-1440 # smallest time range - SDGVM fluxes 
                else: # yearly data
                    start=var.shape[0]-120 # smallest time range - SDGVM fluxes
                # temporal average + difference
                if vn in ['gpp', 'npp', 'ra', 'rh', 'npp_nlim']: # fluxes per sec to total per 120 years
                    var_avg=np.ma.mean(var[start:,:,:], axis=0)*86400*365*120                  
                else:  # for pools we compute mean              
                    var_avg=np.ma.mean(var[start:,:,:], axis=0)    
                #var_avg=np.ma.mean(var[start:,:,:], axis=0)                
                var_diff=var[-1,:,:]-var[start,:,:]
                # masking
                var_masked_avg = np.ma.array(var_avg, mask = model_mask.index_mask)
                var_masked_diff = np.ma.array(var_diff, mask = model_mask.index_mask) 
                # resampling to target grid
                var_resampled_avg=gh.resample_grid (
                    source_coord_mask=model_mask, 
                    target_coord_mask=target_mask, 
                    var=var_masked_avg, 
                    method=method,
                    radius_of_influence=radius_of_influence,
                )
                var_resampled_diff=gh.resample_grid (
                    source_coord_mask=model_mask, 
                    target_coord_mask=target_mask, 
                    var=var_masked_diff, 
                    method=method,
                    radius_of_influence=radius_of_influence,
                )
                # final masking
                var_final_avg = np.ma.array(var_resampled_avg.index_mask, 
                                            mask = target_mask.index_mask)
                var_final_diff = np.ma.array(var_resampled_diff.index_mask, 
                                            mask = target_mask.index_mask)            
                # creating and writing a new NetCDF file               
                s = target_mask.index_mask.shape
                n_lats, n_lons = s
                
                # if vn in ['gpp', 'npp', 'ra', 'rh', 'npp_nlim']: # for fluxes               
                    # new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_sum_res.nc")
                # else: # for pools 
                    # new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_ave_res.nc")
                new_path=dataPath.joinpath(dataPath,experiment_name+vn+"_ave_res.nc")                    
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)
                # creating variables            
                avg = ds_new.createVariable(vn, "float32", ["lat", "lon"])     
                avg[:, :] = var_final_avg                
                diff = ds_new.createVariable(str(vn)+"_diff", "float32", ["lat", "lon"])                 
                diff[:, :] = var_final_diff
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(target_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(target_mask.tr.i2lon, range(n_lons)))        
                # closing NetCDF files
                ds.close()        
                ds_new.close()
            k+=1 # model counter
        print("Done!")

def traceability_iterator_instance(
    mf,  # name of the model folder
    test_args,  # a tuple defined in model-specific test_helpers.py
    delta_t_val,  # model time step
):
    """Wrapper that just needs the folder name and the step
    size.
    """
    mvs_t = gh.mvs(mf)
    dvs_t = test_args.dvs
    cpa_t = test_args.cpa
    epa_t = test_args.epa_opt
    X_0 = gh.msh(mf).numeric_X_0(mvs_t, dvs_t, cpa_t, epa_t)
    func_dict = gh.msh(mf).make_func_dict(mvs_t, dvs_t, cpa_t, epa_t)

    return gh.traceability_iterator(
        X_0,
        func_dict,
        mvs=mvs_t,
        dvs=dvs_t,
        cpa=cpa_t,
        epa=epa_t,
        delta_t_val=delta_t_val,
    )
    
def get_traceable_components(
    model_names,  # dictionary (folder name : model name)
    test_arg_list,  # a list of test_args from all models involved
    delta_t_val,  # model time step
    model_cols,  # dictionary (folder name :color)
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging,  # number of iterator steps over which to average results. 1 for no averaging
    overlap=True,  # compute overlapping timeframe or plot whole duration for all models
):
    if (part < -1) | (part > 1) | (part == 0):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between -1 and 1 excluding 0"
        )
    model_folders = [(k) for k in model_names]
    k = 0
    sum_x = np.array(0)
    sum_x_c = np.array(0)
    sum_x_p = np.array(0)
    sum_u = np.array(0)
    sum_rt = np.array(0)
    all_components = list ()
    for mf in model_folders:
        print ("Computing traceable components for "+mf+"...")
        itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        if overlap == True:
            start_min, stop_max = gh.min_max_index(
                test_arg_list[k],
                delta_t_val,
                *gh.t_min_tmax_overlap(test_arg_list, delta_t_val)
            )
        else:
            start_min, stop_max = gh.min_max_index(
                test_arg_list[k],
                delta_t_val,
                *gh.t_min_tmax_full(test_arg_list, delta_t_val)
            )
        # if we do not want the whole interval but look at a smaller part to observe the dynamics
        if part < 0:
            start, stop = int(stop_max - (stop_max - start_min) * abs(part)), stop_max
        else:
            start, stop = start_min, int(start_min + (stop_max - start_min) * part)
        times = (
            gh.times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            / gh.days_per_year()
        )
        vals = itr[start:stop]
        k += 1       
        
        times_avg = gh.avg_timeline(times, averaging)     
        
        comp_dict = {
            "x": gh.avg_timeline(vals.x, averaging),
            "x_c": gh.avg_timeline(vals.x_c, averaging),
            "x_p": gh.avg_timeline(vals.x_p, averaging),
            "u": gh.avg_timeline(vals.u, averaging),
            "rt": gh.avg_timeline(vals.rt, averaging),
            }
        all_components.append(comp_dict) 
        if sum_x.all()==0:
            sum_x = np.append(sum_x,vals.x)[1:]
            sum_x_c = np.append(sum_x_c,vals.x_c)[1:]
            sum_x_p = np.append(sum_x_p,vals.x_p)[1:]
            sum_u = np.append(sum_u,vals.u)[1:]
            sum_rt = np.append(sum_rt,vals.rt)[1:]
        else:
            sum_x = sum_x + vals.x
            sum_x_c = sum_x_c + vals.x_c
            sum_x_p = sum_x_p + vals.x_p
            sum_u = sum_u + vals.u
            sum_rt = sum_rt + vals.rt   
                
    
    ave_dict = {
            "x": gh.avg_timeline(sum_x / len(model_folders), averaging),
            "x_c": gh.avg_timeline(sum_x_c / len(model_folders), averaging),
            "x_p": gh.avg_timeline(sum_x_p / len(model_folders), averaging),
            "u": gh.avg_timeline(sum_u / len(model_folders), averaging),
            "rt": gh.avg_timeline(sum_rt / len(model_folders), averaging),
            }       
    
    all_components.append(ave_dict)
    all_components.append(times_avg)       
    
    mods=list(model_names.values())
    mods.append("Mean")
    mods.append("Times")
    all_comp_dict = {mods[i]: all_components[i] for i in range(len(mods))}      

    return all_comp_dict 
    
def plot_traceable_component(
    all_comp_dict,
    comp_name,
    model_cols,
    delta=False,
):
    models=list(all_comp_dict.keys())[:-2]
    fig = plt.figure(figsize=(20, 10))
    ax = fig.subplots(1, 1)
    vals_mean=0
    st_dev=0
    if comp_name == "x_x_c":
        ax.plot(
                all_comp_dict["Times"],
                all_comp_dict["Mean"]["x"]
                * 148940000
                * 1000000
                * 0.000000000001,  # convert to global C in Gt
                label="Ensemble mean - X",
                color="black",
                linewidth=3,
        )
        ax.plot(
                all_comp_dict["Times"],
                all_comp_dict["Mean"]["x_c"]
                * 148940000
                * 1000000
                * 0.000000000001,  # convert to global C in Gt
                label="Ensemble mean - X_c",
                color="black",
                linestyle="dashed",
                linewidth=2,
        )
        for m in models:      
            
            ax.plot(
                    all_comp_dict["Times"],
                    all_comp_dict[m]["x"]
                    * 148940000
                    * 1000000
                    * 0.000000000001,  # convert to global C in Gt
                    label=m + " - X",
                    color=model_cols[m],
                    linewidth=1,
                )
            ax.plot(
                    all_comp_dict["Times"],
                    all_comp_dict[m]["x_c"]
                    * 148940000
                    * 1000000
                    * 0.000000000001,  # convert to global C in Gt
                    label=m + " - X_c",
                    color=model_cols[m],
                    linestyle="dashed",
                    linewidth=1,
            )
    else:
        if delta:
            vals_array_mean = all_comp_dict["Mean"][comp_name]
            vals_mean = (vals_array_mean-vals_array_mean[0]) #/ vals_array_mean[0] * 100
        else:
            vals_mean = all_comp_dict["Mean"][comp_name]
                
        ax.plot(
                all_comp_dict["Times"],
                vals_mean,
                label="Ensemble mean",
                color="black",
                linewidth=3, 
                )                
        
        diff_sqr = np.zeros_like(all_comp_dict[models[0]][comp_name])             
        
        for m in models:
            if delta:
                vals_array = all_comp_dict[m][comp_name]
                vals=(vals_array-vals_array[0])#/vals_array[0] * 100
            else:
                vals = all_comp_dict[m][comp_name]  
                      
            diff_sqr=diff_sqr + (vals-vals_mean)**2
            
            ax.plot(
                    all_comp_dict["Times"],
                    vals,
                    label=m,
                    color=model_cols[m],
                    linewidth=1.2, 
                )
        variance = diff_sqr / (len(models)-1)
        st_dev=np.sqrt(variance)
      
        ax.fill_between(
                all_comp_dict["Times"],
                vals_mean-st_dev*2, 
                vals_mean+st_dev*2,
                label="\u00B12$\sigma$ interval",
                color="grey",
                alpha=0.2                
                )                          
    ax.grid()
    if delta:
        ax.set_title("$\Delta$ "+str(comp_name))
    else:
        ax.set_title(comp_name)
    #ax.set_ylabel("Gt C")
    ax.legend(bbox_to_anchor =(1.2, 1))   
    plt.show()
    return(vals_mean, st_dev)
         
def plot_attribution (
    times,
    delta_x_pos,
    delta_x_neg,
    rt_contrib_pos,
    rt_contrib_neg,
    u_contrib_pos,
    u_contrib_neg,
    rt_u_inter_pos,
    rt_u_inter_neg,
    x_p_contrib_pos,
    x_p_contrib_neg,
    percent=False,
):
        x_p_contrib=x_p_contrib_pos-x_p_contrib_neg
        rt_contrib=rt_contrib_pos-rt_contrib_neg
        u_contrib=u_contrib_pos-u_contrib_neg
        rt_u_inter=rt_u_inter_pos-rt_u_inter_neg
        
        fig4=plt.figure(figsize=(15,5))
        axs=fig4.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]})
        
        # contributions timeline
        ax=axs[0]

        # quadratic trends
            
        if abs(np.mean(x_p_contrib_pos))>0.01:
            z = np.polyfit(times,  x_p_contrib_pos, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ X_p positive",
                color="blue",
            ) 
            ax.fill_between(
                    times,
                    x_p_contrib_pos, 
                    p(times),
                    color="blue",
                    alpha=0.2                
                    )         
        if abs(np.mean(x_p_contrib_neg))>0.01:
            z = np.polyfit(times,  x_p_contrib_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ X_p negative",
                color="blue",
                linestyle="dashed",            
            ) 
            ax.fill_between(
                    times,
                    x_p_contrib_neg, 
                    p(times),
                    color="blue",
                    #linewidth=0.1,
                    alpha=0.2                
                    )          
        if abs(np.mean(u_contrib_pos))>0.01:        
            z = np.polyfit(times,  u_contrib_pos, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ NPP positive",
                color="green",
            ) 
            ax.fill_between(
                    times,
                    u_contrib_pos, 
                    p(times),
                    color="green",
                    alpha=0.2                
                    )
        if abs(np.mean(u_contrib_neg))>0.01:        
            z = np.polyfit(times,  u_contrib_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ NPP negative",
                color="green",
                linestyle="dashed",            
            ) 
            ax.fill_between(
                    times,
                    u_contrib_neg, 
                    p(times),
                    color="green",
                    alpha=0.2                
                    )           
        if abs(np.mean(rt_contrib_pos))>0.01:                
            z = np.polyfit(times,  rt_contrib_pos, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ RT positive",
                color="darkorange",
            )              
            ax.fill_between(
                    times,
                    rt_contrib_pos, 
                    p(times),
                    color="darkorange",
                    alpha=0.2                
                    )   
        if abs(np.mean(rt_contrib_neg))>0.01:        
            z = np.polyfit(times,  rt_contrib_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ RT negative",
                color="darkorange",
                linestyle="dashed",
            )              
            ax.fill_between(
                    times,
                    rt_contrib_neg, 
                    p(times),
                    color="darkorange",
                    alpha=0.2                
                    )      
        if abs(np.mean(rt_u_inter_pos))>0.01:          
            z = np.polyfit(times,  rt_u_inter_pos, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$NPP*$\Delta$RT positive",
                color="grey",
            )         
            ax.fill_between(
                    times,
                    rt_u_inter_pos, 
                    p(times),
                    color="grey",
                    alpha=0.2                
                    )  
        if abs(np.mean(rt_u_inter_neg))>0.01:         
            z = np.polyfit(times,  rt_u_inter_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$NPP*$\Delta$RT negative",
                color="grey",
                linestyle="dashed",            
            )         
            ax.fill_between(
                    times,
                    rt_u_inter_neg, 
                    p(times),
                    color="grey",
                    alpha=0.2                
                    )                         
        if abs(np.mean(delta_x_pos+delta_x_neg))>0.01:        
            z = np.polyfit(times,  delta_x_pos+delta_x_neg, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="Net $\Delta$ X",
                color="black",
            )         
            ax.fill_between(
                    times,
                    delta_x_pos+delta_x_neg, 
                    p(times),
                    color="black",
                    alpha=0.2                
                    )                                   
        ax.set_title('Contributions over time')
        ax.grid()                 
                
        # bar charts       
        
        positive_delta_X=sum(delta_x_pos)/len(delta_x_pos)
        negative_delta_X=sum(delta_x_neg)/len(delta_x_neg)
        positive_x_p_contrib=sum(x_p_contrib_pos)/len(x_p_contrib_pos)
        negative_x_p_contrib=sum(x_p_contrib_neg)/len(x_p_contrib_neg)        
        positive_rt_contrib=sum(rt_contrib_pos)/len(rt_contrib_pos)
        negative_rt_contrib=sum(rt_contrib_neg)/len(rt_contrib_neg)
        positive_u_contrib=sum(u_contrib_pos)/len(u_contrib_pos)
        negative_u_contrib=sum(u_contrib_neg)/len(u_contrib_neg)      
        positive_rt_u_inter=sum(rt_u_inter_pos)/len(rt_u_inter_pos)
        negative_rt_u_inter=sum(rt_u_inter_neg)/len(rt_u_inter_neg)        
        
        ax0=axs[1]  
        ax0.set_title('Temporal average of contributions')       
        ax0.axhline(0, color='black', ls='dashed')
        
        if abs(np.mean(positive_delta_X+negative_delta_X))>0.01:
            ax0.bar ('Net $\Delta$ X', positive_delta_X+negative_delta_X, width=0.4, color="black", label='Net $\Delta$ X')        
        ax0.bar ('$\Delta$ RT', positive_rt_contrib, color="darkorange", label='$\Delta$ RT')
        ax0.bar ('$\Delta$ RT', negative_rt_contrib, color="darkorange")        
        ax0.bar ('$\Delta$ NPP', positive_u_contrib, color="green", label='$\Delta$ NPP')
        ax0.bar ('$\Delta$ NPP', negative_u_contrib, color="green")
        ax0.bar ('$\Delta$npp*$\Delta$rt', positive_rt_u_inter, color="lightgrey", label='$\Delta$npp*$\Delta$rt')
        ax0.bar ('$\Delta$npp*$\Delta$rt', negative_rt_u_inter, color="lightgrey")        
        ax0.bar ('$\Delta$ X_p', positive_x_p_contrib, color="blue", label='$\Delta$ X_p')
        ax0.bar ('$\Delta$ X_p', negative_x_p_contrib, color="blue")   
        ax0.legend(bbox_to_anchor =(1,1))  
        ax0.grid() 

        abs_total=x_p_contrib+rt_contrib+u_contrib+rt_u_inter
        percent_x_p=x_p_contrib/abs_total*100
        percent_rt=rt_contrib/abs_total*100
        percent_u=u_contrib/abs_total*100
        percent_inter=rt_u_inter/abs_total*100
        
        ######### Percents
        if percent==True:   
            fig3=plt.figure(figsize=(15,5))
            axs=fig3.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]})         
            # if delta==False: 
            ax=axs[0]      
            
            # % timeline
       
            # quadratic trends
            z = np.polyfit(times,  percent_x_p, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ X_p",
                color="blue",
            ) 
            ax.fill_between(
                    times,
                    percent_x_p, 
                    p(times),
                    color="blue",
                    alpha=0.2                
                    )  
                    
            z = np.polyfit(times,  percent_u, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ NPP",
                color="green",
            ) 
            ax.fill_between(
                    times,
                    percent_u, 
                    p(times),
                    color="green",
                    alpha=0.2                
                    )  

            z = np.polyfit(times,  percent_rt, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$ RT",
                color="darkorange",
            )         
            ax.fill_between(
                    times,
                    percent_rt, 
                    p(times),
                    color="darkorange",
                    alpha=0.2                
                    )  

            z = np.polyfit(times,  percent_inter, 2)
            p = np.poly1d(z)  
            ax.plot(        
                times, 
                p(times),
                label="$\Delta$NPP*$\Delta$RT",
                color="grey",
            ) 
            ax.fill_between(
                    times,
                    percent_inter, 
                    p(times),
                    color="grey",
                    alpha=0.2                
                    )  
            
            ax.set_title('% Contributions over time')
            ax.set_ylabel('%')
            ax.grid()
            
            # pie charts
            ax1=axs[1]
            ax1.set_title('Temporal average % contributions')
              
            if np.mean(percent_inter) > 0.001:
                labels = '$\Delta$ RT', '$\Delta$ NPP', '$\Delta$NPP*$\Delta$RT', '$\Delta$ X_p'
                sizes = [np.mean(percent_rt), np.mean(percent_u), 
                    np.mean(percent_inter), np.mean(percent_x_p)]
                ax1.pie(sizes, autopct='%1.1f%%', 
                    startangle=90, counterclock=False, colors=("darkorange", "green", "lightgrey", "blue"))
            else:
                labels = '$\Delta$ RT','$\Delta$ NPP', '$\Delta$ X_p'
                sizes = [np.mean(percent_rt), np.mean(percent_u), np.mean(percent_x_p)]
                ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                    startangle=90, counterclock=False, colors=("darkorange", "green", "blue"))        
            ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
            ax1.legend(labels, bbox_to_anchor =(1,1))
        
        plt.show()


def plot_attribution_per_model (
    all_comp_dict,
    percent=False,
):
    models=list(all_comp_dict.keys())[:-2]        
    
    for m in models: 
    
        print ('\033[1m'+m) 
        print ('\033[1m'+'Deviation from the mean')
           
        x=all_comp_dict[m]["x"]
        x_c=all_comp_dict[m]["x_c"]
        x_p=all_comp_dict[m]["x_p"]
        u=all_comp_dict[m]["u"]
        rt=all_comp_dict[m]["rt"]   
        
        delta_x=x-all_comp_dict["Mean"]["x"]
        delta_x_c=x_c-all_comp_dict["Mean"]["x_c"]
        delta_x_p=x_p-all_comp_dict["Mean"]["x_p"]
        delta_u=u-all_comp_dict["Mean"]["u"]
        delta_rt=rt-all_comp_dict["Mean"]["rt"] 
        
        # attribution of delta X to delta X_c and delta X_p
        x_c_contrib=delta_x_c
        x_p_contrib=-delta_x_p
         
        # attribution of delta X_c to delta u and delta RT
        rt_contrib=delta_rt*(u-delta_u/2)
        u_contrib=delta_u*(rt-delta_rt/2) 
        rt_u_inter=delta_x_c-rt_contrib-u_contrib           

        # separation of positive and negative contributions
          
        pos_delta_x=delta_x.copy(); pos_delta_x[pos_delta_x<0]=0
        neg_delta_x=delta_x.copy(); neg_delta_x[neg_delta_x>0]=0
 
        pos_cont_rt=rt_contrib.copy(); pos_cont_rt[pos_cont_rt<0]=0
        neg_cont_rt=rt_contrib.copy(); neg_cont_rt[neg_cont_rt>0]=0

        pos_cont_u=u_contrib.copy(); pos_cont_u[pos_cont_u<0]=0
        neg_cont_u=u_contrib.copy(); neg_cont_u[neg_cont_u>0]=0

        pos_cont_rt_u_inter=rt_u_inter.copy(); pos_cont_rt_u_inter[pos_cont_rt_u_inter<0]=0
        neg_cont_rt_u_inter=rt_u_inter.copy(); neg_cont_rt_u_inter[neg_cont_rt_u_inter>0]=0

        pos_cont_x_p=x_p_contrib.copy(); pos_cont_x_p[pos_cont_x_p<0]=0
        neg_cont_x_p=x_p_contrib.copy(); neg_cont_x_p[neg_cont_x_p>0]=0
        
        plot_attribution (
            times=all_comp_dict["Times"],
            delta_x_pos=pos_delta_x,
            delta_x_neg=neg_delta_x,
            rt_contrib_pos=pos_cont_rt,
            rt_contrib_neg=neg_cont_rt,
            u_contrib_pos=pos_cont_u,
            u_contrib_neg=neg_cont_u,
            rt_u_inter_pos=pos_cont_rt_u_inter,
            rt_u_inter_neg=neg_cont_rt_u_inter,
            x_p_contrib_pos=pos_cont_x_p,
            x_p_contrib_neg=neg_cont_x_p,
        )
        
        print ('\033[1m'+'Deviation from other models')
        print("\033[0;0m")
        for m2 in models:       
            if m2!=m:
                print ('Deviation from '+m2)
                delta_x=x-all_comp_dict[m2]["x"]
                delta_x_c=x_c-all_comp_dict[m2]["x_c"]
                delta_x_p=x_p-all_comp_dict[m2]["x_p"]
                delta_u=u-all_comp_dict[m2]["u"]
                delta_rt=rt-all_comp_dict[m2]["rt"]
                
                # attribution of delta X to delta X_c and delta X_p
                x_c_contrib=delta_x_c
                x_p_contrib=-delta_x_p
                 
                # attribution of delta X_c to delta u and delta RT
                rt_contrib=delta_rt*(u-delta_u/2)
                u_contrib=delta_u*(rt-delta_rt/2) 
                rt_u_inter=delta_x_c-rt_contrib-u_contrib           

                # separation of positive and negative contributions                  
                pos_delta_x=delta_x.copy(); pos_delta_x[pos_delta_x<0]=0
                neg_delta_x=delta_x.copy(); neg_delta_x[neg_delta_x>0]=0
         
                pos_cont_rt=rt_contrib.copy(); pos_cont_rt[pos_cont_rt<0]=0
                neg_cont_rt=rt_contrib.copy(); neg_cont_rt[neg_cont_rt>0]=0

                pos_cont_u=u_contrib.copy(); pos_cont_u[pos_cont_u<0]=0
                neg_cont_u=u_contrib.copy(); neg_cont_u[neg_cont_u>0]=0

                pos_cont_rt_u_inter=rt_u_inter.copy(); pos_cont_rt_u_inter[pos_cont_rt_u_inter<0]=0
                neg_cont_rt_u_inter=rt_u_inter.copy(); neg_cont_rt_u_inter[neg_cont_rt_u_inter>0]=0

                pos_cont_x_p=x_p_contrib.copy(); pos_cont_x_p[pos_cont_x_p<0]=0
                neg_cont_x_p=x_p_contrib.copy(); neg_cont_x_p[neg_cont_x_p>0]=0
                
                plot_attribution (
                    times=all_comp_dict["Times"],
                    delta_x_pos=pos_delta_x,
                    delta_x_neg=neg_delta_x,
                    rt_contrib_pos=pos_cont_rt,
                    rt_contrib_neg=neg_cont_rt,
                    u_contrib_pos=pos_cont_u,
                    u_contrib_neg=neg_cont_u,
                    rt_u_inter_pos=pos_cont_rt_u_inter,
                    rt_u_inter_neg=neg_cont_rt_u_inter,
                    x_p_contrib_pos=pos_cont_x_p,
                    x_p_contrib_neg=neg_cont_x_p,
                    percent=percent,
                )

def plot_attribution_sum (
    all_comp_dict,
    percent=False,
    part=1,
):
    models=list(all_comp_dict.keys())[:-2]  

    # if we do not want the whole interval but look at a smaller part to observe the dynamics
    start_min=0
    stop_max=len(all_comp_dict["Times"])
    
    # if we do not want the whole interval but look at a smaller part to observe the dynamics
    if part < 0:
        start, stop = int(stop_max - (stop_max - start_min) * abs(part)), stop_max
    else:
        start, stop = start_min, int(start_min + (stop_max - start_min) * part)
    
    times=all_comp_dict["Times"][start:stop]   
    
    sum_pos_diff_x=0
    sum_pos_cont_rt=0
    sum_pos_cont_u=0
    sum_pos_cont_rt_u_inter=0
    sum_pos_cont_x_p=0
    
    sum_neg_diff_x=0
    sum_neg_cont_rt=0
    sum_neg_cont_u=0
    sum_neg_cont_rt_u_inter=0
    sum_neg_cont_x_p=0  
    
    print ('\033[1m'+'Attribution of summed deviations from the mean for all models ' +
        'to the differences in traceable components')     
    for m in models: 
           
        x=all_comp_dict[m]["x"][start:stop]
        x_c=all_comp_dict[m]["x_c"][start:stop]
        x_p=all_comp_dict[m]["x_p"][start:stop]
        if "gpp" in all_comp_dict[m].keys():
            u=all_comp_dict[m]["gpp"][start:stop]
        else:     
            u=all_comp_dict[m]["u"][start:stop]
        rt=all_comp_dict[m]["rt"][start:stop]
        x_mean=all_comp_dict["Mean"]["x"][start:stop]
        x_c_mean=all_comp_dict["Mean"]["x_c"][start:stop]
        x_p_mean=all_comp_dict["Mean"]["x_p"][start:stop]
        if "gpp" in all_comp_dict["Mean"].keys():
            u_mean=all_comp_dict["Mean"]["gpp"][start:stop]
        else:    
            u_mean=all_comp_dict["Mean"]["u"][start:stop]
        rt_mean=all_comp_dict["Mean"]["rt"][start:stop]           
            
        delta_x=x-x_mean
        delta_x_c=x_c-x_c_mean
        delta_x_p=x_p-x_p_mean
        delta_u=u-u_mean
        delta_rt=rt-rt_mean
        
        # attribution of delta X to delta X_c and delta X_p
        x_c_contrib=delta_x_c
        x_p_contrib=-delta_x_p
         
        # attribution of delta X_c to delta u and delta RT
        rt_contrib=delta_rt*(u-delta_u/2)
        u_contrib=delta_u*(rt-delta_rt/2) 
        rt_u_inter=delta_x_c-rt_contrib-u_contrib         

        # summation of positive and negative contributions separately         
        
        pos_delta_x=delta_x.copy(); pos_delta_x[pos_delta_x<0]=0
        neg_delta_x=delta_x.copy(); neg_delta_x[neg_delta_x>0]=0
        sum_pos_diff_x+=pos_delta_x
        sum_neg_diff_x+=neg_delta_x       

        pos_cont_rt=rt_contrib.copy(); pos_cont_rt[pos_cont_rt<0]=0
        neg_cont_rt=rt_contrib.copy(); neg_cont_rt[neg_cont_rt>0]=0
        sum_pos_cont_rt+=pos_cont_rt
        sum_neg_cont_rt+=neg_cont_rt

        pos_cont_u=u_contrib.copy(); pos_cont_u[pos_cont_u<0]=0
        neg_cont_u=u_contrib.copy(); neg_cont_u[neg_cont_u>0]=0
        sum_pos_cont_u+=pos_cont_u
        sum_neg_cont_u+=neg_cont_u

        pos_cont_rt_u_inter=rt_u_inter.copy(); pos_cont_rt_u_inter[pos_cont_rt_u_inter<0]=0
        neg_cont_rt_u_inter=rt_u_inter.copy(); neg_cont_rt_u_inter[neg_cont_rt_u_inter>0]=0
        sum_pos_cont_rt_u_inter+=pos_cont_rt_u_inter
        sum_neg_cont_rt_u_inter+=neg_cont_rt_u_inter

        pos_cont_x_p=x_p_contrib.copy(); pos_cont_x_p[pos_cont_x_p<0]=0
        neg_cont_x_p=x_p_contrib.copy(); neg_cont_x_p[neg_cont_x_p>0]=0
        sum_pos_cont_x_p+=pos_cont_x_p
        sum_neg_cont_x_p+=neg_cont_x_p 
        
    plot_attribution (
        times=times,
        delta_x_pos=sum_pos_diff_x,
        delta_x_neg=sum_neg_diff_x,
        rt_contrib_pos=sum_pos_cont_rt,
        rt_contrib_neg=sum_neg_cont_rt,
        u_contrib_pos=sum_pos_cont_u,
        u_contrib_neg=sum_neg_cont_u,
        rt_u_inter_pos=sum_pos_cont_rt_u_inter,
        rt_u_inter_neg=sum_neg_cont_rt_u_inter,
        x_p_contrib_pos=sum_pos_cont_x_p,
        x_p_contrib_neg=sum_neg_cont_x_p,
        percent=percent,
    )    

def plot_single_trend(var, times, polynom_order,title):
    fig=plt.figure(figsize=(15,5))
    ax=fig.subplots()
    z = np.polyfit(times, var, polynom_order)
    p = np.poly1d(z)  
    ax.fill_between(
            times,
            var, 
            p(times),
            color="black",
            label="original value",
            alpha=0.2                
            ) 
    ax.plot(        
        times, 
        p(times),
        label="polynomial trend",
        color="black",
        ) 
    ax.grid()
    ax.set_title(title)
    ax.legend()
    plt.show()



    
def get_vars_all_list(model_folders, experiment_names):
    vars_all_list = []
    i=0
    for mf in model_folders:
        print('\033[1m'+mf)
        print ('\033[0m')
        current_var_list = gh.msh(mf).get_global_mean_vars_all(experiment_name=experiment_names[i])
        vars_all_list.append(current_var_list)
        i+=1
    print("Done!")
    return vars_all_list        
    
def get_components_from_output(
    model_names,  # dictionary (folder name : model name)
    vars_all_list,  # a list of test_args from all models involved
    delta_t_val,  # model time step
    #model_cols,  # dictionary (folder name :color)
    part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    averaging=12,  # number of iterator steps over which to average results. 12 - avarage monthly to yearly 
    overlap=True,  # compute overlapping timeframe or plot whole duration for all models
):

    def times_in_days_aD(
        vars_all,  # a tuple defined in model-specific test_helpers.py
        delta_t_val,  # iterator time step
        start_date,
    ):
        n_months = len(vars_all.gpp)*12
        #n_days = month_2_day_index([n_months], start_date)[0]
        n_iter = n_months #int(n_days / delta_t_val)
        #print("n_months: "+ str(n_months))
        #print("n_days: "+ str(n_days))
        #print("n_iter: "+ str(n_iter))
        return np.array(
            tuple(
                (
                    days_since_AD(i, delta_t_val, start_date)
                    for i in np.arange(n_iter)  # days_after_sim_start
                )
            )
        )
    # function to determine overlapping time frames for models simulations
    def t_min_tmax_overlap(
        vars_all_list,  # a list of test_args from all models involved
        delta_t_val,  # iterator time step
        start_date,
    ):
        delta_t_val=1
        td = {
            i: times_in_days_aD(vars_all_list[i], delta_t_val, start_date)
            for i in range(len(vars_all_list))
        }
        t_min = max([t.min() for t in td.values()])
        t_max = min([t.max() for t in td.values()])
        return (t_min, t_max)
    
    def min_max_index(
    vars_all,  # a tuple defined in model-specific test_helpers.py
    delta_t_val,  # iterator time step
    t_min,  # output of t_min_tmax_overlap or t_min_tmax_full
    t_max,  # output of t_min_tmax_overlap or t_min_tmax_full
    start_date,
    ):
        ts = times_in_days_aD(vars_all, delta_t_val, start_date)

        def count(acc, i):
            min_i, max_i = acc
            t = ts[i]
            min_i = min_i + 1 if t < t_min else min_i
            max_i = max_i + 1 if t < t_max else max_i
            return (min_i, max_i)

        return reduce(count, range(len(ts)), (0, 0))

    if (part < -1) | (part > 1) | (part == 0):
        raise Exception(
            "Invalid partitioning in plot_components_combined: use part between -1 and 1 excluding 0"
        )
    model_folders = [(k) for k in model_names]
    k = 0
    sum_cVeg = np.array(0)
    sum_cSoil = np.array(0)
    sum_gpp = np.array(0)    
    sum_npp = np.array(0)
    sum_ra = np.array(0)
    sum_rh = np.array(0)
    sum_x = np.array(0)
    sum_rt = np.array(0)
    sum_x_c = np.array(0)
    sum_x_p = np.array(0)
    sum_nep = np.array(0)    
    sum_rt_soil = np.array(0)
    sum_f_v2s = np.array(0)
    sum_f_v2s_2 = np.array(0)    
    sum_rt_soil = np.array(0)
    sum_rt_soil_2 = np.array(0)    
    sum_rt_veg = np.array(0)
    sum_rt_veg_2 = np.array(0)    
    sum_x_c_soil = np.array(0)
    sum_x_c_veg = np.array(0)
    sum_x_p_veg = np.array(0)
    sum_x_p_soil = np.array(0)
    sum_ra_rate = np.array(0)
    sum_rh_rate = np.array(0)
    sum_delta_cVeg = np.array(0)
    sum_delta_cSoil = np.array(0)
    sum_delta_dist = np.array(0)
    sum_dist = np.array(0)
    all_components = list ()
    
    def annual_delta (data_stream):
        delta=np.zeros_like(data_stream)
        for i in range(len(data_stream)-1):
            delta[i]=data_stream[i+1]-data_stream[i]  
        return delta

    min_len=1000
    for i in range (len(model_folders)):
        if len(vars_all_list[i].cVeg)<min_len: min_len=len(vars_all_list[i].cVeg)
    print(min_len)
    
    for mf in model_folders:
        print ("Getting traceable components for "+mf+"...")
        start_date=gh.msh(mf).start_date()       
        
        stop=len(vars_all_list[k].cVeg)#-20
        start=len(vars_all_list[k].cVeg)-min_len
        
        times=np.array(range(start,stop))+1700
        
        # harmonising model outputs
        start_flux=start
        stop_flux=stop
        #print(len(vars_all_list[k].cVeg)>1000)
        # if len(vars_all_list[k].cVeg)>1000:
            # start_pool=start//12
            # stop_pool=stop//12+1
        # else:
        start_pool=start
        stop_pool=stop       
           
        cVeg=vars_all_list[k].cVeg[start_pool:stop_pool]
        if cVeg.shape[0]>500: cVeg=gh.avg_timeline(cVeg, averaging)
        cSoil=vars_all_list[k].cSoil[start_pool:stop_pool]
        if cSoil.shape[0]>500: cSoil=gh.avg_timeline(cSoil, averaging)
        gpp=vars_all_list[k].gpp[start_flux:stop_flux]
        if gpp.shape[0]>500: gpp=gh.avg_timeline(gpp, averaging)        
        npp=vars_all_list[k].npp[start_flux:stop_flux]
        if npp.shape[0]>500: npp=gh.avg_timeline(npp, averaging)
        rh=vars_all_list[k].rh[start_flux:stop_flux]
        if rh.shape[0]>500: rh=gh.avg_timeline(rh, averaging)
        ra=vars_all_list[k].ra[start_flux:stop_flux]
        if ra.shape[0]>500: ra=gh.avg_timeline(ra, averaging)        
        # print(times.shape)
        # print(cVeg.shape)
        # print(cSoil.shape)
        # print(npp.shape)
        # print(rh.shape)
        # print(start_pool, stop_pool)
        # print(start_flux, stop_flux)        
        
        # correction for ISAM data issue
        if mf=="cj_isam": 
            # print (cVeg[-9])
            # print(cVeg[-8])
            # print(cVeg[-10])
            cVeg[-9]=np.mean((cVeg[-8], cVeg[-10]))
        
        # traditional traceable components 
        n_years=len(cSoil)
        #print(n_years)
        x=cVeg+cSoil
        nep=annual_delta(x)
        delta_cVeg=annual_delta(cVeg)
        delta_cSoil=annual_delta(cSoil)
        rt=x/(rh+ra)
        x_c=rt*gpp
        x_p=x_c-x
        # flux-based                
        dist=gpp-nep-ra-rh
        f_v2s=delta_cSoil+rh
        f_v2s_2=gpp-delta_cVeg-ra-dist
        rt_soil=cSoil/rh
        rt_soil_2=cSoil/(f_v2s-delta_cSoil)
        rt_veg=cVeg/(ra+dist+f_v2s)
        rt_veg_2=cVeg/(gpp-delta_cVeg)
        x_c_soil=f_v2s*rt_soil
        x_c_veg=gpp*rt_veg
        x_p_soil=x_c_soil-cSoil
        x_p_veg=x_c_veg-cVeg
        ra_rate=ra/cVeg
        rh_rate=rh/cSoil
        
        comp_dict = {
            "cVeg": cVeg,
            "cSoil": cSoil,
            "delta_cVeg": delta_cVeg,
            "delta_cSoil": delta_cSoil,
            "gpp":gpp,
            "npp":npp,
            "nep":nep,            
            "rh":rh,
            "ra":ra,
            "f_v2s":f_v2s,
            "f_v2s_2":f_v2s_2,
            "dist":dist,
            "x":x,
            "rt":rt,
            "x_c":x_c,
            "x_p":x_p,
            "rt_veg":rt_veg,
            "rt_veg_2":rt_veg_2,            
            "rt_soil":rt_soil,
            "rt_soil_2":rt_soil_2,
            "x_c_soil":x_c_soil,
            "x_c_veg":x_c_veg,
            #"x_c_veg2":x_c_veg2,
            "x_p_veg":x_p_veg,
            "x_p_soil":x_p_soil,
            #"x_p_veg2":x_p_veg2,
            "ra_rate":ra_rate,
            "rh_rate":rh_rate,
            }  
                                           
        all_components.append(comp_dict)
        if sum_cVeg.all()==0:
            sum_cVeg = np.append(sum_cVeg,cVeg)[1:]
            sum_cSoil = np.append(sum_cSoil,cSoil)[1:]
            sum_gpp = np.append(sum_gpp,gpp)[1:]
            sum_npp = np.append(sum_npp,npp)[1:]
            sum_rh = np.append(sum_rh,rh)[1:]
            sum_ra = np.append(sum_ra,ra)[1:]            
            sum_x = np.append(sum_x,x)[1:]
            sum_rt = np.append(sum_rt,rt)[1:]
            sum_x_c = np.append(sum_x_c,x_c)[1:]
            sum_x_p = np.append(sum_x_p,x_p)[1:]
            #sum_nep = np.append(sum_nep,nep)[1:]  
            sum_f_v2s = np.append(sum_f_v2s,f_v2s)[1:]
            sum_f_v2s_2 = np.append(sum_f_v2s_2,f_v2s_2)[1:]            
            sum_rt_soil = np.append(sum_rt_soil,rt_soil)[1:]
            sum_rt_soil_2 = np.append(sum_rt_soil_2,rt_soil_2)[1:]            
            sum_rt_veg = np.append(sum_rt_veg,rt_veg)[1:]
            sum_rt_veg_2 = np.append(sum_rt_veg_2,rt_veg_2)[1:]            
            sum_x_c_soil = np.append(sum_x_c_soil,x_c_soil)[1:]
            sum_x_c_veg = np.append(sum_x_c_veg,x_c_veg)[1:]
            sum_x_p_veg = np.append(sum_x_p_veg,x_p_veg)[1:]
            sum_x_p_soil = np.append(sum_x_p_soil,x_p_soil)[1:]            
            sum_ra_rate = np.append(sum_ra_rate,ra_rate)[1:]
            sum_rh_rate = np.append(sum_rh_rate,rh_rate)[1:] 
            sum_dist = np.append(sum_dist,dist)[1:] 
            sum_delta_cVeg = np.append(sum_delta_cVeg,delta_cVeg)[1:] 
            sum_delta_cSoil = np.append(sum_delta_cSoil,delta_cSoil)[1:]             
        else:
            sum_cVeg = sum_cVeg + cVeg
            sum_cSoil = sum_cSoil + cSoil
            sum_gpp = sum_gpp + gpp
            sum_npp = sum_npp + npp
            sum_rh = sum_rh + rh
            sum_ra = sum_ra + ra            
            sum_x = sum_x + x
            sum_rt = sum_rt + rt
            sum_x_c = sum_x_c + x_c
            sum_x_p = sum_x_p + x_p
            sum_nep = sum_nep + nep  
            sum_f_v2s = sum_f_v2s + f_v2s
            sum_f_v2s_2 = sum_f_v2s_2 + f_v2s_2            
            sum_rt_soil = sum_rt_soil + rt_soil
            sum_rt_soil_2 = sum_rt_soil_2 + rt_soil_2            
            sum_rt_veg = sum_rt_veg + rt_veg
            sum_rt_veg_2 = sum_rt_veg_2 + rt_veg_2            
            sum_x_c_soil = sum_x_c_soil + x_c_soil
            sum_x_c_veg = sum_x_c_veg + x_c_veg
            sum_x_p_veg = sum_x_p_veg + x_p_veg
            sum_x_p_soil = sum_x_p_soil + x_p_soil
            sum_ra_rate = sum_ra_rate + ra_rate
            sum_rh_rate = sum_rh_rate + rh_rate 
            sum_dist = sum_dist + dist
            sum_delta_cVeg = sum_delta_cVeg + delta_cVeg
            sum_delta_cSoil = sum_delta_cSoil + delta_cSoil
        k += 1
    ave_dict = {
            "cVeg": sum_cVeg / len(model_folders),
            "cSoil": sum_cSoil / len(model_folders), 
            "delta_cVeg": sum_delta_cVeg / len(model_folders),
            "delta_cSoil": sum_delta_cSoil / len(model_folders),
            "nep": sum_nep / len(model_folders),
            "f_v2s": sum_f_v2s / len(model_folders),  
            "f_v2s_2": sum_f_v2s_2 / len(model_folders),            
            "gpp": sum_gpp / len(model_folders), 
            "rh": sum_rh / len(model_folders),
            "ra": sum_ra / len(model_folders),
            "x": sum_x / len(model_folders), 
            "rt": sum_rt / len(model_folders),
            "x_c": sum_x_c / len(model_folders),
            "x_p": sum_x_p / len(model_folders), 
            #"nep": sum_nep,# / len(model_folders), 
            "dist": sum_dist / len(model_folders),
            "rt_veg": sum_rt_veg / len(model_folders),
            "rt_veg_2": sum_rt_veg_2 / len(model_folders),              
            "rt_soil": sum_rt_soil / len(model_folders), 
            "rt_soil_2": sum_rt_soil_2 / len(model_folders),
            "x_c_veg": sum_x_c_veg / len(model_folders),            
            "x_c_soil": sum_x_c_soil / len(model_folders),
            "x_p_veg": sum_x_p_veg / len(model_folders),
            "x_p_soil": sum_x_p_soil / len(model_folders),
            "ra_rate": sum_ra_rate / len(model_folders),
            "rh_rate": sum_rh_rate / len(model_folders),
            }    
    all_components.append(ave_dict) 
    
    times_avg = times
    #times_avg = gh.avg_timeline(times, averaging)
    all_components.append(times_avg)  
    mods=list(model_names.values())
    mods.append("Mean")
    mods.append("Times")
    all_comp_dict = {mods[i]: all_components[i] for i in range(len(mods))}             
    return all_comp_dict       

def add_gridded_vars (
        model_names,
        experiment_names,
        global_mask,       
        ):
    model_folders=[(m) for m in model_names] 
    m_names=list(model_names.values())      
    g_mask=global_mask.index_mask    
    for experiment in experiment_names:
        print('\033[1m'+'. . . Adding variables for '+experiment+' experiment . . .')   
        print('\033[0m')                         
        k=0 # model counter         
        for mf in model_folders:
            print('\033[1m'+model_names[mf]+'\033[0m')
            #print('cSoil_total')
            experiment_name=m_names[k]+"_"+experiment+"_"
            conf_dict = gh.confDict(mf)
            dataPath=Path(conf_dict["dataPath"])              
            file_path = dataPath.joinpath(experiment_name+"cSoil_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_soil_avg=ds.variables["cSoil"][:, :].data
            var_soil_diff=ds.variables["cSoil_diff"][:, :].data
            ds.close()
            if "cLitter" in gh.msh(mf).data_str._fields:           
                file_path = dataPath.joinpath(experiment_name+"cLitter_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_litter_avg=ds.variables["cLitter"][:, :].data
                var_litter_diff=ds.variables["cLitter_diff"][:, :].data                
                ds.close()
                var_soil_total_avg=np.ma.array(var_soil_avg+var_litter_avg, mask=g_mask)
                var_soil_total_diff=np.ma.array(var_soil_diff+var_litter_diff, mask=g_mask)                
            else:
                var_soil_total_avg=np.ma.array(var_soil_avg, mask=g_mask)
                var_soil_total_diff=np.ma.array(var_soil_diff, mask=g_mask)     
            var_soil_total_avg[var_soil_total_avg<0.00001]=0.00001
            #print(np.ma.mean(var_soil_total_avg))                
            #print(np.ma.mean(var_soil_total_diff))                
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"cSoil_total_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var_avg = ds_new.createVariable('cSoil_total', "float32", ["lat", "lon"])
            var_avg[:, :] = var_soil_total_avg
            var_diff = ds_new.createVariable('cSoil_total_diff', "float32", ["lat", "lon"])
            var_diff[:, :] = var_soil_total_diff            
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()                      
            if "npp_nlim" in gh.msh(mf).data_str._fields:  
                #print('npp_nlim')            
                file_path = dataPath.joinpath(experiment_name+"npp_nlim_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_npp_nlim_avg=ds.variables["npp_nlim"][:, :].data
                var_npp_nlim_diff=ds.variables["npp_nlim_diff"][:, :].data                
                ds.close()
                var_npp_avg=np.ma.array(var_npp_nlim_avg, mask=g_mask)
                var_npp_diff=np.ma.array(var_npp_nlim_diff, mask=g_mask)
                #print(np.ma.mean(var_npp_avg))                
                #print(np.ma.mean(var_npp_diff))
                # creating and writing a new NetCDF file 
                s = g_mask.shape
                n_lats, n_lons = s
                new_path=dataPath.joinpath(experiment_name+"npp_ave_res.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)         
                # creating variables                         
                var_avg = ds_new.createVariable('npp', "float32", ["lat", "lon"])
                var_avg[:, :] = var_npp_avg
                var_diff = ds_new.createVariable('npp_diff', "float32", ["lat", "lon"])
                var_diff[:, :] = var_npp_diff            
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
                # closing NetCDF file      
                ds_new.close()                       
            elif not ("npp" in gh.msh(mf).data_str._fields):  
                #print('npp')            
                file_path = dataPath.joinpath(experiment_name+"gpp_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_gpp_avg=ds.variables["gpp"][:, :].data
                var_gpp_diff=ds.variables["gpp_diff"][:, :].data                
                ds.close()
                file_path = dataPath.joinpath(experiment_name+"ra_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_ra_avg=ds.variables["ra"][:, :].data
                var_ra_diff=ds.variables["ra_diff"][:, :].data                
                ds.close()                                              
                var_npp_avg=np.ma.array(var_gpp_avg-var_ra_avg, mask=g_mask)
                var_npp_diff=np.ma.array(var_gpp_diff-var_ra_diff, mask=g_mask)            
                # creating and writing a new NetCDF file 
                s = g_mask.shape
                n_lats, n_lons = s
                new_path=dataPath.joinpath(experiment_name+"npp_ave_res.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)         
                # creating variables                         
                var_avg = ds_new.createVariable('npp', "float32", ["lat", "lon"])
                var_avg[:, :] = var_npp_avg
                var_diff = ds_new.createVariable('npp_diff', "float32", ["lat", "lon"])
                var_diff[:, :] = var_npp_diff            
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
                # closing NetCDF file      
                ds_new.close() 
            if not ("ra" in gh.msh(mf).data_str._fields):  
                #print('ra')             
                file_path = dataPath.joinpath(experiment_name+"gpp_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_gpp_avg=ds.variables["gpp"][:, :].data
                var_gpp_diff=ds.variables["gpp_diff"][:, :].data                
                ds.close()
                file_path = dataPath.joinpath(experiment_name+"npp_ave_res.nc")
                ds = nc.Dataset(str(file_path))
                var_npp_avg=ds.variables["npp"][:, :].data
                var_npp_diff=ds.variables["npp_diff"][:, :].data                
                ds.close()                                              
                var_ra_avg=np.ma.array(var_gpp_avg-var_npp_avg, mask=g_mask)
                var_ra_diff=np.ma.array(var_gpp_diff-var_npp_diff, mask=g_mask)    
                #print(np.ma.mean(var_ra_avg))                
                #print(np.ma.mean(var_ra_diff))                  
                # creating and writing a new NetCDF file 
                s = g_mask.shape
                n_lats, n_lons = s
                new_path=dataPath.joinpath(experiment_name+"ra_ave_res.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)         
                # creating variables                         
                var_avg = ds_new.createVariable('ra', "float32", ["lat", "lon"])
                var_avg[:, :] = var_ra_avg
                var_diff = ds_new.createVariable('ra_diff', "float32", ["lat", "lon"])
                var_diff[:, :] = var_ra_diff            
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
                # closing NetCDF file      
                ds_new.close()
                
            # adding variables for uncertianty propagation
            
            # adding NEP = cVeg_diff + cSoil_total_diff  
            #print('nep')           
            file_path = dataPath.joinpath(experiment_name+"cVeg_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_cVeg_diff_array=ds.variables["cVeg_diff"][:, :].data          
            var_cVeg_diff=np.ma.array(var_cVeg_diff_array, mask=g_mask)
            var_cVeg_array=ds.variables["cVeg"][:, :].data 
            var_cVeg_array[var_cVeg_array<=0.00001]=0.00001           
            var_cVeg=np.ma.array(var_cVeg_array,mask=g_mask)
            #print('cVEG!!! '+star(np.ma.mean(var_cVeg)))
            ds.close()                                        
            var_nep=np.ma.array(var_cVeg_diff+var_soil_total_diff, mask=g_mask)  
            #print(np.ma.mean(var_nep))                                                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"nep_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('nep', "float32", ["lat", "lon"])
            var[:, :] = var_nep   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()
            
            # adding C loss from disturbances: dist = gpp - nep - ra - rh 
            #print('disturbance')           
            file_path = dataPath.joinpath(experiment_name+"gpp_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_gpp_array=ds.variables["gpp"][:, :].data
            var_gpp=np.ma.array(var_gpp_array, mask = g_mask)
            ds.close()                                        
            file_path = dataPath.joinpath(experiment_name+"ra_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_ra_array=ds.variables["ra"][:, :].data 
            var_ra = np.ma.array(var_ra_array, mask=g_mask) 
            var_ra [var_ra<=0.00001]=0.00001            
            ds.close()  
            file_path = dataPath.joinpath(experiment_name+"rh_ave_res.nc")
            ds = nc.Dataset(str(file_path))
            var_rh_array=ds.variables["rh"][:, :].data 
            var_rh = np.ma.array(var_rh_array, mask=g_mask) 
            var_rh [var_rh<=0.00001]=0.00001  
            ds.close()  
            var_dist = np.ma.array(var_gpp - var_nep - var_ra - var_rh, mask = g_mask)
            var_dist[var_dist<0]=0  
            #print(np.ma.mean(var_dist))            
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"dist_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('dist', "float32", ["lat", "lon"])
            var[:, :] = var_dist   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding flux from vegetation to soil (cSoil_diff + rh) or (gpp - cVeg_diff - ra - dist)
            #print ('adding f_v2s')  
            #veg_dist_frac=1            
            var_f_v2s=np.ma.array(var_soil_diff+var_rh, mask = g_mask)
            var_f_v2s[var_f_v2s<=0]=0.00001  
            var_f_v2s_2=np.ma.array(var_gpp-var_cVeg_diff-var_ra-var_dist, mask = g_mask)
            var_f_v2s_2[var_f_v2s_2<=0]=0.00001 
            #print(np.ma.mean(var_f_v2s_1))
            #print(np.ma.mean(var_f_v2s_2))
            #print(np.ma.mean(var_f_v2s_1-var_f_v2s_2))
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"f_v2s_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('f_v2s', "float32", ["lat", "lon"])
            var[:, :] = var_f_v2s  
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding vegetation residence time (RT_veg = cVeg / (f_veg2soil + ra + dist) )
            #print ('adding RT_veg')                      
            var_RT_veg=np.ma.array(var_cVeg / (var_ra + var_dist + var_f_v2s) * 120, mask = g_mask)
            var_RT_veg2=np.ma.array(var_cVeg / (var_gpp - var_cVeg_diff) * 120, mask = g_mask)              
            #print(np.ma.mean(var_RT_veg))            
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"RT_veg_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('RT_veg', "float32", ["lat", "lon"])
            var[:, :] = var_RT_veg   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding soil residence time (RT_soil = cSoil / rh )
            #print ('adding RT_soil')           
            var_RT_soil=np.ma.array(var_soil_total_avg / var_rh * 120, mask = g_mask)            
            #var_RT_soil2=np.ma.array(var_soil_total_avg / (var_f_v2s - var_soil_total_diff) * 120, mask = g_mask)
            var_RT_soil[var_RT_soil>10000]=10000
            #print(np.ma.mean(var_RT_soil))            
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"RT_soil_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('RT_soil', "float32", ["lat", "lon"])
            var[:, :] = var_RT_soil   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()
            
            # adding total ecosystem carbon ( X = cVeg + cSoil_total )
            var_X=np.ma.array(var_cVeg + var_soil_total_avg, mask = g_mask)                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X', "float32", ["lat", "lon"])
            var[:, :] = var_X   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding ecosystem residence time ( RT = X / (ra + rh + dist) )
            var_RT=np.ma.array(var_X / (var_ra + var_rh + var_dist)* 120, mask = g_mask)                     
            var_RT[var_RT>10000]=10000
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"RT_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('RT', "float32", ["lat", "lon"])
            var[:, :] = var_RT   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()
            
            # adding equilibrium carbon storage capacity ( X_c = gpp * RT )
            var_X_c=np.ma.array(var_gpp * var_RT / 120, mask = g_mask)                     
            # print('X_c=gpp*RT/120: '+str(np.ma.mean(var_X_c))+' = '+str(np.ma.mean(var_gpp))+
                # ' * '+str(np.ma.mean(var_RT))+ ' / 120')
            var_X_c[var_X_c<0]=0.000001  
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_c_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_c', "float32", ["lat", "lon"])
            var[:, :] = var_X_c   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding carbon storage potential ( X_p = X_c - X )
            var_X_p=np.ma.array(var_X_c - var_X, mask = g_mask)                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_p_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_p', "float32", ["lat", "lon"])
            var[:, :] = var_X_p   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()
            
            # adding equilibrium carbon storage capacity for Vegetation ( X_c_veg = gpp * RT_veg )
            var_X_c_veg=np.ma.array(var_gpp * var_RT_veg / 120, mask = g_mask)                     
            # print('X_c=gpp*RT/120: '+str(np.ma.mean(var_X_c))+' = '+str(np.ma.mean(var_gpp))+
                # ' * '+str(np.ma.mean(var_RT))+ ' / 120')
            var_X_c_veg[var_X_c_veg<0]=0.000001  
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_c_veg_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_c_veg', "float32", ["lat", "lon"])
            var[:, :] = var_X_c_veg   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding carbon storage potential for Vegetation ( X_p_veg = X_c_veg - C_Veg )
            var_X_p_veg=np.ma.array(var_X_c_veg - var_cVeg, mask = g_mask)                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_p_veg_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_p_veg', "float32", ["lat", "lon"])
            var[:, :] = var_X_p_veg   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()            
            
            # adding equilibrium carbon storage capacity for soil ( X_c_soil = f_v2s * RT_soil )
            var_X_c_soil=np.ma.array(var_f_v2s * var_RT_soil / 120, mask = g_mask)                     
            # print('X_c=gpp*RT/120: '+str(np.ma.mean(var_X_c))+' = '+str(np.ma.mean(var_gpp))+
                # ' * '+str(np.ma.mean(var_RT))+ ' / 120')
            var_X_c_soil[var_X_c_soil<0]=0.000001  
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_c_soil_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_c_soil', "float32", ["lat", "lon"])
            var[:, :] = var_X_c_soil   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()

            # adding carbon storage potential for soiletation ( X_p_soil = X_c_soil - C_soil )
            var_X_p_soil=np.ma.array(var_X_c_soil - var_soil_total_avg, mask = g_mask)                     
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=dataPath.joinpath(experiment_name+"X_p_soil_ave_res.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var = ds_new.createVariable('X_p_soil', "float32", ["lat", "lon"])
            var[:, :] = var_X_p_soil   
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close()             
            
            k+=1 # model counter

            print('cVeg: '+ str(np.ma.mean(var_cVeg)))
            print('cVeg_max: '+ str(np.ma.max(var_cVeg)))            
            print('cVeg_min: '+ str(np.ma.min(var_cVeg)))
            print('cVeg_diff: '+ str(np.ma.mean(var_cVeg_diff)))
            print('cSoil: '+ str(np.ma.mean(var_soil_total_avg)))
            print('cSoil_min: '+ str(np.ma.min(var_soil_total_avg)))            
            #print('cSoil_max: '+ str(np.ma.max(var_soil_total_avg))) 
            #print('cSoil_min: '+ str(np.ma.min(var_soil_total_avg)))             
            print('cSoil_diff: '+ str(np.ma.mean(var_soil_total_diff)))
            print('X: '+str(np.ma.mean(var_X)))    
            print('X_max: '+str(np.ma.max(var_X))) 
            print('X_min: '+str(np.ma.min(var_X)))             
            print('nep: '+str(np.ma.mean(var_nep)))
            print('gpp: '+str(np.ma.mean(var_gpp)))
            #print('gpp_max: '+str(np.ma.max(var_gpp)))
            #print('gpp_min: '+str(np.ma.min(var_gpp)))
            #print('npp: '+str(np.mean(var_npp)))            
            print('ra: '+str(np.ma.mean(var_ra)))
            print('ra_max: '+str(np.ma.max(var_ra)))
            print('ra_min: '+str(np.ma.min(var_ra)))
            print('rh: '+str(np.ma.mean(var_rh)))
            print('rh_max: '+str(np.ma.max(var_rh)))
            print('rh_min: '+str(np.ma.min(var_rh)))            
            print('dist: '+str(np.ma.mean(var_dist)))
            print('dist_max: '+str(np.ma.max(var_dist)))
            print('dist_min: '+str(np.ma.min(var_dist)))            
            print('f_v2s: '+str(np.ma.mean(var_f_v2s)))
            #print('f_v2s_max: '+str(np.ma.max(var_f_v2s)))
            #print('f_v2s_min: '+str(np.ma.min(var_f_v2s)))             
            print('f_v2s2: '+str(np.ma.mean(var_f_v2s_2)))
            
            print('RT_veg: '+str(np.ma.mean(var_RT_veg)))
            print('RT_veg_max: '+str(np.ma.max(var_RT_veg)))
            print('RT_veg_min: '+str(np.ma.min(var_RT_veg)))            
            print('RT_veg2: '+str(np.ma.mean(var_RT_veg2)))
            print('RT_soil: '+str(np.ma.mean(var_RT_soil))) 
            print('RT_soil_max: '+str(np.ma.max(var_RT_soil)))   
            print('RT_soil_min: '+str(np.ma.min(var_RT_soil)))               
            #print('RT_soil2: '+str(np.ma.mean(var_RT_soil2)))               
            #print('dist_veg_error: '+str(np.ma.mean( var_gpp - var_cVeg_diff - var_ra - var_dist - var_f_v2s_1 )))    
            print('RT: '+str(np.ma.mean(var_RT))) 
            print('RT_max: '+str(np.ma.max(var_RT)))  
            print('RT_min: '+str(np.ma.min(var_RT)))  
            print('X_c: '+str(np.ma.mean(var_X_c))) 
            print('X_c_max: '+str(np.ma.max(var_X_c)))  
            print('X_c_min: '+str(np.ma.min(var_X_c))) 
            print('X_p: '+str(np.ma.mean(var_X_p))) 
            print('X_p_max: '+str(np.ma.max(var_X_p)))  
            print('X_p_min: '+str(np.ma.min(var_X_p)))  

            print('X_c_veg '+str(np.ma.max(var_X_c_veg)))  
            print('X_p_veg '+str(np.ma.max(var_X_p_veg)))  
            print('X_c_soil '+str(np.ma.max(var_X_c_soil)))  
            print('X_p_soil '+str(np.ma.max(var_X_p_soil)))              
            
    print('\033[1m'+'Done!')
    
def uncertainty_grids (
        model_names,
        experiment_names,
        global_mask,
        output_path,
        ):
    model_folders=[(m) for m in model_names]
    m_names=list(model_names.values())
    g_mask=global_mask.index_mask
    for experiment in experiment_names:
        print('\033[1m'+'. . . Computing uncertainty for '+experiment+' experiment . . .')
        print('\033[0m')
                  
        print('computing uncertainties...')
        for vn in ['cVeg', 'cSoil_total', 'gpp','ra','rh','f_v2s','nep',
            'RT_veg', 'RT_soil', 'dist', 'X', 'RT', 'X_c', 'X_p', 'X_c_veg', 'X_p_veg', 'X_c_soil', 'X_p_soil']:
            #print(vn)
            var_sum_zero=np.zeros_like(g_mask)
            var_sum=np.ma.array(var_sum_zero,mask = g_mask)
            var_diff_sqr_zero=np.zeros_like(g_mask)
            var_diff_sqr=np.ma.array(var_diff_sqr_zero,mask = g_mask)
            var_abs_diff_zero=np.zeros_like(g_mask)
            var_abs_diff=np.ma.array(var_abs_diff_zero,mask = g_mask)            
            delta_sum_zero=np.zeros_like(g_mask)
            delta_sum=np.ma.array(delta_sum_zero,mask = g_mask)
            delta_diff_sqr_zero=np.zeros_like(g_mask)
            delta_diff_sqr=np.ma.array(delta_diff_sqr_zero,mask = g_mask)   
            delta_abs_diff_zero=np.zeros_like(g_mask)
            delta_abs_diff=np.ma.array(delta_abs_diff_zero,mask = g_mask)            
            # computing ensemble mean       
            k=0 # model counter        
            for mf in model_folders:
                experiment_name=m_names[k]+"_"+experiment+"_"
                conf_dict = gh.confDict(mf)
                dataPath=Path(conf_dict["dataPath"])       
                # if vn in ['gpp','ra','rh','cVeg','cSoil_total']: # for fluxes we use sum
                    # file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc")                
                # else: # for pools we use mean
                    # file_path = dataPath.joinpath(experiment_name+vn+"_res.nc")
                file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc") 
                ds = nc.Dataset(str(file_path))
                var=ds.variables[vn][:, :].data
                var_sum=var_sum+var
                if vn=="cVeg" or vn=="cSoil_total":
                    delta=ds.variables[vn+"_diff"][:, :].data
                    delta_sum=delta_sum+delta
                k+=1 # model counter
                ds.close() 
            mean=var_sum/len(model_folders) 
            delta_mean=delta_sum/len(model_folders)
            # computing uncertainty measures: standard deviation and average deviation                   
            k=0 # model counter        
            for mf in model_folders:
                experiment_name=m_names[k]+"_"+experiment+"_"
                conf_dict = gh.confDict(mf)
                dataPath=Path(conf_dict["dataPath"])             
                # if vn in ['gpp','ra','rh','cVeg','cSoil_total']: # for fluxes we use sum    
                    # file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc")                
                # else:
                    # file_path = dataPath.joinpath(experiment_name+vn+"_res.nc")    
                file_path = dataPath.joinpath(experiment_name+vn+"_ave_res.nc")    
                ds = nc.Dataset(str(file_path))
                var=ds.variables[vn][:, :].data           
                var_diff_sqr = var_diff_sqr + (var-mean)**2
                var_abs_diff = var_abs_diff + np.abs(var-mean)
                if vn=="cVeg" or vn=="cSoil_total":
                    delta=ds.variables[vn+"_diff"][:, :].data
                    delta_diff_sqr=delta_diff_sqr+(delta-delta_mean)**2 
                    delta_abs_diff=delta_abs_diff+np.abs(delta-delta_mean)
                k+=1 # model counter 
                ds.close()             
            variance = var_diff_sqr / (len(model_folders)-1)
            st_dev=np.sqrt(variance) 
            ave_dev=var_abs_diff / len(model_folders)
            # final masking
            var_mean_final = np.ma.array(mean,mask = g_mask)
            var_sd_final = np.ma.array(st_dev,mask = g_mask)
            var_avd_final = np.ma.array(ave_dev,mask = g_mask)
            # creating and writing a new NetCDF file 
            s = g_mask.shape
            n_lats, n_lons = s
            new_path=Path(output_path).joinpath(experiment+"_"+vn+"_uncertainty.nc")
            ds_new = nc.Dataset(str(new_path), "w", persist=True)
            # creating dimensions
            lat = ds_new.createDimension("lat", size=n_lats)
            lon = ds_new.createDimension("lon", size=n_lons)         
            # creating variables                         
            var_mean = ds_new.createVariable(vn+'_mean', "float32", ["lat", "lon"])
            var_mean[:, :] = var_mean_final
            var_sd = ds_new.createVariable(vn+'_sd', "float32", ["lat", "lon"])
            var_sd[:, :] = var_sd_final            
            var_avd = ds_new.createVariable(vn+'_avd', "float32", ["lat", "lon"])
            var_avd[:, :] = var_avd_final
            var_avd_relative = ds_new.createVariable(vn+'_avd_relative', "float32", ["lat", "lon"])
            relative_uncertainty = var_avd_final / abs(var_mean_final) *100
            relative_uncertainty [relative_uncertainty>300]=300
            var_avd_relative[:, :] = relative_uncertainty              
            lats = ds_new.createVariable("lat", "float32", ["lat"])
            lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
            lons = ds_new.createVariable("lon", "float32", ["lon"])
            lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
            # closing NetCDF file      
            ds_new.close() 
            print(vn+': '+str(np.ma.mean(var_mean_final)))
            if vn=="cVeg" or vn=="cSoil_total":
                delta_variance = delta_diff_sqr / (len(model_folders)-1)
                delta_st_dev=np.sqrt(delta_variance)
                delta_ave_dev=delta_abs_diff / len(model_folders)
                # final masking
                var_mean_final = np.ma.array(delta_mean,mask = g_mask)
                var_sd_final = np.ma.array(delta_st_dev,mask = g_mask)
                var_avd_final = np.ma.array(delta_st_dev,mask = g_mask)                
                # creating and writing a new NetCDF file 
                s = g_mask.shape
                n_lats, n_lons = s
                new_path=Path(output_path).joinpath(experiment+"_"+vn+"_diff_uncertainty.nc")
                ds_new = nc.Dataset(str(new_path), "w", persist=True)
                # creating dimensions
                lat = ds_new.createDimension("lat", size=n_lats)
                lon = ds_new.createDimension("lon", size=n_lons)         
                # creating variables                         
                var_mean = ds_new.createVariable(vn+'_diff_mean', "float32", ["lat", "lon"])
                var_mean[:, :] = var_mean_final
                var_sd = ds_new.createVariable(vn+'_diff_sd', "float32", ["lat", "lon"])
                var_sd[:, :] = var_sd_final            
                var_avd = ds_new.createVariable(vn+'_diff_avd', "float32", ["lat", "lon"])
                var_avd[:, :] = var_avd_final
                var_avd_relative = ds_new.createVariable(vn+'_diff_avd_relative', "float32", ["lat", "lon"])
                var_avd_relative[:, :] = var_avd_final / var_mean_final               
                lats = ds_new.createVariable("lat", "float32", ["lat"])
                lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))
                lons = ds_new.createVariable("lon", "float32", ["lon"])
                lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))               
                # closing NetCDF file      
                ds_new.close()                 
            k=0 # model counter 
        
    print('\033[1m'+'Done!')    

def grid_attribution (
        experiment_names,
        global_mask,    
        data_path,    
        ):
    g_mask=global_mask.index_mask
    for experiment in experiment_names:
        print('\033[1m'+experiment+'\033[0m')  
        # attribution of uncertainty in X to GPP, RT and X_p 
        file_path = Path(data_path).joinpath(experiment+"_X_c_uncertainty.nc")
        ds = nc.Dataset(str(file_path))        
        X_c_mean_array=ds.variables["X_c_mean"][:, :].data
        X_c_mean=np.ma.array(X_c_mean_array, mask=g_mask)                
        X_c_avd_array=ds.variables["X_c_avd"][:, :].data
        X_c_avd=np.ma.array(X_c_avd_array, mask=g_mask)
        ds.close()        
        file_path = Path(data_path).joinpath(experiment+"_gpp_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        gpp_avd_array=ds.variables["gpp_avd"][:, :].data/120
        gpp_avd=np.ma.array(gpp_avd_array, mask=g_mask)
        gpp_mean_array=ds.variables["gpp_mean"][:, :].data/120
        gpp_mean=np.ma.array(gpp_mean_array, mask=g_mask)
        ds.close() 
        file_path = Path(data_path).joinpath(experiment+"_RT_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        rt_avd_array=ds.variables["RT_avd"][:, :].data
        rt_avd=np.ma.array(rt_avd_array, mask=g_mask)
        rt_mean_array=ds.variables["RT_mean"][:, :].data  
        rt_mean=np.ma.array(rt_mean_array, mask=g_mask)
        ds.close() 
        file_path = Path(data_path).joinpath(experiment+"_X_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_avd_array=ds.variables["X_avd"][:, :].data
        X_avd=np.ma.array(X_avd_array, mask=g_mask)
        X_mean_array=ds.variables["X_mean"][:, :].data  
        X_mean=np.ma.array(X_mean_array, mask=g_mask)
        ds.close()
        file_path = Path(data_path).joinpath(experiment+"_X_p_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_p_avd_array=ds.variables["X_p_avd"][:, :].data
        X_p_avd=np.ma.array(X_p_avd_array, mask=g_mask)
        X_p_mean_array=ds.variables["X_p_mean"][:, :].data  
        X_p_mean=np.ma.array(X_p_mean_array, mask=g_mask)
        ds.close()
        
        gpp_contrib_ecosystem = np.ma.array(gpp_avd * rt_mean, mask=g_mask)
        rt_contrib = np.ma.array(rt_avd * gpp_mean, mask=g_mask)
        #gpp_rt_contrib = np.ma.array(X_c_avd - (gpp_contrib + rt_contrib), mask=g_mask)
        X_c=np.ma.array(X_c_avd, mask=g_mask)
        
        # print('RT_contrib: '+str(np.ma.mean(rt_contrib))) 
        # print('gpp_contrib: '+str(np.ma.mean(gpp_contrib_ecosystem)))         
        # print('X_c_avd: '+str(np.ma.mean(X_c_avd))) 
        # print('X_p_contrib: '+str(np.ma.mean(X_p_avd)))
       
        # print('gpp: '+str(np.ma.mean(gpp_mean)))   
        # print('RT: '+str(np.ma.mean(rt_mean))) 
        # print('X_c: '+str(np.ma.mean(X_c_mean))) 

        # print('X_c reconstr: '+str(np.ma.mean(rt_mean*gpp_mean)))
        # print('X_p: '+str(np.ma.mean(X_p_mean)))        
        # print('X: '+str(np.ma.mean(X_mean)))         
        # print('X_p reconstr: '+str(np.ma.mean(X_c_mean - X_mean)))   

        file_path = Path(data_path).joinpath(experiment+"_RT_veg_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        RT_veg_avd_array=ds.variables["RT_veg_avd"][:, :].data
        RT_veg_mean_array=ds.variables["RT_veg_mean"][:, :].data
        RT_veg_avd=np.ma.array(RT_veg_avd_array, mask=g_mask)
        RT_veg_mean=np.ma.array(RT_veg_mean_array, mask=g_mask)        
        ds.close()          
        file_path = Path(data_path).joinpath(experiment+"_RT_soil_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        RT_soil_avd_array=ds.variables["RT_soil_avd"][:, :].data
        RT_soil_mean_array=ds.variables["RT_soil_mean"][:, :].data
        RT_soil_avd=np.ma.array(RT_soil_avd_array, mask=g_mask)
        RT_soil_mean=np.ma.array(RT_soil_mean_array, mask=g_mask)         
        ds.close()          
        file_path = Path(data_path).joinpath(experiment+"_f_v2s_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        f_v2s_avd_array=ds.variables["f_v2s_avd"][:, :].data/120
        f_v2s_mean_array=ds.variables["f_v2s_mean"][:, :].data/120
        f_v2s_avd=np.ma.array(f_v2s_avd_array, mask=g_mask)
        f_v2s_mean=np.ma.array(f_v2s_mean_array, mask=g_mask)         
        ds.close()          
        
        gpp_contrib = np.ma.array(gpp_avd * RT_veg_mean, mask=g_mask)
        RT_veg_contrib = np.ma.array(RT_veg_avd * gpp_mean, mask=g_mask)
        X_c_veg_mean = np.ma.array(gpp_mean*RT_veg_mean, mask=g_mask)    
     
        
        f_v2s_contrib = np.ma.array(f_v2s_avd * RT_soil_mean, mask=g_mask)
        RT_soil_contrib = np.ma.array(RT_soil_avd * f_v2s_mean, mask=g_mask)
        X_c_soil_mean = np.ma.array(f_v2s_mean*RT_soil_mean, mask=g_mask)           


        file_path = Path(data_path).joinpath(experiment+"_cVeg_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        cVeg_avd_array=ds.variables["cVeg_avd"][:, :].data
        cVeg_mean_array=ds.variables["cVeg_mean"][:, :].data
        cVeg_avd=np.ma.array(cVeg_avd_array, mask=g_mask)
        cVeg_mean=np.ma.array(cVeg_mean_array, mask=g_mask)        
        ds.close()          
        file_path = Path(data_path).joinpath(experiment+"_cSoil_total_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        cSoil_avd_array=ds.variables["cSoil_total_avd"][:, :].data
        cSoil_mean_array=ds.variables["cSoil_total_mean"][:, :].data
        cSoil_avd=np.ma.array(cSoil_avd_array, mask=g_mask)
        cSoil_mean=np.ma.array(cSoil_mean_array, mask=g_mask)         
        ds.close() 
        file_path = Path(data_path).joinpath(experiment+"_X_c_veg_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_c_veg_avd_array=ds.variables["X_c_veg_avd"][:, :].data
        X_c_veg_mean_array=ds.variables["X_c_veg_mean"][:, :].data
        X_c_veg_avd=np.ma.array(X_c_veg_avd_array, mask=g_mask)
        X_c_veg_mean=np.ma.array(X_c_veg_mean_array, mask=g_mask)        
        ds.close()
        file_path = Path(data_path).joinpath(experiment+"_X_p_veg_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_p_veg_avd_array=ds.variables["X_p_veg_avd"][:, :].data
        X_p_veg_mean_array=ds.variables["X_p_veg_mean"][:, :].data
        X_p_veg_avd=np.ma.array(X_p_veg_avd_array, mask=g_mask)
        X_p_veg_mean=np.ma.array(X_p_veg_mean_array, mask=g_mask)        
        ds.close() 
        file_path = Path(data_path).joinpath(experiment+"_X_c_soil_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_c_soil_avd_array=ds.variables["X_c_soil_avd"][:, :].data
        X_c_soil_mean_array=ds.variables["X_c_soil_mean"][:, :].data
        X_c_soil_avd=np.ma.array(X_c_soil_avd_array, mask=g_mask)
        X_c_soil_mean=np.ma.array(X_c_soil_mean_array, mask=g_mask)        
        ds.close()
        file_path = Path(data_path).joinpath(experiment+"_X_p_soil_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        X_p_soil_avd_array=ds.variables["X_p_soil_avd"][:, :].data
        X_p_soil_mean_array=ds.variables["X_p_soil_mean"][:, :].data
        X_p_soil_avd=np.ma.array(X_p_soil_avd_array, mask=g_mask)
        X_p_soil_mean=np.ma.array(X_p_soil_mean_array, mask=g_mask)        
        ds.close()          
        
        #X_p_veg_contrib = np.ma.array(RT_veg_contrib+gpp_contrib-cVeg_avd)
        #X_p_veg_mean = np.ma.array(RT_veg_mean+gpp_mean-cVeg_mean)   
        
        #X_p_soil_contrib = np.ma.array(RT_soil_contrib+f_v2s_contrib-cSoil_avd)
        #X_p_soil_mean = np.ma.array(RT_soil_mean+f_v2s_mean-cSoil_mean)
                
        # print('veg:')
        # print('gpp_contrib: '+str(np.ma.mean(gpp_contrib)))         
        # print('RT_contrib: '+str(np.ma.mean(RT_veg_contrib))) 
        # print('cVeg_avd: '+str(np.ma.mean(cVeg_avd))) 
        # print('X_p_veg_contrib: '+str(np.ma.mean(X_p_veg_avd))) 
        # print('gpp_mean: '+str(np.ma.mean(gpp_mean)))         
        # print('RT_mean: '+str(np.ma.mean(RT_veg_mean)))         
        # print('X_c_mean: '+str(np.ma.mean(X_c_veg_mean))) 
        # print('cVeg_mean:'+str(np.ma.mean(cVeg_mean)))  
        # print('X_p_veg_mean: '+str(np.ma.mean(X_p_veg_mean)))         
            
        # print('soil:')
        # print('f_v2s_contrib: '+str(np.ma.mean(f_v2s_contrib)))   
        # print('RT_contrib: '+str(np.ma.mean(RT_soil_contrib))) 
        # print('cSoil_avd: '+str(np.ma.mean(cSoil_avd)))         
        # print('X_p_soil_contrib: '+str(np.ma.mean(X_p_soil_avd)))         
        # print('f_v2s_mean: '+str(np.ma.mean(f_v2s_mean)))   
        # print('RT_mean: '+str(np.ma.mean(RT_soil_mean)))         
        # print('X_c_mean: '+str(np.ma.mean(X_c_soil_mean)))  
        # print('cSoil_mean:'+str(np.ma.mean(cSoil_mean))) 
        # print('X_p_soil_mean: '+str(np.ma.mean(X_p_soil_mean)))  
        
        
        gcm=global_mask

        print("computing means, this may take some minutes...")
        
        # data_str = namedtuple(
        # 'data_str',
        # ["cVeg", "cSoil_total","cVeg_diff", "cSoil_total_diff","nep", "gpp", "ra", "rh",  "dist", "f_v2s", 
        # "X", "X_c", "X_p","RT", "RT_veg", "RT_soil"]
        # )  
                
        def global_mean_var(
            lats: np.ma.core.MaskedArray,
            lons: np.ma.core.MaskedArray,
            #mask: np.array,
            var: nc._netCDF4.Variable,
        ) -> np.array:
            #lats=np.ma.array(gcm.lats),
            #lons=np.ma.array(gcm.lons),
            mask=gcm.index_mask,
            weight_mat = np.ma.array(gh.get_weight_mat(lats, lons), mask=mask)
            wms = weight_mat.sum()
            res = (weight_mat * var[:, :]).sum() / wms #np.zeros(n_t)
            return res
        
        
        
        rt_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=rt_contrib)
        print("rt_contrib_global: "+str(rt_contrib_global))        
        gpp_contrib_eco_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=gpp_contrib_ecosystem)
        print("gpp_contrib_eco_global: "+str(gpp_contrib_eco_global)) 
        X_p_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=X_p_avd)
        print("X_p_contrib_global: "+str(X_p_contrib_global))         


        cVeg_avd_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=cVeg_avd)
        print("cVeg_avd_global: "+str(cVeg_avd_global)) 

        cSoil_avd_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=cSoil_avd)
        print("cSoil_total_avd_global: "+str(cSoil_avd_global)) 

        rt_veg_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=RT_veg_contrib)
        print("rt_veg_contrib_global: "+str(rt_veg_contrib_global))        
        gpp_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=gpp_contrib)
        print("gpp_contrib_global: "+str(gpp_contrib_global)) 
        X_p_veg_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=X_p_veg_avd)
        print("X_p_veg_contrib_global: "+str(X_p_veg_contrib_global))

        rt_soil_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=RT_soil_contrib)
        print("rt_soil_contrib_global: "+str(rt_soil_contrib_global))        
        f_v2s_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=f_v2s_contrib)
        print("f_v2s_contrib_global: "+str(f_v2s_contrib_global)) 
        X_p_soil_contrib_global=global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=X_p_soil_avd)
        print("rt_soil_contrib_global: "+str(X_p_soil_contrib_global))        

        total_eco_grid=rt_contrib+gpp_contrib_ecosystem+X_p_avd
        percent_rt_grid=rt_contrib/total_eco_grid
        percent_gpp_eco_grid=gpp_contrib_ecosystem/total_eco_grid
        percent_X_p_grid=X_p_avd/total_eco_grid   
        
        total_veg_soil_grid=cVeg_avd+cSoil_avd
        percent_veg_grid=cVeg_avd/total_veg_soil_grid
        percent_soil_grid=cSoil_avd/total_veg_soil_grid

        percent_rt = global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=percent_rt_grid)
        percent_gpp_eco = global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=percent_gpp_eco_grid)
        percent_X_p = global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,               
                var=percent_X_p_grid)
                
        percent_veg = global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,
                var=percent_veg_grid)
        percent_soil = global_mean_var(
                lats=gcm.lats,
                lons=gcm.lons,
                #mask=gcm.index_mask,                                
                var=percent_soil_grid)                
        #total_eco=rt_contrib_global+gpp_contrib_eco_global+X_p_contrib_global
        #percent_rt=rt_contrib_global/total_eco
        #percent_gpp_eco=gpp_contrib_eco_global/total_eco
        #percent_X_p=X_p_contrib_global/total_eco        
               
        total_C=cVeg_avd_global+cSoil_avd_global
        cVeg_avd_part=cVeg_avd_global/total_C
        cSoil_avd_part=cSoil_avd_global/total_C
        print('total_C: '+str(total_C))
        print('cVeg_avd_part: '+str(cVeg_avd_part))
        print('cSoil_avd_part: '+str(cSoil_avd_part))
        
        # total_veg_soil=(rt_veg_contrib_global+gpp_contrib_global+X_p_veg_contrib_global+
                            # rt_soil_contrib_global+f_v2s_contrib_global+X_p_soil_contrib_global)        
        
        total_veg=rt_veg_contrib_global+gpp_contrib_global+X_p_veg_contrib_global
        total_soil=rt_soil_contrib_global+f_v2s_contrib_global+X_p_soil_contrib_global
        
        percent_rt_veg=rt_veg_contrib_global/total_veg*100*cVeg_avd_part
        percent_gpp=gpp_contrib_global/total_veg*100*cVeg_avd_part
        percent_X_p_veg=X_p_veg_contrib_global/total_veg*100*cVeg_avd_part      
        percent_rt_soil=rt_soil_contrib_global/total_soil*100*cSoil_avd_part
        percent_f_v2s=f_v2s_contrib_global/total_soil*100*cSoil_avd_part
        percent_X_p_soil=X_p_soil_contrib_global/total_soil*100*cSoil_avd_part       
        
        # pie charts
        
        fig=plt.figure(figsize=(7,7))
        ax1=fig.subplots()#axs[1]
        ax1.set_title('Global average % attribution of uncertainty in C storage')
                        
        #labels = '$\Delta$ GPP', '$\Delta$ RA', '$\Delta$NPP*$\Delta$ RH'#, '$\Delta$ Dist'
        labels = 'GPP', 'RT', 'X_p'#, '$\Delta$ Dist'
        sizes = [percent_gpp_eco, percent_rt, percent_X_p]
        ax1.pie(sizes, autopct='%1.1f%%', 
            startangle=90, counterclock=False, colors=("green", "darkorange", "blue"))#, "purple"))
   
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
        ax1.legend(labels, bbox_to_anchor =(1,1))
    
        plt.show() 

        # pie charts
        
        fig=plt.figure(figsize=(7,7))
        ax1=fig.subplots()#axs[1]
        ax1.set_title('Global average % attribution of uncertainty in C storage: veg + soil')
                        
        #labels = '$\Delta$ GPP', '$\Delta$ RA', '$\Delta$NPP*$\Delta$ RH'#, '$\Delta$ Dist'
        labels = 'GPP', 'RT_veg', 'X_p_veg', 'f_v2s', 'RT_soil', 'X_p_soil'#, '$\Delta$ Dist'
        sizes = [percent_gpp, percent_rt_veg, percent_X_p_veg, percent_f_v2s, percent_rt_soil, percent_X_p_soil ]
        ax1.pie(sizes, autopct='%1.1f%%', 
            startangle=90, counterclock=False, colors=("green", "darkorange", "blue", "cyan", "brown", "purple"))#, "purple"))
   
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
        plt.rcParams.update({'font.size': 20})
        ax1.legend(labels, bbox_to_anchor =(1,1))
    
        plt.show()
        
        # RGB map of attribution to GPP, RT and X_p
        rgba_map (
            mask = gcm.index_mask,
            grid_red = percent_rt_grid, 
            grid_green = percent_gpp_eco_grid, 
            grid_blue = percent_X_p_grid, 
            grid_alpha = X_avd / np.ma.max(X_avd)*10, 
            labels = ["GPP 80% - RT 20%",
                    #"GPP 75% - RT 25%",
                    "GPP 70% - RT 30%",
                    #"GPP 65% - RT 35%",
                    "GPP 60% - RT 40%",
                    #"GPP 55% - RT 45%",
                    "GPP 50% - RT 50%",
                    #"GPP 45% - RT 55%",
                    "GPP 40% - RT 60%",
                    #"GPP 35% - RT 65%",
                    "GPP 30% - RT 70%",
                    #"GPP 35% - RT 75%",
                    "GPP 20% - RT 80%"],
            title = "Spatial attribution of uncetainty in X to uncertainty in GPP (green) and RT (red)",         
        )
        
        # creating and writing a new NetCDF file 
        s = g_mask.shape
        n_lats, n_lons = s
        new_path=Path(data_path).joinpath("gpp_vs_RT_ecosystem.nc")
        ds_new = nc.Dataset(str(new_path), "w", persist=True)
        # creating dimensions
        lat = ds_new.createDimension("lat", size=n_lats)
        lon = ds_new.createDimension("lon", size=n_lons)         
        # creating variables                         
        var = ds_new.createVariable('percent_rt_uncertainty', "float32", ["lat", "lon"])
        var[:, :] = percent_rt_grid 
        var2 = ds_new.createVariable('percent_gpp_uncertainty', "float32", ["lat", "lon"])
        var2[:, :] = percent_gpp_eco_grid 
        var3 = ds_new.createVariable('percent_X_p_uncertainty', "float32", ["lat", "lon"])
        var3[:, :] = percent_X_p_grid
        var4 = ds_new.createVariable('X_uncertainty', "float32", ["lat", "lon"])
        var4[:, :] = X_avd                
        lats = ds_new.createVariable("lat", "float32", ["lat"])
        lats[:] = list(map(global_mask.tr.i2lat, range(n_lats)))      
        lons = ds_new.createVariable("lon", "float32", ["lon"])
        lons[:] = list(map(global_mask.tr.i2lon, range(n_lons)))
       
        # closing NetCDF file      
        ds_new.close()    
                
        file_path = Path(data_path).joinpath(experiment+"_ra_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        ra_avd=ds.variables["ra_avd"][:, :].data
        ds.close()         
        file_path = Path(data_path).joinpath(experiment+"_rh_uncertainty.nc")
        ds = nc.Dataset(str(file_path))
        rh_avd=ds.variables["rh_avd"][:, :].data
        ds.close() 

def rgba_map (mask, grid_red, grid_green, grid_blue, grid_alpha, labels, title = "RGB map"):

        #height = len(gcm.lats) # can be any value
        #width = len(gcm.lons) # can be any value

        height = grid_red.shape[0]
        width = grid_red.shape[1]

        array = np.zeros((height, width, 4))
        array[:,:,0] = grid_red # red
        array[:,:,1] = grid_green # green
        array[:,:,2] = grid_blue # blue
        array[:,:,3] = grid_alpha # alpha
        #array[:,:,3] = np.zeros_like(X_p_2d) # blue
        
        rgb_map = array.copy()
        #rgb_map[rgb_map>1]=0
        
        for i in range(height): 
            rgb_map[i,:,:] = array[height-i-1,:,:]
        mask_map = mask.copy()
        for i in range(height): 
            mask_map[i,:] = mask[height-i-1,:]   
        
        mask_map=(mask_map-1)*(-1)
        mask_4d_map=np.zeros_like(array)
        mask_4d_map[:,:,0] = mask_map # red
        mask_4d_map[:,:,1] = mask_map # green
        mask_4d_map[:,:,2] = mask_map # blue
        mask_4d_map[:,:,3] = mask_map # alpha
        fig=plt.figure(figsize=(200,100))
        #plt.rcParams.update({'font.size': 150})
        axs=fig.subplots(1,2, gridspec_kw={'width_ratios': [15, 1]})        
        ax=axs[0]
        ax.axis("off")
        ax.imshow(rgb_map*mask_4d_map)
        
        ax.set_title(title, fontsize=160)
        ax_legend=axs[1]
        green=np.arange(0.8,0.2,-0.05)
        red=np.arange(0.2,0.8,0.05)        
        blue=np.zeros(13)
        #alpha=np.arange(0.2,0.8,0.05)
        legend = np.zeros((13, 1, 3))
        legend [:,:,0]=red.reshape(13,1)
        legend [:,:,1]=green.reshape(13,1)
        legend [:,:,2]=blue.reshape(13,1)
        #plt.rcParams.update({'font.size': 8})
        ax_legend.set_title(" ", fontsize=140)
        #ax3.axis("off")
        ax_legend.get_xaxis().set_visible(False)
        ticks=np.array(range(0,13,round(13/len(labels))))
        ax_legend.set_yticks (ticks)
        ax_legend.set_yticklabels(labels,fontsize=120)

        ax_legend.imshow(legend)
        plt.show()
        # save RGB pdf
        fig.savefig("rgb_map.pdf")       
        plt.rcParams.update({'font.size': 15})
        return ()
       
      
def get_global_mean_uncertainty(dataPath,  
                            experiment_name, # string, e.g. "S2"
                            data_str, # named tuple
                            avd,
                            #lat_var,
                            #lon_var,
                            ):
    gcm=gh.globalMask(file_name="common_mask_all_models.nc")
    print("computing means, this may take some minutes...")

    def global_mean_var(
        lats: np.ma.core.MaskedArray,
        lons: np.ma.core.MaskedArray,
        mask: np.array,
        var: nc._netCDF4.Variable,
    ) -> np.array:
        """As the signature shows this function expects a netCDF4.Variable
        This is basically metadata which allows us to compute the maean even
        if the whole array would not fit into memory.

        ds = nc.Dataset("example.nc")
        var=ds.variables['npp'] #->type(var)=netCDF4._netCDF4.Variable

        the mask array is used to block out extra pixels that are not
        masked in var
        """
        weight_mat = np.ma.array(gh.get_weight_mat(lats, lons), mask=mask)

        # to compute the sum of weights we add only those weights that
        # do not correspond to an unmasked grid cell
        wms = weight_mat.sum()
        res = (weight_mat * var[:, :]).sum() / wms 
        
        return res
    
    if avd: suffix = '_avd' 
    else: suffix = '_mean'
    
    def compute_and_cache_global_mean(vn):
        path = Path(dataPath).joinpath(experiment_name+'_'+vn+'_uncertainty.nc')
        #if vn=="npp_nlim": path=dataPath.joinpath(gh.msh(model_folder).nc_file_name("npp", experiment_name=experiment_name))
        print(path)
        ds = nc.Dataset(str(path))
        vs=ds.variables
        lats= vs['lat'].__array__()
        lons= vs['lon'].__array__()
        var=ds.variables[vn+suffix]
        # check if we have a cached version (which is much faster)
        #gm_path = dataPath.joinpath(nc_global_mean_file_name(experiment_name=experiment_name))
        
        #model_mask = gh.msh(model_folder).spatial_mask(dataPath=Path(conf_dict["dataPath"]))
        #combined_mask = combine_masks ([model_mask,gcm])
        gm=global_mean_var(
                lats,
                lons,
                #combined_mask.index_mask,
                gcm.index_mask,
                var      
        )     
        return gm #* 86400 if vn in ["gpp", "npp", "npp_nlim", "rh", "ra"] else gm
        
    output=data_str(*map(compute_and_cache_global_mean, data_str._fields)) 
    return (output) 


# def plot_attribution_output (
    # all_comp_dict,
    # percent=False,
    # part=1,
# ):
    # models=list(all_comp_dict.keys())[:-2]  

    # # if we do not want the whole interval but look at a smaller part to observe the dynamics
    # start_min=0
    # stop_max=len(all_comp_dict["Times"])
    
    # # if we do not want the whole interval but look at a smaller part to observe the dynamics
    # if part < 0:
        # start, stop = int(stop_max - (stop_max - start_min) * abs(part)), stop_max
    # else:
        # start, stop = start_min, int(start_min + (stop_max - start_min) * part)
    
    # times=all_comp_dict["Times"][start:stop]   
    
    # sum_pos_diff_x=0
    # sum_pos_cont_rt=0
    # sum_pos_cont_u=0
    # sum_pos_cont_rt_u_inter=0
    # sum_pos_cont_x_p=0
    
    # sum_neg_diff_x=0
    # sum_neg_cont_rt=0
    # sum_neg_cont_u=0
    # sum_neg_cont_rt_u_inter=0
    # sum_neg_cont_x_p=0  
    
    # print ('\033[1m'+'Attribution of summed deviations from the mean for all models ' +
        # 'to the differences in traceable components')     
    # for m in models: 
           
        # x=all_comp_dict[m]["x"][start:stop]
        # x_c=all_comp_dict[m]["x_c"][start:stop]
        # x_p=all_comp_dict[m]["x_p"][start:stop]
        # u=all_comp_dict[m]["gpp"][start:stop]
        # rt=all_comp_dict[m]["rt"][start:stop]
        # x_mean=all_comp_dict["Mean"]["x"][start:stop]
        # x_c_mean=all_comp_dict["Mean"]["x_c"][start:stop]
        # x_p_mean=all_comp_dict["Mean"]["x_p"][start:stop]
        # u_mean=all_comp_dict["Mean"]["gpp"][start:stop]
        # rt_mean=all_comp_dict["Mean"]["rt"][start:stop]           
            
        # delta_x=x-x_mean
        # delta_x_c=x_c-x_c_mean
        # delta_x_p=x_p-x_p_mean
        # delta_u=u-u_mean
        # delta_rt=rt-rt_mean
        
        # # attribution of delta X to delta X_c and delta X_p
        # x_c_contrib=delta_x_c
        # x_p_contrib=-delta_x_p
         
        # # attribution of delta X_c to delta u and delta RT
        # rt_contrib=delta_rt*(u-delta_u/2)
        # u_contrib=delta_u*(rt-delta_rt/2) 
        # rt_u_inter=delta_x_c-rt_contrib-u_contrib         

        # # summation of positive and negative contributions separately         
        
        # pos_delta_x=delta_x.copy(); pos_delta_x[pos_delta_x<0]=0
        # neg_delta_x=delta_x.copy(); neg_delta_x[neg_delta_x>0]=0
        # sum_pos_diff_x+=pos_delta_x
        # sum_neg_diff_x+=neg_delta_x       

        # pos_cont_rt=rt_contrib.copy(); pos_cont_rt[pos_cont_rt<0]=0
        # neg_cont_rt=rt_contrib.copy(); neg_cont_rt[neg_cont_rt>0]=0
        # sum_pos_cont_rt+=pos_cont_rt
        # sum_neg_cont_rt+=neg_cont_rt

        # pos_cont_u=u_contrib.copy(); pos_cont_u[pos_cont_u<0]=0
        # neg_cont_u=u_contrib.copy(); neg_cont_u[neg_cont_u>0]=0
        # sum_pos_cont_u+=pos_cont_u
        # sum_neg_cont_u+=neg_cont_u

        # pos_cont_rt_u_inter=rt_u_inter.copy(); pos_cont_rt_u_inter[pos_cont_rt_u_inter<0]=0
        # neg_cont_rt_u_inter=rt_u_inter.copy(); neg_cont_rt_u_inter[neg_cont_rt_u_inter>0]=0
        # sum_pos_cont_rt_u_inter+=pos_cont_rt_u_inter
        # sum_neg_cont_rt_u_inter+=neg_cont_rt_u_inter

        # pos_cont_x_p=x_p_contrib.copy(); pos_cont_x_p[pos_cont_x_p<0]=0
        # neg_cont_x_p=x_p_contrib.copy(); neg_cont_x_p[neg_cont_x_p>0]=0
        # sum_pos_cont_x_p+=pos_cont_x_p
        # sum_neg_cont_x_p+=neg_cont_x_p 
        
    # plot_attribution (
        # times=times,
        # delta_x_pos=sum_pos_diff_x,
        # delta_x_neg=sum_neg_diff_x,
        # rt_contrib_pos=sum_pos_cont_rt,
        # rt_contrib_neg=sum_neg_cont_rt,
        # u_contrib_pos=sum_pos_cont_u,
        # u_contrib_neg=sum_neg_cont_u,
        # rt_u_inter_pos=sum_pos_cont_rt_u_inter,
        # rt_u_inter_neg=sum_neg_cont_rt_u_inter,
        # x_p_contrib_pos=sum_pos_cont_x_p,
        # x_p_contrib_neg=sum_neg_cont_x_p,
        # percent=percent,        
    # )    
    
    # print ('\033[1m'+'Attribution of summed deviations from the mean for all models ' +
        # 'to the differences in fluxes')     

    # sum_pos_diff_nep=0
    # sum_neg_diff_nep=0
    # sum_gpp_pos=0
    # sum_gpp_neg=0        
    # sum_rh_pos=0
    # sum_rh_neg=0
    # sum_ra_pos=0
    # sum_ra_neg=0
    # sum_dist_pos=0
    # sum_dist_neg=0     
    # sum_cVeg_contrib_pos=0
    # sum_cVeg_contrib_neg=0
    # sum_ra_rate_contrib_pos=0
    # sum_ra_rate_contrib_neg=0
    # sum_ra_cVeg_inter_pos=0
    # sum_ra_cVeg_inter_neg=0      
    # sum_cSoil_contrib_pos=0
    # sum_cSoil_contrib_neg=0   
    # sum_rh_rate_contrib_pos=0
    # sum_rh_rate_contrib_neg=0   
    # sum_rh_cSoil_inter_pos=0
    # sum_rh_cSoil_inter_neg=0        
    # sum_f_v2s_pos=0
    # sum_f_v2s_neg=0    

    # for m in models: 
           
        # nep=all_comp_dict[m]["nep"][start:stop]
        # gpp=all_comp_dict[m]["gpp"][start:stop]        
        # ra=all_comp_dict[m]["ra"][start:stop]
        # rh=all_comp_dict[m]["rh"][start:stop]
        # dist=all_comp_dict[m]["dist"][start:stop]
        # cVeg=all_comp_dict[m]["cVeg"][start:stop]
        # cSoil=all_comp_dict[m]["cSoil"][start:stop]
        # ra_rate=all_comp_dict[m]["ra_rate"][start:stop]
        # rh_rate=all_comp_dict[m]["rh_rate"][start:stop] 
        # f_v2s=all_comp_dict[m]["f_v2s"][start:stop] 
         
        # nep_mean=all_comp_dict["Mean"]["nep"][start:stop]
        # gpp_mean=all_comp_dict["Mean"]["gpp"][start:stop]        
        # ra_mean=all_comp_dict["Mean"]["ra"][start:stop]
        # rh_mean=all_comp_dict["Mean"]["rh"][start:stop]
        # dist_mean=all_comp_dict["Mean"]["dist"][start:stop]        
        # cVeg_mean=all_comp_dict["Mean"]["cVeg"][start:stop]
        # cSoil_mean=all_comp_dict["Mean"]["cSoil"][start:stop]
        # ra_rate_mean=all_comp_dict["Mean"]["ra_rate"][start:stop]
        # rh_rate_mean=all_comp_dict["Mean"]["rh_rate"][start:stop] 
        # f_v2s_mean=all_comp_dict["Mean"]["f_v2s"][start:stop]         
        
        # delta_nep=nep-nep_mean
        # delta_gpp=gpp-gpp_mean
        # delta_ra_contrib=ra-ra_mean
        # delta_rh_contrib=rh-rh_mean
        # delta_dist_contrib=dist-dist_mean        
        # delta_cVeg=cVeg-cVeg_mean
        # delta_cSoil=cSoil-cSoil_mean
        # delta_ra_rate=ra_rate-ra_rate_mean
        # delta_rh_rate=rh_rate-rh_rate_mean
        # delta_f_v2s_contrib=f_v2s-f_v2s_mean
        
        # # attribution of delta nep to delta gpp, ra, rh, disturbance
        # #x_c_contrib=delta_x_c
        # #x_p_contrib=-delta_x_p
         
        # # attribution of product of pools and rates
        # delta_cVeg_contrib=delta_cVeg*(ra_rate-delta_ra_rate/2)
        # delta_ra_rate_contrib=delta_ra_rate*(cVeg-delta_cVeg/2) 
        # delta_ra_cVeg_inter=delta_ra_contrib-delta_ra_rate_contrib-delta_cVeg_contrib         

        # delta_cSoil_contrib=delta_cSoil*(rh_rate-delta_rh_rate/2)
        # delta_rh_rate_contrib=delta_rh_rate*(cSoil-delta_cSoil/2) 
        # delta_rh_cSoil_inter=delta_rh_contrib-delta_rh_rate_contrib-delta_cSoil_contrib  

        # # summation of positive and negative contributions separately         
        
        # pos_delta_nep=delta_nep.copy(); pos_delta_nep[pos_delta_nep<0]=0
        # neg_delta_nep=delta_nep.copy(); neg_delta_nep[neg_delta_nep>0]=0
        # sum_pos_diff_nep+=pos_delta_nep
        # sum_neg_diff_nep+=neg_delta_nep     

        # pos_delta_gpp=delta_gpp.copy(); pos_delta_gpp[pos_delta_gpp<0]=0
        # neg_delta_gpp=delta_gpp.copy(); neg_delta_gpp[neg_delta_gpp>0]=0
        # sum_gpp_pos+=pos_delta_gpp
        # sum_gpp_neg+=neg_delta_gpp  
        
        # pos_delta_ra=delta_ra_contrib.copy(); pos_delta_ra[pos_delta_ra<0]=0
        # neg_delta_ra=delta_ra_contrib.copy(); neg_delta_ra[neg_delta_ra>0]=0
        # sum_ra_pos+=pos_delta_ra
        # sum_ra_neg+=neg_delta_ra    

        # pos_delta_rh=delta_rh_contrib.copy(); pos_delta_rh[pos_delta_rh<0]=0
        # neg_delta_rh=delta_rh_contrib.copy(); neg_delta_rh[neg_delta_rh>0]=0
        # sum_rh_pos+=pos_delta_rh
        # sum_rh_neg+=neg_delta_rh    

        # pos_delta_dist=delta_dist_contrib.copy(); pos_delta_dist[pos_delta_dist<0]=0
        # neg_delta_dist=delta_dist_contrib.copy(); neg_delta_dist[neg_delta_dist>0]=0
        # sum_dist_pos+=pos_delta_dist
        # sum_dist_neg+=neg_delta_dist    

        # pos_delta_cVeg_contrib=delta_cVeg_contrib.copy(); pos_delta_cVeg_contrib[pos_delta_cVeg_contrib<0]=0
        # neg_delta_cVeg_contrib=delta_cVeg_contrib.copy(); neg_delta_cVeg_contrib[neg_delta_cVeg_contrib>0]=0
        # sum_cVeg_contrib_pos+=pos_delta_cVeg_contrib
        # sum_cVeg_contrib_neg+=neg_delta_cVeg_contrib 

        # pos_delta_cSoil_contrib=delta_cSoil_contrib.copy(); pos_delta_cSoil_contrib[pos_delta_cSoil_contrib<0]=0
        # neg_delta_cSoil_contrib=delta_cSoil_contrib.copy(); neg_delta_cSoil_contrib[neg_delta_cSoil_contrib>0]=0
        # sum_cSoil_contrib_pos+=pos_delta_cSoil_contrib
        # sum_cSoil_contrib_neg+=neg_delta_cSoil_contrib 
        
        # pos_delta_ra_rate_contrib=delta_ra_rate_contrib.copy(); pos_delta_ra_rate_contrib[pos_delta_ra_rate_contrib<0]=0
        # neg_delta_ra_rate_contrib=delta_ra_rate_contrib.copy(); neg_delta_ra_rate_contrib[neg_delta_ra_rate_contrib>0]=0
        # sum_ra_rate_contrib_pos+=pos_delta_ra_rate_contrib
        # sum_ra_rate_contrib_neg+=neg_delta_ra_rate_contrib  
        
        # pos_delta_ra_cVeg_inter=delta_ra_cVeg_inter.copy(); pos_delta_ra_cVeg_inter[pos_delta_ra_cVeg_inter<0]=0
        # neg_delta_ra_cVeg_inter=delta_ra_cVeg_inter.copy(); neg_delta_ra_cVeg_inter[neg_delta_ra_cVeg_inter>0]=0
        # sum_ra_cVeg_inter_pos+=pos_delta_ra_cVeg_inter
        # sum_ra_cVeg_inter_neg+=neg_delta_ra_cVeg_inter          

        # pos_delta_rh_rate_contrib=delta_rh_rate_contrib.copy(); pos_delta_rh_rate_contrib[pos_delta_rh_rate_contrib<0]=0
        # neg_delta_rh_rate_contrib=delta_rh_rate_contrib.copy(); neg_delta_rh_rate_contrib[neg_delta_rh_rate_contrib>0]=0
        # sum_rh_rate_contrib_pos+=pos_delta_rh_rate_contrib
        # sum_rh_rate_contrib_neg+=neg_delta_rh_rate_contrib    

        # pos_delta_rh_cSoil_inter=delta_rh_cSoil_inter.copy(); pos_delta_rh_cSoil_inter[pos_delta_rh_cSoil_inter<0]=0
        # neg_delta_rh_cSoil_inter=delta_rh_cSoil_inter.copy(); neg_delta_rh_cSoil_inter[neg_delta_rh_cSoil_inter>0]=0
        # sum_rh_cSoil_inter_pos+=pos_delta_rh_cSoil_inter
        # sum_rh_cSoil_inter_neg+=neg_delta_rh_cSoil_inter         

        # pos_delta_f_v2s=delta_f_v2s_contrib.copy(); pos_delta_f_v2s[pos_delta_f_v2s<0]=0
        # neg_delta_f_v2s=delta_f_v2s_contrib.copy(); neg_delta_f_v2s[neg_delta_f_v2s>0]=0
        # sum_f_v2s_pos+=pos_delta_f_v2s
        # sum_f_v2s_neg+=neg_delta_f_v2s         
        
        
    # plot_attribution_2 (
        # times=times,
        # delta_nep_pos=sum_pos_diff_nep,
        # delta_nep_neg=sum_neg_diff_nep,
        # gpp_pos=sum_gpp_pos,
        # gpp_neg=sum_gpp_neg,        
        # rh_pos=sum_rh_pos,
        # rh_neg=sum_rh_neg,
        # ra_pos=sum_ra_pos,
        # ra_neg=sum_ra_neg,
        # dist_pos=sum_dist_pos,
        # dist_neg=sum_dist_neg,     
        # cVeg_contrib_pos=sum_cVeg_contrib_pos,
        # cVeg_contrib_neg=sum_cVeg_contrib_neg,
        # ra_rate_contrib_pos=sum_ra_rate_contrib_pos,
        # ra_rate_contrib_neg=sum_ra_rate_contrib_neg,
        # ra_cVeg_inter_pos=sum_ra_cVeg_inter_pos,
        # ra_cVeg_inter_neg=sum_ra_cVeg_inter_neg,      
        # cSoil_contrib_pos=sum_cSoil_contrib_pos,
        # cSoil_contrib_neg=sum_cSoil_contrib_neg,   
        # rh_rate_contrib_pos=sum_rh_rate_contrib_pos,
        # rh_rate_contrib_neg=sum_rh_rate_contrib_neg,   
        # rh_cSoil_inter_pos=sum_rh_cSoil_inter_pos,
        # rh_cSoil_inter_neg=sum_rh_cSoil_inter_neg,        
        # f_v2s_pos=sum_f_v2s_pos,
        # f_v2s_neg=sum_f_v2s_neg,        
        # percent=percent,        
    # )      
    
    
# def plot_attribution_2 (
        # times,
        # delta_nep_pos,
        # delta_nep_neg,
        # gpp_pos,
        # gpp_neg,        
        # rh_pos,
        # rh_neg,
        # ra_pos,
        # ra_neg,
        # dist_pos,
        # dist_neg,     
        # cVeg_contrib_pos,
        # cVeg_contrib_neg,
        # ra_rate_contrib_pos,
        # ra_rate_contrib_neg,
        # ra_cVeg_inter_pos,
        # ra_cVeg_inter_neg,      
        # cSoil_contrib_pos,
        # cSoil_contrib_neg,   
        # rh_rate_contrib_pos,
        # rh_rate_contrib_neg,   
        # rh_cSoil_inter_pos,
        # rh_cSoil_inter_neg, 
        # f_v2s_pos,
        # f_v2s_neg,
        # percent,  
# ):
        
        # gpp_contrib=gpp_pos-gpp_neg 
        # rh_contrib=rh_pos-rh_neg
        # ra_contrib=ra_pos-ra_neg
        # dist_contrib=dist_pos-dist_neg
        # cVeg_contrib=cVeg_contrib_pos-cVeg_contrib_neg
        # ra_rate_contrib=ra_rate_contrib_pos-ra_rate_contrib_neg
        # cSoil_contrib=cSoil_contrib_pos-cSoil_contrib_neg
        # rh_rate_contrib=rh_rate_contrib_pos-rh_rate_contrib_neg
        # ra_cVeg_inter_contrib=ra_cVeg_inter_pos-ra_cVeg_inter_neg
        # rh_cSoil_inter_contrib=rh_cSoil_inter_pos-rh_cSoil_inter_neg
        # f_v2s_contrib=f_v2s_pos-f_v2s_neg        
       
        # fig5=plt.figure(figsize=(15,5))
        # axs=fig5.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]})
        
        # # contributions timeline
        # ax=axs[0]

        # # quadratic trends
            

        # z = np.polyfit(times,  gpp_pos, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # times, 
            # p(times),
            # label="$\Delta$ GPP positive",
            # color="green",
        # ) 
        # ax.fill_between(
                # times,
                # gpp_pos, 
                # p(times),
                # color="green",
                # alpha=0.2                
                # )         

        # z = np.polyfit(times,  gpp_neg, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # times, 
            # p(times),
            # label="$\Delta$ GPP negative",
            # color="green",
            # linestyle="dashed",            
        # ) 
        # ax.fill_between(
                # times,
                # gpp_neg, 
                # p(times),
                # color="green",
                # #linewidth=0.1,
                # alpha=0.2                
                # )          
   
        # z = np.polyfit(times,  ra_pos, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # times, 
            # p(times),
            # label="$\Delta$ RA positive",
            # color="blue",
        # ) 
        # ax.fill_between(
                # times,
                # ra_pos, 
                # p(times),
                # color="blue",
                # alpha=0.2                
                # )
   
        # z = np.polyfit(times,  ra_neg, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # times, 
            # p(times),
            # label="$\Delta$ RA negative",
            # color="blue",
            # linestyle="dashed",            
        # ) 
        # ax.fill_between(
                # times,
                # ra_neg, 
                # p(times),
                # color="blue",
                # alpha=0.2                
                # )   
   
        # z = np.polyfit(times,  rh_pos, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # times, 
            # p(times),
            # label="$\Delta$ rh positive",
            # color="red",
        # ) 
        # ax.fill_between(
                # times,
                # rh_pos, 
                # p(times),
                # color="red",
                # alpha=0.2                
                # )
   
        # z = np.polyfit(times,  rh_neg, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # times, 
            # p(times),
            # label="$\Delta$ rh negative",
            # color="red",
            # linestyle="dashed",            
        # ) 
        # ax.fill_between(
                # times,
                # rh_neg, 
                # p(times),
                # color="red",
                # alpha=0.2                
                # )  
                    
        # # if abs(np.mean(dist_pos))>0.01:        
            # # z = np.polyfit(times,  dist_pos, 2)
            # # p = np.poly1d(z)  
            # # ax.plot(        
                # # times, 
                # # p(times),
                # # label="$\Delta$ dist positive",
                # # color="purple",
            # # ) 
            # # ax.fill_between(
                    # # times,
                    # # dist_pos, 
                    # # p(times),
                    # # color="purple",
                    # # alpha=0.2                
                    # # )
        # # if abs(np.mean(dist_neg))>0.01:        
            # # z = np.polyfit(times,  dist_neg, 2)
            # # p = np.poly1d(z)  
            # # ax.plot(        
                # # times, 
                # # p(times),
                # # label="$\Delta$ dist negative",
                # # color="purple",
                # # linestyle="dashed",            
            # # ) 
            # # ax.fill_between(
                    # # times,
                    # # dist_neg, 
                    # # p(times),
                    # # color="purple",
                    # # alpha=0.2                
                    # # )                 

        # # bar charts       

        # positive_delta_nep=sum(delta_nep_pos)/len(delta_nep_pos)
        # negative_delta_nep=sum(delta_nep_neg)/len(delta_nep_neg)
        # positive_gpp_contrib=sum(gpp_pos)/len(gpp_pos)
        # negative_gpp_contrib=sum(gpp_neg)/len(gpp_neg)        
        # positive_ra_contrib=sum(ra_pos)/len(ra_pos)
        # negative_ra_contrib=sum(ra_neg)/len(ra_neg)
        # positive_rh_contrib=sum(rh_pos)/len(rh_pos)
        # negative_rh_contrib=sum(rh_neg)/len(rh_neg)
        # positive_dist_contrib=sum(dist_pos)/len(dist_pos)
        # negative_dist_contrib=sum(dist_neg)/len(dist_neg)        
        
        # ax0=axs[1]  
        # ax0.set_title('Temporal average of contributions')       
        # ax0.axhline(0, color='black', ls='dashed')
        
        # if abs(np.mean(positive_delta_nep+negative_delta_nep))>0.01:
            # ax0.bar ('Net $\Delta$ nep', positive_delta_nep+negative_delta_nep, width=0.4, color="black", label='Net $\Delta$ NEP')        
        # ax0.bar ('$\Delta$ gpp', positive_gpp_contrib, color="green", label='$\Delta$ gpp')
        # ax0.bar ('$\Delta$ gpp', negative_gpp_contrib, color="green")        
        # ax0.bar ('$\Delta$ ra', positive_ra_contrib, color="blue", label='$\Delta$ ra')
        # ax0.bar ('$\Delta$ ra', negative_ra_contrib, color="blue")
        # ax0.bar ('$\Delta$ rh', positive_rh_contrib, color="red", label='$\Delta$ rh')
        # ax0.bar ('$\Delta$ rh', negative_rh_contrib, color="red")        
        # # ax0.bar ('$\Delta$ dist', positive_dist_contrib, color="purple", label='$\Delta$ dist')
        # # ax0.bar ('$\Delta$ dist', negative_dist_contrib, color="purple")   
        # ax0.legend(bbox_to_anchor =(1,1))  
        # ax0.grid() 

        # abs_total=gpp_contrib+rh_contrib+ra_contrib+dist_contrib
                
        # percent_gpp=gpp_contrib/abs_total*100
        # percent_ra=ra_contrib/abs_total*100
        # percent_rh=rh_contrib/abs_total*100
        # percent_dist=dist_contrib/abs_total*100
        
        # ######### Percents
        # if percent==True:   
            # fig3=plt.figure(figsize=(15,5))
            # axs=fig3.subplots(1,2, gridspec_kw={'width_ratios': [2, 1]})         
            # # if delta==False: 
            # ax=axs[0]      
            
            # # % timeline
       
            # # quadratic trends
            # z = np.polyfit(times,  percent_gpp, 2)
            # p = np.poly1d(z)  
            # ax.plot(        
                # times, 
                # p(times),
                # label="$\Delta$ gpp",
                # color="green",
            # ) 
            # ax.fill_between(
                    # times,
                    # percent_gpp, 
                    # p(times),
                    # color="green",
                    # alpha=0.2                
                    # )  
                    
            # z = np.polyfit(times,  percent_ra, 2)
            # p = np.poly1d(z)  
            # ax.plot(        
                # times, 
                # p(times),
                # label="$\Delta$ NPP",
                # color="blue",
            # ) 
            # ax.fill_between(
                    # times,
                    # percent_ra, 
                    # p(times),
                    # color="blue",
                    # alpha=0.2                
                    # )  

            # z = np.polyfit(times,  percent_rh, 2)
            # p = np.poly1d(z)  
            # ax.plot(        
                # times, 
                # p(times),
                # label="$\Delta$ rh",
                # color="red",
            # )         
            # ax.fill_between(
                    # times,
                    # percent_rh, 
                    # p(times),
                    # color="red",
                    # alpha=0.2                
                    # )  

            # # z = np.polyfit(times,  percent_dist, 2)
            # # p = np.poly1d(z)  
            # # ax.plot(        
                # # times, 
                # # p(times),
                # # label="$\Delta$NPP*$\Delta$RT",
                # # color="purple",
            # # ) 
            # # ax.fill_between(
                    # # times,
                    # # percent_dist, 
                    # # p(times),
                    # # color="purple",
                    # # alpha=0.2                
                    # # )  
            
            # ax.set_title('% Contributions over time')
            # ax.set_ylabel('%')
            # ax.grid()
            
            # # pie charts
            # ax1=axs[1]
            # ax1.set_title('Temporal average % contributions')
              

            # labels = '$\Delta$ GPP', '$\Delta$ RA', '$\Delta$NPP*$\Delta$ RH'#, '$\Delta$ Dist'
            # sizes = [np.mean(percent_gpp), np.mean(percent_ra), 
                # np.mean(percent_rh)]#, np.mean(percent_dist)]
            # ax1.pie(sizes, autopct='%1.1f%%', 
                # startangle=90, counterclock=False, colors=("green", "blue", "red"))#, "purple"))
       
            # ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
            # ax1.legend(labels, bbox_to_anchor =(1,1))
        
        # plt.show()       