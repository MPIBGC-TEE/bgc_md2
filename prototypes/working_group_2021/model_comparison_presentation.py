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

# # Matrix Model Intercomparison Project: Traceability Analysis
# <p>We use the <a href="https://github.com/MPIBGC-TEE/bgc_md2">bgc_md2</a> package infrastructure to compare two or more models using traceability analysis approach (<a href="https://doi.org/10.1111/gcb.12172">Xia et al. 2013</a>, <a href="https://doi.org/10.5194/bg-14-145-2017">Luo et al. 2017</a>, <a href="https://doi.org/10.1002/2017MS001004">Jiang et al. 2017</a>) modified for transient simulations of the TRENDY model intercomparison project. </p>
# <p>The biogeochemical models from the TRENDY project have been reconstructed from literature using matrix approach (<a href="https://doi.org/10.1029/2002GB001923">Luo et al. 2003</a>, <a href=" https://doi.org/10.1111/gcb.12766">2015</a>, <a href="https://doi.org/10.5194/gmd-5-1045-2012">Sierra et al. 2012</a>) and <a href="https://www.sympy.org/en/index.html">sympy</a> python package in the <a href="https://github.com/MPIBGC-TEE/bgc_md2">bgc_md2</a> database. Currently the models are simplified to be driven by NPP and thus omitting explicit simulation of leaf processes and autotrophic respiration. Data assimilation (<a href=" https://doi.org/10.1890/09-1275.1">Luo et al. 2011</a>, <a href="https://doi.org/10.1002/2015GB005239">2015</a>) was used to optimize parameters of reconstructed models to fit TRENDY output.</p>
# <p>Currently our analysis includes the following steps:
# <ol>
# <li>We start by comparing model outputs - <em><strong>C storage</strong></em> ($X$) over time. Then we compute and compare traceable components to investigate sources of discrepancy between model predictions of C storage, and thus sources of uncertainty in our understanding of global C dynamics. </li>
# <li>We compute <em><strong>C storage capacity </strong></em> $X_{C}$ for each model. $X_{C}$ (as a function of time) represents the maximum C storage for the system under current conditions. 
# It is a theoretically predicted steady state of the autonomous system that we would obtain if we froze our original non-autonomous system at time point $t_{freeze}$ keeping the conditions (C input and environmental factors) constant after $t_{freeze}$. While $X_{C}$ is not a correct prediction for the solution $X$ of the original time dependent system, it is rather a momentary characteristic of the system that depends on conditions at each point in time.  
# $X_{C}$ is an attractor of the frozen autonomous system, although it is not attracting the real solution  $X$. In a non-autonomous system with variable C input and environmental factors (such as temperature and moisture)  <em><strong>ecosystem</strong></em> $X$ (if we consider ecosystem as a whole - as a <em>"surrogate" 1-pool system</em>, not pool-wise) is constantly chasing $X_{C}$, but cannot consistently reach it due to variable conditions: when $X_{C}$ &gt $X$, $X$ is increasing, and when $X_{C}$ &lt $X$, $X$ is decreasing.</li>
# <li>The difference between $X_{C}$ and $X$ at each point in time is called  <em><strong>C storage potential </strong></em>($X_{P}$): it shows how far the current C storage of the system is from the theoretical steady state (in a positive or negative direction).</li> 
# <li>$X_{C}$ of each model depends on <em><strong>C input</strong></em> (in our case - $NPP$) and <em><strong>Equilibrium Residence Time </strong></em>($RT$). We compare $NPP$ and $RT$ for each model and attribute the contribution of $NPP$ and $RT$ to the discrepancy of $X_{C}$.</li>
# <li>$RT$ represents the time which it would take on average for C particle to exit the system after it entered it, if the system was at the equilibrium. In a transient simulation, where the system is not at the steady state, $RT$ becomes a dynamic model characteristic that can be interpreted as a measure of the inverse C turnover rate of the system at each point in time. $RT$ depends on the model structure and on environmental factors. To single out environmental effects we determine <em><strong>temperature</strong></em> and <em><strong>moisture sensitivity</strong></em> of $RT$, and compare it between models. We also compare $RT$ at fixed temperatures and moisture conditions (gauge conditions) including minimum, maximum and mean conditions across the period of the simulation.</li>
# <li> Based on the traceable components: $NPP$, $RT$, <em><strong>temperature</strong></em> and <em><strong>moisture sensitivity</strong></em> of $RT$, we can make conclusions regarding the inherent similarity/dissimilarity between the models and single out mechanisms that contribute most to the discrepancy of C predictions. 
# </ol>
# <p>The analysis is currently performed for total global C storage, including vegetation and soil parts. <br> Next steps to expand the analysis may include the following: </p>
# <ul>
# <li>Explore temperature and moisture sensitivities of major fluxes (e.g. litter decomposition, soil respiration);</li> 
# <li>Separetely compare vegetation and soil components for each model;</li> 
# <li>Investigate differences over different biomes;</li>
# <li>Include additional diagnostic variables (e.g. transit time, carbon use efficiency, etc.);</li>
# <li>Expand models by explicitly including autotrophic processes (start from GPP as C input);</li>
# <li>Include sensitivity to more environmental factors and anthropogenic disturbances</li>
# </ul>
# <p>The short description of methodology for deriving traceable components is given below. </p>

from IPython.display import Markdown, display
display(Markdown("TracebilityText.md"))

# ### Loading required packages  and functions

# %load_ext autoreload
# %autoreload 2
import matplotlib.pyplot as plt
import numpy as np
from functools import lru_cache
import general_helpers as gh
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    StateVariableTuple
)

# ### Selecting models to compare

# define models to compare as a dictionary (folder name : model name)
model_names={
    "yz_jules": "JULES",
    "kv_visit2": "VISIT",
    "jon_yib": "YIBs",
    "kv_ft_dlem": "DLEM",
    #"Aneesh_SDGVM":"SDGVM",
    "cj_isam": "ISAM",
    "bian_ibis2":"IBIS",
    "ORCHIDEE-V2":"OCN",
}

# ### Exploring model structures

# bgc_md2 automatically generates a flow diagram, input allocation vector and compartmental matrix for each model
comparison_table=gh.model_table(model_names)
comparison_table

# ### Loading TRENDY data and model parameters

# define same step size for each model (in days)
delta_t_val=30

# +
# load data and parameters
model_folders=[(k) for k in model_names]
test_arg_list=gh.get_test_arg_list(model_folders)

# fixme mm 8-12: 
# it would make sense to create a dictionary indexed by the model name 
# so it can be used for model comparisons by name like this one
test_args_dictionary={mf: gh.test_args(mf) for mf in model_folders}
# -

# ### Checking data assimilation quality for all models

# +
# plt.rcParams.update({'font.size': 32})

# +
# import matplotlib.lines as mlines

# for j, mf in enumerate(model_folders):
#     print ('\033[1m'+"Matrix version vs original Trendy output: "+ mf)
#     print ('\033[0m')
#     fig = plt.figure(figsize=(50,10))
#     axs=fig.subplots(1, len(gh.msh(mf).Observables._fields))
    
#     mvs=test_arg_list[j].mvs
#     dvs=test_arg_list[j].dvs
#     cpa=test_arg_list[j].cpa
#     epa_opt=test_arg_list[j].epa_opt
    
#     param2res_sym = gh.msh(mf).make_param2res_sym(mvs,cpa,dvs)
#     out_simu=param2res_sym(epa_opt)._asdict()
#     obs=test_arg_list[j].svs._asdict()
#     print ("Amount of variance explined: ")
#     for i,f in enumerate(gh.msh(mf).Observables._fields):
#         resid=out_simu[f]-obs[f]
#         mean_obs = obs[f].mean()
#         mean_centered_obs = obs[f] - mean_obs
#         AVE=1 - np.sum( resid**2) / np.sum( mean_centered_obs**2 )
#         print(f+ " : " + str(round(AVE,3)) )
#     for i,f in enumerate(gh.msh(mf).Observables._fields):
#         axs[i].scatter(out_simu[f], obs[f], c='black')
#         line = mlines.Line2D([0, 1], [0, 1], color='red')
#         transform = axs[i].transAxes
#         line.set_transform(transform)
#         axs[i].add_line(line)
#         axs[i].set_title(f)
#         axs[i].set_xlabel('Matrix output')
#         axs[i].set_ylabel('Original Trendy output')
#         axs[i].grid()
#     plt.show()
# -


# ### Plots of traceable components

# +
# plt.rcParams.update({'font.size': 18})
# var_names={
#     "x": "X and X_c",
#     #"x_p": "X_p",
#     "u": "C input",
#     "rt": "Residence time",
# }
# gh.plot_components_combined(model_names=model_names,
#                         test_arg_list=test_arg_list,   
#                         var_names=var_names,
#                         delta_t_val=delta_t_val,
#                         model_cols=model_cols,
#                         part=1,
#                         averaging=12*30//delta_t_val # yearly averaging
#                        )

# +
# gh.plot_x_xc(model_names=model_names,
#              test_arg_list=test_arg_list,
#              delta_t_val=delta_t_val, 
#              model_cols=model_cols,
#              part=1,
#              averaging=12,
#              overlap=True
#              )

# +
# gh.plot_normalized_x(model_names=model_names,
#                      delta_t_val=delta_t_val,
#                      test_arg_list=test_arg_list,
#                      model_cols=model_cols,
#                      part=1,
#                      averaging=12*30//delta_t_val,
#                      overlap=True
#                      )

# +
# gh.plot_normalized_xc(model_names=model_names,
#                       test_arg_list=test_arg_list,
#                       delta_t_val=delta_t_val, 
#                       model_cols=model_cols,
#                       part=1,
#                       averaging=12*30*5//delta_t_val,
#                       overlap=True
#                      )

# +
# gh.plot_xp(model_names=model_names,
#            test_arg_list=test_arg_list,
#            delta_t_val=delta_t_val, 
#            model_cols=model_cols,
#            part=1,
#            averaging=12*30*5//delta_t_val,
#            overlap=True
#           )

# +
# gh.plot_u(model_names=model_names,
#           test_arg_list=test_arg_list,
#           delta_t_val=delta_t_val, 
#           model_cols=model_cols,
#           part=1,
#           averaging=30*12*5//delta_t_val,
#           overlap=True
#          )

# +
# gh.plot_normalized_u(model_names=model_names,
#                      test_arg_list=test_arg_list,
#                      delta_t_val=delta_t_val, 
#                      model_cols=model_cols,
#                      part=1,
#                      averaging=12*30*5//delta_t_val,
#                      overlap=True
#                      )

# +
# gh.plot_rt(model_names=model_names,
#            test_arg_list=test_arg_list,
#            delta_t_val=delta_t_val, 
#            model_cols=model_cols,
#            part=1,
#            averaging=12*30//delta_t_val,
#            overlap=True
#           )

# +
# gh.plot_normalized_rt(model_names=model_names,
#                       test_arg_list=test_arg_list,
#                      delta_t_val=delta_t_val, 
#                      model_cols=model_cols,
#                      part=1,
#                      averaging=12*30*5//delta_t_val,
#                      overlap=True
#                      )
# -

# ## Contribution of Residense Time and C Input to the Differences in C Storage Capacity

# +
# # mm test for two files,
# # fixme mm 8-12: 
# # test_args_dictionary={mf: gh.test_args(mf) for mf in model_folders}

# mf_1="yz_jules"
# mf_2="kv_visit2"
# gh.plot_attribution_X_c(
#     mf_1=mf_1,
#     mf_2=mf_2,
#     ta_1=test_args_dictionary[mf_1],
#     ta_2=test_args_dictionary[mf_2],
#     delta_t_val=delta_t_val,
#     part=1
# )


# +
# count_rt_weighted=0
# count_u_weighted=0
# count_combined_weighted=0
# count_delta_x_c=0
# for i in range(len(model_folders)-1):
#     j=i
#     while j<len(model_folders)-1:
#         j+=1
#         mf_1=model_folders[i]
#         mf_2=model_folders[j]
#         ta_1=test_arg_list[i]
#         ta_2=test_arg_list[j]
#         print("Attribution of difference in C storage capacity between "+model_names[mf_1]+" and "+model_names[mf_2])
#         rt,u,combined,delta_x_c=gh.plot_attribution_X_c(mf_1=mf_1, mf_2=mf_2, ta_1=ta_1,ta_2=ta_2, delta_t_val=delta_t_val, part=1)
#         count_rt_weighted=count_rt_weighted+rt*abs(delta_x_c)
#         count_u_weighted=count_u_weighted+u*abs(delta_x_c)
#         count_combined_weighted=count_combined_weighted+combined*abs(delta_x_c)
#         count_delta_x_c=count_delta_x_c+abs(delta_x_c)

# +
# overall_rt_perc=count_rt_weighted/count_delta_x_c
# overall_u_perc=count_u_weighted/count_delta_x_c
# overall_combined_perc=count_combined_weighted/count_delta_x_c

# # Pie chart, where the slices will be ordered and plotted counter-clockwise:
# fig1=plt.figure(figsize=(8,8))
# ax1 = fig1.subplots()
# if overall_combined_perc > 0.001:
#     ax1.pie([overall_rt_perc, overall_combined_perc, overall_u_perc], 
#             labels=('$\Delta$ RT', ' $\Delta$ u * $\Delta$ RT', '$\Delta$ u'), 
#             autopct='%1.1f%%',
#             startangle=90, counterclock=False, colors=("darkorange", "lightgrey", "green"))
# else: 
#     ax1.pie([overall_rt_perc, overall_u_perc], 
#             labels= ('$\Delta$ RT', '$\Delta$ u'), 
#             autopct='%1.1f%%',
#             startangle=90, counterclock=False, colors=("darkorange", "green"))    
# ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
# ax1.set_title('Average Contribution of $\Delta$ Residense Time (RT) and $\Delta$ C Input (u) Across All Pairs of Models')
# plt.show()
# -
plt.rcParams.update({'font.size': 16})

model_cols={
    "JULES": "blue",
    "VISIT": "orange",
    "YIBs": "green",
    "DLEM": "red",
    "SDGVM":"yellow",
    "ISAM": "purple",
    "IBIS":"magenta",
    "OCN":"teal",
}

all_comp_dict= gh.get_traceable_components(model_names=model_names,
             test_arg_list=test_arg_list,
             delta_t_val=delta_t_val, 
             model_cols=model_cols,
             part=1,
             averaging=12*30//delta_t_val, # yearly averaging
             overlap=True
             )

gh.plot_traceable_component(
    all_comp_dict,
    "x_x_c",
    model_cols,
    #delta=True,
)

gh.plot_traceable_component(
    all_comp_dict,
    "x",
    model_cols,
    delta=True,
)

gh.plot_traceable_component(
    all_comp_dict,
    "x_c",
    model_cols,
    delta=True,
)

gh.plot_traceable_component(
    all_comp_dict,
    "x_p",
    model_cols,
    #delta=True,
)

gh.plot_traceable_component(
    all_comp_dict,
    "u",
    model_cols,
    delta=True,
)

gh.plot_traceable_component(
    all_comp_dict,
    "rt",
    model_cols,
    delta=True,
)

x=gh.plot_attribution(
    all_comp_dict=all_comp_dict,
    uncertainty=True,
)

x

delta_x=x[0]
rt_contrib=x[1]
u_contrib=x[2]
rt_u_inter=x[3]
x_p_contrib=x[4]
abs_total=abs(x_p_contrib)+abs(rt_contrib)+abs(u_contrib)+abs(rt_u_inter)
percent_x_p=abs(x_p_contrib)/abs_total
percent_rt=abs(rt_contrib)/abs_total
percent_u=abs(u_contrib)/abs_total
percent_inter=abs(rt_u_inter)/abs_total

        ####### Figure 3 
        print ('\033[1m'+'Attribution of sum of absolute differences from the multi-model mean ' +
            'to the absolute differences of traceable components')         
        fig3=plt.figure(figsize=(15,15))
        axs=fig3.subplots(2,2)
        
        # contributions timeline
        ax=axs[0,0]

        # quadratic trends
        z = np.polyfit(all_comp_dict["Times"],  x_p_contrib, 2)
        p = np.poly1d(z)  
        ax.plot(        
            all_comp_dict["Times"], 
            p(all_comp_dict["Times"]),
            label="$\Delta$ X_p",
            color="blue",
        ) 
        ax.fill_between(
                all_comp_dict["Times"],
                x_p_contrib, 
                p(all_comp_dict["Times"]),
                #label="\u00B12$\sigma$ confidence interval",
                color="blue",
                #linewidth=0.1,
                alpha=0.2                
                )         

        z = np.polyfit(all_comp_dict["Times"],  u_contrib, 2)
        p = np.poly1d(z)  
        ax.plot(        
            all_comp_dict["Times"], 
            p(all_comp_dict["Times"]),
            label="$\Delta$ X_p",
            color="green",
        ) 
        ax.fill_between(
                all_comp_dict["Times"],
                u_contrib, 
                p(all_comp_dict["Times"]),
                #label="\u00B12$\sigma$ confidence interval",
                color="green",
                #linewidth=0.1,
                alpha=0.2                
                )   
                
        z = np.polyfit(all_comp_dict["Times"],  rt_contrib, 2)
        p = np.poly1d(z)  
        ax.plot(        
            all_comp_dict["Times"], 
            p(all_comp_dict["Times"]),
            label="$\Delta$ X_p",
            color="darkorange",
        )              
        ax.fill_between(
                all_comp_dict["Times"],
                rt_contrib, 
                p(all_comp_dict["Times"]),
                #label="\u00B12$\sigma$ confidence interval",
                color="darkorange",
                #linewidth=0.1,
                alpha=0.2                
                )         

        z = np.polyfit(all_comp_dict["Times"],  rt_u_inter, 2)
        p = np.poly1d(z)  
        ax.plot(        
            all_comp_dict["Times"], 
            p(all_comp_dict["Times"]),
            label="$\Delta$ X_p",
            color="grey",
        )         
        ax.fill_between(
                all_comp_dict["Times"],
                rt_u_inter, 
                p(all_comp_dict["Times"]),
                #label="\u00B12$\sigma$ confidence interval",
                color="grey",
                #linewidth=0.1,
                alpha=0.2                
                )          
                 
        
        #ax.legend()
        ax.set_title('Contributions over time')
        #ax.set_ylabel('%')
        ax.grid()
        
        # % timeline
        ax=axs[1,0]
        ax.plot(          
            all_comp_dict["Times"], 
            percent_rt*100,
            label="contribution of $\Delta$ R_t",
            color="darkorange",
            linewidth=0.8,
            alpha=0.5,
        )
        ax.plot(        
            all_comp_dict["Times"], 
            percent_u*100,
            label="contribution of $\Delta$ u",
            color="green",
            linewidth=0.8,
            alpha=0.5,
        )     
        
        # ax.plot(        
            # all_comp_dict["Times"], 
            # delta_x,
            # label="$\Delta$ X",
            # color="black",
        # )
        ax.plot(        
            all_comp_dict["Times"], 
            percent_x_p*100,
            label="$\Delta$ X_p",
            color="blue",
            linewidth=0.8,
            alpha=0.5,
        )        
        if np.mean(percent_inter) > 0.001:
            ax.plot(            
                all_comp_dict["Times"], 
                percent_inter*100,
                label="contribution of interaction terms",
                color="grey",
                linewidth=0.5,
                #alpha=0.1,
            )      
        # quadratic trends
        z = np.polyfit(all_comp_dict["Times"],  percent_x_p*100, 2)
        p = np.poly1d(z)  
        ax.plot(        
            all_comp_dict["Times"], 
            p(all_comp_dict["Times"]),
            label="$\Delta$ X_p",
            color="blue",
        ) 
        ax.fill_between(
                all_comp_dict["Times"],
                percent_x_p*100, 
                p(all_comp_dict["Times"]),
                #label="\u00B12$\sigma$ confidence interval",
                color="blue",
                #linewidth=0.1,
                alpha=0.2                
                )  
                
        z = np.polyfit(all_comp_dict["Times"],  percent_u*100, 2)
        p = np.poly1d(z)  
        ax.plot(        
            all_comp_dict["Times"], 
            p(all_comp_dict["Times"]),
            label="$\Delta$ X_p",
            color="green",
        ) 
        ax.fill_between(
                all_comp_dict["Times"],
                percent_u*100, 
                p(all_comp_dict["Times"]),
                #label="\u00B12$\sigma$ confidence interval",
                color="green",
                #linewidth=0.1,
                alpha=0.2                
                )  

        z = np.polyfit(all_comp_dict["Times"],  percent_rt*100, 2)
        p = np.poly1d(z)  
        ax.plot(        
            all_comp_dict["Times"], 
            p(all_comp_dict["Times"]),
            label="$\Delta$ X_p",
            color="darkorange",
        )         
        ax.fill_between(
                all_comp_dict["Times"],
                percent_rt*100, 
                p(all_comp_dict["Times"]),
                #label="\u00B12$\sigma$ confidence interval",
                color="darkorange",
                #linewidth=0.1,
                alpha=0.2                
                )  

        z = np.polyfit(all_comp_dict["Times"],  percent_inter*100, 2)
        p = np.poly1d(z)  
        ax.plot(        
            all_comp_dict["Times"], 
            p(all_comp_dict["Times"]),
            label="$\Delta$ X_p",
            color="grey",
        ) 
        ax.fill_between(
                all_comp_dict["Times"],
                percent_inter*100, 
                p(all_comp_dict["Times"]),
                #label="\u00B12$\sigma$ confidence interval",
                color="grey",
                #linewidth=0.1,
                alpha=0.2                
                )  
        
        #ax.legend()
        ax.set_title('% Contributions over time')
        ax.set_ylabel('%')
        ax.grid()

        # bar charts
        
        positive_delta_X=sum(delta_x[delta_x>0])/len(delta_x)
        negative_delta_X=sum(delta_x[delta_x<0])/len(delta_x)
        #positive_x_c_contrib=sum(x_c_contrib[x_c_contrib>0])/len(x_c_contrib)
        #negative_x_c_contrib=sum(x_c_contrib[x_c_contrib<0])/len(x_c_contrib)
        positive_x_p_contrib=sum(x_p_contrib[x_p_contrib>0])/len(x_p_contrib)
        negative_x_p_contrib=sum(x_p_contrib[x_p_contrib<0])/len(x_p_contrib)        
        positive_rt_contrib=sum(rt_contrib[rt_contrib>0])/len(rt_contrib)
        negative_rt_contrib=sum(rt_contrib[rt_contrib<0])/len(rt_contrib)
        positive_u_contrib=sum(u_contrib[u_contrib>0])/len(u_contrib)
        negative_u_contrib=sum(u_contrib[u_contrib<0])/len(u_contrib)      
        positive_rt_u_inter=sum(rt_u_inter[rt_u_inter>0])/len(rt_u_inter)
        negative_rt_u_inter=sum(rt_u_inter[rt_u_inter<0])/len(rt_u_inter)        
        
        ax0=axs[0,1]  
        ax0.set_title('Average contributions')       
        ax0.axhline(0, color='black', ls='dashed')
        
        # ax0.bar ('$\Delta$ X', positive_delta_X, color="black")
        # ax0.bar ('$\Delta$ X', negative_delta_X, color="black")
        ax0.bar ('$\Delta$ RT', positive_rt_contrib, color="darkorange")
        ax0.bar ('$\Delta$ RT', negative_rt_contrib, color="darkorange")        
        ax0.bar ('$\Delta$ NPP', positive_u_contrib, color="green")
        ax0.bar ('$\Delta$ NPP', negative_u_contrib, color="green")
        ax0.bar ('$\Delta$npp*$\Delta$rt', positive_rt_u_inter, color="lightgrey")
        ax0.bar ('$\Delta$npp*$\Delta$rt', negative_rt_u_inter, color="lightgrey")        
        ax0.bar ('$\Delta$ X_p', positive_x_p_contrib, color="blue")
        ax0.bar ('$\Delta$ X_p', negative_x_p_contrib, color="blue")   
        
        # if abs(np.mean(combined_contrib)) > 0.001:
           # ax0.bar ('$\Delta$ U * $\Delta$ RT', np.mean(combined_contrib), color="lightgrey")   
        ax0.grid() 
        
        # pie charts

        ax1=axs[1,1]
        ax1.set_title('% Contributions')
          
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
        ax1.legend(labels, bbox_to_anchor =(1.42, 2.22))
        plt.show()


