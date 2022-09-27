 
# def plot_attribution_X_c (mf_1, mf_2, ta_1, ta_2, delta_t_val, part):
    # if part<0 | part >1: 
        # raise Exception('Invalid partitioning in plot_diff: use part between 0 and 1')
    # test_arg_list=[ta_1,ta_2]

    # itr_1=traceability_iterator_instance(mf_1,ta_1,delta_t_val)
    # itr_2=traceability_iterator_instance(mf_2,ta_2,delta_t_val)

    # start_min_1,stop_max_1=min_max_index(ta_1,delta_t_val,*t_min_tmax_overlap(test_arg_list,delta_t_val))
    # start_min_2,stop_max_2=min_max_index(ta_2,delta_t_val,*t_min_tmax_overlap(test_arg_list,delta_t_val))

    # # if we do not want the whole interval but look at a smaller part to observe the dynamics
    # start_1,stop_1 = int(stop_max_1-(stop_max_1-start_min_1)*part), stop_max_1
    # start_2,stop_2 = int(stop_max_2-(stop_max_2-start_min_2)*part), stop_max_2
    # times_1=times_in_days_aD(ta_1,delta_t_val)[start_1:stop_1]/days_per_year()
    # times_2=times_in_days_aD(ta_2,delta_t_val)[start_2:stop_2]/days_per_year()
    # vals_1=itr_1[start_1:stop_1]
    # vals_2=itr_2[start_2:stop_2]

    # x_c_1=interp1d(times_1,vals_1.x_c)
    # x_c_2=interp1d(times_2,vals_2.x_c)
    # u_1=interp1d(times_1,vals_1.u)
    # u_2=interp1d(times_2,vals_2.u)
    # rt_1=interp1d(times_1,vals_1.rt)
    # rt_2=interp1d(times_2,vals_2.rt)

    # # common plot times
    # start=max(times_1.min(),times_2.min())
    # stop=min(times_1.max(),times_2.max())
    # nstep=min(len(times_1),len(times_2))
    # times=np.linspace(start,stop,nstep)
    
    # # values for plots
    # delta_u=u_1(times)-u_2(times)
    # delta_rt=rt_1(times)-rt_2(times)
    # delta_x_c=x_c_1(times)-x_c_2(times)

    # rt_contrib=delta_rt*(u_1(times)-delta_u/2)
    # u_contrib=delta_u*(rt_1(times)-delta_rt/2)
    # combined_contrib=delta_x_c-rt_contrib-u_contrib

    # percent_rt=abs(rt_contrib)/(abs(rt_contrib)+abs(u_contrib)+abs(combined_contrib))*100
    # percent_u=abs(u_contrib)/(abs(rt_contrib)+abs(u_contrib)+abs(combined_contrib))*100
    # percent_combined=abs(combined_contrib)/(abs(rt_contrib)+abs(u_contrib)+abs(combined_contrib))*100
    
    # averaging=12*30//delta_t_val # yearly averages
    # fig1=plt.figure(figsize=(17,8))
    # ax=fig1.subplots()
    # ax.plot(
        # #times[60*12:60*12+36],                    
        # #rt_contrib[60*12:60*12+36],            
        # avg_timeline(times, averaging), 
        # avg_timeline(percent_rt, averaging),
        # label="contribution of $\Delta$ RT",
        # color="darkorange",
    # )
    # ax.plot(
        # #times[60*12:60*12+36],                    
        # #u_contrib[60*12:60*12+36],            
        # avg_timeline(times, averaging), 
        # avg_timeline(percent_u, averaging),
        # label="contribution of $\Delta$ u",
        # color="green",
    # )
    # if np.mean(percent_combined) > 0.001:
        # ax.plot(
            # #times[60*12:60*12+36],                    
            # #u_contrib[60*12:60*12+36],            
            # avg_timeline(times, averaging),
            # avg_timeline(percent_combined, averaging),
            # label="contribution of $\Delta$ u * $\Delta$ RT",
            # color="grey",
        # )
    
    # ax.legend()
    # ax.set_title('Contribution of $\Delta$ Residense Time (RT) and $\Delta$ C Input (u) Over Time')
    # ax.set_ylabel('%')
    # ax.grid()
    
    # fig=plt.figure(figsize=(17,8))
    # axs=fig.subplots(1,2)
    
    # # bar charts
    # ax0=axs[0]  
    # ax0.set_title('Average Contribution of RT and u to the difference in X_c')
    # #labels = '$\Delta$ RT', '$\Delta$ U * $\Delta$ RT', '$\Delta$ U'
    # #sizes = [np.mean(percent_rt), np.mean(percent_combined), np.mean(percent_u)]
    # #bar_labels = '$\Delta$ X_c', '$\Delta$ RT', '$\Delta$ U', '$\Delta$ U * $\Delta$ RT'
    # #bar_sizes=[round(np.mean(rt_contrib),0), np.mean(u_contrib), np.mean(combined_contrib)]
    
    # ax0.axhline(0, color='black', ls='dashed')
    # ax0.bar ('$\Delta$ X_c', np.mean(delta_x_c), color="blue")
    # ax0.bar ('$\Delta$ RT', np.mean(rt_contrib), color="darkorange")
    # ax0.bar ('$\Delta$ u', np.mean(u_contrib), color="green")
    # if abs(np.mean(combined_contrib)) > 0.001:
        # ax0.bar ('$\Delta$ U * $\Delta$ RT', np.mean(combined_contrib), color="lightgrey")   
    # ax0.grid()  
    # # pie charts
    # ax1=axs[1]
    # ax1.set_title('Average Contribution of RT and u')
      
    # if np.mean(percent_combined) > 0.001:
        # labels = '$\Delta$ RT', '$\Delta$ U * $\Delta$ RT', '$\Delta$ U'
        # sizes = [np.mean(percent_rt), np.mean(percent_combined), np.mean(percent_u)]
        # ax1.pie(sizes, labels=labels, autopct='%1.1f%%', 
            # startangle=90, counterclock=False, colors=("darkorange", "lightgrey", "green"))
    # else:
        # labels = '$\Delta$ RT', '$\Delta$ U'
        # sizes = [np.mean(percent_rt), np.mean(percent_u)]
        # ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
            # startangle=90, counterclock=False, colors=("darkorange", "green"))        
    # ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    # plt.show()
    
    # return (np.mean(percent_rt), np.mean(percent_u), np.mean(percent_combined), np.mean(delta_x_c)) 



               # ####### Figure 1  
        # print ('\033[1m'+'Attribution of X to X_c and X_p (differences from the mean)')        
        # fig1=plt.figure(figsize=(17,8))
        # axs=fig1.subplots(1,2)
        # ax=axs[0]
                       
        # ax.plot(          
            # all_comp_dict["Times"], 
            # x_c_contrib,
            # label="contribution of $\Delta$ X_c",
            # color="red",
        # )
        # ax.plot(        
            # all_comp_dict["Times"], 
            # x_p_contrib,
            # label="contribution of $\Delta$ X_p",
            # color="blue",
        # )
        # ax.plot(        
            # all_comp_dict["Times"], 
            # delta_x,
            # label="$\Delta$ X",
            # color="grey",
        # )
  
        # ax.axhline(0, color='grey', ls='dashed')
        # #ax.legend()
        # ax.set_title('Contributions over time')
        # #ax.set_ylabel('%')
        # ax.grid()
 
         # # bar charts
        # ax0=axs[1]  
        # ax0.set_title('Average contributions')
        # positive_delta_X=sum(delta_x[delta_x>0])/len(delta_x)
        # negative_delta_X=sum(delta_x[delta_x<0])/len(delta_x)
        # positive_x_c_contrib=sum(x_c_contrib[x_c_contrib>0])/len(x_c_contrib)
        # negative_x_c_contrib=sum(x_c_contrib[x_c_contrib<0])/len(x_c_contrib)
        # positive_x_p_contrib=sum(x_p_contrib[x_p_contrib>0])/len(x_p_contrib)
        # negative_x_p_contrib=sum(x_p_contrib[x_p_contrib<0])/len(x_p_contrib)
        
        # ax0.axhline(0, color='grey', ls='dashed')
        
        # ax0.bar ('$\Delta$ X', positive_delta_X, color="grey")
        # ax0.bar ('$\Delta$ X', negative_delta_X, color="grey")
        # ax0.bar ('$\Delta$ X_c', positive_x_c_contrib, color="red")
        # ax0.bar ('$\Delta$ X_c', negative_x_c_contrib, color="red")        
        # ax0.bar ('$\Delta$ X_p', positive_x_p_contrib, color="blue")
        # ax0.bar ('$\Delta$ X_p', negative_x_p_contrib, color="blue")
  
        # ax0.grid()  
        # plt.show()
        
        # # # pie charts

        # # ax1=axs[2]
        # # ax1.set_title('Average contributions in %')
          
        # # if np.mean(percent_inter_1) > 0.001:
            # # labels = '$\Delta$ X_c', 'Residual', '$\Delta$ X_p'
            # # sizes = [np.mean(percent_x_c), np.mean(percent_inter_1), np.mean(percent_x_p)]
            # # ax1.pie(sizes, labels=labels, autopct='%1.1f%%', 
                # # startangle=90, counterclock=False, colors=("blue", "purple", "red"))
        # # else:
            # # labels = '$\Delta$ X_c', '$\Delta$ X_p'
            # # sizes = [np.mean(percent_x_c), np.mean(percent_x_p)]
            # # ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                # # startangle=90, counterclock=False, colors=("blue", "red"))        
        # # ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.        
        # # ax1.legend(bbox_to_anchor =(1.3, 1))
        # # plt.show()
        
        # print("X: "+str(positive_delta_X-negative_delta_X)+
            # "; X_c: "+str(positive_x_c_contrib-negative_x_c_contrib) + 
            # "; X_p: "+str(positive_x_p_contrib-negative_x_p_contrib)  
        # )  
        # print("sum check: "+str(
            # (positive_delta_X-negative_delta_X)-
            # (positive_x_c_contrib-negative_x_c_contrib)+
            # (positive_x_p_contrib-negative_x_p_contrib)
            # )
        # )
        
        # ####### Figure 2
        # print ('\033[1m'+'Attribution of X_c to NPP and RT (differences from the mean)') 
        # fig2=plt.figure(figsize=(17,8))
        # axs=fig2.subplots(1,2)
        # ax=axs[0]
        # ax.plot(          
            # all_comp_dict["Times"], 
            # rt_contrib,
            # label="contribution of $\Delta$ RT",
            # color="darkorange",
        # )
        # ax.plot(        
            # all_comp_dict["Times"], 
            # u_contrib,
            # label="contribution of $\Delta$ NPP",
            # color="green",
        # )
        # ax.plot(        
            # all_comp_dict["Times"], 
            # delta_x_c,
            # label="$\Delta$ X_c",
            # color="red",
        # )        
        # if np.mean(percent_rt_u_inter) > 0.001:
            # ax.plot(            
                # all_comp_dict["Times"], 
                # #percent_rt_u_inter,
                # rt_u_inter,
                # label="contribution of $\Delta$NPP*$\Delta$RT",
                # color="purple",
            # )        
        # #ax.legend()
        # ax.axhline(0, color='grey', ls='dashed')
        # ax.set_title('Contributions over time')
        # #ax.set_ylabel('%')
        # ax.grid()  
        
        # # bar charts
        # ax0=axs[1]  
        # ax0.set_title('Average contributions')
        # positive_delta_X_c=sum(delta_x_c[delta_x_c>0])/len(delta_x_c)
        # negative_delta_X_c=sum(delta_x_c[delta_x_c<0])/len(delta_x_c)
        # positive_rt_contrib=sum(rt_contrib[rt_contrib>0])/len(rt_contrib)
        # negative_rt_contrib=sum(rt_contrib[rt_contrib<0])/len(rt_contrib)
        # positive_u_contrib=sum(u_contrib[u_contrib>0])/len(u_contrib)
        # negative_u_contrib=sum(u_contrib[u_contrib<0])/len(u_contrib)      
        # positive_rt_u_inter=sum(rt_u_inter[rt_u_inter>0])/len(rt_u_inter)
        # negative_rt_u_inter=sum(rt_u_inter[rt_u_inter<0])/len(rt_u_inter)
        
        # ax0.axhline(0, color='grey', ls='dashed')
        
        # ax0.bar ('$\Delta$ X_c', positive_delta_X_c, color="red")
        # ax0.bar ('$\Delta$ X_c', negative_delta_X_c, color="red")
        # ax0.bar ('$\Delta$ RT', positive_rt_contrib, color="darkorange")
        # ax0.bar ('$\Delta$ RT', negative_rt_contrib, color="darkorange")        
        # ax0.bar ('$\Delta$ u', positive_u_contrib, color="green")
        # ax0.bar ('$\Delta$ u', negative_u_contrib, color="green")
        # ax0.bar ('$\Delta$NPP*$\Delta$RT', positive_rt_u_inter, color="purple")
        # ax0.bar ('$\Delta$NPP*$\Delta$RT', negative_rt_u_inter, color="purple")        
        
        # ax0.grid() 
        # plt.show()
        
        # # # pie charts

        # # ax1=axs[2]
        # # ax1.set_title('Average Contributions in %')
          
        # # if np.mean(percent_rt_u_inter) > 0.001:
            # # labels = '$\Delta$ RT', '$\Delta$NPP*$\Delta$RT', '$\Delta$ NPP'
            # # sizes = [np.mean(percent_rt), np.mean(percent_rt_u_inter), np.mean(percent_u)]
            # # ax1.pie(sizes, labels=labels, autopct='%1.1f%%', 
                # # startangle=90, counterclock=False, colors=("darkorange", "purple", "green"))
        # # else:
            # # labels = '$\Delta$ RT', '$\Delta$ NPP'
            # # sizes = [np.mean(percent_x_c), np.mean(percent_x_p)]
            # # ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                # # startangle=90, counterclock=False, colors=("darkorange", "green"))        
        # # ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle. 
        # # ax1.legend(bbox_to_anchor =(1.3, 1))
        # # plt.show()

        # print("X_c: "+str(positive_delta_X_c-negative_delta_X_c)+
            # "; RT: "+str(positive_rt_contrib-negative_rt_contrib)+
            # "; NPP: "+str(positive_u_contrib-negative_u_contrib)+
            # "; rt_u_inter: "+ str(positive_rt_u_inter-negative_rt_u_inter)
        # )   
        # print("sum check: "+str(
            # (positive_delta_X_c-negative_delta_X_c)-
            # (positive_rt_contrib-negative_rt_contrib)-
            # (positive_u_contrib-negative_u_contrib)-
            # (positive_rt_u_inter-negative_rt_u_inter)
            # )
        # )        
        
        # ####### Figure 3
        # print ('\033[1m'+'Attribution of X_p to NEP and RT (differences from the mean)') 
        # fig2=plt.figure(figsize=(17,8))
        # axs=fig2.subplots(1,2)
        # ax=axs[0]
        # ax.plot(          
            # all_comp_dict["Times"], 
            # rt_contrib2,
            # label="contribution of $\Delta$ RT",
            # color="darkorange",
        # )
        # ax.plot(        
            # all_comp_dict["Times"], 
            # nep_contrib,
            # label="contribution of $\Delta$ NEP",
            # color="teal",
        # )
        # ax.plot(        
            # all_comp_dict["Times"], 
            # delta_x_p,
            # label="$\Delta$ X_p",
            # color="red",
        # )        
        # if np.mean(rt_nep_inter) > 0.001:
            # ax.plot(            
                # all_comp_dict["Times"], 
                # rt_nep_inter,
                # label="contribution of $\Delta$NEP*$\Delta$RT",
                # color="purple",
            # )        
        # #ax.legend()
        # ax.axhline(0, color='grey', ls='dashed')
        # ax.set_title('Contributions over time')
        # #ax.set_ylabel('%')
        # ax.grid()  
        
        # # bar charts
        # ax0=axs[1]  
        # ax0.set_title('Average contributions')
        # positive_delta_X_p=sum(delta_x_p[delta_x_p>0])/len(delta_x_p)
        # negative_delta_X_p=sum(delta_x_p[delta_x_p<0])/len(delta_x_p)
        # positive_rt_contrib2=sum(rt_contrib2[rt_contrib2>0])/len(rt_contrib2)
        # negative_rt_contrib2=sum(rt_contrib2[rt_contrib2<0])/len(rt_contrib2)
        # positive_nep_contrib=sum(nep_contrib[nep_contrib>0])/len(nep_contrib)
        # negative_nep_contrib=sum(nep_contrib[nep_contrib<0])/len(nep_contrib)      
        # positive_rt_nep_inter=sum(rt_nep_inter[rt_nep_inter>0])/len(rt_nep_inter)
        # negative_rt_nep_inter=sum(rt_nep_inter[rt_nep_inter<0])/len(rt_nep_inter)
        # ax0.axhline(0, color='grey', ls='dashed')
        
        # ax0.bar ('$\Delta$ X_p', positive_delta_X_p, color="red")
        # ax0.bar ('$\Delta$ X_p', negative_delta_X_p, color="red")
        # ax0.bar ('$\Delta$ RT', positive_rt_contrib2, color="darkorange")
        # ax0.bar ('$\Delta$ RT', negative_rt_contrib2, color="darkorange")        
        # ax0.bar ('$\Delta$ NEP', positive_nep_contrib, color="teal")
        # ax0.bar ('$\Delta$ NEP', negative_nep_contrib, color="teal")
        # ax0.bar ('$\Delta$NEP*$\Delta$RT', positive_rt_nep_inter, color="purple")
        # ax0.bar ('$\Delta$NEP*$\Delta$RT', negative_rt_nep_inter, color="purple")        
        
        # ax0.grid() 
        # plt.show()
        
        # # # pie charts

        # # ax1=axs[2]
        # # ax1.set_title('Average Contributions in %')
          
        # # if np.mean(percent_rt_u_inter) > 0.001:
            # # labels = '$\Delta$ RT', '$\Delta$NPP*$\Delta$RT', '$\Delta$ NPP'
            # # sizes = [np.mean(percent_rt), np.mean(percent_rt_u_inter), np.mean(percent_u)]
            # # ax1.pie(sizes, labels=labels, autopct='%1.1f%%', 
                # # startangle=90, counterclock=False, colors=("darkorange", "purple", "green"))
        # # else:
            # # labels = '$\Delta$ RT', '$\Delta$ NPP'
            # # sizes = [np.mean(percent_x_c), np.mean(percent_x_p)]
            # # ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                # # startangle=90, counterclock=False, colors=("darkorange", "green"))        
        # # ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle. 
        # # ax1.legend(bbox_to_anchor =(1.3, 1))
        # # plt.show()

        # print("X_c: "+str(positive_delta_X_c-negative_delta_X_c)+
            # "; RT: "+str(positive_rt_contrib-negative_rt_contrib)+
            # "; NPP: "+str(positive_u_contrib-negative_u_contrib)+
            # "; rt_u_inter: "+ str(positive_rt_u_inter-negative_rt_u_inter)
        # )   
        # print("sum check: "+str(
            # (positive_delta_X_c-negative_delta_X_c)-
            # (positive_rt_contrib-negative_rt_contrib)-
            # (positive_u_contrib-negative_u_contrib)-
            # (positive_rt_u_inter-negative_rt_u_inter)
            # )
        # )                         







        
               
        # # % timeline
        # ax=axs[1,0]
        # ax.plot(          
            # all_comp_dict["Times"], 
            # percent_rt*100,
            # label="contribution of $\Delta$ R_t",
            # color="darkorange",
            # linewidth=0.8,
            # alpha=0.5,
        # )
        # ax.plot(        
            # all_comp_dict["Times"], 
            # percent_u*100,
            # label="contribution of $\Delta$ u",
            # color="green",
            # linewidth=0.8,
            # alpha=0.5,
        # )     
        
        # # ax.plot(        
            # # all_comp_dict["Times"], 
            # # delta_x,
            # # label="$\Delta$ X",
            # # color="black",
        # # )
        # ax.plot(        
            # all_comp_dict["Times"], 
            # percent_x_p*100,
            # label="$\Delta$ X_p",
            # color="blue",
            # linewidth=0.8,
            # alpha=0.5,
        # )        
        # if np.mean(percent_inter) > 0.001:
            # ax.plot(            
                # all_comp_dict["Times"], 
                # percent_inter*100,
                # label="contribution of interaction terms",
                # color="lightgrey",
                # linewidth=0.5,
                # #alpha=0.1,
            # )      
        # # quadratic trends
        # z = np.polyfit(all_comp_dict["Times"],  percent_x_p*100, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_p",
            # color="blue",
        # ) 
        # ax.fill_between(
                # all_comp_dict["Times"],
                # percent_x_p*100, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="blue",
                # #linewidth=0.1,
                # alpha=0.2                
                # )  
                
        # z = np.polyfit(all_comp_dict["Times"],  percent_u*100, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_p",
            # color="green",
        # ) 
        # ax.fill_between(
                # all_comp_dict["Times"],
                # percent_u*100, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="green",
                # #linewidth=0.1,
                # alpha=0.2                
                # )  

        # z = np.polyfit(all_comp_dict["Times"],  percent_rt*100, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_p",
            # color="darkorange",
        # )         
        # ax.fill_between(
                # all_comp_dict["Times"],
                # percent_rt*100, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="darkorange",
                # #linewidth=0.1,
                # alpha=0.2                
                # )  

        # z = np.polyfit(all_comp_dict["Times"],  percent_inter*100, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_p",
            # color="lightgrey",
        # ) 
        # ax.fill_between(
                # all_comp_dict["Times"],
                # percent_inter*100, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="lightgrey",
                # #linewidth=0.1,
                # alpha=0.2                
                # )  
                
        
        # #ax.legend()
        # ax.set_title('% Contributions over time')
        # ax.set_ylabel('%')
        # ax.grid()
        
        
        

        
            
    # if part<0 | part >1: 
        # raise Exception('Invalid partitioning in plot_diff: use part between 0 and 1')
    # test_arg_list=[ta_1,ta_2]

    # itr_1=traceability_iterator_instance(mf_1,ta_1,delta_t_val)
    # itr_2=traceability_iterator_instance(mf_2,ta_2,delta_t_val)

    # start_min_1,stop_max_1=min_max_index(ta_1,delta_t_val,*t_min_tmax_overlap(test_arg_list,delta_t_val))
    # start_min_2,stop_max_2=min_max_index(ta_2,delta_t_val,*t_min_tmax_overlap(test_arg_list,delta_t_val))

    # # if we do not want the whole interval but look at a smaller part to observe the dynamics
    # start_1,stop_1 = int(stop_max_1-(stop_max_1-start_min_1)*part), stop_max_1
    # start_2,stop_2 = int(stop_max_2-(stop_max_2-start_min_2)*part), stop_max_2
    # times_1=times_in_days_aD(ta_1,delta_t_val)[start_1:stop_1]/days_per_year()
    # times_2=times_in_days_aD(ta_2,delta_t_val)[start_2:stop_2]/days_per_year()
    # vals_1=itr_1[start_1:stop_1]
    # vals_2=itr_2[start_2:stop_2]

    # x_c_1=interp1d(times_1,vals_1.x_c)
    # x_c_2=interp1d(times_2,vals_2.x_c)
    # u_1=interp1d(times_1,vals_1.u)
    # u_2=interp1d(times_2,vals_2.u)
    # rt_1=interp1d(times_1,vals_1.rt)
    # rt_2=interp1d(times_2,vals_2.rt)

    # # common plot times
    # start=max(times_1.min(),times_2.min())
    # stop=min(times_1.max(),times_2.max())
    # nstep=min(len(times_1),len(times_2))
    # times=np.linspace(start,stop,nstep)
    
    # # values for plots
    # delta_u=u_1(times)-u_2(times)
    # delta_rt=rt_1(times)-rt_2(times)
    # delta_x_c=x_c_1(times)-x_c_2(times)

    # rt_contrib=delta_rt*(u_1(times)-delta_u/2)
    # u_contrib=delta_u*(rt_1(times)-delta_rt/2)
    # combined_contrib=delta_x_c-rt_contrib-u_contrib

    # percent_rt=abs(rt_contrib)/(abs(rt_contrib)+abs(u_contrib)+abs(combined_contrib))*100
    # percent_u=abs(u_contrib)/(abs(rt_contrib)+abs(u_contrib)+abs(combined_contrib))*100
    # percent_combined=abs(combined_contrib)/(abs(rt_contrib)+abs(u_contrib)+abs(combined_contrib))*100
    
    # averaging=12*30//delta_t_val # yearly averages
    # fig1=plt.figure(figsize=(17,8))
    # ax=fig1.subplots()
    # ax.plot(
        # #times[60*12:60*12+36],                    
        # #rt_contrib[60*12:60*12+36],            
        # avg_timeline(times, averaging), 
        # avg_timeline(percent_rt, averaging),
        # label="contribution of $\Delta$ RT",
        # color="darkorange",
    # )
    # ax.plot(
        # #times[60*12:60*12+36],                    
        # #u_contrib[60*12:60*12+36],            
        # avg_timeline(times, averaging), 
        # avg_timeline(percent_u, averaging),
        # label="contribution of $\Delta$ u",
        # color="green",
    # )
    # if np.mean(percent_combined) > 0.001:
        # ax.plot(
            # #times[60*12:60*12+36],                    
            # #u_contrib[60*12:60*12+36],            
            # avg_timeline(times, averaging),
            # avg_timeline(percent_combined, averaging),
            # label="contribution of $\Delta$ u * $\Delta$ RT",
            # color="lightgrey",
        # )
    
    # ax.legend()
    # ax.set_title('Contribution of $\Delta$ Residense Time (RT) and $\Delta$ C Input (u) Over Time')
    # ax.set_ylabel('%')
    # ax.grid()
    
    # fig=plt.figure(figsize=(17,8))
    # axs=fig.subplots(1,2)
    
    # # bar charts
    # ax0=axs[0]  
    # ax0.set_title('Average Contribution of RT and u to the difference in X_c')
    # #labels = '$\Delta$ RT', '$\Delta$ U * $\Delta$ RT', '$\Delta$ U'
    # #sizes = [np.mean(percent_rt), np.mean(percent_combined), np.mean(percent_u)]
    # #bar_labels = '$\Delta$ X_c', '$\Delta$ RT', '$\Delta$ U', '$\Delta$ U * $\Delta$ RT'
    # #bar_sizes=[round(np.mean(rt_contrib),0), np.mean(u_contrib), np.mean(combined_contrib)]
    
    # ax0.axhline(0, color='black', ls='dashed')
    # ax0.bar ('$\Delta$ X_c', np.mean(delta_x_c), color="blue")
    # ax0.bar ('$\Delta$ RT', np.mean(rt_contrib), color="darkorange")
    # ax0.bar ('$\Delta$ u', np.mean(u_contrib), color="green")
    # if abs(np.mean(combined_contrib)) > 0.001:
        # ax0.bar ('$\Delta$ U * $\Delta$ RT', np.mean(combined_contrib), color="lightgrey")   
    # ax0.grid()  
    # # pie charts
    # ax1=axs[1]
    # ax1.set_title('Average Contribution of RT and u')
      
    # if np.mean(percent_combined) > 0.001:
        # labels = '$\Delta$ RT', '$\Delta$ U * $\Delta$ RT', '$\Delta$ U'
        # sizes = [np.mean(percent_rt), np.mean(percent_combined), np.mean(percent_u)]
        # ax1.pie(sizes, labels=labels, autopct='%1.1f%%', 
            # startangle=90, counterclock=False, colors=("darkorange", "lightgrey", "green"))
    # else:
        # labels = '$\Delta$ RT', '$\Delta$ U'
        # sizes = [np.mean(percent_rt), np.mean(percent_u)]
        # ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
            # startangle=90, counterclock=False, colors=("darkorange", "green"))        
    # ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    # plt.show()
    
    # return (np.mean(percent_rt), np.mean(percent_u), np.mean(percent_combined), np.mean(delta_x_c)) 

 
################################## Deprecated ########################################

# def difference_attribution (
    # x_1, x_c_1, x_p_1, u_1, rt_1, # traceble components of 1st data stream
    # x_2, x_c_2, x_p_2, u_2, rt_2, # traceble components of 2nd data stream
    # percent=False # results in real values or percents
    # ):
    # x=x_1
    # x_c=x_c_1
    # x_p=x_p_1
    # u=u_1
    # rt=rt_1

    # # differences
    # delta_x=x-x_2
    # delta_x_c=x_c-x_c_2
    # delta_x_p=x_p-x_p_2
    # delta_u=u-u_2
    # delta_rt=rt-rt_2
    
    # # contributions in the difference
    # x_c_contrib=delta_x_c 
    # x_p_contrib=-delta_x_p
    # rt_contrib=delta_rt*(u-delta_u/2)
    # u_contrib=delta_u*(rt-delta_rt/2) 
    # rt_u_inter=delta_x_c-rt_contrib-u_contrib#delta_u*delta_rt     

    # # contributions in percent
    # abs_total=abs(x_p_contrib)+abs(rt_contrib)+abs(u_contrib)+abs(rt_u_inter)
    # percent_x_p=abs(x_p_contrib)/abs_total*100
    # percent_x_c=abs(x_c_contrib)/abs_total*100
    # percent_rt=abs(rt_contrib)/abs_total*100
    # percent_u=abs(u_contrib)/abs_total*100
    # percent_inter=abs(rt_u_inter)/abs_total*100
    
    # # output    
    # attribution_dict = {
        # "x": 100 if percent else delta_x,
        # "x_c": percent_x_c if percent else x_c_contrib,
        # "x_p": percent_x_p if percent else x_p_contrib,
        # "rt": percent_rt if percent else rt_contrib,
        # "u": percent_u if percent else u_contrib,
        # "rt_u_inter": percent_inter if percent else rt_u_inter,
    # }       
    # return (attribution_dict) 

# def plot_attribution_all (
    # all_comp_dict,
    # #uncertainty,
    # #attribution_dict
# ):
    # models=list(all_comp_dict.keys())[:-2]   
    
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
    
    # for m in models: 
           
        # x=all_comp_dict[m]["x"]
        # x_c=all_comp_dict[m]["x_c"]
        # x_p=all_comp_dict[m]["x_p"]
        # u=all_comp_dict[m]["u"]
        # rt=all_comp_dict[m]["rt"]   
        # #nep=all_comp_dict[m]["u"]-all_comp_dict[m]["x"]/all_comp_dict[m]["rt"]  
        
        # delta_x=x-all_comp_dict["Mean"]["x"]
        # delta_x_c=x_c-all_comp_dict["Mean"]["x_c"]
        # delta_x_p=x_p-all_comp_dict["Mean"]["x_p"]
        # delta_u=u-all_comp_dict["Mean"]["u"]
        # delta_rt=rt-all_comp_dict["Mean"]["rt"]
        # #delta_nep=nep-(all_comp_dict["Mean"]["u"]-all_comp_dict["Mean"]["x"]/all_comp_dict["Mean"]["rt"]) 
        
        # # attribution of delta X to delta X_c and delta X_p
        # x_c_contrib=delta_x_c
        # x_p_contrib=-delta_x_p
         
        # # attribution of delta X_c to delta u and delta RT
        # rt_contrib=delta_rt*(u-delta_u/2)
        # u_contrib=delta_u*(rt-delta_rt/2) 
        # rt_u_inter=delta_x_c-rt_contrib-u_contrib#delta_u*delta_rt  

        # # # attribution of delta X_p to delta NEP and delta RT
        # # rt_contrib2=delta_rt*(nep-delta_nep/2)
        # # nep_contrib=delta_nep*(rt-delta_rt/2) 
        # # rt_nep_inter=delta_x_p-rt_contrib2-nep_contrib#delta_nep*delta_rt 

        # # attribution in percents
        # abs_total=abs(x_p_contrib)+abs(rt_contrib)+abs(u_contrib)+abs(rt_u_inter)
        # percent_x_p=abs(x_p_contrib)/abs_total*100
        # percent_rt=abs(rt_contrib)/abs_total*100
        # percent_u=abs(u_contrib)/abs_total*100
        # percent_inter=abs(rt_u_inter)/abs_total*100             
        
        # print ('\033[1m'+m)
                
        # ####### Figure 3 
        # print ('\033[1m'+'Attribution of X differences from the multi-model mean ' +
            # 'to the differences of traceable components')         
        # fig3=plt.figure(figsize=(15,15))
        # axs=fig3.subplots(2,2)
        
        # # contributions timeline
        # ax=axs[0,0]
        # # ax.plot(          
            # # all_comp_dict["Times"], 
            # # rt_contrib,
            # # label="contribution of $\Delta$ R_t",
            # # color="darkorange",
            # # linewidth=0.8,
            # # alpha=0.5,            
        # # )
        # # ax.plot(        
            # # all_comp_dict["Times"], 
            # # u_contrib,
            # # label="contribution of $\Delta$ u",
            # # color="green",
            # # linewidth=0.8,
            # # alpha=0.5,            
        # # )     
        # # ax.plot(        
            # # all_comp_dict["Times"], 
            # # delta_x,
            # # label="$\Delta$ X",
            # # color="black",
        # # )
        # # ax.plot(        
            # # all_comp_dict["Times"], 
            # # x_p_contrib,
            # # label="$\Delta$ X_p",
            # # color="blue",
            # # linewidth=0.8,
            # # alpha=0.5,            
        # # )        
        # # if np.mean(percent_inter) > 0.001:
            # # ax.plot(            
                # # all_comp_dict["Times"], 
                # # #percent_inter,
                # # rt_u_inter,
                # # label="contribution of interaction terms",
                # # color="lightgrey",
                # # linewidth=0.8,
                # # alpha=0.5,                
            # # )    
        # # quadratic trends
        # z = np.polyfit(all_comp_dict["Times"],  x_p_contrib, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_p",
            # color="blue",
        # ) 
        # ax.fill_between(
                # all_comp_dict["Times"],
                # x_p_contrib, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="blue",
                # #linewidth=0.1,
                # alpha=0.2                
                # )         

        # z = np.polyfit(all_comp_dict["Times"],  u_contrib, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_p",
            # color="green",
        # ) 
        # ax.fill_between(
                # all_comp_dict["Times"],
                # u_contrib, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="green",
                # #linewidth=0.1,
                # alpha=0.2                
                # )   
                
        # z = np.polyfit(all_comp_dict["Times"],  rt_contrib, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_p",
            # color="darkorange",
        # )              
        # ax.fill_between(
                # all_comp_dict["Times"],
                # rt_contrib, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="darkorange",
                # #linewidth=0.1,
                # alpha=0.2                
                # )         

        # z = np.polyfit(all_comp_dict["Times"],  rt_u_inter, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_p",
            # color="lightgrey",
        # )         
        # ax.fill_between(
                # all_comp_dict["Times"],
                # rt_u_inter, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="lightgrey",
                # #linewidth=0.1,
                # alpha=0.2                
                # )   

        # z = np.polyfit(all_comp_dict["Times"],  delta_x, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_p",
            # color="black",
        # )         
        # ax.fill_between(
                # all_comp_dict["Times"],
                # delta_x, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="black",
                # #linewidth=0.1,
                # alpha=0.2                
                # )                 
        # z = np.polyfit(all_comp_dict["Times"],  x_c_contrib, 2)
        # p = np.poly1d(z)  
        # ax.plot(        
            # all_comp_dict["Times"], 
            # p(all_comp_dict["Times"]),
            # label="$\Delta$ X_c",
            # color="red",
        # )         
        # ax.fill_between(
                # all_comp_dict["Times"],
                # x_c_contrib, 
                # p(all_comp_dict["Times"]),
                # #label="\u00B12$\sigma$ confidence interval",
                # color="red",
                # #linewidth=0.1,
                # alpha=0.2                
                # )                  

        # #ax.legend()
        # ax.set_title('Contributions over time')
        # #ax.set_ylabel('%')
        # ax.grid()        

        # # bar charts
        
        # positive_delta_X=sum(delta_x[delta_x>0])/len(delta_x)
        # negative_delta_X=sum(delta_x[delta_x<0])/len(delta_x)
        # positive_x_c_contrib=sum(x_c_contrib[x_c_contrib>0])/len(x_c_contrib)
        # negative_x_c_contrib=sum(x_c_contrib[x_c_contrib<0])/len(x_c_contrib)
        # positive_x_p_contrib=sum(x_p_contrib[x_p_contrib>0])/len(x_p_contrib)
        # negative_x_p_contrib=sum(x_p_contrib[x_p_contrib<0])/len(x_p_contrib)        
        # positive_rt_contrib=sum(rt_contrib[rt_contrib>0])/len(rt_contrib)
        # negative_rt_contrib=sum(rt_contrib[rt_contrib<0])/len(rt_contrib)
        # positive_u_contrib=sum(u_contrib[u_contrib>0])/len(u_contrib)
        # negative_u_contrib=sum(u_contrib[u_contrib<0])/len(u_contrib)      
        # positive_rt_u_inter=sum(rt_u_inter[rt_u_inter>0])/len(rt_u_inter)
        # negative_rt_u_inter=sum(rt_u_inter[rt_u_inter<0])/len(rt_u_inter)         
        
        # ax0=axs[0,1]  
        # ax0.set_title('Average contributions')       
        # ax0.axhline(0, color='black', ls='dashed')
        
        # ax0.bar ('$\Delta$ X', positive_delta_X, color="black")
        # ax0.bar ('$\Delta$ X', negative_delta_X, color="black")
        # ax0.bar ('$\Delta$ X_p', positive_x_p_contrib, color="blue")
        # ax0.bar ('$\Delta$ X_p', negative_x_p_contrib, color="blue")
        # ax0.bar ('$\Delta$ X_c', positive_x_c_contrib, color="red")
        # ax0.bar ('$\Delta$ X_c', negative_x_c_contrib, color="red") 
        # ax0.bar ('$\Delta$ RT', positive_rt_contrib, color="darkorange")
        # ax0.bar ('$\Delta$ RT', negative_rt_contrib, color="darkorange")        
        # ax0.bar ('$\Delta$ NPP', positive_u_contrib, color="green")
        # ax0.bar ('$\Delta$ NPP', negative_u_contrib, color="green")
        # ax0.bar ('$\Delta$npp*$\Delta$rt', positive_rt_u_inter, color="lightgrey")
        # ax0.bar ('$\Delta$npp*$\Delta$rt', negative_rt_u_inter, color="lightgrey")        
  
        
        # # if abs(np.mean(combined_contrib)) > 0.001:
           # # ax0.bar ('$\Delta$ U * $\Delta$ RT', np.mean(combined_contrib), color="lightgrey")   
        # ax0.grid() 

        # # pie chart incorrect
        # ax=axs[1,0]
        
        # x_c_in_x=abs(x_c_contrib)/(abs(x_c_contrib)+abs(x_p_contrib))   
        # x_p_in_x=abs(x_p_contrib)/(abs(x_c_contrib)+abs(x_p_contrib))   
        
        # rt_in_xc=abs(rt_contrib)/(abs(rt_contrib)+abs(u_contrib)+abs(rt_u_inter))  
        # u_in_xc=abs(u_contrib)/(abs(rt_contrib)+abs(u_contrib)+abs(rt_u_inter))  
        # inter_in_xc=abs(rt_u_inter)/(abs(rt_contrib)+abs(u_contrib)+abs(rt_u_inter))  
        # print(np.mean(rt_in_xc))
        # print(np.mean(u_in_xc))
        # print(np.mean(inter_in_xc))
        # #print(np.mean(x_p_in_x))

        # percent_rt1=rt_in_xc*x_c_in_x*100
        # percent_u1=u_in_xc*x_c_in_x*100
        # percent_inter1=inter_in_xc*x_c_in_x*100 
        # percent_x_p1=x_p_in_x*100        

        # print(np.mean(percent_rt1))  
        # print(np.mean(percent_u1))
        # print(np.mean(percent_x_p))
        
        # ax.set_title('% Contributions - if applied hierarchically')
          
        # if np.mean(percent_inter) > 0.001:
            # labels = '$\Delta$ RT', '$\Delta$ NPP', '$\Delta$NPP*$\Delta$RT', '$\Delta$ X_p'
            # sizes = [np.mean(percent_rt1), np.mean(percent_u1), 
                # np.mean(percent_inter1), np.mean(percent_x_p1)]
            # ax.pie(sizes, autopct='%1.1f%%', 
                # startangle=90, counterclock=False, colors=("darkorange", "green", "lightgrey", "blue"))
        # else:
            # labels = '$\Delta$ RT','$\Delta$ NPP', '$\Delta$ X_p'
            # sizes = [np.mean(percent_rt1), np.mean(percent_u1), np.mean(percent_x_p1)]
            # ax.pie(sizes, labels=labels, autopct='%1.1f%%',
                # startangle=90, counterclock=False, colors=("darkorange", "green", "blue"))        
        # ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
        # #ax.legend(labels, bbox_to_anchor =(1.42, 2.22)) 
        
        # # pie charts

        # ax1=axs[1,1]
        # ax1.set_title('% Contributions')
          
        # if np.mean(percent_inter) > 0.001:
            # labels = '$\Delta$ RT', '$\Delta$ NPP', '$\Delta$NPP*$\Delta$RT', '$\Delta$ X_p'
            # sizes = [np.mean(percent_rt), np.mean(percent_u), 
                # np.mean(percent_inter), np.mean(percent_x_p)]
            # ax1.pie(sizes, autopct='%1.1f%%', 
                # startangle=90, counterclock=False, colors=("darkorange", "green", "lightgrey", "blue"))
        # else:
            # labels = '$\Delta$ RT','$\Delta$ NPP', '$\Delta$ X_p'
            # sizes = [np.mean(percent_rt), np.mean(percent_u), np.mean(percent_x_p)]
            # ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                # startangle=90, counterclock=False, colors=("darkorange", "green", "blue"))        
        # ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.       
        # ax1.legend(labels, bbox_to_anchor =(1.42, 2.22))
        # plt.show()
 
        # # print("X: "+str(positive_delta_X-negative_delta_X)+
            # # "; RT: "+str(positive_rt_contrib-negative_rt_contrib)+
            # # "; NPP: "+str(positive_u_contrib-negative_u_contrib)+
            # # "; rt_u_inter: "+ str(positive_rt_u_inter-negative_rt_u_inter)+
            # # "; X_p: "+str(positive_x_p_contrib-negative_x_p_contrib)
        # # )   
        # # print("sum check: "+str(
            # # (positive_delta_X-negative_delta_X)-
            # # (positive_rt_contrib-negative_rt_contrib)-
            # # (positive_u_contrib-negative_u_contrib)-
            # # (positive_rt_u_inter-negative_rt_u_inter)-
            # # (positive_x_p_contrib-negative_x_p_contrib)
            # # )
        # # )  
        # plt.show()
        
        # # summation of positive and negative contributions separately
        # pos_delta_x=delta_x.copy(); pos_delta_x[pos_delta_x<0]=0
        # neg_delta_x=delta_x.copy(); neg_delta_x[neg_delta_x>0]=0
        # sum_pos_diff_x+=pos_delta_x
        # sum_neg_diff_x+=neg_delta_x
        
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
        
        
    # return (sum_pos_diff_x, sum_neg_diff_x,
        # sum_pos_cont_rt, sum_neg_cont_rt,
        # sum_pos_cont_u, sum_neg_cont_u,
        # sum_pos_cont_rt_u_inter, sum_neg_cont_rt_u_inter,
        # sum_pos_cont_x_p, sum_neg_cont_x_p,
        # )


# def plot_attribution_model_vs_mean (
    # all_comp_dict,
# ):
    # models=list(all_comp_dict.keys())[:-2]        
    
    # for m in models: 
    
        # print ('\033[1m'+m)    
           
        # x=all_comp_dict[m]["x"]
        # x_c=all_comp_dict[m]["x_c"]
        # x_p=all_comp_dict[m]["x_p"]
        # u=all_comp_dict[m]["u"]
        # rt=all_comp_dict[m]["rt"]   
        # #nep=all_comp_dict[m]["u"]-all_comp_dict[m]["x"]/all_comp_dict[m]["rt"]  
        
        # delta_x=x-all_comp_dict["Mean"]["x"]
        # delta_x_c=x_c-all_comp_dict["Mean"]["x_c"]
        # delta_x_p=x_p-all_comp_dict["Mean"]["x_p"]
        # delta_u=u-all_comp_dict["Mean"]["u"]
        # delta_rt=rt-all_comp_dict["Mean"]["rt"]
        # #delta_nep=nep-(all_comp_dict["Mean"]["u"]-all_comp_dict["Mean"]["x"]/all_comp_dict["Mean"]["rt"]) 
        
        # # attribution of delta X to delta X_c and delta X_p
        # x_c_contrib=delta_x_c
        # x_p_contrib=-delta_x_p
         
        # # attribution of delta X_c to delta u and delta RT
        # rt_contrib=delta_rt*(u-delta_u/2)
        # u_contrib=delta_u*(rt-delta_rt/2) 
        # rt_u_inter=delta_x_c-rt_contrib-u_contrib           

        # # separation of positive and negative contributions
          
        # pos_delta_x=delta_x.copy(); pos_delta_x[pos_delta_x<0]=0
        # neg_delta_x=delta_x.copy(); neg_delta_x[neg_delta_x>0]=0
 
        # pos_cont_rt=rt_contrib.copy(); pos_cont_rt[pos_cont_rt<0]=0
        # neg_cont_rt=rt_contrib.copy(); neg_cont_rt[neg_cont_rt>0]=0

        # pos_cont_u=u_contrib.copy(); pos_cont_u[pos_cont_u<0]=0
        # neg_cont_u=u_contrib.copy(); neg_cont_u[neg_cont_u>0]=0

        # pos_cont_rt_u_inter=rt_u_inter.copy(); pos_cont_rt_u_inter[pos_cont_rt_u_inter<0]=0
        # neg_cont_rt_u_inter=rt_u_inter.copy(); neg_cont_rt_u_inter[neg_cont_rt_u_inter>0]=0

        # pos_cont_x_p=x_p_contrib.copy(); pos_cont_x_p[pos_cont_x_p<0]=0
        # neg_cont_x_p=x_p_contrib.copy(); neg_cont_x_p[neg_cont_x_p>0]=0
        
        # plot_attribution (
            # times=all_comp_dict["Times"],
            # delta_x_pos=pos_delta_x,
            # delta_x_neg=neg_delta_x,
            # rt_contrib_pos=pos_cont_rt,
            # rt_contrib_neg=neg_cont_rt,
            # u_contrib_pos=pos_cont_u,
            # u_contrib_neg=neg_cont_u,
            # rt_u_inter_pos=pos_cont_rt_u_inter,
            # rt_u_inter_neg=neg_cont_rt_u_inter,
            # x_p_contrib_pos=pos_cont_x_p,
            # x_p_contrib_neg=neg_cont_x_p,
        # )  

 
# def plot_components_combined(
    # model_names,  # dictionary (folder name : model name)
    # test_arg_list,  # a list of test_args from all models involved
    # var_names,  # dictionary (trace_tuple name : descriptive name)
    # delta_t_val,  # model time step
    # model_cols,  # dictionary (folder name :color)
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
# ):
    # if (part < 0) | (part > 1):
        # raise Exception(
            # "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        # )
    # model_folders = [(k) for k in model_names]
    # variables = [(k) for k in var_names]
    # n = len(variables)

    # fig = plt.figure(figsize=(17, n * 8))
    # axs = fig.subplots(n, 1)

    # for i, name in enumerate(variables):
        # k = 0
        # for mf in model_folders:
            # itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_overlap(test_arg_list, delta_t_val)
            # )
            # # if we do not want the whole interval but look at a smaller part to observe the dynamics
            # # start,stop = start_min, int(start_min+(stop_max-start_min)*part)
            # start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
            # times = (
                # times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
                # / days_per_year()
            # )
            # print("Plotting " + str(name) + " for " + str(mf) + " model")
            # vals = itr[start:stop]
            # ax = axs[i]
            # times_for_plot = avg_timeline(times, averaging)
            # vals_for_plot = avg_timeline(vals.__getattribute__(name), averaging)
            # if name == "x":
                # vals_for_plot = (
                    # vals_for_plot * 148940000 * 1000000 * 0.000000000001
                # )  # convert to global C in Pg
            # ax.plot(
                # times_for_plot,
                # vals_for_plot,
                # label=model_names[mf] + " - " + name,
                # color=model_cols[mf],
            # )
            # if name == "x":  # we plot x together with x_c
                # ax.plot(
                    # times_for_plot,
                    # avg_timeline(vals.x_c, averaging)
                    # * 148940000
                    # * 1000000
                    # * 0.000000000001,
                    # label=model_names[mf] + " - x_c",
                    # color=model_cols[mf],
                    # linestyle="dashed",
                # )
            # if name == "x_p":  # 0 line for X_p
                # ax.plot(
                    # times_for_plot,
                    # np.zeros_like(avg_timeline(times, averaging)),
                    # color="black",
                    # linestyle="dotted",
                    # alpha=0.5,
                # )
            # k += 1
        # ax.legend()
        # ax.set_title(var_names[name])
        # ax.grid()


# def plot_x_xc(
    # model_names,  # dictionary (folder name : model name)
    # test_arg_list,  # a list of test_args from all models involved
    # delta_t_val,  # model time step
    # model_cols,  # dictionary (folder name :color)
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
    # overlap=True,  # compute overlapping timeframe or plot whole duration for all models
# ):
    # if (part < 0) | (part > 1):
        # raise Exception(
            # "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        # )
    # model_folders = [(k) for k in model_names]
    # fig = plt.figure(figsize=(17, 8))
    # ax = fig.subplots(1, 1)
    # k = 0
    # for mf in model_folders:
        # itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        # if overlap == True:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_overlap(test_arg_list, delta_t_val)
            # )
        # else:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_full(test_arg_list, delta_t_val)
            # )
        # # if we do not want the whole interval but look at a smaller part to observe the dynamics
        # start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        # times = (
            # times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            # / days_per_year()
        # )
        # vals = itr[start:stop]
        # # print("vals.x")
        # # print(vals.x[0:240])
        # ax.plot(
            # avg_timeline(times, averaging),
            # avg_timeline(vals.x, averaging)
            # * 148940000
            # * 1000000
            # * 0.000000000001,  # convert to global C in Gt
            # label=model_names[mf] + " - X",
            # color=model_cols[mf],
        # )
        # ax.plot(
            # avg_timeline(times, averaging),
            # avg_timeline(vals.x_c, averaging)
            # * 148940000
            # * 1000000
            # * 0.000000000001,  # convert to global C in Gt
            # label=model_names[mf] + " - X_c",
            # color=model_cols[mf],
            # linestyle="dashed",
        # )
        # k += 1
    # ax.legend()
    # ax.set_title("Total Carbon (X) and Carbon Storage Capacity (X_c)")
    # ax.set_ylabel("Gt C")
    # ax.grid()


# # change of X since the start of simulation
# def plot_normalized_x(
    # model_names,  # dictionary (folder name : model name)
    # test_arg_list,  # a list of test_args from all models involved
    # delta_t_val,  # model time step
    # model_cols,  # dictionary (folder name :color)
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
    # overlap=True,  # compute overlapping timeframe or plot whole duration for all models
# ):
    # if (part < 0) | (part > 1):
        # raise Exception(
            # "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        # )
    # model_folders = [(k) for k in model_names]

    # fig = plt.figure(figsize=(17, 8))
    # ax = fig.subplots(1, 1)
    # k = 0
    # for mf in model_folders:
        # itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        # if overlap == True:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_overlap(test_arg_list, delta_t_val)
            # )
        # else:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_full(test_arg_list, delta_t_val)
            # )
        # # if we do not want the whole interval but look at a smaller part to observe the dynamics
        # start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        # times = (
            # times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            # / days_per_year()
        # )
        # vals = itr[start:stop]
        # vals_for_plot = (
            # avg_timeline(vals.x, averaging) * 148940000 * 1000000 * 0.000000000001,
        # )  # convert to global C in Gt
        # vals_array = vals_for_plot[0]
        # vals_for_plot_norm = (
            # (vals_array - vals_array[0]) / vals_array[0] * 100
        # )  # convert to % change
        # ax.plot(
            # avg_timeline(times, averaging),
            # vals_for_plot_norm,
            # label=model_names[mf],
            # color=model_cols[mf],
        # )
        # k += 1
    # ax.legend()
    # ax.set_title("Total Carbon (X) change since the start of the simulation")
    # ax.set_ylabel("% change")
    # ax.grid()


# def plot_normalized_xc(
    # model_names,  # dictionary (folder name : model name)
    # test_arg_list,  # a list of test_args from all models involved
    # delta_t_val,  # model time step
    # model_cols,  # dictionary (folder name :color)
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
    # overlap=True,  # compute overlapping timeframe or plot whole duration for all models
# ):
    # if (part < 0) | (part > 1):
        # raise Exception(
            # "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        # )
    # model_folders = [(k) for k in model_names]

    # fig = plt.figure(figsize=(17, 8))
    # ax = fig.subplots(1, 1)
    # k = 0
    # for mf in model_folders:
        # itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        # if overlap == True:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_overlap(test_arg_list, delta_t_val)
            # )
        # else:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_full(test_arg_list, delta_t_val)
            # )
        # # if we do not want the whole interval but look at a smaller part to observe the dynamics
        # start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        # times = (
            # times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            # / days_per_year()
        # )
        # vals = itr[start:stop]
        # vals_for_plot = (
            # avg_timeline(vals.x_c, averaging) * 148940000 * 1000000 * 0.000000000001,
        # )  # convert to global C in Gt
        # vals_array = vals_for_plot[0]
        # vals_for_plot_norm = (
            # (vals_array - vals_array[0]) / vals_array[0] * 100
        # )  # convert to % change
        # ax.plot(
            # avg_timeline(times, averaging),
            # vals_for_plot_norm,
            # label=model_names[mf],
            # color=model_cols[mf],
        # )
        # k += 1
    # ax.legend()
    # ax.set_title("Carbon Storage Capacity (X) change since the start of the simulation")
    # ax.set_ylabel("% change")
    # ax.grid()


# def plot_xp(
    # model_names,  # dictionary (folder name : model name)
    # test_arg_list,  # a list of test_args from all models involved
    # delta_t_val,  # model time step
    # model_cols,  # dictionary (folder name :color)
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
    # overlap=True,  # compute overlapping timeframe or plot whole duration for all models
# ):
    # if (part < 0) | (part > 1):
        # raise Exception(
            # "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        # )
    # model_folders = [(k) for k in model_names]

    # fig = plt.figure(figsize=(17, 8))
    # ax = fig.subplots(1, 1)
    # k = 0
    # for mf in model_folders:
        # itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        # if overlap == True:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_overlap(test_arg_list, delta_t_val)
            # )
        # else:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_full(test_arg_list, delta_t_val)
            # )
        # # if we do not want the whole interval but look at a smaller part to observe the dynamics
        # start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        # times = (
            # times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            # / days_per_year()
        # )
        # vals = itr[start:stop]
        # ax.plot(
            # avg_timeline(times, averaging),
            # avg_timeline(vals.x_p, averaging)
            # * 148940000
            # * 1000000
            # * 0.000000000001,  # convert to global C in Gt
            # label=model_names[mf],
            # color=model_cols[mf],
        # )
        # k += 1
    # ax.legend()
    # ax.set_title("Carbon Storage Potential (X_p)")
    # ax.set_ylabel("Gt C")
    # ax.grid()


# def plot_u(
    # model_names,  # dictionary (folder name : model name)
    # test_arg_list,  # a list of test_args from all models involved
    # delta_t_val,  # model time step
    # model_cols,  # dictionary (folder name :color)
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
    # overlap=True,  # compute overlapping timeframe or plot whole duration for all models
# ):
    # if (part < 0) | (part > 1):
        # raise Exception(
            # "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        # )
    # model_folders = [(k) for k in model_names]

    # fig = plt.figure(figsize=(17, 8))
    # ax = fig.subplots(1, 1)
    # k = 0
    # for mf in model_folders:
        # itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        # if overlap == True:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_overlap(test_arg_list, delta_t_val)
            # )
        # else:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_full(test_arg_list, delta_t_val)
            # )
        # # if we do not want the whole interval but look at a smaller part to observe the dynamics
        # start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        # times = (
            # times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            # / days_per_year()
        # )
        # vals = itr[start:stop]
        # ax.plot(
            # avg_timeline(times, averaging),
            # avg_timeline(vals.u, averaging)
            # * 148940000
            # * 1000000
            # * 0.000000000001,  # convert to global C in Gt
            # label=model_names[mf],
            # color=model_cols[mf],
        # )
        # k += 1
    # ax.legend()
    # ax.set_title("Carbon Input (NPP)")
    # ax.set_ylabel("Gt C / day")
    # ax.grid()


# # u change since the start of simulation
# def plot_normalized_u(
    # model_names,  # dictionary (folder name : model name)
    # test_arg_list,  # a list of test_args from all models involved
    # delta_t_val,  # model time step
    # model_cols,  # dictionary (folder name :color)
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
    # overlap=True,  # compute overlapping timeframe or plot whole duration for all models
# ):
    # if (part < 0) | (part > 1):
        # raise Exception(
            # "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        # )
    # model_folders = [(k) for k in model_names]

    # fig = plt.figure(figsize=(17, 8))
    # ax = fig.subplots(1, 1)
    # k = 0
    # for mf in model_folders:
        # itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        # if overlap == True:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_overlap(test_arg_list, delta_t_val)
            # )
        # else:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_full(test_arg_list, delta_t_val)
            # )
        # # if we do not want the whole interval but look at a smaller part to observe the dynamics
        # start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        # times = (
            # times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            # / days_per_year()
        # )
        # vals = itr[start:stop]
        # vals_for_plot = (
            # avg_timeline(vals.u, averaging) * 148940000 * 1000000 * 0.000000000001,
        # )  # convert to global C in Gt
        # vals_array = vals_for_plot[0]
        # vals_for_plot_norm = (
            # (vals_array - vals_array[0]) / vals_array[0] * 100
        # )  # convert to % change
        # ax.plot(
            # avg_timeline(times, averaging),
            # vals_for_plot_norm,
            # label=model_names[mf],
            # color=model_cols[mf],
        # )
        # k += 1
    # ax.legend()
    # ax.set_title("Carbon Input (NPP) change since the start of the simulation")
    # ax.set_ylabel("% change")
    # ax.grid()


# def plot_rt(
    # model_names,  # dictionary (folder name : model name)
    # test_arg_list,  # a list of test_args from all models involved
    # delta_t_val,  # model time step
    # model_cols,  # dictionary (folder name :color)
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
    # overlap=True,  # compute overlapping timeframe or plot whole duration for all models
# ):
    # if (part < 0) | (part > 1):
        # raise Exception(
            # "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        # )
    # model_folders = [(k) for k in model_names]

    # fig = plt.figure(figsize=(17, 8))
    # ax = fig.subplots(1, 1)
    # k = 0
    # for mf in model_folders:
        # itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        # if overlap == True:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_overlap(test_arg_list, delta_t_val)
            # )
        # else:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_full(test_arg_list, delta_t_val)
            # )
        # # if we do not want the whole interval but look at a smaller part to observe the dynamics
        # start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        # times = (
            # times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            # / days_per_year()
        # )
        # vals = itr[start:stop]
        # ax.plot(
            # avg_timeline(times, averaging),
            # avg_timeline(vals.rt, averaging) / days_per_year(),  # convert to years
            # label=model_names[mf],
            # color=model_cols[mf],
        # )
        # k += 1
    # ax.legend()
    # ax.set_title("Equilibrium Residense Time (RT)")
    # ax.set_ylabel("Years")
    # ax.grid()


# def plot_normalized_rt(
    # model_names,  # dictionary (folder name : model name)
    # test_arg_list,  # a list of test_args from all models involved
    # delta_t_val,  # model time step
    # model_cols,  # dictionary (folder name :color)
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
    # overlap=True,  # compute overlapping timeframe or plot whole duration for all models
# ):
    # if (part < 0) | (part > 1):
        # raise Exception(
            # "Invalid partitioning in plot_components_combined: use part between 0 and 1"
        # )
    # model_folders = [(k) for k in model_names]

    # fig = plt.figure(figsize=(17, 8))
    # ax = fig.subplots(1, 1)
    # k = 0
    # for mf in model_folders:
        # itr = traceability_iterator_instance(mf, test_arg_list[k], delta_t_val)
        # if overlap == True:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_overlap(test_arg_list, delta_t_val)
            # )
        # else:
            # start_min, stop_max = min_max_index(
                # test_arg_list[k],
                # delta_t_val,
                # *t_min_tmax_full(test_arg_list, delta_t_val)
            # )
        # # if we do not want the whole interval but look at a smaller part to observe the dynamics
        # start, stop = int(stop_max - (stop_max - start_min) * part), stop_max
        # times = (
            # times_in_days_aD(test_arg_list[k], delta_t_val)[start:stop]
            # / days_per_year()
        # )
        # vals = itr[start:stop]
        # vals_for_plot = (
            # avg_timeline(vals.rt, averaging) / days_per_year()
        # )  # convert to years
        # vals_array = vals_for_plot
        # vals_for_plot_norm = (
            # (vals_array - vals_array[0]) / vals_array[0] * 100
        # )  # convert to % change
        # ax.plot(
            # avg_timeline(times, averaging),
            # vals_for_plot_norm,
            # label=model_names[mf],
            # color=model_cols[mf],
        # )
        # k += 1
    # ax.legend()
    # ax.set_title(
        # "Equilibrium Residense Time (RT) change since the start of the simulation"
    # )
    # ax.set_ylabel("% change")
    # ax.grid()


# from scipy.interpolate import interp1d, splprep

# function to compute a difference between traceable companents of 2 models
# Since the two models do not necessarily share the same point in time and not
# even the same stepsize or number of steps we compute interpolating functions
# to make them comparable


# def timeline_diff(
    # name,  # name of variable in the tuple)
    # mf_1,  # model folder (1st model)
    # times_1,  # array of time steps (1st model)
    # vals_1,  # tuple of values for each time step (1st model)
    # mf_2,  # model folder (2nd model)
    # times_2,  # array of time steps (2nd model)
    # vals_2,  # tuple of values for each time step (2nd model)
# ):
    # f1 = interp1d(times_1, vals_1.__getattribute__(name))
    # f2 = interp1d(times_2, vals_2.__getattribute__(name))
    # # chose the interval covered by both to avoid extrapolation
    # start = max(times_1.min(), times_2.min())
    # stop = min(times_1.max(), times_2.max())
    # nstep = min(len(times_1), len(times_2))
    # times = np.linspace(start, stop, nstep)

    # diff = f1(times) - f2(times)
    # # diff=vals_1.__getattribute__(name)-vals_2.__getattribute__(name) # if time step is same, this should work instead of interpolation
    # return (diff, times)


# # plotting differences between traceable companents of 2 models
# def plot_diff(
    # model_names,  # dictionary (folder name : model name)
    # var_names,  # dictionary (trace_tuple name : descriptive name)
    # delta_t_val,  # model time step
    # part,  # 0<part<1 to plot only a part of the whole timeling, e.g. 1 (whole timeline) or 0.1 (10%)
    # averaging,  # number of iterator steps over which to average results. 1 for no averaging
# ):
    # if part < 0 | part > 1:
        # raise Exception("Invalid partitioning in plot_diff: use part between 0 and 1")
    # model_folders = [(k) for k in model_names]
    # variables = [(k) for k in var_names]

    # n = int(len(variables) * (len(model_folders) * (len(model_folders) - 1) / 2))
    # # yr=int(365/delta_t_val)

    # fig = plt.figure(figsize=(17, n * 8))
    # axs = fig.subplots(n, 1)
    # plot_number = 0
    # for i, name in enumerate(variables):
        # for j, mf_1 in enumerate(model_folders[0 : len(model_folders) - 1]):
            # for mf_2 in model_folders[j + 1 : len(model_folders)]:
                # start_min_1, stop_max_1 = min_max_index(
                    # mf_1, delta_t_val, *t_min_tmax_overlap([mf_1, mf_2], delta_t_val)
                # )
                # # if we do not want the whole interval but look at a smaller part to observe the dynamics
                # start_1, stop_1 = start_min_1, int(
                    # start_min_1 + (stop_max_1 - start_min_1) * part
                # )
                # itr_1 = traceability_iterator_instance(mf_1, delta_t_val)
                # # parts_1=partitions(start_1,stop_1,yr)
                # times_1 = (
                    # times_in_days_aD(mf_1, delta_t_val)[start_1:stop_1]
                    # / days_per_year()
                # )
                # vals_1 = itr_1[start_1:stop_1]

                # start_min_2, stop_max_2 = min_max_index(
                    # mf_2, delta_t_val, *t_min_tmax_overlap([mf_2, mf_2], delta_t_val)
                # )
                # # if we do not want the whole interval but look at a smaller part to observe the dynamics
                # start_2, stop_2 = start_min_2, int(
                    # start_min_2 + (stop_max_2 - start_min_2) * part
                # )
                # itr_2 = traceability_iterator_instance(mf_2, delta_t_val)
                # # parts_2=partitions(start_2,stop_2,yr)
                # times_2 = (
                    # times_in_days_aD(mf_2, delta_t_val)[start_2:stop_2]
                    # / days_per_year()
                # )
                # vals_2 = itr_1[start_1:stop_1]

                # print("Plotting " + str(plot_number + 1) + " out of " + str(n))
                # diff, times = timeline_diff(
                    # name, mf_1, times_1, vals_1, mf_2, times_2, vals_2
                # )
                # ax = axs[plot_number]
                # ax.plot(times, diff, color="black")
                # ax.set_title(
                    # "Delta {0} for {1}-{2}".format(
                        # name, model_names[mf_1], model_names[mf_2]
                    # )
                # )
                # plot_number += 1