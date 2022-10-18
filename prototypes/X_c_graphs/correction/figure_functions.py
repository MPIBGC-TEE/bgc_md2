import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb
import plot_functions as pf



def combined_timelines_and_2d_phase_space(
    size,
    intervals,
    rectangle_colors,
    linestyles,
    color_dict,
    label_dict,
    state_vector,
    x_dim,
    y_dim,
    z_dim,
    ss,
    t_eval,
    Is,
    Xis,
    Xs,
    X_cs,
    X_dots
):
    fig = plt.figure(
        figsize=size,
        constrained_layout=True
    )
    #outer_gridspec = fig.add_gridspec(
    #    2, #5,
    #    1,
    #    #width_ratios=(1, 1,1),
    #    height_ratios=(1, 1),#, 1, 2),#,8),
    #    #left=0.1,
    #    #right=0.9,
    #    #bottom=0.1,
    #    #top=0.9,
    #    #wspace=0.05,
    #    #hspace=0.05
    #)
    #timelines=fig.add_subfigure(outer_gridspec[0, 0])
    #phase_portraits = fig.add_subfigure(outer_gridspec[1, 0])
    timelines, phase_portraits = fig.subfigures(
        2, 1,
        wspace=0.05,
        hspace=0.05
    )
    # facecolor='0.75'
    facecolor='1'  # 1 white , 0 black

    np.random.seed(19680808)

    timelines.set_facecolor(facecolor)
    timelines.suptitle('Timelines', fontsize='x-large')
    axs_tl = timelines.subplots(3, 1)#, sharey=True)
    ax=axs_tl[0]
    ax.set_title("Inputs")
    ax.set_xlabel("time")
    ax.set_ylabel("")
    for dim in (x_dim,y_dim):
        ax.plot(
            t_eval,
            Is[:,dim],
            color=color_dict["I"],
            label=f"$I_{{{state_vector[dim]}}}$",
            linestyle=linestyles[dim]
        )
    ax.legend()
    ns=[1]
    pf.plot_timerectangles(ax, intervals[ns], rectangle_colors[ns], t_eval, ss )
    
    ######################################## 
    ax=axs_tl[1]
    ax.set_title(f"{label_dict['Xi']}")
    for dim in (x_dim,y_dim):
        ax.plot(
            t_eval,
            Xis[:,dim],
            color=color_dict["Xi"],
            label=f"$\\xi_{{{state_vector[dim]}}}$",
            linestyle=linestyles[dim]
        )
    ax.legend()
    ns=[2]
    pf.plot_timerectangles(ax, intervals[ns], rectangle_colors[ns], t_eval, ss )
    
    y_min, y_max=ax.get_ybound()
    ax=axs_tl[2]
    ######################################## 
    ax.set_title(f"{label_dict['X']}")
    ax.set_xlabel("time")
    ax.set_ylabel("")
    for dim in (x_dim,y_dim):
        ax.plot(
            t_eval,
            Xs[:,dim],
            color=color_dict["X"],
            label=state_vector[dim],
            linestyle=linestyles[dim]
        )
        ax.plot(
            t_eval,
            X_cs[:,dim],
            color=color_dict["X_c"],
            label=f"$\\mathbf{{X_c}}_{{{state_vector[dim]}}}$",
            linestyle=linestyles[dim]
        )
    
    ax.legend()
    ns = [0,1,2]
    pf.plot_timerectangles(ax, intervals[ns], rectangle_colors[ns], t_eval, ss )

    #######################################################################
    axs_pp = phase_portraits.subplots(1, 3)
    for i,iv in enumerate(intervals[ns]):
        ax = axs_pp[i]
        pf.plot_vector2d_single(
            ax,
            iv[0],
            iv[1],
            ss,
            color_dict,
            label_dict,
            rectangle_colors[i],
            state_vector,
            Xs,
            X_cs,
            X_dots
        )
        ax.set_title(f"timeinterval {i}.")
        ax.legend()
        if i ==0:
        # if True:
            ax.set_xlabel(f"$x={state_vector[x_dim]}$")
            ax.set_ylabel(f"$y={state_vector[y_dim]}$")
    
    phase_portraits.set_facecolor(facecolor)
    #phase_portraits.colorbar(pc, shrink=0.6, ax=axsRight)
    phase_portraits.suptitle('Phase Portraits', fontsize='x-large')
    
    fig.suptitle('Figure suptitle', fontsize='xx-large')
    fig.suptitle(
        f'Geometric Relationship between {label_dict["X"]} and {label_dict["X_c"]} for changing Input and {label_dict["Xi"]}',
        fontsize="xx-large"
    )
    
    return fig

def plot_vector_3d_X_Xc_time_z(
    size,
    intervals,
    rectangle_colors,
    color_dict,
    label_dict,
    state_vector,
    x_dim,
    y_dim,
    z_dim,
    ss,
    t_eval,
    Is,
    Xis,
    Xs,
    X_cs,
    X_dots,
    X_ps
):
    f = plt.figure( figsize=size)  
    ax=f.add_subplot(111, projection="3d", box_aspect=[1,1,np.sqrt(2)])
    plot_params = {"elev": 20, "azim": -110}
    ax.view_init(**plot_params)
    pf.plot_time_z(
        ax,
        color_dict,
        label_dict,
        state_vector,
        x_dim,
        y_dim,
        ss,
        t_eval,
        Xs,
        X_cs,
        X_ps,
        X_dots,
    )
    pf.plot_timerectangles_3d(ax, intervals, rectangle_colors, t_eval, ss)
    return f


def plot_vector_3d_X_Xc_time_x(
    size,
    intervals,
    rectangle_colors,
    color_dict,
    label_dict,
    state_vector,
    x_dim,
    y_dim,
    z_dim,
    ss,
    t_eval,
    Is,
    Xis,
    Xs,
    X_cs,
    X_dots,
    X_ps
):
    f = plt.figure( figsize=size)  
    ax=f.add_subplot(111, projection="3d", box_aspect=[1,1,np.sqrt(2)])
    plot_params = {"elev": 20, "azim": -120}
    ax.view_init(**plot_params)
    pf.plot_time_x(
        ax,
        color_dict,
        label_dict,
        state_vector,
        y_dim,
        z_dim,
        ss,
        t_eval,
        Xs,
        X_cs,
        X_ps,
        X_dots,
    )
    return f


def plot_vector_3d_I_Xi_X_Xc_time_z(
    size,
    intervals,
    rectangle_colors,
    color_dict,
    label_dict,
    state_vector,
    x_dim,
    y_dim,
    z_dim,
    ss,
    t_eval,
    Is,
    Xis,
    Xs,
    X_cs,
    X_dots,
    X_ps
):
    # figure 3d Input Xi and X
    f = plt.figure( figsize=(16,24))
    gs = f.add_gridspec(
       1,
       3,
       width_ratios=(1, 1, 1),
       #height_ratios=(1, 1, 1),
       left=0.1,
       right=0.9,
       #bottom=0.1,
       #top=0.9,
       wspace=0.01,
       hspace=0.01
    )
    plot_params = {"elev": 20, "azim": 120}
    ax = f.add_subplot(
       gs[0, 0],
       projection="3d",
       box_aspect=(1, 1, 4),
    )
    ax.view_init(**plot_params)
    #pf.ticks_test(ax)I
    pf.plot_Inputs_time_z(
       ax,
       color_dict,
       label_dict,
       state_vector,
       x_dim,
       y_dim,
       ss,
       t_eval,
       Is
    )
    ax = f.add_subplot(
       gs[0, 1],
       projection="3d",
       box_aspect=(1, 1, 4),
    )
    ax.view_init(**plot_params)
    pf.plot_xi_time_z(
       ax,
       color_dict,
       label_dict,
       state_vector,
       x_dim,
       y_dim,
       ss,
       t_eval,
       Xis
    )
    ax = f.add_subplot(
       gs[0, 2],
       projection="3d",
       box_aspect=(1, 1, 4),
    )
    ax.view_init(**plot_params)
    # plot_vectors(mvs,param_dict,X_0,t_eval)
    pf.plot_time_z(
       ax,
       color_dict,
       label_dict,
       state_vector,
       x_dim,
       y_dim,
       ss,
       t_eval,
       Xs,
       X_cs,
       X_ps,
       X_dots,
    )
    return f

#################################################################################
def plot_vector_3d_I_Xi_X_Xc_time_x(
    size,
    intervals,
    rectangle_colors,
    color_dict,
    label_dict,
    state_vector,
    x_dim,
    y_dim,
    z_dim,
    ss,
    t_eval,
    Is,
    Xis,
    Xs,
    X_cs,
    X_dots,
    X_ps
):
    f = plt.figure(figsize=(8,12), tight_layout=True)
    # Add a gridspec with three rows and two columns and a ratio of 1 to 4 between
    gs = f.add_gridspec(
       3,
       1,
       #width_ratios=(4, 1),
       height_ratios=(1, 1, 1),
       #left=0.1,
       #right=0.9,
       #bottom=0.1,
       #top=0.9,
       #wspace=0.05,
       #hspace=0.05
    )
    # Create the Axes.
    ax = f.add_subplot(
       gs[0, 0],
       projection="3d",
       box_aspect=(4, 1, 1),
    )
    plot_params = {"elev": 20, "azim": -120}
    ax.view_init(**plot_params)
    pf.plot_Inputs_time_x(
       ax,
       color_dict,
       label_dict,
       state_vector,
       y_dim,
       z_dim,
       ss,
       t_eval,
       Is,
    )
    ax.legend()
    
    ax= f.add_subplot(
       gs[1, 0],
       projection="3d",
       box_aspect=(4, 1, 1),
    )
    ax.view_init(**plot_params)
    pf.plot_time_x(
       ax,
       color_dict,
       label_dict,
       state_vector,
       y_dim,
       z_dim,
       ss,
       t_eval,
       Xs,
       X_cs,
       X_ps,
       X_dots,
    )
    return f
    
