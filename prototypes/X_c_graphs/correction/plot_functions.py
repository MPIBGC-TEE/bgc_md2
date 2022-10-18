import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

###############################################################################
def plot_vector2d_single(
    ax,
    start,
    stop,
    #n: int,  # number of subdivisions of the curve
    # p:  int,  # subdivision to be plotted
    ss: int,  # stepsize (number of points before next arrow is drawn)
    color_dict: dict,  # color dictionary
    label_dict: dict,  # tex for labels
    rectangle_color,
    state_vector,
    Xs: np.ndarray,
    X_cs: np.ndarray,
    X_dots: np.ndarray,
):
    def plot_segment(ax, xs, ys, key, avoid_label=False):
        if avoid_label:
            ax.plot(xs, ys, color=color_dict[key])
        else:
            ax.plot(xs, ys, color=color_dict[key], label=label_dict[key])
    
    mutation_scale = 15
    linewidth = 3
    ax.axis("equal")
    x_dim = 0
    y_dim = 1
    xs = Xs[:, x_dim]
    ys = Xs[:, y_dim]

    x_cs = X_cs[:, x_dim]
    y_cs = X_cs[:, y_dim]

    x_dots = X_dots[:, x_dim]
    y_dots = X_dots[:, y_dim]

    x_ccol = color_dict["X_c"]
    sl = slice(start * ss, (start + 1) * ss)
    nt = len(Xs[:, 0])
    s = slice(0, nt, ss)
    # plot the first interval between arrows wiht labels for the legend
    key="X"
    ax.plot(
        xs[sl],
        ys[sl],
        color=color_dict[key],
        label=label_dict[key],
        linewidth=linewidth
    )
    
    key="X_c"
    ax.plot(
        x_cs[sl],
        y_cs[sl],
        color=color_dict[key],
        label=label_dict[key],
        linewidth=linewidth
    )

    arrow_base = (xs[s][start], ys[s][start])
    arrow_head_xp = (x_cs[s][start], y_cs[s][start])
    arrow = mpatches.FancyArrowPatch(
        arrow_base,
        (xs[s][start] + x_dots[s][start], ys[s][start] + y_dots[s][start]),
        mutation_scale=mutation_scale,
        color=color_dict["X_dot"],
        label=label_dict["X_dot"],
    )
    ax.add_patch(arrow)
    arrow = mpatches.FancyArrowPatch(
        # (xs[s][i], ys[s][i]),
        arrow_base,
        arrow_head_xp,
        mutation_scale=mutation_scale,
        color=color_dict["X_p"],
        label=label_dict["X_p"],
    )
    ax.add_patch(arrow)
    for i in range(start + 1, stop):
        sl = slice(i * ss, (i + 1) * ss)
        #  trajectory for X
        arrow_base = (xs[s][i], ys[s][i])
        arrow_head_xp = (x_cs[s][i], y_cs[s][i])
        key = "X"
        ax.plot(
            xs[sl],
            ys[sl],
            color=color_dict[key],
            linewidth=linewidth
        )
        ax.plot(*arrow_base, "*", color=color_dict["X"])
        # trajectory for X_c
        key = "X_c"
        ax.plot(
            x_cs[sl],
            y_cs[sl],
            color=color_dict[key],
            linewidth=linewidth
        )
        ax.plot(
            x_cs[s][i],
            y_cs[s][i],
            "*",
            color=x_ccol,
            # label=label_dict["X_c"]
        )
        arrow = mpatches.FancyArrowPatch(
            # (xs[s][i], ys[s][i]),
            arrow_base,
            arrow_head_xp,
            mutation_scale=mutation_scale,
            color=color_dict["X_p"],
            # label=label_dict["X_p"]
        )
        ax.add_patch(arrow)
        arrow = mpatches.FancyArrowPatch(
            arrow_base,
            (xs[s][i] + x_dots[s][i], ys[s][i] + y_dots[s][i]),
            mutation_scale=mutation_scale,
            color=color_dict["X_dot"],
            # label=label_dict["X_dot"]
        )
        ax.add_patch(arrow)



    #from IPython import embed;embed()
    x_min, x_max = ax.get_xbound()
    y_min, y_max = ax.get_ybound()
    #poly_xs = [x_min, x_min, x_max, x_max, x_min]
    #poly_ys = [y_min, y_max, y_max, y_min, y_min]
    poly_xs = [x_min, x_min, x_min, x_min, x_min]
    poly_ys = [y_min, y_min, y_min, y_min, y_min]
    ax.fill(
        poly_xs, poly_ys,
        color=rectangle_color,
        alpha=.1,
        label="time interval"
    )

###############################################################################
def plot_vector2d(
    f: plt.figure,
    n: int,  # number of subdivisions of the curve
    p_range: range,  # subdivisions to be plotted
    ss: int,  # stepsize (number of points before next arrow is drawn)
    color_dict: dict,  # color dictionary
    label_dict: dict,  # tex for labels
    state_vector,
    Xs: np.ndarray,
    X_cs: np.ndarray,
    X_dots: np.ndarray,
):

    nt = len(Xs[:, 0])
    s = slice(0, nt, ss)
    l = len(Xs[s, 0])  # number of arrows
    start = int(p * l / n)
    stop = int(min((p + 1) * l / n, l))


    for p in p_range:
        ax = f.add_subplot(len(p_range), 1, p + 1 - p_range[0])
        ax.axis("equal")
        plot_vector2d_single(
            ax,
            start,
            stop,
            ss,
            color_dict,
            label_dict,
            state_vector,
            Xs,
            X_cs,
            X_dots
        )


###############################################################################
def plot_time_z(
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
):
    nt = len(t_eval)
    s = slice(0, nt, ss)
    # scale_factor = np.max(t_eval)/np.max(X_ps)/4
    scale_factor = 1
    # arrowhaeds are plotted data space and look funny if the aspect ratios
    # between the axes are very different (not the box aspect)

    #ax.invert_xaxis()
    #ax.invert_yaxis()
    # trajectory for X
    xs = scale_factor * Xs[:, x_dim]
    ys = scale_factor * Xs[:, y_dim]
    ax.plot3D(xs, ys, t_eval, label=label_dict["X"], color=color_dict["X"])
    ax.scatter(xs[s], ys[s], t_eval[s], s=10, color=color_dict["X"])

    x_cs = scale_factor * X_cs[:, x_dim]
    y_cs = scale_factor * X_cs[:, y_dim]
    ax.plot3D(x_cs, y_cs, t_eval, label=label_dict["X_c"], color=color_dict["X_c"])
    ax.scatter(x_cs[s], y_cs[s], t_eval[s], s=10, color=color_dict["X_c"])

    x_dots = scale_factor * X_dots[:, x_dim]
    y_dots = scale_factor * X_dots[:, y_dim]
    #projs = [np.dot(X_dots[i, :], X_ps[i, :]) for i in range(len(t_eval))]
    #print(min(projs))

    ax.quiver(
        xs[s],
        ys[s],
        t_eval[s],
        x_dots[s],
        y_dots[s],
        0,
        length=1,
        color=color_dict["X_dot"],
        label=label_dict["X_dot"],
    )

    x_ps = scale_factor * X_ps[:, x_dim]
    y_ps = scale_factor * X_ps[:, y_dim]

    ax.quiver(
        xs[s],
        ys[s],
        t_eval[s],
        x_ps[s],
        y_ps[s],
        0,
        length=1,
        color="green",
        label=label_dict["X_p"],
    )
    #ax.set_xticks([])
    #ax.set_xticklabels([])
    #ax.set_yticks([])
    #ax.set_yticklabels([])
    ax.set_zlabel("time")
    ax.set_xlabel(f"$x={state_vector[x_dim]}$")
    ax.set_ylabel(f"$y={state_vector[y_dim]}$")
    ax.legend()


def plot_time_x(
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
):
    nt = len(t_eval)
    s = slice(0, nt, ss)
    #ax.invert_yaxis()
    #ax.invert_zaxis()
    # trajectory for X
    zs = Xs[:, z_dim]
    ys = Xs[:, y_dim]
    ax.plot3D(t_eval, ys, zs, label=label_dict["X"], color=color_dict["X"])
    ax.scatter(t_eval[s], ys[s], zs[s], s=10, color=color_dict["X"])

    z_cs = X_cs[:, z_dim]
    y_cs = X_cs[:, y_dim]
    ax.plot3D(t_eval, y_cs, z_cs, label=label_dict["X_c"], color=color_dict["X_c"])
    ax.scatter(t_eval[s], y_cs[s], z_cs[s], s=10, color=color_dict["X_c"])

    z_dots = X_dots[:, z_dim]
    y_dots = X_dots[:, y_dim]
    projs = [np.dot(X_dots[i, :], X_ps[i, :]) for i in range(len(t_eval))]
    print(min(projs))

    ax.quiver(
        t_eval[s],
        ys[s],
        zs[s],
        0,
        y_dots[s],
        z_dots[s],
        length=1,
        color="red",
        label=label_dict["X_dot"],
    )

    z_ps = X_ps[:, z_dim]
    y_ps = X_ps[:, y_dim]

    ax.quiver(
        t_eval[s],
        ys[s],
        zs[s],
        0,
        y_ps[s],
        z_ps[s],
        length=1,
        color="green",
        label=label_dict["X_p"],
    )
    #ax.legend()

    # since the arrows are drawn in data space we adjust the range
    vals=np.concatenate([zs,ys,z_cs,y_cs]) 
    minv=min(vals)
    maxv=max(vals)
    ax.set_zlim( (minv,maxv))
    ax.set_ylim( (minv,maxv))
    
    ax.set_xlabel("time")
    ax.set_zlabel(f"$z={state_vector[z_dim]}$")
    ax.set_ylabel(f"$y={state_vector[y_dim]}$")


def plot_Inputs_time_x(
    ax, color_dict, label_dict, state_vector, y_dim, z_dim, ss, t_eval, Is
):
    nt = len(t_eval)
    s = slice(0, nt, ss)
    zs = Is[:, z_dim, 0]
    ys = Is[:, y_dim, 0]
    #ax.invert_yaxis()
    #ax.invert_zaxis()
    ax.plot3D(t_eval, ys, zs, label=label_dict["I"], color=color_dict["I"])
    ax.scatter(
        t_eval[s], ys[s], zs[s],
        s=10,
        color=color_dict["I"]
    )
    ax.quiver(
        t_eval[s], 0, 0,
        0, ys[s], zs[s],
        length=1,
        color=color_dict["I"]
        #label=label_dict["X_dot"],
    )
    #ax.set_zlim(0, max(ys))
    ax.set_zlabel(f"$z=I_{{{state_vector[z_dim]}}}$")
    ax.set_ylabel(f"$y=I_{{{state_vector[y_dim]}}}$")
    #ax.set_zticks = [min(zs), 0 ,max(zs)]
    #ax.set_ztickilabelss = ["a","b",",c"]


def plot_Inputs_time_z(
    ax, color_dict, label_dict, state_vector, x_dim, y_dim, ss, t_eval, Is
):
    nt = len(t_eval)
    s = slice(0, nt, ss)

    # arrowhaeds are plotted data space and look funny if the aspect ratios
    # between the axes are very different (not the box aspect)
    scale_factor = np.max(t_eval)/(np.max(Is)-np.min(Is))/8
    xs = scale_factor * Is[:, x_dim, 0]
    ys = scale_factor * Is[:, y_dim, 0]
    #ax.invert_xaxis()
    #ax.invert_yaxis()
    vals=np.concatenate([xs,ys])
    min_v, max_v = min(vals), max(vals)
    

    #ax.set_xlim(min_v, max_v)
    #ax.set_ylim(min_v, max_v)
    
    ax.set_xticks([0])
    #ax.set_xticklabels(["0", "x_max"])
    ax.set_yticks([0])
    #ax.set_yticklabels(["0", "y_max"])
    ax.plot3D(
        xs,
        ys,
        t_eval,
        label=label_dict["I"],
        color=color_dict["I"]
    )
    ax.scatter(
        xs[s], ys[s], t_eval[s],
        s=10,
        color=color_dict["I"],
        #label=label_dict["I"],
    )
    ax.quiver(
        0, 0, t_eval[s], 
        xs[s], ys[s], 0, 
        length=1,
        color=color_dict["I"],
        #label=label_dict["I"],
    )
    ax.set_zlabel("time")
    ax.set_xlabel(f"$x=I_{{{state_vector[x_dim]}}}$")
    ax.set_ylabel(f"$y=I_{{{state_vector[y_dim]}}}$")
    ax.legend()


def plot_xi_time_z(
    ax, color_dict, label_dict, state_vector, x_dim, y_dim, ss, t_eval, Xis
):
    nt = len(t_eval)
    s = slice(0, nt, ss)

    # arrowhaeds are plotted data space and look funny if the aspect ratios
    # between the axes are very different (not the box aspect)
    scale_factor = np.max(t_eval)/(np.max(Xis)-np.min(Xis))/8
    xs = scale_factor * Xis[:, x_dim]
    ys = scale_factor * Xis[:, y_dim]
    #ax.invert_xaxis()
    #ax.invert_yaxis()
    vals=np.concatenate([xs,ys])
    min_v, max_v = min(vals), max(vals)
    

    #ax.set_xlim(min_v, max_v)
    #ax.set_ylim(min_v, max_v)
    
    ax.set_xticks([0])
    #ax.set_xticklabels(["0", "x_max"])
    ax.set_yticks([0])
    #ax.set_yticklabels(["0", "y_max"])
    ax.plot3D(
        xs,
        ys,
        t_eval,
        label=label_dict["Xi"],
        color=color_dict["Xi"]
    )
    ax.scatter(
        xs[s], ys[s], t_eval[s],
        s=10,
        color=color_dict["Xi"],
    )
    ax.quiver(
        0, 0, t_eval[s], 
        xs[s], ys[s], 0, 
        length=1,
        color=color_dict["Xi"],
    )
    ax.set_zlabel("time")
    ax.set_xlabel(f"$x=xi_{{{state_vector[x_dim]}}}$")
    ax.set_ylabel(f"$y=xi_{{{state_vector[y_dim]}}}$")
    ax.legend()


def plot_timerectangles(ax, intervals, colors, t_eval, ss):
    y_min, y_max = ax.get_ybound()
    for i, _ in enumerate(intervals):
        t_min, t_max = t_eval[intervals[i, 0] * ss] ,t_eval[intervals[i, 1]*ss]
        poly_ts = [t_min, t_min, t_max, t_max, t_min]
        poly_ys = [y_min, y_max, y_max, y_min, y_min]
        ax.fill(
            poly_ts, poly_ys,
            color=colors[i],
            alpha=.1
        )


def plot_timerectangles_3d(ax, intervals, colors, t_eval, ss):
    x_min, x_max = ax.get_xbound()
    y_min, y_max = ax.get_ybound()
    alpha=0.1
    for i, _ in enumerate(intervals):
        t_min, t_max = t_eval[intervals[i, 0] * ss] ,t_eval[intervals[i, 1]*ss]
        poly_ts = [t_min, t_min, t_max, t_max, t_min]
        poly_xs = [x_min, x_max, x_max, x_min, x_min]
        poly_ys = [y_max, y_max, y_max, y_max, y_max]
        vertices = [list(zip(poly_xs, poly_ys, poly_ts))]

        poly1 = Poly3DCollection(vertices, color=colors[i], alpha=alpha)
        ax.add_collection3d(poly1)
        poly_ts = [t_min, t_min, t_max, t_max, t_min]
        poly_xs = [x_max, x_max, x_max, x_max, x_max]
        poly_ys = [y_min, y_max, y_max, y_min, y_min]
        vertices = [list(zip(poly_xs, poly_ys, poly_ts))]
        poly2 = Poly3DCollection(vertices, color=colors[i], alpha=alpha)
        ax.add_collection3d(poly2)

        poly_ts = [t_max, t_max, t_max, t_max, t_max]
        poly_xs = [x_min, x_max, x_max, x_min, x_min]
        poly_ys = [y_min, y_min, y_max, y_max, y_min]
        vertices = [list(zip(poly_xs, poly_ys, poly_ts))]
        poly2 = Poly3DCollection(vertices, color=colors[i], alpha=alpha)
        ax.add_collection3d(poly2)

def ticks_test(ax):
    theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
    z = np.linspace(-2, 2, 100)
    r = z**2 + 1
    x = r * np.sin(theta)
    y = r * np.cos(theta)
    ax.plot(x, y, z, label='parametric curve')
    ax.plot(x, y, z*2, label='test')
    ax.legend()

    ax.set_xlabel('$X$', fontsize=20)
    ax.set_ylabel('$Y$')
    ax.yaxis._axinfo['label']['space_factor'] = 3.0
    # set z ticks and labels
    ax.set_zticks([-2, 0, 2])
    ax.set_zticklabels(["a", 0, "b"])
    ax.set_xticks([min(x), 0, max(x)])
    ax.set_xticklabels(["x_min", 0, "x_max"])
    ax.set_yticks([min(y), 0, max(y)])
    ax.set_yticklabels(["y_min", 0, "y_max"])
    # change fontsize
    for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(10)
    # disable auto rotation
    ax.zaxis.set_rotate_label(False) 
    ax.set_zlabel('$\gamma$', fontsize=30, rotation = 0)
