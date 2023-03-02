# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
from sympy import diff, Matrix, Function, var, sin, cos, lambdify, exp, pi
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

var("t x")
k = (sin(t+pi/2) + 2)/2
# k=3+t-t
# u=cos(t)+2
# u=exp(t/10)+2
u_min=5
u_amp=5
u = (sin(t+pi)+1)*u_amp+u_min
x_dot = u - k * x
x_c = 1 / k * u
x_p = x_c - x
x_of_t = Function("x")(t)
x_dot_of_t = x_dot.subs({x: x_of_t})
x_dot_of_t
# +

x_dot_num = lambdify(args=(t, x), expr=x_dot, modules=["numpy"])
# x_dot_num(0,1)
n_t = 400
t_start = 0
t_end = 4 * np.pi
t_step = (t_end - t_start)/n_t
times = np.arange(t_start, t_end, t_step)
sol = solve_ivp(
    x_dot_num, 
    t_span=[t_start, t_end], 
    #t_eval=times, 
    y0=[20],
    dense_output=True
)
xs = sol.y.transpose()

u_num = lambdify(args=t, expr=u, modules=["numpy"])
k_num = lambdify(args=t, expr=k, modules=["numpy"])
x_num = lambda t: sol.sol(t).reshape(-1)
x_c_num = lambdify(args=t, expr=x_c, modules=["numpy"])
nt = len(times)
rnt = range(nt)
n_fr = 8
# x_dots= make_time_line(x_dot)
# x_ps= make_time_line(x_p)

c_map = plt.colormaps['magma']
cols = c_map.resampled(len(times)).colors
small=1
large=20
transparent=0.2
def plot_cooling_trajectory(
    ax,
    f_num,
    **kwargs,
):
    ax.set_title="u(t)"
    vals = f_num(times)
    ax.scatter(
        x=times,
        y=vals,
        c=cols,
        **kwargs,
    )
    last_x = times[-1]
    last_y = vals[-1]
    last_c = cols[-1]
    ax.scatter(
        last_x,
        last_y,
        color=last_c,
        edgecolors="black"
    )

def plot_trajectory_forward(
    ax,
    f_num,
    i_fr_1,
    **kwargs
):
    #b_times_1 = times[0:i_fr_1]

    f_times_1 = times[i_fr_1:]
    c_fr_1 = cols[i_fr_1]
    #ax.scatter(
    #    x=b_times_1, 
    #    y=np.ones_like(b_times_1)*f_num(t_fr_1),
    #    color=c_fr_1,
    #)
    ax.scatter(
        x=f_times_1, 
        y=f_num(f_times_1),
        color=c_fr_1,
        **kwargs
    )


fig = plt.figure(figsize=(20, 20))
axs = fig.subplots(3, sharex=True)
#############################################################
ax=axs[0]
ax.set_title("Input $u(t)$,and frozen input $u_{f_1}=u(t_{f_1})$")
plot_cooling_trajectory(
    ax,
    u_num,
    s=large,
    alpha=1
)
i_fr_1 = int(n_t/4)*1
i_fr_2 = int(n_t/4)*2
t_fr_1 = times[i_fr_1]
t_fr_2 = times[i_fr_2]
plot_trajectory_forward(
    ax,
    lambda t: u_num(t_fr_1)*np.ones_like(t),
    i_fr_1,
    s=large,
    alpha=transparent
)
ax.axvline(x=t_fr_1,lw=0.5,color="black")
#############################################################
ax=axs[1]
plot_cooling_trajectory(
    ax, 
    k_num,
    s=large,
    alpha=1
)
plot_trajectory_forward(
    ax,
    lambda t: k_num(t_fr_1)*np.ones_like(t),
    i_fr_1,
    s=large,
    alpha=transparent
)
ax.axvline(x=t_fr_1,lw=0.5,color="black")
#############################################################
ax=axs[2]
ax.set_title("X,X_c, X_f1, X_f2...")
plot_cooling_trajectory(
    ax,
    x_num,
    s=large,
    alpha=1,
    label="x"
)
plot_cooling_trajectory(
    ax,
    x_c_num,
    s=large,
    alpha=transparent,
    label="x_c"
)
sol_fr_1 = solve_ivp(
    lambda t, x: u_num(t_fr_1) - k_num(t_fr_1) * x,
    t_span=[t_fr_1, t_end], 
    y0=x_num(t_fr_1),
    dense_output=True
)
plot_trajectory_forward(
    ax, 
    lambda t: sol_fr_1.sol(t).reshape(-1), 
    i_fr_1,
    s=large,
    alpha=transparent,
    label="X_{fr_1}(t)"
)
ax.axvline(x=t_fr_1,lw=0.5,color="black")
ax.axhline(y=x_c_num(t_fr_1),lw=0.5,color="black")
sol_fr_2 = solve_ivp(
    lambda t, x: u_num(t_fr_2) - k_num(t_fr_2) * x,
    t_span=[t_fr_2, t_end], 
    y0=x_num(t_fr_2),
    dense_output=True
)
plot_trajectory_forward(
    ax, 
    lambda t: sol_fr_2.sol(t).reshape(-1), 
    i_fr_2,
    s=large,
    alpha=transparent,
    label="X_{fr_2}(t)"
)
ax.axvline(x=t_fr_2,lw=0.5,color="black")
ax.axhline(y=x_c_num(t_fr_2),lw=0.5,color="black")
# ax.plot(times,xs,label="x")
# -

#ax.legend()
fig.savefig("x_c.pdf")
