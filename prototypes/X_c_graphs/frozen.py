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
from sympy import diff,Matrix,Function,var,sin,cos,lambdify,exp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

var("t x")
k=sin(t)+2
#k=3+t-t
#u=cos(t)+2
#u=exp(t/10)+2
u=cos(t/2)+2
x_dot=u-k*x
x_c=1/k*u
x_p=x_c-x
x_of_t=Function("x")(t)
x_dot_of_t=x_dot.subs({x: x_of_t})
x_dot_of_t
# +
u_num=lambdify(args=t,expr=u,modules=['numpy'])
k_num=lambdify(args=t,expr=k,modules=['numpy'])

x_dot_num=lambdify(args=(t,x),expr=x_dot,modules=['numpy'])
#x_dot_num(0,1)
t_start=0
t_end=1*np.pi
times=np.linspace(t_start,t_end,1000)
sol = solve_ivp(
    x_dot_num,
    t_span=[t_start,t_end],
    t_eval=times,
    y0=[2]
)
xs=sol.y.transpose()

nt=len(times)
rnt=range(nt)
n_freeze=4
ts_freeze=[t_start+(t_end-t_start)*1.0*i_freeze/n_freeze for i_freeze in range(n_freeze)] 
#x_dots= make_time_line(x_dot)
#x_ps= make_time_line(x_p)

fig=plt.figure(figsize=(20,20))
axs=fig.subplots(3,sharex=True)
ax=axs[0]
name="u"
expr=globals()[name]
ax.plot(
    times,
    u_num(times),
    label=name
)
def plot_const(u_f,t_f,ax):
    u_f = u_num(t_f)
    n_p=100
    ts = np.linspace(t_f,t_end,n_p)
    ax.plot(ts,np.ones_like(ts)*u_f)
    
for i in range(len(ts_freeze)):
    plot_frozen(u_num,ts_freeze[i],ax)
#ax=axs[1]
#name="k"
#expr=globals()[name]
#ax.plot(
#    times,
#    make_time_line(expr),
#    label=name
#)
#ax=axs[2]
#for name in ["x_c"]:
#    expr=globals()[name]
#    ax.plot(
#        times,
#        make_time_line(expr),
#        label=name
#    )
#ax.plot(times,xs,label="x")
ax.legend()
# -






