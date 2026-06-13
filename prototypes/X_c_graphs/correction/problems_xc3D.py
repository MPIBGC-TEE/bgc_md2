# %matplotlib inline
# %load_ext autoreload
# %autoreload 2

# ## paragraph in Luo2017 paper 3.2 Direction and rate of C storage change at a given time
from typing import Tuple
from IPython.display import display_pretty, display
from sympy import (
    var,
	Function,
	sin,
	cos,
	diag,
	simplify,
	sympify,
	lambdify,
	Piecewise,
	pi
)
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np
from scipy.integrate import solve_ivp
from sympy import latex
import CompartmentalSystems.helpers_reservoir as hr
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers
from bgc_md2.display_helpers import mass_balance_equation
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)

import helpers as h 
import plot_functions  as pf 
import figure_functions  as ff 
# define the system symbolically

# all the symbols you are going to use in the equations
var(
    """
    C_veg 
    C_soil
    r_veg2soil_0
    r_soil2veg
    r_veg2out
    r_soil2out
    I_veg_0
"""
)
xi_veg=Function("xi_veg")
I_veg=Function("I_veg")
t = TimeSymbol("t")  # the symbol used for time since it has a special role

# now we create intervals of different behaviour
# by defining the functions piecewise we start with a frozen system (I(t)=I_0=const. B(t)=B_0=const.)
mvs_1 = CMTVS(
    {
        StateVariableTuple((C_veg, C_soil)),  # the pool names in your preferred order
        t,
        InFluxesBySymbol({C_veg: I_veg(t), C_soil: 0}),
        OutFluxesBySymbol(
            {
                C_veg: xi_veg(t) * r_veg2out * C_veg,
                C_soil: r_soil2out * C_soil
            }
        ),
        InternalFluxesBySymbol(
            {
                (C_veg, C_soil): xi_veg(t) * r_veg2soil_0 * C_veg,
                #(C_soil, C_veg): r_soil2veg * C_soil
            }
        ),
    },
    bgc_md2_computers(),
)

display(mass_balance_equation(mvs_1))
#mvs_1.provided_mvar_values

#write it to a file for the figure caption in the paper
from helpers import mass_balance_equation
from pathlib import Path
lsnp=Path("latex_snips")
fp=Path("figures")
with open(lsnp.joinpath("Massbalance.tex"),"w") as f:
    f.write(mass_balance_equation(mvs_1))


import bgc_md2.helper as bh 
bh.compartmental_graph(mvs_1)

M = mvs_1.get_CompartmentalMatrix()
Tau = - M.inverse_LU()
with open(lsnp.joinpath("Tau.tex"),"w") as f:
    f.write(latex(Tau))
simplify(Tau)
In=mvs_1.get_InputTuple()
X_csym=Tau * In
with open(lsnp.joinpath("X_csym.tex"),"w") as f:
    f.write(latex(X_csym))
Tau*In

calm = 2 * pi
crazy = calm
delta = pi
# we let the input loose first
t_thaw_I_0 = calm + 1 / 2 * pi
# we let the input loose first
t_freeze_I_0 = t_thaw_I_0 + crazy
# t_thaw_I_1=t_freeze_I_0+calm
# t_freeze_I_1=t_thaw_I_1+crazy

t_thaw_xi_0 = t_freeze_I_0 + calm
t_freeze_xi_0 = t_thaw_xi_0 + crazy
# r_veg2soil=r_veg2soil*Piecewise((1,True))
#os_i = 1.1  # >1 (otherwise negative influx)
os_i = 2.1  # >1 (otherwise negative influx)
os_xi = 8  # >1 (otherwise negative influx)
assert(os_i > 1)
assert(os_xi > 1)
spec_dict={
    xi_veg(t) : Piecewise(
        (1, t < t_thaw_xi_0),
        ((sin(t) + os_xi) / (os_xi + 1), t < t_freeze_xi_0),
        # ((sin(t-delta)+os_xi)*r_veg2soil_0,t<t_freeze_xi_0),
        (1, True),
    ),
    I_veg(t) : Piecewise(
        (I_veg_0, t < t_thaw_I_0),
        ((sin(t) + os_i) / (os_i + 1) * I_veg_0, t < t_freeze_I_0),
        # (I_veg_0,t<t_thaw_I_1),
        # ((sin(t)+os_i)/(os_i+1)*I_veg_0,t<t_freeze_I_1),
        (I_veg_0, True),
    ),
}    
mvs_2 =   CMTVS(
    { v.subs(spec_dict) for v in mvs_1.provided_mvar_values },
    bgc_md2_computers()
)    
# start to query the model description..
param_dict = {
    r_veg2soil_0: 0.6,
    #r_veg2soil_0: 0.75,
    r_soil2veg: 0.075,
    r_veg2out: 0.1,
    #r_soil2out: 1.25,
    r_soil2out: 0.125,
    I_veg_0: 5.5,
}
func_dict = {}

t_start = 0
# t_end=15*np.pi
t_end = 5 * float(calm.subs({pi: np.pi}))

#X_0 = np.array([1, 0.125])*5  # .reshape((2,1))
#X_0 = np.array([0, 0])  # .reshape((2,1))
#X_0 = np.array([15,0])  # .reshape((2,1))
X_0 = np.array([15,3.5])  # .reshape((2,1))

t_eval = np.linspace(
    start=t_start,
    stop=t_end,
    # num=1000
    num=5000,
)

Xs, X_cs, X_ps, X_dots, Is, Bs = h.compute_trajectories(mvs_2, param_dict, func_dict,  X_0, t_eval)
state_vector = mvs_2.get_StateVariableTuple()
Xi_vegs = lambdify(t,xi_veg(t).subs(spec_dict))(t_eval)
Xi_d = {C_veg:Xi_vegs,C_soil:np.ones_like(Xi_vegs)}
Xis = np.stack([Xi_d[pn] for pn in state_vector],axis=1) # same order as for state 
# plots are either x,y or y,z 
x_dim = 0
y_dim = 1
z_dim = 0
ss = 100
color_dict = {
    "X": "black",
    "X_c": "orange",
    "X_p": "green",
    "X_dot": "red",
    "I": "blue",
    "Xi": "purple",
}
label_dict = {
    "X": "$\mathbf{X}$",
    "X_c": "$\mathbf{X_c}$",
    "X_p": "$\mathbf{X_p}$",
    "X_dot": "$\mathbf{\dot{X}}$",
    "I": "$\\mathbf{I}$",
    "Xi": "$\\mathbf{\\xi}$",
}
rectangle_colors=np.array( ["red","green","blue","purple"])
cm_per_inch=2.54
A4 = tuple((cm/cm_per_inch for cm in (21,29.7)))
A3 = tuple((cm/cm_per_inch for cm in (29.7,42)))
intervals = np.array([
    (0, 7), 
    #(13,15), 
    (16, 19), 
    (37, 42)
])

linestyles = ["solid","dashed"]
#f1=ff.combined_timelines_and_2d_phase_space(
f1=h.fig_save(
    ff,
    "combined_timelines_and_2d_phase_space",
    fp,
    size=A3,
    intervals=intervals,
    rectangle_colors=rectangle_colors,
    linestyles=linestyles,
    color_dict=color_dict,
    label_dict=label_dict,
    state_vector=state_vector,
    x_dim=x_dim,
    y_dim=y_dim,
    z_dim=z_dim,
    ss=ss,
    Is=Is,
    Xis=Xis,
    t_eval=t_eval,
    Xs=Xs,
    X_cs=X_cs,
    X_dots=X_dots
)


# +
#f2=ff.plot_vector_3d_X_Xc_time_z(
f2=h.fig_save(
    ff,
    "plot_vector_3d_X_Xc_time_z",
    fp,
    size=A3,
    intervals=intervals,
    rectangle_colors=rectangle_colors,
    color_dict=color_dict,
    label_dict=label_dict,
    state_vector=state_vector,
    x_dim=x_dim,
    y_dim=y_dim,
    z_dim=z_dim,
    ss=ss,
    Is=Is,
    Xis=Xis,
    t_eval=t_eval,
    Xs=Xs,
    X_cs=X_cs,
    X_dots=X_dots,
    X_ps=X_ps
)

#f3=ff.plot_vector_3d_I_Xi_X_Xc_time_z(
#    size=A3,
#    intervals=intervals,
#    rectangle_colors=rectangle_colors,
#    color_dict=color_dict,
#    label_dict=label_dict,
#    state_vector=state_vector,
#    x_dim=x_dim,
#    y_dim=y_dim,
#    z_dim=z_dim,
#    ss=ss,
#    Is=Is,
#    Xis=Xis,
#    t_eval=t_eval,
#    Xs=Xs,
#    X_cs=X_cs,
#    X_dots=X_dots,
#    X_ps=X_ps
#)
#f4=ff.plot_vector_3d_I_Xi_X_Xc_time_x(
#    size=A3,
#    intervals=intervals,
#    rectangle_colors=rectangle_colors,
#    color_dict=color_dict,
#    label_dict=label_dict,
#    state_vector=state_vector,
#    x_dim=x_dim,
#    y_dim=y_dim,
#    z_dim=z_dim,
#    ss=ss,
#    Is=Is,
#    Xis=Xis,
#    t_eval=t_eval,
#    Xs=Xs,
#    X_cs=X_cs,
#    X_dots=X_dots,
#    X_ps=X_ps
#)
#################################################################################
# -


