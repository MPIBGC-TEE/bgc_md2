from pathlib import Path
from sympy import lambdify
from scipy.integrate import solve_ivp
from importlib import import_module
import numpy as np
import matplotlib.pyplot as plt

def solve_and_plot(mod_name):
    name=f"{Path(__file__).stem}_{mod_name}"
    mod=import_module(mod_name)
    sol_Monod=solve_ivp(
        fun=lambdify(
            (
                mod.t,
                mod.state_var_tuple
            ),
            mod.rhs_sym.subs(mod.par_dict),"numpy"),
        t_span=(0,500),
        y0=mod.state_var_tuple.subs(mod.start_value_dict),
        dense_output=True
    ).sol
    
    plot_times=np.linspace(sol_Monod.t_min,sol_Monod.t_max,100)
    fig=plt.figure()
    ax=fig.add_subplot()
    values=sol_Monod(plot_times)
    for i in range(values.shape[0]):
        ax.plot(
            plot_times,
            values[i,:],
            label=str(mod.state_var_tuple[i])
        ) 
    ax.legend()
    ax.set_title(label=name)
    fig_path=Path("figures") 
    if not fig_path.exists():
        fig_path.mkdir()
    fig.savefig(fig_path.joinpath(f"{name}.pdf"))

solve_and_plot("P2TherMic_dormancy_2_Monod")
solve_and_plot("P2TherMic_dormancy_2_MTS")
