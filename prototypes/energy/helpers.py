from pathlib import Path
from sympy import lambdify
from scipy.integrate import solve_ivp
from importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import inspect
import torch
# sbi
from sbi import analysis as analysis
from sbi import utils as utils
from sbi.inference import SNPE, simulate_for_sbi
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)


def solve_and_plot(mod_name):
    func_name=inspect.currentframe().f_code.co_name
    name=f"{func_name}_{mod_name}"
    mod=import_module(mod_name)
    # print(mod.rhs_sym.subs(mod.par_dict).free_symbols)
    sol_Monod=solve_ivp(
        fun=lambdify(
            (
                mod.t,
                mod.state_var_tuple
            ),
            mod.rhs_sym.subs(mod.par_dict),"numpy"),
        t_span=(0,250),
        y0=mod.state_var_tuple.subs(mod.start_value_dict),
        dense_output=True
    ).sol
    
    plot_times=np.linspace(sol_Monod.t_min,sol_Monod.t_max,100)
    fig=plt.figure()
    ax=fig.add_subplot()
    for i in range(values.shape[0]):
        ax.plot(
            plot_times,
            sol_Monod(plot_times),
            label=str(mod.state_var_tuple[i])
        ) 
    ax.legend(loc="upper right")
    ax.set_title(label=name)
    
    fig_path=Path("figures") 
    if not fig_path.exists():
        fig_path.mkdir()
    fig.savefig(fig_path.joinpath(f"{name}.pdf"))


def fit(mod_name):
    func_name=inspect.currentframe().f_code.co_name
    name=f"{func_name}_{mod_name}"
    mod=import_module(mod_name)
    sv=mod.mod.state_var_tuple
    #res=da_mod.simulation_wrapper(mod.prior_min)
    prior, num_parameters, prior_returns_numpy = process_prior(mod.prior)
    simulation_wrapper = process_simulator(mod.simulation_wrapper, prior, prior_returns_numpy)
    check_sbi_inputs(simulation_wrapper, prior)

    # Create inference object. Here, NPE is used.
    inference = SNPE(prior=prior)
    
    # generate simulations and pass to the inference object

    theta, x = simulate_for_sbi(
            simulation_wrapper, 
            proposal=prior,
            num_simulations=2000,
            num_workers=8
    )
    inference = inference.append_simulations(
        theta, 
        x
    )
    
    # train the density estimator and build the posterior
    density_estimator = inference.train()
    posterior = inference.build_posterior(density_estimator)
    default=mod.ep0.reshape(1,len(mod.ep0))
    default_obs= simulation_wrapper(default)
    fig=plt.figure()
    axs=fig.subplots(len(sv))
    values = default_obs.reshape(
        len(sv),
        -1
    )
    cmap=colormaps['tab20']
    color_dict={
        str(sv[i]): cmap.colors[i]
        for i in range(len(sv))
    }
    for i in range(values.shape[0]):
        ax=axs[i]
        ax.plot(
            mod.simulation_times(values),
            values[i,:],
            label=str(mod.mod.state_var_tuple[i]),
            color=color_dict[str(sv[i])]
        ) 
        ax.set_title(str(sv[i]))
    

    n_samples=3000
    samples=posterior.sample(
        (n_samples,),
        default_obs
    )
    sample_obs=simulation_wrapper(samples).reshape(
        n_samples,
        len(sv),
        -1
    )
    for n_s in range(n_samples):
        values = sample_obs[n_s].reshape(
            len(mod.mod.state_var_tuple),
            -1
        )
        for i in range(values.shape[0]):
            ax=axs[i]
            ax.plot(
                mod.simulation_times(values),
                values[i,:],
                color=color_dict[ str(sv[i])],
                alpha=0.01
            ) 
            ax.set_title(str(sv[i]))
    
    fig_path=Path("figures") 
    if not fig_path.exists():
        fig_path.mkdir()
    fig.savefig(fig_path.joinpath(f"{name}.pdf"))
    fig2,axs2 = analysis.pairplot(
            samples,           
            points=default,
            #limits=list(zip(mod.prior_min,mod.prior_max)),
            points_color="r",
            labels=mod.estimated_parameters._fields

    )
    fig2.savefig(fig_path.joinpath(f"{name}_pairs.pdf"))
    #from IPython import embed;embed()

