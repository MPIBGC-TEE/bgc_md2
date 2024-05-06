# visualization
from pathlib import Path
from sympy import lambdify
from scipy.integrate import solve_ivp
from importlib import import_module
import matplotlib.pyplot as plt
import numpy as np
import torch

from sbi import analysis as analysis

# sbi
from sbi import utils as utils
from sbi.inference import SNPE, simulate_for_sbi
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)
import sys

#sys.path.insert(0,"DA")
sys.path.insert(0,".")

def fit(mod_name):
    name=f"{Path(__file__).stem}_{mod_name}"
    mod=import_module(mod_name)
    #res=da_mod.simulation_wrapper(mod.prior_min)
    prior, num_parameters, prior_returns_numpy = process_prior(mod.prior)
    simulation_wrapper = process_simulator(mod.simulation_wrapper, prior, prior_returns_numpy)
    check_sbi_inputs(simulation_wrapper, prior)

    # Create inference object. Here, NPE is used.
    inference = SNPE(prior=prior)
    
    # generate simulations and pass to the inference object

    theta, x = simulate_for_sbi(simulation_wrapper, proposal=prior,
                                 num_simulations=100000, num_workers=16)
    inference = inference.append_simulations(theta, x)
    
    # train the density estimator and build the posterior
    density_estimator = inference.train()
    posterior = inference.build_posterior(density_estimator)
    from IPython import embed;embed()

for mod_name in [
        (
            "DA_P2TherMic_dormancy_2_Monod"
        ),
        #"P2TherMic_dormancy_2_MTS",
        #"P2TherMic_dormancy_3_Monod",
        #"P2TherMic_dormancy_3_MTS"
    ]:    
    fit(mod_name)
