# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Use the common infrastructure to inspect the model
#
# This illustrative notebook shows how to use the factored out model description in ```source.py``` the model specific functions in ```model_specific_helpers.py``` and the general infrastructure in ```general_helpers.py``` to inspect the model.
#
# ## Preconditions
# This is the next step after the ```createModel.py``` notebook where everything was defined in one place, the notebook.
# In order to be able to use the symbolic descriptions and functions in other code (e.g tests, scripts or other  notebooks)
# we have disassembled it, moving the functions into the seperate file ```model_specific_helpers.py``` and the model description to ```source.py```.
#
# ## Applications
# 1. Inspect the structure of the model with symbolic tools
# 1. Run the model forward with a guessed  parameter set
# 1. Optimize the parameter set using data assimilation
# 1. Use the optimized paramerts to run the tracability analysis.
#


# +

# adjust the output to full width
from IPython.display import HTML

# display(HTML("<style>.container { width:100% !important; }</style>"))

# make changes to imported files immidiately available
# avoiding the need to reload (in most cases)
# %load_ext autoreload
# %autoreload 2

from pathlib import Path
import json
from sympy import Symbol, Function
import numpy as np
import matplotlib.pyplot as plt
from ComputabilityGraphs.CMTVS import CMTVS
import CompartmentalSystems.helpers_reservoir as hr
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c
import bgc_md2.helper as h
import bgc_md2.display_helpers as dh

# imports from new files
import sys

sys.path.insert(0, '..')
from source import mvs
import model_specific_helpers_2 as msh
from general_helpers import day_2_month_index, make_B_u_funcs_2
import general_helpers as gh

# -

# we can also print the whole mass balance equation
dh.mass_balance_equation(mvs)

# we can also plot a picture
h.compartmental_graph(mvs)

# +
with Path('config.json').open(mode='r') as f:
    conf_dict = json.load(f)

# msh.download_my_TRENDY_output(conf_dict)

# + codehighlighter=[[5, 6], [23, 33], [5, 6], [23, 33]]
#     # Read NetCDF data  ********************************************************
svs, dvs = msh.get_global_mean_vars(dataPath=Path(conf_dict["dataPath"]))
print(dvs.npp.shape)
print(svs.cVeg.shape)
svs_0 = msh.Observables(*map(lambda v: v[0], svs))
dvs_0 = msh.Drivers(*map(lambda v: v[0], dvs))
print('svs_0=', svs_0)
print('dvs_0=', dvs_0)
# -
# To be able to run the model forward we not only have to replace parameter symbols by values but symbolic functions by normal python functions.
# In our case the functions for NPP and ξ have to be provided. NPP_fun will interpolate the NPP for the day in question from the data. Which we have to load.
# We will later store these functions in  `model_specific_helpers.py` which resides in the same folder as this notebook. You will have to adapt them to your data set.

# +
# ## Data assimilation
# Until now we have used only the initial values of the observations.
# The next step is to decide which parameters we want to consider fixed and which to be estimated.
# This distinction helps, to keep the to create a function which only takes the estimated parameters and thus can be used by a generalized mcmc as will become clear.
#
# We can change which parameters we fix and which we estimate later or can have several approaches for the same symbolic model.
# The distinction is not model inherent but just a reflection of our choice for data assimilation.
# The more parameter values we can find out from the literature the fewer values we have to estimate.

# +
cpa = msh.Constants(
    cVeg_0=svs_0.cVeg,
    cLitter_0=svs_0.cLitter,
    cSoil_0=svs_0.cSoil,
    npp_0=dvs.npp[0],  # * 86400,   # kg/m2/s kg/m2/day
#    xi_0=dvs.xi[0],
    rh_0=svs_0.rh,  # * 86400,   # kg/m2/s kg/m2/day
    ra_0=svs_0.ra,  # * 86400,   # kg/m2/s kg/m2/day
    # r_C_root_litter_2_C_soil_slow=3.48692403486924e-5,
    # r_C_root_litter_2_C_soil_passive=1.74346201743462e-5,
    number_of_months=len(svs.rh)
    # number_of_months=120 # for testing and tuning mcmc
    # number_of_months=3840 # for testing and tuning mcmc
)

# ### Finding better start values for the data assimilation
# You don't have to do this. It's a heuristic approach to find a better starting position.
epa_0 = msh.EstimatedParameters(
    beta_wood1=0.229,
    beta_wood2=0.27112,
    beta_leaf=0.2032,
    beta_root=0.1303,
    #   beta_fruit=0.1,
    T_0=1.865,
    E=3.246,
    KM=10.2,
    # r_C_leaf_rh=0,
    # r_C_wood_rh=0,
    # r_C_root_rh=0,
    r_C_litter1_rh=0.0012934,
    r_C_litter2_rh=0.59747,
    r_C_litter3_rh=0.01977,
    r_C_litter4_rh=0.0006258,
    r_C_litter5_rh=0.0000016211,
    r_C_litter6_rh=0.002435,
    r_C_som1_rh=0.00229726,
    r_C_som2_rh=0.0010939,
    r_C_som3_rh=0.0000455408,
    r_C_som4_rh=0.00028037,

    r_C_wood1_2_C_wood3=0.0107,
    r_C_wood1_2_C_litter1=0.066445,
    r_C_wood2_2_C_wood4=0.00063076,
    r_C_wood2_2_C_litter2=0.00168526,
    r_C_wood3_2_C_litter1=0.000007645231,
    r_C_wood4_2_C_litter2=0.000209176,
    r_C_leaf_2_C_litter3=0.0020595,
    r_C_leaf_2_C_litter5=0.0000250875,
    r_C_root_2_C_litter4=0.0000076896,
    r_C_root_2_C_litter6=0.000703836,
    r_C_fruit_2_C_litter3=0.00026986,
    r_C_fruit_2_C_litter5=0.00000463,
    r_C_litter1_2_C_som1=0.000001348,
    r_C_litter1_2_C_som2=0.0000192111,
    r_C_litter2_2_C_som2=0.00058924,
    r_C_litter2_2_C_som3=0.000645606,
    r_C_litter3_2_C_som1=0.000183481,
    r_C_litter3_2_C_som3=0.0001882826,
    r_C_litter4_2_C_som1=0.0001623,
    r_C_litter4_2_C_som2=0.0007858,
    r_C_litter5_2_C_som1=0.0000005513,
    r_C_litter6_2_C_som2=0.0829,
    r_C_som1_2_C_som3=0.000414588,
    r_C_som2_2_C_som3=0.005084,
    r_C_som2_2_C_som4=0.0000151975,
    r_C_som3_2_C_som2=0.0000509621,
    r_C_som3_2_C_som4=0.00000000396756,
    r_C_som4_2_C_som2=0.0004856444,
    r_C_som4_2_C_som3=0.00,

    C_wood1_0=0.1675,
    C_wood2_0=0.4386777,
    C_wood3_0=3.2784,
    C_wood4_0=0.0293,
    C_leaf_0=0.283,
    C_root_0=0.039,
    C_litter1_0=0.1177,
    C_litter2_0=0.0903,
    C_litter3_0=0.0243,
    C_litter4_0=0.01146,
    C_litter5_0=1.9707,
    C_som1_0=0.0479,
    C_som2_0=0.68511,
    C_som3_0=6.8025,
)
# -

from sympy import Symbol

par_dict = {
    Symbol(k): v for k, v in
    {
        "beta_wood1": 0.25,
        "beta_wood2": 0.1,
        "beta_leaf": 0.4,
        "beta_root": 0.15,
        #       "beta_fruit": 0.1,
        "T_0": 2,
        "E": 4,
        "KM": 10,
        "r_C_litter1_rh": 0.0000000000275,
        "r_C_litter2_rh": 0.000000000228,
        "r_C_litter3_rh": 0.000000000143,
        "r_C_litter4_rh": 0.0000000000441,
        "r_C_litter5_rh": 0.000000000119,
        "r_C_litter6_rh": 0.000000000177,
        "r_C_som1_rh": 0.000000000223,
        "r_C_som2_rh": 0.000000000116,
        "r_C_som3_rh": 0.00000000281,
        "r_C_som4_rh": 0.00000000476,
        "r_C_wood1_2_C_wood3": 0.00000000134,
        "r_C_wood1_2_C_litter1": 0.0000234,
        "r_C_wood2_2_C_wood4": 0.0000148,
        "r_C_wood2_2_C_litter2": 0.00000154,
        "r_C_wood3_2_C_litter1": 0.0000169,
        "r_C_wood4_2_C_litter2": 0.0000000253,
        "r_C_leaf_2_C_litter3": 0.000000719,
        "r_C_leaf_2_C_litter5": 0.00000945,
        "r_C_root_2_C_litter4": 0.00000743,
        "r_C_root_2_C_litter6": 0.00000212,
        "r_C_fruit_2_C_litter3": 0.0000159,
        "r_C_fruit_2_C_litter5": 0.000000728,
        "r_C_litter1_2_C_som1": 0.00001,
        "r_C_litter1_2_C_som2": 0.00000143,
        "r_C_litter2_2_C_som2": 0.0000185,
        "r_C_litter2_2_C_som3": 0.00000635,
        "r_C_litter3_2_C_som1": 0.00000117,
        "r_C_litter3_2_C_som3": 0.00000000659,
        "r_C_litter4_2_C_som1": 0.00000446,
        "r_C_litter4_2_C_som2": 0.00000270,
        "r_C_litter5_2_C_som1": 0.0000288,
        "r_C_litter6_2_C_som2": 0.0000000747,
        "r_C_som1_2_C_som3": 0.00000363,

        "r_C_som2_2_C_som3": 0.0074,
        "r_C_som2_2_C_som4": 0.01278,
        "r_C_som3_2_C_som2": 0.02027,
        "r_C_som3_2_C_som4": 0.008349,
        "r_C_som4_2_C_som2": 0.018376,
        "r_C_som4_2_C_som3": 0.0,
    }.items()
}

# +
import numpy as np

# svs_0=msh.Observables(*map(lambda v: v[0],svs))

X_0 = np.array((
    svs_0.cVeg / 7,
    svs_0.cVeg / 7,
    svs_0.cVeg / 7,
    svs_0.cVeg / 7,
    svs_0.cVeg / 7,
    svs_0.cVeg / 7,
    svs_0.cVeg / 7,
    svs_0.cLitter / 6,
    svs_0.cLitter / 6,
    svs_0.cLitter / 6,
    svs_0.cLitter / 6,
    svs_0.cLitter / 6,
    svs_0.cLitter / 6,
    svs_0.cSoil / 4,
    svs_0.cSoil / 4,
    svs_0.cSoil / 4,
    svs_0.cSoil / 4,

))  # .reshape(9,)

# +
import sys

sys.path.insert(0, '..')  # necessary to import general_helpers
import general_helpers as gh
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2


def make_steady_state_iterator_sym(
        mvs,
        V_init,
        par_dict,
        func_dict
):
    B_func, u_func = gh.make_B_u_funcs_2(mvs, par_dict, func_dict)

    def f(it, tup):
        X, _, _ = tup
        b = u_func(it, X)
        B = B_func(it, X)
        return (X, b, B)

    return TimeStepIterator2(
        initial_values=V_init,
        f=f)


# calculate steady state
func_dict = msh.make_func_dict(svs, dvs)
B_func, u_func = gh.make_B_u_funcs_2(mvs, par_dict, func_dict)

# +
it_sym = make_steady_state_iterator_sym(
    mvs,
    V_init=(X_0, u_func(0, X_0), B_func(0, X_0)),
    par_dict=par_dict,
    func_dict=func_dict
)

Bs = []
bs = []
max_day = (cpa.number_of_months - 1) * 30
print(max_day, gh.day_2_month_index(max_day))

for i in range(max_day):
    res = it_sym.__next__()
    bs.append(res[1])
    Bs.append(res[2])
B_mean = np.stack(Bs).mean(axis=0)
b_mean = np.stack(bs).mean(axis=0)
B_mean, b_mean
np.linalg.inv(Bs[0])

# calculate pseudo steady state
X_ss = np.linalg.solve(B_mean, (-b_mean))

steady_state_dict = {str(name): X_ss[i, 0] for i, name in enumerate(mvs.get_StateVariableTuple())}
# -

# svs
# cpa
print('Steadt_state_dict:', steady_state_dict)

# create a start parameter tuple for the mcmc. The order has to be the same as when you created the namedtupl3
# If you don't you get a "TypeError".

# -


# +
# svs_0.cLitter

# +
# now test it

import matplotlib.pyplot as plt
from general_helpers import plot_solutions

param2res_sym = msh.make_param2res_sym(mvs, cpa, dvs)

print(type(param2res_sym))

obs_0 = param2res_sym(epa_0)
# obs=np.column_stack([ np.array(v) for v in svs])
# obs=np.column_stack((np.repeat(svs.cVeg, 12),np.repeat(svs.cLitter, 12),np.repeat(svs.cSoil, 12),svs.rh,svs.ra))
# xs.shape

# +

# obs=obs[0:cpa.number_of_months,:] #cut
# obs[:,3:4]=obs[:,3:4]
# n=cpa.number_of_months

# convert to yearly output if necessary
# obs_yr=np.zeros(int(cpa.number_of_months/12)*obs.shape[1]).reshape([int(cpa.number_of_months/12),obs.shape[1]])
# for i in range(obs.shape[1]):
#    obs_yr[:,i]=gh.monthly_to_yearly(obs[:,i])
# obs=obs_yr
# n=int(cpa.number_of_months/12)

out_simu_d = obs_0._asdict()
obs_d = svs._asdict()

print("Forward run with initial parameters (blue) vs TRENDY output (orange)")
obs_d['cVeg'], svs.cVeg

fig = plt.figure(figsize=(10, 50))
axs = fig.subplots(5, 1)

for ind, f in enumerate(msh.Observables._fields):
    val_sim = out_simu_d[f]
    val_obs = obs_d[f]
    axs[ind].plot(range(len(val_obs)), val_obs, label=f + "_obs")
    axs[ind].plot(range(len(val_sim)), val_sim, label=f + "_sim")
    axs[ind].legend()

# fig = plt.figure(figsize=(12, 4), dpi=80)
# plot_solutions(
#         fig,
#         times=range(n),
#         var_names=msh.Observables._fields,
#         tup=(xs,obs)
#         #tup=(obs,)
# )
fig.savefig('solutions.pdf')

# +

epa_min = msh.EstimatedParameters(
    beta_wood1=0.01,
    beta_wood2=0.01,
    beta_leaf=0.01,
    beta_root=0.01,
    #   beta_fruit=0.1,
    T_0=1,
    E=1,
    KM=1,
    # r_C_leaf_rh=0,
    # r_C_wood_rh=0,
    # r_C_root_rh=0,

    r_C_litter1_rh=0.0012934/100,
    r_C_litter2_rh=0.59747/100,
    r_C_litter3_rh=0.01977/100,
    r_C_litter4_rh=0.0006258/100,
    r_C_litter5_rh=0.0000016211/100,
    r_C_litter6_rh=0.002435/100,
    r_C_som1_rh=0.00229726/100,
    r_C_som2_rh=0.0010939/100,
    r_C_som3_rh=0.0000455408/100,
    r_C_som4_rh=0.00028037/100,

    r_C_wood1_2_C_wood3=0.0107/100,
    r_C_wood1_2_C_litter1=0.066445/100,
    r_C_wood2_2_C_wood4=0.00063076/100,
    r_C_wood2_2_C_litter2=0.00168526/100,
    r_C_wood3_2_C_litter1=0.000007645231/100,
    r_C_wood4_2_C_litter2=0.000209176/100,
    r_C_leaf_2_C_litter3=0.0020595/100,
    r_C_leaf_2_C_litter5=0.0000250875/100,
    r_C_root_2_C_litter4=0.0000076896/100,
    r_C_root_2_C_litter6=0.000703836/100,
    r_C_fruit_2_C_litter3=0.00026986/100,
    r_C_fruit_2_C_litter5=0.00000463/100,
    r_C_litter1_2_C_som1=0.000001348/100,
    r_C_litter1_2_C_som2=0.0000192111/100,
    r_C_litter2_2_C_som2=0.00058924/100,
    r_C_litter2_2_C_som3=0.000645606/100,
    r_C_litter3_2_C_som1=0.000183481/100,
    r_C_litter3_2_C_som3=0.0001882826/100,
    r_C_litter4_2_C_som1=0.0001623/100,
    r_C_litter4_2_C_som2=0.0007858/100,
    r_C_litter5_2_C_som1=0.0000005513/100,
    r_C_litter6_2_C_som2=0.0829/100,
    r_C_som1_2_C_som3=0.000414588/100,
    r_C_som2_2_C_som3=0.005084/100,
    r_C_som2_2_C_som4=0.0000151975/100,
    r_C_som3_2_C_som2=0.0000509621/100,
    r_C_som3_2_C_som4=0.00000000396756/100,
    r_C_som4_2_C_som2=0.0004856444/100,
    r_C_som4_2_C_som3=0.00,


    C_wood1_0=0.016,
    C_wood2_0=0.0072,
    C_wood3_0=1,
    C_wood4_0=0.001,
    C_leaf_0=0.002,
    C_root_0=0.0013,
    C_litter1_0=0.0005,
    C_litter2_0=0.0015,
    C_litter3_0=0.0001,
    C_litter4_0=0.00002,
    C_litter5_0=0.1,
    C_som1_0=0.0001,
    C_som2_0=0.0001,
    C_som3_0=1,
)

epa_max = msh.EstimatedParameters(
    beta_wood1=0.99,
    beta_wood2=0.99,
    beta_leaf=0.99,
    beta_root=0.99,
    #   beta_fruit=0.1,
    T_0=10,
    E=10,
    KM=100,

    r_C_litter1_rh=0.0012934*100,
    r_C_litter2_rh=0.59747*100,
    r_C_litter3_rh=0.01977*100,
    r_C_litter4_rh=0.0006258*100,
    r_C_litter5_rh=0.0000016211*100,
    r_C_litter6_rh=0.002435*100,
    r_C_som1_rh=0.00229726*100,
    r_C_som2_rh=0.0010939*100,
    r_C_som3_rh=0.0000455408*100,
    r_C_som4_rh=0.00028037*100,

    r_C_wood1_2_C_wood3=0.0107*100,
    r_C_wood1_2_C_litter1=0.066445*100,
    r_C_wood2_2_C_wood4=0.00063076*100,
    r_C_wood2_2_C_litter2=0.00168526*100,
    r_C_wood3_2_C_litter1=0.000007645231*100,
    r_C_wood4_2_C_litter2=0.000209176*100,
    r_C_leaf_2_C_litter3=0.0020595*100,
    r_C_leaf_2_C_litter5=0.0000250875*100,
    r_C_root_2_C_litter4=0.0000076896*100,
    r_C_root_2_C_litter6=0.000703836*100,
    r_C_fruit_2_C_litter3=0.00026986*100,
    r_C_fruit_2_C_litter5=0.00000463*100,
    r_C_litter1_2_C_som1=0.000001348*100,
    r_C_litter1_2_C_som2=0.0000192111*100,
    r_C_litter2_2_C_som2=0.00058924*100,
    r_C_litter2_2_C_som3=0.000645606*100,
    r_C_litter3_2_C_som1=0.000183481*100,
    r_C_litter3_2_C_som3=0.0001882826*100,
    r_C_litter4_2_C_som1=0.0001623*100,
    r_C_litter4_2_C_som2=0.0007858*100,
    r_C_litter5_2_C_som1=0.0000005513*100,
    r_C_litter6_2_C_som2=0.0829*100,
    r_C_som1_2_C_som3=0.000414588*100,
    r_C_som2_2_C_som3=0.005084*100,
    r_C_som2_2_C_som4=0.0000151975*100,
    r_C_som3_2_C_som2=0.0000509621*100,
    r_C_som3_2_C_som4=0.00000000396756*100,
    r_C_som4_2_C_som2=0.0004856444*100,
    r_C_som4_2_C_som3=0.00,

    C_wood1_0=5,
    C_wood2_0=5,
    C_wood3_0=10,
    C_wood4_0=5,
    C_leaf_0=5,
    C_root_0=5,
    C_litter1_0=1,
    C_litter2_0=1,
    C_litter3_0=1,
    C_litter4_0=1,
    C_litter5_0=17,
    C_som1_0=10,
    C_som2_0=10,
    C_som3_0=20,
)

# +
# import test_helpers as th
# ta=th.make_test_args(conf_dict,msh,mvs)
# -

# ### mcmc to optimize parameters
#

# +

from general_helpers import autostep_mcmc, make_param_filter_func, make_feng_cost_func

isQualified = make_param_filter_func(epa_max, epa_min)
param2res = msh.make_param2res_sym(mvs, cpa, dvs)

print("Starting data assimilation...")
# Autostep MCMC: with uniform proposer modifying its step every 100 iterations depending on acceptance rate
C_autostep, J_autostep = autostep_mcmc(
    initial_parameters=epa_0,
    filter_func=isQualified,
    param2res=param2res,
    # costfunction=make_feng_cost_func(obs),
    costfunction=msh.make_weighted_cost_func(svs),
    nsimu=1500,  # for testing and tuning mcmc
    # nsimu=20000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=40,  # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=10,  # default value | increase value to reduce initial step size
    K=2  # default value | increase value to reduce acceptance of higher cost functions
)
print("Data assimilation finished!")

# +
# optimized parameter set (lowest cost function)
par_opt = np.min(
    C_autostep[:, np.where(J_autostep[1] == np.min(J_autostep[1]))].reshape(len(msh.EstimatedParameters._fields), 1),
    axis=1)
epa_opt = msh.EstimatedParameters(*par_opt)
mod_opt = param2res(epa_opt)
DA = mod_opt._asdict()
obs_d = svs._asdict()
print("Forward run with optimized parameters (orange) vs TRENDY output (blue)")
fig = plt.figure(figsize=(4, 10), dpi=80)
axs_DA = fig.subplots(5, 1)

for a, b in enumerate(msh.Observables._fields):
    val_DA = DA[b]
    val_obs = obs_d[b]
    axs_DA[a].plot(range(len(val_obs)), val_obs, label=b + "_obs")
    axs_DA[a].plot(range(len(val_DA)), val_DA, label=b + "_DA")
    axs_DA[a].legend()

fig.savefig('solution_DA.pdf')

fig.savefig('solutions_opt.pdf')

# save the parameters and cost function values for postprocessing
outputPath = Path(conf_dict["dataPath"])  # save output to data directory (or change it)

import pandas as pd

pd.DataFrame(C_autostep).to_csv(outputPath.joinpath('OCN-V2_da_aa.csv'), sep=',')
pd.DataFrame(J_autostep).to_csv(outputPath.joinpath('OCN-V2_j_aa.csv'), sep=',')
pd.DataFrame(epa_opt).to_csv(outputPath.joinpath('OCN-V2_optimized_pars.csv'), sep=',')
pd.DataFrame(mod_opt).to_csv(outputPath.joinpath('OCN-V2_optimized_solutions.csv'), sep=',')
# -

print("Optimized parameters: ", epa_opt)

# +
# Traceability Analysis
it_sym_trace = msh.make_traceability_iterator(mvs, dvs, cpa, epa_opt)
ns = 1500  # 10*360
StartVectorTrace = gh.make_StartVectorTrace(mvs)
nv = len(StartVectorTrace._fields)
res_trace = np.zeros((ns, nv))
for i in range(ns):
    res_trace[i, :] = it_sym_trace.__next__().reshape(nv)
# res_trace

import matplotlib.pyplot as plt

n = len(mvs.get_StateVariableTuple())
fig = plt.figure(figsize=(20, (n + 1) * 10), dpi=80)
axs = fig.subplots(n + 1, 2)
days = list(range(ns))

for i in range(n):
    ax = axs[i, 0]
    #  the solution
    pos = i
    ax.plot(
        days,
        res_trace[:, i],
        label=StartVectorTrace._fields[pos],
        color='blue'
    )
    # X_p
    pos = i + n
    ax.plot(
        days,
        res_trace[:, pos],
        label=StartVectorTrace._fields[pos],
        color='red'
    )
    # X_c
    pos = i + 2 * n
    ax.plot(
        days,
        res_trace[:, pos],
        label=StartVectorTrace._fields[pos],
        color='yellow'
    )
    ax.legend()

    ax = axs[i, 1]
    # RT
    pos = i + 3 * n
    ax.plot(
        days,
        res_trace[:, pos],
        label=StartVectorTrace._fields[pos],
        color='black'
    )
    ax.legend()

axs[n, 0].plot(
    days,
    [msh.make_npp_func(dvs)(d) for d in days],
    label='NPP',
    color='green'
)
axs[n, 0].legend()

fig.savefig('solution_Traceability.pdf')
# -
