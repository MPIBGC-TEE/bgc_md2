# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
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
    #number_of_months=len(svs.rh)
    # number_of_months=120 # for testing and tuning mcmc
    number_of_months=3840 # for testing and tuning mcmc
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
    r_C_litter2_rh=0.0059747,
    r_C_litter3_rh=0.001977,
    r_C_litter4_rh=0.006258,
    r_C_litter5_rh=0.00016211,
    r_C_litter6_rh=0.002435,
    r_C_som1_rh=0.000229726,
    r_C_som2_rh=0.00010939,
    r_C_som3_rh=0.0000455408,
    r_C_som4_rh=0.00028037,

#     r_C_wood1_2_C_wood3=0.0107,
#     r_C_wood1_2_C_litter1=0.066445,
#     r_C_wood2_2_C_wood4=0.00063076,
#     r_C_wood2_2_C_litter2=0.00168526,
#     r_C_wood3_2_C_litter1=0.000007645231,
#     r_C_wood4_2_C_litter2=0.000209176,
#     r_C_leaf_2_C_litter3=0.0020595,
#     r_C_leaf_2_C_litter5=0.0000250875,
#     r_C_root_2_C_litter4=0.0000076896,
#     r_C_root_2_C_litter6=0.000703836,
#     r_C_fruit_2_C_litter3=0.00026986,
#     r_C_fruit_2_C_litter5=0.00000463,
#     r_C_litter1_2_C_som1=0.000001348,
#     r_C_litter1_2_C_som2=0.0000192111,
#     r_C_litter2_2_C_som2=0.00058924,
#     r_C_litter2_2_C_som3=0.000645606,
#     r_C_litter3_2_C_som1=0.000183481,
#     r_C_litter3_2_C_som3=0.0001882826,
#     r_C_litter4_2_C_som1=0.0001623,
#     r_C_litter4_2_C_som2=0.0007858,
#     r_C_litter5_2_C_som1=0.0000005513,
#     r_C_litter6_2_C_som2=0.0829,
#     r_C_som1_2_C_som3=0.000414588,
#     r_C_som2_2_C_som3=0.005084,
#     r_C_som2_2_C_som4=0.0000151975,
#     r_C_som3_2_C_som2=0.0000509621,
#     r_C_som3_2_C_som4=0.00000000396756,
#     r_C_som4_2_C_som2=0.0004856444,
#     r_C_som4_2_C_som3=0.00,
    
    r_C_wood1_2_C_wood3=0.00107,
    r_C_wood1_2_C_litter1=0.0066445,
    r_C_wood2_2_C_wood4=0.0063076,
    r_C_wood2_2_C_litter2=0.000168526,
    r_C_wood3_2_C_litter1=0.007645231,
    r_C_wood4_2_C_litter2=0.000209176,
    r_C_leaf_2_C_litter3=0.0020595,
    r_C_leaf_2_C_litter5=0.000250875,
    r_C_root_2_C_litter4=0.000076896,
    r_C_root_2_C_litter6=0.000703836,
    r_C_fruit_2_C_litter3=0.0026986,
    r_C_fruit_2_C_litter5=0.0000463,
    r_C_litter1_2_C_som1=0.00001348,
    r_C_litter1_2_C_som2=0.000192111,
    r_C_litter2_2_C_som2=0.0058924,
    r_C_litter2_2_C_som3=0.00645606,
    r_C_litter3_2_C_som1=0.00183481,
    r_C_litter3_2_C_som3=0.001882826,
    r_C_litter4_2_C_som1=0.001623,
    r_C_litter4_2_C_som2=0.007858,
    r_C_litter5_2_C_som1=0.00005513,
    r_C_litter6_2_C_som2=0.0829,
    r_C_som1_2_C_som3=0.000414588,
    r_C_som2_2_C_som3=0.005084,
    r_C_som2_2_C_som4=0.0000151975,
    r_C_som3_2_C_som2=0.0000509621,
    r_C_som3_2_C_som4=0.0000000396756,
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

# +
# par_dict = {
#     Symbol(k): v for k, v in
#     {
#         "beta_wood1": 0.22877154801496366,
#         "beta_wood2": 0.27105886616734787, 
#         "beta_leaf": 0.22107830601036774,
#         "beta_root": 0.15714387838179705,
#         #beta_fruit": 0.1,
#         "T_0": 2,
#         "E": 4,
#         "KM": 10,
#         "r_C_litter1_rh": 0.00010671175898322236, 
#         "r_C_litter2_rh": 0.010178710800936974, 
#         "r_C_litter3_rh": 0.015203061301603886, 
#         "r_C_litter4_rh": 0.016702123445176017, 
#         "r_C_litter5_rh": 0.00023514296836874357, 
#         "r_C_litter6_rh": 0.008289055384865513,
#         "r_C_som1_rh": 0.0007171384519944395, 
#         "r_C_som2_rh": 0.00022301039347960097, 
#         "r_C_som3_rh": 0.0002825022995321478, 
#         "r_C_som4_rh": 0.0001162776599638064, 
#         "r_C_wood1_2_C_wood3": 0.004928919484099116, 
#         "r_C_wood1_2_C_litter1": 0.006256684737690185, 
#         "r_C_wood2_2_C_wood4": 0.00217953422350438, 
#         "r_C_wood2_2_C_litter2": 0.0018073911322216475, 
#         "r_C_wood3_2_C_litter1": 0.0384553289605736, 
#         "r_C_wood4_2_C_litter2": 5.431884592671132e-05, 
#         "r_C_leaf_2_C_litter3": 0.006284885870497801, 
#         "r_C_leaf_2_C_litter5": 0.00037806917561193106, 
#         "r_C_root_2_C_litter4": 0.00015681794955980385, 
#         "r_C_root_2_C_litter6": 0.00046519784183239924, 
#         "r_C_fruit_2_C_litter3": 0.004753571144341833, 
#         "r_C_fruit_2_C_litter5": 0.00014812848285782678, 
#         "r_C_litter1_2_C_som1": 5.025438958467177e-05, 
#         "r_C_litter1_2_C_som2": 0.00044648860679763177, 
#         "r_C_litter2_2_C_som2": 0.0156547203735249, 
#         "r_C_litter2_2_C_som3": 0.011965105209883626, 
#         "r_C_litter3_2_C_som1": 0.0014996323616631663, 
#         "r_C_litter3_2_C_som3": 0.009082763807079477, 
#         "r_C_litter4_2_C_som1": 0.008895958091190646, 
#         "r_C_litter4_2_C_som2": 0.070351520657772, 
#         "r_C_litter5_2_C_som1": 0.00026601269284842814, 
#         "r_C_litter6_2_C_som2": 0.06559236914167413, 
#         "r_C_som1_2_C_som3": 0.0034632100878440554,
#         "r_C_som2_2_C_som3": 0.02752553514851257, 
#         "r_C_som2_2_C_som4": 9.401293235880684e-05, 
#         "r_C_som3_2_C_som2": 0.00024267005147667766, 
#         "r_C_som3_2_C_som4": 6.8351229369685e-08, 
#         "r_C_som4_2_C_som2": 0.0013613403610407717, 
#         "r_C_som4_2_C_som3": 0.0, 
#     }.items()
# }

# +
# # updated version
# epa_0=msh.EstimatedParameters(
#     beta_wood1=0.22877154801496366, 
#     beta_wood2=0.27105886616734787, 
#     beta_leaf=0.22107830601036774, 
#     beta_root=0.15714387838179705, 
#     T_0=1.9567731838557805, 
#     E=3.503132360776974, 
#     KM=5.75856954136544, 
#     r_C_wood1_2_C_wood3=0.004928919484099116, 
#     r_C_wood1_2_C_litter1=0.006256684737690185,
#     r_C_wood2_2_C_wood4=0.0217953422350438, 
#     r_C_wood2_2_C_litter2=0.0018073911322216475, 
#     r_C_wood3_2_C_litter1=0.0384553289605736, 
#     r_C_wood4_2_C_litter2=5.431884592671132e-05, 
#     r_C_leaf_2_C_litter3=0.006284885870497801, 
#     r_C_leaf_2_C_litter5=0.00037806917561193106, 
#     r_C_root_2_C_litter4=0.00015681794955980385, 
#     r_C_root_2_C_litter6=0.00046519784183239924, 
#     r_C_fruit_2_C_litter3=0.004753571144341833, 
#     r_C_fruit_2_C_litter5=0.00014812848285782678, 
#     r_C_litter1_2_C_som1=5.025438958467177e-05, 
#     r_C_litter1_2_C_som2=0.00044648860679763177, 
#     r_C_litter2_2_C_som2=0.0156547203735249, 
#     r_C_litter2_2_C_som3=0.011965105209883626, 
#     r_C_litter3_2_C_som1=0.0014996323616631663, 
#     r_C_litter3_2_C_som3=0.009082763807079477, 
#     r_C_litter4_2_C_som1=0.008895958091190646, 
#     r_C_litter4_2_C_som2=0.070351520657772, 
#     r_C_litter5_2_C_som1=0.00026601269284842814, 
#     r_C_litter6_2_C_som2=0.06559236914167413, 
#     r_C_som1_2_C_som3=0.0034632100878440554, 
#     r_C_som2_2_C_som3=0.02752553514851257, 
#     r_C_som2_2_C_som4=9.401293235880684e-05, 
#     r_C_som3_2_C_som2=0.00024267005147667766, 
#     r_C_som3_2_C_som4=6.8351229369685e-08, 
#     r_C_som4_2_C_som2=0.0013613403610407717, 
#     r_C_som4_2_C_som3=0.0, 
#     r_C_litter1_rh=0.00010671175898322236, 
#     r_C_litter2_rh=0.010178710800936974, 
#     r_C_litter3_rh=0.015203061301603886, 
#     r_C_litter4_rh=0.016702123445176017, 
#     r_C_litter5_rh=0.00023514296836874357, 
#     r_C_litter6_rh=0.008289055384865513, 
#     r_C_som1_rh=0.0007171384519944395, 
#     r_C_som2_rh=0.00022301039347960097, 
#     r_C_som3_rh=0.0002825022995321478, 
#     r_C_som4_rh=0.0001162776599638064, 
#     C_wood1_0=0.14945959558401622, 
#     C_wood2_0=0.09999039520811101, 
#     C_wood3_0=3.4279199917419625, 
#     C_wood4_0=0.09057172587032938, 
#     C_leaf_0=0.39570986570575556, 
#     C_root_0=0.16841302725748428, 
#     C_litter1_0=0.12614309296839382, 
#     C_litter2_0=0.17328448144838912, 
#     C_litter3_0=0.004934392639447546, 
#     C_litter4_0=0.028918928860438016, 
#     C_litter5_0=1.9534741480016722, 
#     C_som1_0=0.22620867532497674, 
#     C_som2_0=0.6855397390527145, 
#     C_som3_0=6.542657926539519
# )

epa_0=msh.EstimatedParameters(
    beta_wood1=0.2183040322253002, 
    beta_wood2=0.26863096288181115, 
    beta_leaf=0.21943807491798112, 
    beta_root=0.16040533474423047, 
    T_0=1.950717867245791, 
    E=3.5293391839615964, 
    KM=5.159397430100239, 
    r_C_wood1_2_C_wood3=0.006981115219153814, 
    r_C_wood1_2_C_litter1=0.008361941501650308, 
    r_C_wood2_2_C_wood4=0.02473652272277666, 
    r_C_wood2_2_C_litter2=0.0012921776077177455, 
    r_C_wood3_2_C_litter1=0.021703768513993506, 
    r_C_wood4_2_C_litter2=0.00010215013792353313, 
    r_C_leaf_2_C_litter3=0.012906426484573155, 
    r_C_leaf_2_C_litter5=0.0001362984730925845, 
    r_C_root_2_C_litter4=7.511838182173387e-05, 
    r_C_root_2_C_litter6=0.000139216665393605, 
    r_C_fruit_2_C_litter3=0.005327641304999624, 
    r_C_fruit_2_C_litter5=0.0001564841236483191, 
    r_C_litter1_2_C_som1=8.065742707223306e-05, 
    r_C_litter1_2_C_som2=0.0001945680038686853, 
    r_C_litter2_2_C_som2=0.021855567720853678, 
    r_C_litter2_2_C_som3=0.008862640704241595, 
    r_C_litter3_2_C_som1=0.00044753314434899067, 
    r_C_litter3_2_C_som3=0.002920218828837676, 
    r_C_litter4_2_C_som1=0.0013627372863071966, 
    r_C_litter4_2_C_som2=0.08667028867036282, 
    r_C_litter5_2_C_som1=0.0003058160346834833, 
    r_C_litter6_2_C_som2=0.03198987317695222, 
    r_C_som1_2_C_som3=0.004199397889363635, 
    r_C_som2_2_C_som3=0.02386802237580712, 
    r_C_som2_2_C_som4=4.926778642156811e-05, 
    r_C_som3_2_C_som2=2.808847983166608e-05, 
    r_C_som3_2_C_som4=3.166786733299793e-08, 
    r_C_som4_2_C_som2=0.0030270107698128577, 
    r_C_som4_2_C_som3=0.0, 
    r_C_litter1_rh=8.999935134317453e-05, 
    r_C_litter2_rh=0.037467962743264034, 
    r_C_litter3_rh=0.06307149375148956, 
    r_C_litter4_rh=0.021679193567117776, 
    r_C_litter5_rh=0.0004665369868242155, 
    r_C_litter6_rh=0.015307799551151899, 
    r_C_som1_rh=0.00016380762353055268, 
    r_C_som2_rh=0.00014083989056868948, 
    r_C_som3_rh=0.00018008218714862852, 
    r_C_som4_rh=3.1578556336141246e-05, 
    C_wood1_0=0.12249144169921905, 
    C_wood2_0=0.07827830341457838, 
    C_wood3_0=3.378036354758859, 
    C_wood4_0=0.08151041816435299, 
    C_leaf_0=0.3994811118842396, 
    C_root_0=0.14114710798737556, 
    C_litter1_0=0.11245652183364362, 
    C_litter2_0=0.17008440913158368, 
    C_litter3_0=0.019144268016656297, 
    C_litter4_0=0.06666712778703499, 
    C_litter5_0=1.9444305717528556, 
    C_som1_0=0.2729065826369784, 
    C_som2_0=0.6556404511921593, 
    C_som3_0=6.4814853538118955
)

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
# -

apa = {**cpa._asdict(), **epa_0._asdict()}
par_dict = {
    Symbol(k): v for k, v in apa.items() if Symbol(k) in par_dict
}
par_dict

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
func_dict = msh.make_func_dict(svs, dvs, cpa, epa_0)
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


svs_0.cVeg

# +
# deriving initail pool sizes proportional to steady state 
C_Veg_steady=steady_state_dict.get("C_wood1")+steady_state_dict.get("C_wood2")+steady_state_dict.get("C_wood3")
+steady_state_dict.get("C_wood4")+steady_state_dict.get("C_leaf")+steady_state_dict.get("C_root")
+steady_state_dict.get("C_fruit")

C_wood1_0_st=svs_0.cVeg*steady_state_dict.get("C_wood1")/C_Veg_steady
C_wood2_0_st=svs_0.cVeg*steady_state_dict.get("C_wood2")/C_Veg_steady
C_wood3_0_st=svs_0.cVeg*steady_state_dict.get("C_wood3")/C_Veg_steady
C_wood4_0_st=svs_0.cVeg*steady_state_dict.get("C_wood4")/C_Veg_steady
C_leaf_0_st=svs_0.cVeg*steady_state_dict.get("C_leaf")/C_Veg_steady
C_root_0_st=svs_0.cVeg*steady_state_dict.get("C_root")/C_Veg_steady

C_Litter_steady=steady_state_dict.get("C_litter1")+steady_state_dict.get("C_litter2")+steady_state_dict.get("C_litter3")
+steady_state_dict.get("C_litter4")+steady_state_dict.get("C_litter5")+steady_state_dict.get("C_litter6")

C_litter1_0_st=svs_0.cLitter*steady_state_dict.get("C_litter1")/C_Litter_steady          
C_litter2_0_st=svs_0.cLitter*steady_state_dict.get("C_litter2")/C_Litter_steady             
C_litter3_0_st=svs_0.cLitter*steady_state_dict.get("C_litter3")/C_Litter_steady           
C_litter4_0_st=svs_0.cLitter*steady_state_dict.get("C_litter4")/C_Litter_steady           
C_litter5_0_st=svs_0.cLitter*steady_state_dict.get("C_litter5")/C_Litter_steady              

C_Soil_steady=steady_state_dict.get("C_som1")+steady_state_dict.get("C_som2")
+steady_state_dict.get("C_som3")+steady_state_dict.get("C_som4")
         
C_som1_0_st=svs_0.cSoil*steady_state_dict.get("C_som1")/C_Soil_steady             
C_som2_0_st=svs_0.cSoil*steady_state_dict.get("C_som2")/C_Soil_steady           
C_som3_0_st=svs_0.cSoil*steady_state_dict.get("C_som3")/C_Soil_steady  
# -

# this includes derive initial pool sizes
epa_1=msh.EstimatedParameters(
    beta_wood1=epa_0.beta_wood1, 
    beta_wood2=epa_0.beta_wood2, 
    beta_leaf=epa_0.beta_leaf, 
    beta_root=epa_0.beta_root, 
    T_0=epa_0.T_0, 
    E=epa_0.E, 
    KM=epa_0.KM, 
    r_C_wood1_2_C_wood3=epa_0.r_C_wood1_2_C_wood3, 
    r_C_wood1_2_C_litter1=epa_0.r_C_wood1_2_C_litter1,
    r_C_wood2_2_C_wood4=epa_0.r_C_wood2_2_C_wood4, 
    r_C_wood2_2_C_litter2=epa_0.r_C_wood2_2_C_litter2, 
    r_C_wood3_2_C_litter1=epa_0.r_C_wood3_2_C_litter1, 
    r_C_wood4_2_C_litter2=epa_0.r_C_wood4_2_C_litter2, 
    r_C_leaf_2_C_litter3=epa_0.r_C_leaf_2_C_litter3, 
    r_C_leaf_2_C_litter5=epa_0.r_C_leaf_2_C_litter5, 
    r_C_root_2_C_litter4=epa_0.r_C_root_2_C_litter4, 
    r_C_root_2_C_litter6=epa_0.r_C_root_2_C_litter6, 
    r_C_fruit_2_C_litter3=epa_0.r_C_fruit_2_C_litter3, 
    r_C_fruit_2_C_litter5=epa_0.r_C_fruit_2_C_litter5, 
    r_C_litter1_2_C_som1=epa_0.r_C_litter1_2_C_som1, 
    r_C_litter1_2_C_som2=epa_0.r_C_litter1_2_C_som2, 
    r_C_litter2_2_C_som2=epa_0.r_C_litter2_2_C_som2, 
    r_C_litter2_2_C_som3=epa_0.r_C_litter2_2_C_som3, 
    r_C_litter3_2_C_som1=epa_0.r_C_litter3_2_C_som1, 
    r_C_litter3_2_C_som3=epa_0.r_C_litter3_2_C_som3, 
    r_C_litter4_2_C_som1=epa_0.r_C_litter4_2_C_som1, 
    r_C_litter4_2_C_som2=epa_0.r_C_litter4_2_C_som2, 
    r_C_litter5_2_C_som1=epa_0.r_C_litter5_2_C_som1, 
    r_C_litter6_2_C_som2=epa_0.r_C_litter6_2_C_som2, 
    r_C_som1_2_C_som3=epa_0.r_C_som1_2_C_som3, 
    r_C_som2_2_C_som3=epa_0.r_C_som2_2_C_som3, 
    r_C_som2_2_C_som4=epa_0.r_C_som2_2_C_som4, 
    r_C_som3_2_C_som2=epa_0.r_C_som3_2_C_som2, 
    r_C_som3_2_C_som4=epa_0.r_C_som3_2_C_som4, 
    r_C_som4_2_C_som2=epa_0.r_C_som4_2_C_som2, 
    r_C_som4_2_C_som3=epa_0.r_C_som4_2_C_som3, 
    r_C_litter1_rh=epa_0.r_C_litter1_rh, 
    r_C_litter2_rh=epa_0.r_C_litter2_rh, 
    r_C_litter3_rh=epa_0.r_C_litter3_rh, 
    r_C_litter4_rh=epa_0.r_C_litter4_rh, 
    r_C_litter5_rh=epa_0.r_C_litter5_rh, 
    r_C_litter6_rh=epa_0.r_C_litter6_rh, 
    r_C_som1_rh=epa_0.r_C_som1_rh, 
    r_C_som2_rh=epa_0.r_C_som2_rh, 
    r_C_som3_rh=epa_0.r_C_som3_rh, 
    r_C_som4_rh=epa_0.r_C_som4_rh, 
    C_wood1_0=C_wood1_0_st, 
    C_wood2_0=C_wood2_0_st, 
    C_wood3_0=C_wood3_0_st, 
    C_wood4_0=C_wood4_0_st, 
    C_leaf_0=C_leaf_0_st, 
    C_root_0=C_root_0_st, 
    C_litter1_0=C_litter1_0_st, 
    C_litter2_0=C_litter2_0_st, 
    C_litter3_0=C_litter3_0_st, 
    C_litter4_0=C_litter4_0_st, 
    C_litter5_0=C_litter5_0_st, 
    C_som1_0=C_som1_0_st, 
    C_som2_0=C_som2_0_st, 
    C_som3_0=C_som3_0_st,
)

# +
# now test it

import matplotlib.pyplot as plt
from general_helpers import plot_solutions

param2res_sym = msh.make_param2res_sym(mvs, cpa, dvs)

#print(type(param2res_sym))

obs_0 = param2res_sym(epa_0)
# obs=np.column_stack([ np.array(v) for v in svs])
# obs=np.column_stack((np.repeat(svs.cVeg, 12),np.repeat(svs.cLitter, 12),np.repeat(svs.cSoil, 12),svs.rh,svs.ra))
# xs.shape
obs_0[0].shape

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

    r_C_litter1_rh=epa_0.r_C_litter1_rh/100,
    r_C_litter2_rh=epa_0.r_C_litter2_rh/100,
    r_C_litter3_rh=epa_0.r_C_litter3_rh/100,
    r_C_litter4_rh=epa_0.r_C_litter4_rh/100,
    r_C_litter5_rh=epa_0.r_C_litter5_rh/100,
    r_C_litter6_rh=epa_0.r_C_litter6_rh/100,
    r_C_som1_rh=epa_0.r_C_som1_rh/100,
    r_C_som2_rh=epa_0.r_C_som2_rh/100,
    r_C_som3_rh=epa_0.r_C_som3_rh/100,
    r_C_som4_rh=epa_0.r_C_som4_rh/100,

    r_C_wood1_2_C_wood3=epa_0.r_C_wood1_2_C_wood3/100,
    r_C_wood1_2_C_litter1=epa_0.r_C_wood1_2_C_litter1/100,
    r_C_wood2_2_C_wood4=epa_0.r_C_wood2_2_C_wood4/100,
    r_C_wood2_2_C_litter2=epa_0.r_C_wood2_2_C_litter2/100,
    r_C_wood3_2_C_litter1=epa_0.r_C_wood3_2_C_litter1/100,
    r_C_wood4_2_C_litter2=epa_0.r_C_wood4_2_C_litter2/100,
    r_C_leaf_2_C_litter3=epa_0.r_C_leaf_2_C_litter3/100,
    r_C_leaf_2_C_litter5=epa_0.r_C_leaf_2_C_litter5/100,
    r_C_root_2_C_litter4=epa_0.r_C_root_2_C_litter4/100,
    r_C_root_2_C_litter6=epa_0.r_C_root_2_C_litter6/100,
    r_C_fruit_2_C_litter3=epa_0.r_C_fruit_2_C_litter3/100,
    r_C_fruit_2_C_litter5=epa_0.r_C_fruit_2_C_litter5/100,
    r_C_litter1_2_C_som1=epa_0.r_C_litter1_2_C_som1/100,
    r_C_litter1_2_C_som2=epa_0.r_C_litter1_2_C_som2/100,
    r_C_litter2_2_C_som2=epa_0.r_C_litter2_2_C_som2/100,
    r_C_litter2_2_C_som3=epa_0.r_C_litter2_2_C_som3/100,
    r_C_litter3_2_C_som1=epa_0.r_C_litter3_2_C_som1/100,
    r_C_litter3_2_C_som3=epa_0.r_C_litter3_2_C_som3/100,
    r_C_litter4_2_C_som1=epa_0.r_C_litter4_2_C_som1/100,
    r_C_litter4_2_C_som2=epa_0.r_C_litter4_2_C_som2/100,
    r_C_litter5_2_C_som1=epa_0.r_C_litter5_2_C_som1/100,
    r_C_litter6_2_C_som2=epa_0.r_C_litter6_2_C_som2/100,
    r_C_som1_2_C_som3=epa_0.r_C_som1_2_C_som3/100,
    r_C_som2_2_C_som3=epa_0.r_C_som2_2_C_som3/100,
    r_C_som2_2_C_som4=epa_0.r_C_som2_2_C_som4/100,
    r_C_som3_2_C_som2=epa_0.r_C_som3_2_C_som2/100,
    r_C_som3_2_C_som4=epa_0.r_C_som3_2_C_som4/100,
    r_C_som4_2_C_som2=epa_0.r_C_som4_2_C_som2/100,
    r_C_som4_2_C_som3=0.00,

#     C_wood1_0=0.016,
#     C_wood2_0=0.0072,
#     C_wood3_0=1,
#     C_wood4_0=0.001,
#     C_leaf_0=0.002,
#     C_root_0=0.0013,
#     C_litter1_0=0.0005,
#     C_litter2_0=0.0015,
#     C_litter3_0=0.0001,
#     C_litter4_0=0.00002,
#     C_litter5_0=0.1,
#     C_som1_0=0.0001,
#     C_som2_0=0.0001,
#     C_som3_0=1,
    
    C_wood1_0=0.00001,
    C_wood2_0=0.00001,
    C_wood3_0=0.00001,
    C_wood4_0=0.00001,
    C_leaf_0=0.00001,
    C_root_0=0.00001,
    C_litter1_0=0.00001,
    C_litter2_0=0.00001,
    C_litter3_0=0.00001,
    C_litter4_0=0.00001,
    C_litter5_0=0.00001,
    C_som1_0=0.00001,
    C_som2_0=0.00001,
    C_som3_0=0.00001,
    
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

    r_C_litter1_rh=epa_0.r_C_litter1_rh*100,
    r_C_litter2_rh=epa_0.r_C_litter2_rh*100,
    r_C_litter3_rh=epa_0.r_C_litter3_rh*100,
    r_C_litter4_rh=epa_0.r_C_litter4_rh*100,
    r_C_litter5_rh=epa_0.r_C_litter5_rh*100,
    r_C_litter6_rh=epa_0.r_C_litter6_rh*100,
    r_C_som1_rh=epa_0.r_C_som1_rh*100,
    r_C_som2_rh=epa_0.r_C_som2_rh*100,
    r_C_som3_rh=epa_0.r_C_som3_rh*100,
    r_C_som4_rh=epa_0.r_C_som4_rh*100,

    r_C_wood1_2_C_wood3=epa_0.r_C_wood1_2_C_wood3*100,
    r_C_wood1_2_C_litter1=epa_0.r_C_wood1_2_C_litter1*100,
    r_C_wood2_2_C_wood4=epa_0.r_C_wood2_2_C_wood4*100,
    r_C_wood2_2_C_litter2=epa_0.r_C_wood2_2_C_litter2*100,
    r_C_wood3_2_C_litter1=epa_0.r_C_wood3_2_C_litter1*100,
    r_C_wood4_2_C_litter2=epa_0.r_C_wood4_2_C_litter2*100,
    r_C_leaf_2_C_litter3=epa_0.r_C_leaf_2_C_litter3*100,
    r_C_leaf_2_C_litter5=epa_0.r_C_leaf_2_C_litter5*100,
    r_C_root_2_C_litter4=epa_0.r_C_root_2_C_litter4*100,
    r_C_root_2_C_litter6=epa_0.r_C_root_2_C_litter6*100,
    r_C_fruit_2_C_litter3=epa_0.r_C_fruit_2_C_litter3*100,
    r_C_fruit_2_C_litter5=epa_0.r_C_fruit_2_C_litter5*100,
    r_C_litter1_2_C_som1=epa_0.r_C_litter1_2_C_som1*100,
    r_C_litter1_2_C_som2=epa_0.r_C_litter1_2_C_som2*100,
    r_C_litter2_2_C_som2=epa_0.r_C_litter2_2_C_som2*100,
    r_C_litter2_2_C_som3=epa_0.r_C_litter2_2_C_som3*100,
    r_C_litter3_2_C_som1=epa_0.r_C_litter3_2_C_som1*100,
    r_C_litter3_2_C_som3=epa_0.r_C_litter3_2_C_som3*100,
    r_C_litter4_2_C_som1=epa_0.r_C_litter4_2_C_som1*100,
    r_C_litter4_2_C_som2=epa_0.r_C_litter4_2_C_som2*100,
    r_C_litter5_2_C_som1=epa_0.r_C_litter5_2_C_som1*100,
    r_C_litter6_2_C_som2=epa_0.r_C_litter6_2_C_som2*100,
    r_C_som1_2_C_som3=epa_0.r_C_som1_2_C_som3*100,
    r_C_som2_2_C_som3=epa_0.r_C_som2_2_C_som3*100,
    r_C_som2_2_C_som4=epa_0.r_C_som2_2_C_som4*100,
    r_C_som3_2_C_som2=epa_0.r_C_som3_2_C_som2*100,
    r_C_som3_2_C_som4=epa_0.r_C_som3_2_C_som4*100,
    r_C_som4_2_C_som2=epa_0.r_C_som4_2_C_som2*100,
    r_C_som4_2_C_som3=0.00,

    C_wood1_0=svs_0.cVeg,
    C_wood2_0=svs_0.cVeg,
    C_wood3_0=svs_0.cVeg,
    C_wood4_0=svs_0.cVeg,
    C_leaf_0=svs_0.cVeg,
    C_root_0=svs_0.cVeg,
    C_litter1_0=svs_0.cLitter,
    C_litter2_0=svs_0.cLitter,
    C_litter3_0=svs_0.cLitter,
    C_litter4_0=svs_0.cLitter,
    C_litter5_0=svs_0.cLitter,
    C_som1_0=svs_0.cSoil,
    C_som2_0=svs_0.cSoil,
    C_som3_0=svs_0.cSoil,
)
# -

costfunction=msh.make_weighted_cost_func(svs)
print(costfunction(obs_0))
print(np.array(epa_0)-np.array(epa_min))
print(np.array(epa_max)-np.array(epa_0))

# +
# import test_helpers as th
# ta=th.make_test_args(conf_dict,msh,mvs)
# -

# ### mcmc to optimize parameters
#

svs.rh.shape[0]

cpa.number_of_months

# adding additional year to rh and ra (due to wrong files used)
# !!! Remove as soon as right files are used !!!
if svs.rh.shape[0]<cpa.number_of_months:
    last12rh=svs.rh[-12:]
    last12ra=svs.ra[-12:]
    rh=np.ma.append(svs.rh, last12rh, axis=None)
    ra=np.ma.append(svs.ra, last12ra, axis=None)
    #svs
    svs1=msh.Observables(
        cVeg=svs.cVeg,
        cLitter=svs.cLitter,
        cSoil=svs.cSoil,
        rh=rh,
        ra=ra,
        )
    svs=svs1
svs.rh.shape[0]

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
    #nsimu=200,  # for testing and tuning mcmc
    nsimu=5000,
    c_max=np.array(epa_max),
    c_min=np.array(epa_min),
    acceptance_rate=10,  # default value | target acceptance rate in %
    chunk_size=100,  # default value | number of iterations to calculate current acceptance ratio and update step size
    D_init=2,  # increase value to reduce initial step size
    K=1.5  # increase value to reduce acceptance of higher cost functions
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
