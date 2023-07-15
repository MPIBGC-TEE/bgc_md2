# +
# adjust the output to full width
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# make changes to imported files immidiately available 
# avoiding the need to reload (in most cases)
# #%load_ext autoreload
# #%autoreload 2
# -

# This illustrative notebook shows how to create a representation for a new model.
#
# For the sake of simplicity we assume here that we start with a description of pools and fluxes.
# This is the form the models are usually described in the literature.
# We will see that the framework can derive a matrix representation automatically. 
#
# It would also be possible to start from the matrices or even
# mix the two approaches. 
# We will point to some more advanced examples where you can see this in action. If you get familiar with the framework you will find many more combinations than we can cover. 
#
# ## Inspect a minimal model
#
# We will start with an extremely simple model.
# whose original you can find in 
# ```bash
# bgc_md2/src/bgc_md2/models/testVectorFree/source.py
# ```
# or via the following code:

import inspect
import bgc_md2.models.testVectorFree.source 
#print(inspect.getsource(bgc_md2.models.testVectorFree.source)) #uncomment to see the code

# # copy its contents into a new cell. and start (we can later save it under a new name) 
# It looks a bit like this.

# +
from sympy import var, Symbol, Function 
from ComputabilityGraphs.CMTVS import CMTVS
from ComputabilityGraphs.helpers import module_computers
#from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
from importlib import import_module

# Make a small dictionary for the variables we will use
sym_dict={
    "leaf": "vegegation leaf pool",
    "wood": "vegetation wood pool content",
    "I_wood": "Influx into vegetation wood pool",
    "k_wood_o": "out flux rate of wood pool",
    "k_leaf_2_wood": "constant internal flux rate from leaf to wood", 
    "k_wood_2_leaf": "constant internal flux rate from wood to leaf", 
}
# Make symbols from  the strings that we can later use in expressions  
# leaf, wood,...
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)
# Note:
# The use of exec to execute python code strings seems a bit funny 
# The normal way woud be: 
# leaf = Symbol("leaf")
# wood = Symbol("wood")
# ...
# we do it this way, because we wantthe dictionary for documentation later...

# We will also use some symbolic functions ("symbols" with an argument) 
func_dict={
    "I_leaf": "Influx into vegetation leaf pool",
    "k_leaf_o": "out flux rate of leaf pool",
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)
# Note:
# Again we want the dictionary for later

t=TimeSymbol("t")

mvs = CMTVS(
    {
        StateVariableTuple((leaf, wood)),
        t,
        InFluxesBySymbol({leaf: I_leaf(t), wood: I_wood}),
        OutFluxesBySymbol({leaf: k_leaf_o(t) * leaf, wood: k_wood_o * wood}),
        InternalFluxesBySymbol({(leaf, wood): k_leaf_2_wood * leaf, (wood, leaf): k_wood_2_leaf * wood}),
    },

    computers=module_computers(bgc_md2.resolve.computers)
)
# -

# The last statement in the code defines a variable `mvs` which is 
# an instance of CMTVS which stands for `C`onnected`M`ulti`T`ype`V`ariable`S`et".
# It contains information in two forms. 
# 1. Variables of certain types (like InFluxesBySymbol)
# 2. Computers, Functions that "connect" these Variables and to other results we did not specify but which can be computed.
#    
#
#
# Tasks:
# To see what it can do with this information add a new cell an type `mvs.` and press the `tab` key `->`. This will show you the available methods, in other words what can be computed from the provided information.
#
# I wanted to see the compartmental matrix and the inputs.

mvs.get_CompartmentalMatrix()


mvs.get_CompartmentalMatrix()

mvs.get_InputTuple()

# we can also print the whole mass balance equation
import bgc_md2.display_helpers as dh
dh.mass_balance_equation(mvs)

# we can also plot a picture
import bgc_md2.helper as h
h.compartmental_graph(mvs)

# ## Extend the minimal model
#
# ### Add more state variables and fluxes for soil pools
# Task:
# Add new pools for litter (slit), coarse woody debris (scwd) and soil organic matter (som)
# and fluxes from the vegetation pools to the new pools.
#
# To avoid (at the moment) unintelligibleerror messages from the depth of the package, go in small steps:
# - first add only a new symbol (e.g. for a new pool), 
# - then add the pool to StateVariableTuple
# - then add a new symbol for a fluxrate or a symbolic function  
# - then add the flux using the rates 
# - repeat with the next pool,flux
#
# The result will look some what like this:

# +
# Make a small dictionary for the variables we will use
sym_dict={
    "leaf": "vegegation leaf pool",
    "wood": "vegetation wood pool content",
    "lit":"soil litter pool",
    "som": "soil organic matter",
    "slowsom": "slowly decomposing soil organic matter",
    "cwd": "soil coarse woody debris",
    "k_wood_o": "out flux rate of wood pool",
    "k_som_o": "out flux rate of som pool",
    "k_slowsom_o ": "out flux rate of slowsom pool",
    "k_leaf_2_wood": "constant internal flux rate from leaf to wood", 
    "k_wood_2_leaf": "constant internal flux rate from wood to leaf", 
    "k_leaf_2_lit": "flux rate from leaf to litter",
    "k_wood_2_cwd": "flux rate from wood to coarse woody debris",
    "k_lit_2_som": "flux rate from litter to som",
    "k_cwd_2_som": "flux rate from coarse woody debris to som",
    "k_som_2_slowsom": "flux rate from coarse woody debris to som",
}
# Make symbols from  the strings that we can later use in expressions  
# leaf, wood,...
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# We will also use some symbolic functions ("symbols" with an argument) 
func_dict={
    "I_leaf": "Influx into vegetation leaf pool",
    "k_leaf_o": "out flux rate of leaf pool",
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")

mvs = CMTVS(
    {
        StateVariableTuple((leaf, wood, lit, cwd, som, slowsom)),
        t,
        InFluxesBySymbol({leaf: I_leaf(t)}),
        OutFluxesBySymbol({
            leaf: k_leaf_o(t) * leaf,
            wood: k_wood_o * wood,
            som: k_som_o * som,
            slowsom: k_slowsom_o * slowsom,
        }),
        InternalFluxesBySymbol({
            (leaf, wood): k_leaf_2_wood * leaf, 
            (wood, leaf): k_wood_2_leaf * wood,
            (leaf, lit): k_leaf_2_lit * leaf,
            (wood, cwd): k_wood_2_cwd * wood,
            (lit, som): k_lit_2_som *lit,
            (cwd, som): k_cwd_2_som *cwd,
            (som, slowsom): k_som_2_slowsom *som,
        }),
    },

    computers=module_computers(bgc_md2.resolve.computers)
)
# -

h.compartmental_graph(mvs)

# # Task: run the model

# We first look at the computable types (find that nothing numeric is among them)
mvs.computable_mvar_types()

# We then look at all types that could possibly computed and find one that looks as if it 
# sounds like a solution.
import ComputabilityGraphs.helpers as cgh
computers=module_computers(bgc_md2.resolve.computers)
types=cgh.all_mvars(computers)
types

# NumericSolutionArray sounds promising 
# We then use the graph library to find out what the missing information is and 
# finally produce it
mvs.jupyter_widget(
    root_type=bgc_md2.resolve.mvars.NumericSolutionArray
)

mvs.get_StateVariableTuple()

# +
# We look at the ligth red dots. These are the missing bits
# In this case
# - NumericParameterization,
# - NumericSimulationTimes,
# - NumericStartValueDict/ or NumericStartValueA
#  (We could also have also chosen some computable values (shortcutting the computation) )

# For a real model these values are essential, and hard to get
par_dict={
    Symbol(k):v for k,v in {
        "k_wood_o": 0.1,
        "k_som_o": 0.1,
        "k_slowsom_o ": 0.1,
        "k_leaf_2_wood": 0.1, 
        "k_wood_2_leaf": 0.1, 
        "k_leaf_2_lit": 0.1,
        "k_wood_2_cwd": 0.1,
        "k_lit_2_som": 0.1,
        "k_cwd_2_som": 0.1,
        "k_som_2_slowsom": 0.1,
    }.items()     
}
start_value_dict={
    Symbol(k):v for k,v in {
        "leaf": 1,
        "wood": 2,
        "lit": 3,
        "cwd": 4, 
        "som": 1, 
        "slowsom": 10,
    }.items()     
    
}
# in real life this would almost always be functions interpolated from data
# technically important is only the fact that they are normal python functions
# with the write signature (the sympolic functions are functions of time, so these functions
# will be called with the time argument in the ode solver)
func_dict={
    "I_leaf": lambda t: np.sin(t)+1,
    "k_leaf_o": lambda t: np.sin(t+np.pi)+1
}
from bgc_md2.resolve.mvars import ( 
    NumericParameterization,
    NumericSimulationTimes,
    NumericStartValueDict
)
import numpy as np
times=np.linspace(0,10,100)
mvs=mvs.update(
    {
        NumericStartValueDict(start_value_dict),
        NumericParameterization(
            par_dict=par_dict,
            func_dict=func_dict
        ),
        NumericSimulationTimes(times),
    }
)
# now we can compute more results including the one we want
mvs.computable_mvar_types()
# -

for v in mvs.provided_mvar_values:
    print(type(v))


sol_arr = mvs.get_NumericSolutionArray()
import matplotlib.pyplot as plt
sv=mvs.get_StateVariableTuple()
n_pools=len(sv)
fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
axs = fig1.subplots(n_pools, 1)
for i in range(n_pools):
    ax = axs[i]
    ax.plot(times, sol_arr[:, i], label=str(sv[i]))
    ax.legend()
    ax.set_title(f"{sv[i]} solution")

# ### Task compute the mean backward transit time 
# We proceed in the same way we look at all possilbe results (see list above) and look for a promising name. We then check what missing bits we need to provide to obtain this result. 

# This time NumericMeanBackwardTransitTimeSolutionArray sounds promising 
# We then use the graph library to find out what the missing information is and 
# finally produce it
mvs.jupyter_widget(
    root_type=bgc_md2.resolve.mvars.NumericMeanBackwardTransitTimeSolution
)

# This time the missing bit is NumericStartMeanAgeTuple
from bgc_md2.resolve.mvars import NumericStartMeanAgeTuple
# NumericStartMeanAgeVector? #uncomment to look up the help text
mvs=mvs.update( { NumericStartMeanAgeTuple([1,2,3,3,3,3]) } )
mvs.computable_mvar_types()

# now we can get it but also gain the meanages ...
mean_btts2 = mvs.get_NumericMeanBackwardTransitTimeSolution()
# a single line (one value per timestep)
mean_btts2.shape

mvs.provided_mvar_types

# ### Check
# In the previous cells we just invented start values for the pool contents and their mean ages independently.
# Actually there is a common source for both: the start age density function (mass/age as function of age), that we will need later when we compute the full time dependent age distribution.
# If we integrate this density over all ages we get the pool contents.If we compute it's mean we get the mean ages.
# So we could have chosen such an age density (for every pool) and derived compatible  zero'th and first moments.
#
# But the age density is still arbitraty, which creates possibility for contradiction with our assumptions about how we got to the starting position. 
# While one can imagine a system (with the same number of pools as ours) that exibits our arbitrarily chosen age density , we cannot produce a history of *our* system, that would lead to this combination. 
# Totally different rates and functions might have been necessary to create our starting conditions.
# If we want start values that are consistent with a certain kind of pre start history we have to make this assumption concrete first and consider the start condition as a consequence. Possible examples for such assumed pre start histories are:
# 1. no history at all, all pools start empty with zero content of age zero.
# 1. start empty but run for some time (we have to know (or make up) the time dependent function for this  time period)
# 1. freeze the time dependent part and compute an equilibrium if possible (This is easy for linear systems (if we freeze all time dependent parts) but not necessarily possible for nonlinear systems.
# For all three options we provide tooling but, it's important to be aware of the choice. 
#

# We choose option 3 here which relies on the fact that we deal with a linear system and is therefore not
# completly automated in the framework ( no mvs.get_Fixpoint )
# We would have to create special Types for linear systems and only for those offer a computation method:
# To avoid accidents we keep the fixpoint computation in user responsibility for now. 
import CompartmentalSystems.start_distributions as sd
t0=times[0]
srm = mvs.get_SmoothReservoirModel()
a_dens_function, X_fix = sd.start_age_distributions_from_steady_state(
    srm, 
    t0=t0, 
    parameter_dict= par_dict, 
    func_set=func_dict, 
    #x0=smr.start_values#x0=X_0
)
start_mean_age_vec = sd.start_age_moments_from_steady_state(
    srm,
    t0=t0,
    parameter_dict=par_dict,
    func_set=func_dict,
    max_order=1,
    x0=X_fix
).reshape(-1)
from bgc_md2.resolve.mvars import NumericStartValueArray
mvs=mvs.update(
    {
        NumericStartMeanAgeTuple(start_mean_age_vec),
        NumericStartValueArray(X_fix),
    } 
)
start_mean_age_vec

# ### Task: Compute the subsystem backward transit times for the vegetation and soil part of the model

# The list of all mvars showed two interesting variables:
# - NumericVegetationCarbonMeanBackwardTransitTimeSolution
# - NumericVegetationSoilMeanBackwardTransitTimeSolution
# lets find out what we have to provide to be able to compute this.
mvs.jupyter_widget(
    root_type=bgc_md2.resolve.mvars.NumericVegetationCarbonMeanBackwardTransitTimeSolution
)
# Looking at the graph we see that the requirements could be satisfied by providing
# 1.)
#   VegetationCarbonStateVariableTuple 
#   NumericVegetationCarbonStartMeanAgeTuple
#   (or alternatively by the
# 
# 2.)
# VegetationCarbonStateVariableTuple and
# StartConditionMaker (which we will come to later)
# 
# (A similar Graph for the soil part shows symetric requirements)

# +
from bgc_md2.resolve.mvars import VegetationCarbonStateVariableTuple,SoilCarbonStateVariableTuple
# lets add them, uncomment the next lines for help if you want to see the help
# #VegetationCarbonStateVariableTuple?
# #SoilCarbonStateVariableTuple?
mvs = mvs.update(
    {
        VegetationCarbonStateVariableTuple((wood,leaf)),
        SoilCarbonStateVariableTuple((som,slowsom))
    }
)

# +
# Note that this is not a complete division into veg and soil (we left  out the litter pool )
# to look at the veg and soil parts we use some not very well documented tools. 
# (documenatation is improving)
part_dict =  {
    frozenset(mvs.get_VegetationCarbonStateVariableTuple()):'green',
    frozenset(mvs.get_SoilCarbonStateVariableTuple()):'brown',
}
import CompartmentalSystems.helpers_reservoir as hr
hr.igraph_part_plot(
    mvs.get_StateVariableTuple(),
    mvs.get_InFluxesBySymbol(),#in_fluxes,
    mvs.get_InternalFluxesBySymbol(),#internal_fluxes,
    mvs.get_OutFluxesBySymbol(),#out_fluxes,
    part_dict,
)

#mvs.computable_mvar_types()

# +
#mvs.provided_mvar_types,mvs.computable_mvar_types() 
# -

# As the plot shows there are no DIRECT fluxes from Vegetation to Soil (in our devision where litter is not part of the soil)
mvs.get_AggregatedVegetation2SoilCarbonFlux()
# But we can now also compute a lot of interesting numeric results

# +
# The consistency of the subsystem start mean ages with the system start mean ages also has to be guaranteed
# by computing them according to the same pre start past. 
soil_smr=mvs.get_SoilCarbonSmoothModelRun()
soil_start_mean_age_vec = sd.start_age_moments_from_steady_state(
    soil_smr.model,
    t0=np.min(times),
    parameter_dict=soil_smr.parameter_dict,
    func_set=soil_smr.func_set,
    max_order=1,
).reshape(-1)
veg_smr=mvs.get_VegetationCarbonSmoothModelRun()

veg_start_mean_age_vec = sd.start_age_moments_from_steady_state(
    veg_smr.model,
    t0=np.min(times),
    parameter_dict=veg_smr.parameter_dict,
    func_set=veg_smr.func_set,
    max_order=1,
).reshape(-1)
from bgc_md2.resolve.mvars import NumericSoilCarbonStartMeanAgeTuple, NumericVegetationCarbonStartMeanAgeTuple
mvs = mvs.update(
    {
        NumericSoilCarbonStartMeanAgeTuple(soil_start_mean_age_vec),
        NumericVegetationCarbonStartMeanAgeTuple(veg_start_mean_age_vec),
    }
)
#soil_start_mean_age_vec,veg_start_mean_age_vec


# +
# we could now compute the desired properties
soil_mean_btts2 = mvs.get_NumericSoilCarbonMeanBackwardTransitTimeSolution()
mean_btts2 = mvs.get_NumericMeanBackwardTransitTimeSolution()
veg_mean_btts2 = mvs.get_NumericVegetationCarbonMeanBackwardTransitTimeSolution()

fig2 = plt.figure(figsize=(10, 10))
axs2 = fig2.subplots(1, 1)
ax = axs2
ax.plot(times, mean_btts2, color="blue", label="mean btts2")
ax.plot(times, veg_mean_btts2, color="green", label="veg_mean btts2")
ax.plot(times, soil_mean_btts2, color="brown", label="soil_mean btts2")
ax.legend()

#fig2.savefig(Path(mf).joinpath("btts.pdf"))

# +
# there is however some duplication in the above code. Actually the way how the start_values for the subsystems
# are computed is similar and we could factor this out into a function that computes the start values
# consistently for the subsystems.
# The framework allows to formulate such a function and will apply it automatically to the main system
# and the subsystems. Since we know that our system is linear we can formulate a function that 
# computes an equilibrium and derives the equilibrium pool contents, mean ages and the age_density function
# we could also have formulated another function that consistenly initializes empty pools or a spinup
# generated values
def scm(
        npsrm #: NumericParameterizedSmoothReservoirModel or a subclass 
    ):
    srm=npsrm.srm # the Smoth reservoir model (or a subclass)
    nupa=npsrm.parameterization
    par_dict=nupa.par_dict
    func_dict=nupa.func_dict
    a_dens_function, X_fix = sd.start_age_distributions_from_steady_state(
        srm, 
        t0=t0, 
        parameter_dict= par_dict, 
        func_set=func_dict, 
        #x0=smr.start_values#x0=X_0
    )
    start_mean_age_vec = sd.start_age_moments_from_steady_state(
        srm,
        t0=t0,
        parameter_dict=par_dict,
        func_set=func_dict,
        max_order=1,
        x0=X_fix
    ).reshape(-1)
    return X_fix,start_mean_age_vec, a_dens_function

# this function could actually replace several of the previously provided values
# we will demonstrate this by first removing the concrete values from our mvs object:
mvs=mvs.remove(
    [
        NumericVegetationCarbonStartMeanAgeTuple,
        NumericSoilCarbonStartMeanAgeTuple,
        NumericStartMeanAgeTuple,
        NumericStartValueDict
    ]
)    
# which removes our ability to compute the backward transittimes (and even the solution)
mvs.provided_mvar_types, mvs.computable_mvar_types()
# -

from bgc_md2.resolve.mvars import StartConditionMaker
mvs=mvs.update(
    {
        StartConditionMaker(scm)
    }
)

npsrm=mvs.get_NumericParameterizedSmoothReservoirModel()

# +
veg_npsrm=mvs.get_NumericParameterizedVegetationCarbonSmoothReservoirModel()

soir_smr=mvs.get_SoilCarbonSmoothModelRun()
soil_npsrm=mvs.get_NumericParameterizedSoilCarbonSmoothReservoirModel()
func=mvs.get_StartConditionMaker()
# -

soil_smr.func_set.keys()

func(soil_npsrm)
#soil_npsrm.parameterization.func_dict.keys()

# +
soil_mean_btts2 = mvs.get_NumericSoilCarbonMeanBackwardTransitTimeSolution()
mean_btts2 = mvs.get_NumericMeanBackwardTransitTimeSolution()
veg_mean_btts2 = mvs.get_NumericVegetationCarbonMeanBackwardTransitTimeSolution()

fig2 = plt.figure(figsize=(10, 10))
axs2 = fig2.subplots(1, 1)
ax = axs2
ax.plot(times, mean_btts2, color="blue", label="mean btts2")
ax.plot(times, veg_mean_btts2, color="green", label="veg_mean btts2")
ax.plot(times, soil_mean_btts2, color="brown", label="soil_mean btts2")
ax.legend()

#fig2.savefig(Path(mf).joinpath("btts.pdf"))
# -

mean_btts2 = mvs.get_NumericMeanBackwardTransitTimeSolution()
veg_mean_btts2 = mvs.get_NumericVegetationCarbonMeanBackwardTransitTimeSolution()

# +
fig2 = plt.figure(figsize=(10, 10))
axs2 = fig2.subplots(1, 1)
ax = axs2
ax.plot(times, mean_btts2, color="blue", label="mean btts2")
ax.plot(times, veg_mean_btts2, color="green", label="veg_mean btts2")
ax.plot(times, soil_mean_btts2, color="brown", label="soil_mean btts2")
ax.legend()

#fig2.savefig(Path(mf).joinpath("btts.pdf"))
# -


# ### Task: Change the model parameters so that ALL material leaving has to pass through the soil sub system first. 
# - The easiest way to achieve this is to set all output rates for the veg pools to zero (parameters and one function). 
# - Presently the mean backward transit time of the soil subsystem is longer than the mean for the whole system. Is this possible under the new circumstances (of all input coming throug the veg subsystem and all output leaving from the soil subsystem?
# - Hint: Go back and change the parameters, ( to avoid caching issues with the jupyter notebook restart the kernel and rerun. You should see qualitative changes in the plot above.
#

# ### Task: Compute the mean ages for the pools (for the whole system and the sub systems.) 
# - The transit times is the age of the material leaving (the time since it entered the SYSTEM or SUBSYSTEM), so it depends on this age and the outflux. The age of material in pools with large outfluxes will dominate. 
# - In the next plot we look at the system/subsystem ages. Note that the system ages for the soil subsystem are smaller that the system ages for the whole system, because the soil system clock starts ticking when the soil subsystem is entered (after traversing the veg pools...)
# - The veg pools have the same age with respect to the whole system and the subsystem (no system pools to pass before entering the veg subsystem

# +
m_a_arr2 = mvs.get_NumericMeanAgeSolutionArray()
veg_m_a_arr2 = mvs.get_NumericVegetationCarbonMeanAgeSolutionArray()
soil_m_a_arr2 = mvs.get_NumericVegetationCarbonMeanAgeSolutionArray()
vcsv=mvs.get_VegetationCarbonStateVariableTuple()
scsv=mvs.get_SoilCarbonStateVariableTuple()

fig1 = plt.figure(figsize=(2 * 10, n_pools * 10))
axs = fig1.subplots(n_pools)
#start_ind=200
start_ind=0
n_steps=len(times)
for i in range(n_pools):
    
    ax = axs[i]
    #ax.plot(ad_times[start_ind:n_steps], m_a_arr[start_ind:n_steps, i]/dpy, label="1")
    ax.plot(times[start_ind:n_steps], m_a_arr2[start_ind:n_steps, i], label="sys")
    if sv[i] in vcsv:
        j=list(vcsv).index(sv[i])
        ax.plot(times[start_ind:n_steps], veg_m_a_arr2[start_ind:n_steps,j], "+", color="green",label="veg")
    if sv[i] in scsv:                                                         
        j=list(scsv).index(sv[i])                                             
        ax.plot(times[start_ind:n_steps], soil_m_a_arr2[start_ind:n_steps,j], color="brown", label="soil")
    ax.legend()
    ax.set_title(f"{sv[i]} mean_age")
# -

mvs.get_NumericStartMeanAgeTuple()

# ### Task: Compute the age and transit time densities 
# This task is not fully integrated into the framework in the sense that there are functions that compute the desired result directly, so the graph library cannot show us the way.
# This is the typical situation for the development of the package:
# 1. Somebody asks for a new feature, say the backward transsit time density of the vegetation subsystem 
# 1. We implement it, using and potetially extending `CompartmentalSystems` and `LAPM` (using minimal ingredients that the framework can already provide. (The less we need the easier the new result will be available for a potential user later, in other words the less information they will have to provide.)
# 1. We check very carefully if the result depends on any special assumptions. If it does, we 
#    resist the temptation to make it available via the framework, which will be percieved
#    by many users as a black box and should therefore only claim results it can ALWAYS  deliver.
# 1. We give the new result a unique name and create a type (i.e `NewResult`) for it in `bgc_md2.resolve.mvars` (usually a pretty empty python class derived from something similar, but possibly a sanitizing `__init__` method to catch some user error early.)
# 1. We transform the implementation into a function in `bgc_md2.resolve.computers` with a type signature reflecting the minimal ingredients and make it return the new type.
# 1. Hope for the best add the ingredients to a new CMTVS object (let's call it  mvs as in this notebook) and check mvs.computable_mvars(). If it doesn't show up, check the signature and types of the ingredients until it does. Then call `mvs.get_NewResult()` which the framework will have automatically added as a method to `mvs`.  
# 1. Read and address error messages...
#
# At the moment this tutorial shows only the first two points...

# construct a function p that takes an age array "ages" as argument
# and gives back a three-dimensional ndarray (ages x times x pools)
# from the a array-valued function of a single age a_dens_function
#srm = mvs.get_SmoothReservoirModel()
veg_smr = mvs.get_VegetationCarbonSmoothModelRun()
func= mvs.get_StartConditionMaker()
veg_npsrm = mvs.get_NumericParameterizedVegetationCarbonSmoothReservoirModel()
X_fix,veg_start_mean_age_vec,veg_a_dens_function=func(veg_npsrm)

veg_p = veg_smr.pool_age_densities_func(veg_a_dens_function)

veg_ages = np.linspace(
    0,
    (np.array(veg_start_mean_age_vec, dtype=float).reshape(-1)).max() * 50,
    1000
)  
veg_age_densities = veg_p(veg_ages)

from pathlib import Path
from plotly.offline import plot,iplot
for n in range(veg_smr.nr_pools):
    max_ind = np.argmin(veg_ages < start_mean_age_vec[n] * 2)
    fig = veg_smr.plot_3d_density_plotly(
        "age distribution pool {0}".format(sv[n]),
        veg_age_densities[0:max_ind, :, n],
        veg_ages[0:max_ind],
    )
    # plot the computed start age density for t0 on top
    fig.add_scatter3d(
        x=np.array([-t0 for a in veg_ages]),
        y=np.array([a for a in veg_ages]),
        z=np.array([a_dens_function(a)[n] for a in veg_ages]),
        mode="lines",
        line=dict(color="#FF0000", width=15),
    )
    veg_smr.add_line_to_density_plot_plotly(
        fig,
        data=veg_m_a_arr2[:, n],
        color="#FF0000",
        name="mean age",
        time_stride=1,
        on_surface=True,
        bottom=True,
        legend_on_surface=True,
        legend_bottom=False,
    )
    iplot(fig)
    #plot(
    #    fig,
    #    filename=str(
    #        Path(mf).joinpath(f"age_distribution_{sv[n]}.html")
    #    ),
    #    auto_open=False,
    #)

veg_btt_dens = veg_smr.backward_transit_time_density(veg_age_densities)
fig_btt = veg_smr.plot_3d_density_plotly(
    "backward_transit_time_density_steady_state",
    veg_btt_dens,
    ages,
    y_label="transit time",
)
veg_smr.add_line_to_density_plot_plotly(
    fig_btt,
    data=veg_mean_btts2,
    color="#FF0000",
    name="mean age",
    time_stride=1,
    on_surface=True,
    bottom=True,
    legend_on_surface=True,
    legend_bottom=False,
)
iplot(
    fig_btt
)
#plot(
#    fig_btt,
#    filename=str(
#        Path(mf).joinpath("btt_distribution.html")
#    ),
#    auto_open=False,
#)


