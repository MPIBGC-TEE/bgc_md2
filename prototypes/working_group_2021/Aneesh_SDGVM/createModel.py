# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

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
# We will start with an extremely simple model, that you can find in 
# ```bash
# bgc_md2/src/bgc_md2/models/testVectorFree/source.py
# ```
# Copy its contents into a new cell. (we will later save it under a new name) 

# +
from sympy import var 
from ComputabilityGraphs.CMTVS import CMTVS
#from bgc_md2.helper import bgc_md2_computers
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

var("I_vl I_vw vl vw k_vl k_vw")

mvs = CMTVS(
    {
        StateVariableTuple((vl, vw)),
        TimeSymbol("t"),
        InFluxesBySymbol({vl: I_vl, vw: I_vw}),
        OutFluxesBySymbol({vl: k_vl * vl, vw: k_vw * vw}),
        InternalFluxesBySymbol({(vl, vw): k_vl * vl, (vw, vl): k_vw * vw}),
    },

    computers=module_computers(bgc_c)
)
# -

# The last statement in the code defines a variable `mvs` which is 
# an instance of CMTVS which stands for `C`onnected`M`ulti`T`ype`V`ariable`S`et".
# It contains information in two forms. 
# 1. Variables of certain type (like InFluxesBySymbol)
# 2. a set of functions that "connect" these Variables and to other results we did not specify 
#    explicitly.
#
# Taks:
# To see what it can do with this information add a new cell an type `mvs.` and press the `tab` key `->`. This will show you the available methods, in other words what can be computed from the provided information.
#
# I wanted to see the compartmental the pools, the matrix and the inputs.

mvs.get_StateVariableTuple()

mvs.get_CompartmentalMatrix()

mvs.get_InputTuple()

# we can also print the whole mass balance equation
import bgc_md2.display_helpers as dh
dh.mass_balance_equation(mvs)

# we can also plot a picture
import bgc_md2.helper as h
h.compartmental_graph(mvs)

# +
# adjust the output to full width
from IPython.display import HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

# make changes to imported files immidiately available 
# avoiding the need to reload (in most cases)
# %load_ext autoreload
# %autoreload 2
# -

# ## Extend the minimal model
#
# ### add more state variables and add a small description
#
# In a first step I will 
# - add new symbols to desricbe pools for the cable model (from Alisons code) but
#   leave the old ones there since they are presently used in the fluxes.
#   We will replace them one bye one later 
#   but this way we are not forced to replace all the fluxes at
#   once.It's always good to be able to make small steps...
# - use a different way to declare symbols 
#

# +
from sympy import var, exp
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

#var("""
#C_leaf C_wood C_root C_abvstrlit C_abvmetlit C_belowstrlit C_belowmetlit C_surface_microbe C_soil_microbe _slowsom C_passsom
#k_C_leaf k_C_wood k_C_root k_C_abvmetlit k_C_abvstrlit k_C_belowstrlit k_C_belowmetlit 
#k_C_surface_microbe k_C_soil_microbe k_C_slowsom k_C_passsom 
#beta1 beta2 beta3
#leached ls_aboveground 
#ls_belowground f41 f51 f52 f63 f73 f84 f104 f85 f96 f106 f97 f910 f911 f108 f119 f109 f1110 
#Ca clay silt_clay leachedwtr30""")
# we organize the new symbols a bit better by putting them in a dictionary
# together with some description that we can use later to display some metainformation
sym_dict = {
        'C_leaf': 'content of leaf pool',
        'C_wood': 'content of wood pool',
        'C_root': 'content of root pool',
        'C_abvstrlit': 'content of aboveground structural litter pool',
        'C_abvmetlit': 'content of aboveground metabolic litter pool',
        'C_belowstrlit': 'content of belowground structural litter pool',
        'C_belowmetlit': 'content of belowground metabolic litter pool',
        'C_surface_microbe': 'content of surface microbe pool',
        'C_soil_microbe': 'content of soil microbial pool',
        'C_slowsom': 'content of the slow soil pool',
        'C_passsom': 'content of the passive soil pool',
        'beta1': 'NPP partitioning to leaf', 
        'beta2': 'NPP partitioning to wood',
        'beta3': 'NPP partitioning to root',
        'k_C_leaf': 'Turn over rate of leaf pool',
        'k_C_wood':'Turn over rate of wood pool',
        'k_C_root': 'Turn over rate of root pool',
        'k_C_abvmetlit':'Turn over rate of aboveground metabolic litter pool',
        'k_C_abvstrlit' :'Turn over rate of aboveground structural litter pool',
        'k_C_belowstrlit' :'Turn over rate of belowground structural litter pool',
        'k_C_belowmetlit':'Turn over rate of belowground metabolic litter pool',
        'k_C_surface_microbe' :'Turn over rate of surface microbe pool',
        'k_C_soil_microbe' :'Turn over rate of soil microbial pool',
        'k_C_slowsom' :'Turn over rate of slow soil pool',
        'k_C_passsom':'Turn over rate of passive soil pool',
        'ls_aboveground': 'Fraction of structural C that is lignin in aboveground',
        'ls_belowground': 'Fraction of structural C that is lignin in belowground',
        'f_leaf2abvstrlit': 'Transfer coefficient of C from leaf to aboveground structural litter',
        'f_leaf2abvmetlit':'Transfer coefficient of C from leaf to aboveground metabolic litter',
        'f_wood2abvstrlit':'Transfer coefficient of C from wood to aboveground structural litter',
        'f_wood2abvmetlit':'Transfer coefficient of C from wood to aboveground metabolic litter',
        'f_root2belowstrlit':'Transfer coefficient of C from root to belowground structural litter',
        'f_root2belowmetlit':'Transfer coefficient of C from root to belowground metabolic litter',
        'f_abvstrlit2surface_microbe':'Transfer coefficient of C from aboveground structural litter to surface microbe',
        'f_abvstrlit2slowsom':'Transfer coefficient of C from aboveground structural litter to slow soil',
        'f_abvmetlit2surface_microbe':'Transfer coefficient of C from aboveground metabolic litter to surface microbe',
        'f_belowstrlit2soil_microbe':'Transfer coefficient of C from belowground structural litter to soil microbe',
        'f_belowstrlit2slowsom':'Transfer coefficient of C from belowground structural litter to slow soil',
        'f_belowmetlit2soil_microbe':'Transfer coefficient of C from belowground metabolic litter to soil microbe',
        'f_slowsom2soil_microbe':'Transfer coefficient of C from soil soil to soil microbe',
        'f_passsom2soil_microbe':'Transfer coefficient of C from passive soil to soil microbe',
        'f_surface_microbe2slowsom':'Transfer coefficient of C from surface microbe to slow soil',
        'f_soil_microbe2passsom':'Transfer coefficient of C from soil microbe to passive soil',
        'f_soil_microbe2slowsom':'Transfer coefficient of C from soil microbe to slow soil',
        'f_slowsom2passsom':'Transfer coefficient of C from slow soil to passive soil',
        'Ca': 'Need to find (unknown)',
        'clay': 'Clay content (fraction)',
        'silt_clay': 'Silt plus clay content (fraction)',
        'leachedwtr30': 'H20 leached below 30cm (cm/m)'
    
}
# for the moment only use the 
var(list(sym_dict.keys()))

ls_aboveground = 0.12
ls_belowground = 0.35

f_leaf2abvstrlit = 0.91986 + 0.00324* Ca
f_leaf2abvmetlit = 1 - f_leaf2abvstrlit
f_wood2abvstrlit = f_leaf2abvstrlit
f_wood2abvmetlit = 1 - f_wood2abvstrlit
f_root2belowstrlit = f_leaf2abvstrlit
f_root2belowmetlit = 1 - f_root2belowstrlit
f_abvstrlit2surface_microbe = (1 - exp(-3*ls_aboveground))*0.4
f_abvstrlit2slowsom = exp(-3*ls_aboveground) * 0.7
f_abvmetlit2surface_microbe = 0.4
f_belowstrlit2soil_microbe = exp(1 - (-3*ls_belowground))*0.45
f_belowstrlit2slowsom = exp(-3*ls_aboveground) * 0.7
f_belowmetlit2soil_microbe = 0.45
f_slowsom2soil_microbe = 0.447 + 0.009*clay
f_passsom2soil_microbe = 0.45
f_surface_microbe2slowsom = 0.4
f_soil_microbe2passsom = 0.003 + 0.032 * clay
f_soil_microbe2slowsom = 1 - f_soil_microbe2passsom - (leachedwtr30)/18 * (0.01 + 0.04* (1- silt_clay)) - 0.85 - 0.68*silt_clay
f_slowsom2passsom = 0.003 + 0.009*clay


mvs = CMTVS(
    {
        StateVariableTuple((C_leaf, C_wood, C_root, C_abvstrlit, C_abvmetlit, C_belowstrlit, C_belowmetlit, C_surface_microbe, C_soil_microbe, C_slowsom, C_passsom)),
        TimeSymbol("t"),
        InFluxesBySymbol({C_leaf: beta1, C_wood: beta2, C_root: beta3}),
        OutFluxesBySymbol({C_soil_microbe: leached}),
        InternalFluxesBySymbol({(C_leaf, C_abvstrlit): k_C_leaf *f_leaf2abvstrlit* C_leaf, 
                                (C_leaf, C_abvmetlit): k_C_leaf *f_leaf2abvmetlit* C_leaf, 
                                (C_wood, C_abvstrlit): k_C_wood *f_wood2abvstrlit* C_wood, 
                                (C_wood, C_abvmetlit): k_C_wood *f_wood2abvmetlit* C_wood, 
                                (C_root, C_belowstrlit): k_C_root * f_root2belowstrlit* C_root, 
                                (C_root, C_belowmetlit): k_C_root * f_root2belowmetlit * C_root, 
                                (C_abvstrlit , C_surface_microbe ): k_C_abvstrlit *f_abvstrlit2surface_microbe* C_abvstrlit, 
                                (C_abvstrlit , C_slowsom ): k_C_abvstrlit*f_abvstrlit2slowsom*C_abvstrlit,
                                (C_abvmetlit, C_surface_microbe ): k_C_abvmetlit *f_abvmetlit2surface_microbe* C_abvmetlit, 
                                (C_belowstrlit, C_soil_microbe): k_C_belowstrlit * f_belowstrlit2soil_microbe*C_belowstrlit, 
                                (C_belowmetlit , C_soil_microbe): k_C_belowmetlit  * f_belowmetlit2soil_microbe*C_belowmetlit, 
                                (C_belowstrlit, C_slowsom): k_C_belowstrlit *f_belowstrlit2slowsom* C_belowstrlit, 
                                (C_surface_microbe , C_slowsom): k_C_surface_microbe*f_surface_microbe2slowsom*C_surface_microbe, 
                                (C_soil_microbe, C_slowsom): k_C_soil_microbe *f_soil_microbe2slowsom*C_soil_microbe, 
                                (C_slowsom, C_soil_microbe): k_C_slowsom *f_slowsom2soil_microbe*C_slowsom, 
                                (C_soil_microbe, C_passsom): k_C_soil_microbe*f_soil_microbe2passsom*C_soil_microbe, 
                                (C_slowsom , C_passsom): k_C_slowsom*f_slowsom2passsom*C_slowsom,
                               (C_passsom, C_soil_microbe): k_C_passsom * f_passsom2soil_microbe * C_passsom})
    },

    computers=module_computers(bgc_c)
)
# -

# Nothing has changed in the model description but be have some more symbols to work with.
# We could type `C_leaf` somewhere without an error since it is now known as a variable.
# In the next step we replace all occurences of `vl` by `C_leaf` 

# we can also plot a picture
import bgc_md2.helper as h
h.compartmental_graph(mvs)

mvs.get_StateVariableTuple()

mvs.get_CompartmentalMatrix().free_symbols


mvs.get_CompartmentalMatrix()


mvs.get_InputTuple()

# we can also plot a picture
import bgc_md2.helper as h
h.compartmental_graph(mvs)

import sys
import numpy as np
sys.path.insert(0,'..') # necessary to import general_helpers
from general_helpers import make_B_u_funcs
from model_specific_helpers import UnEstimatedParameters,EstimatedParameters,ModelParameters,StateVariables
def construct_V0(
        cpa :UnEstimatedParameters,
        epa :EstimatedParameters,
    ) -> np.ndarray:
    """Construct the initial values for the forward simulation
    from constant and evaluated parameters
    param: cpa : constant parameeters
    param: epa : estimated parameters 
    """
    # to make sure that we are in the right order we use the 
    # StateVariables namedtuple 
    X_0 = StateVariables( 
        C_leaf = epa.C_leaf_0,
        C_wood = cpa.C_Veg_0 - epa.C_leaf_0 - cpa.C_root_0,
        C_root=cpa.C_root_0,
        C_abvstructural_lit = epa.C_abvstrlit_0_frac*cpa.C_litter_0,
        C_abvmetabolic_lit = epa.C_abvmetlit_0_frac*cpa.C_litter_0,
        C_belowstructual__lit = epa.C_blwstrlit_0_frac*cpa.C_litter_0,
        C_belowmetabolic_lit = cpa.C_litter_0 - epa.C_abvstrlit_0_frac*cpa.C_litter_0 - epa.C_abvmetlit_0_frac*cpa.C_litter_0 - epa.C_blwstrlit_0_frac*cpa.C_litter_0,
        C_surface_microbe = epa.C_surfacemic_0_frac*cpa.C_soil_0,
        C_soil_microbe = epa.C_soilmic_0_frac*cpa.C_soil_0,
        C_slow_soil = epa.C_slow_0_frac*cpa.C_soil_0,
        C_passive_soil = cpa.C_soil_0 - epa.C_soilmic_0_frac*cpa.C_soil_0 - epa.C_slow_0_frac*cpa.C_soil_0 - epa.C_surfacemic_0_frac*cpa.C_soil_0
    )
    # add the respiration start value to the tuple
    cpa = UnEstimatedParameters(
    C_Veg_0 = cveg[0],
    C_root_0 = croot[0],
    C_litter_0 = clitter[0],
    npp = npp,
    number_of_months = tot_len,
    C_soil_0 = csoil[0],
    rh_0=rh[0]
    )
    V_0 = (*X_0,cpa.rh_0)
    return np.array(V_0).reshape(11,1)  



# +
from pathlib import Path
import json
from model_specific_helpers import get_example_site_vars

with Path('config.json').open(mode='r') as f:
    conf_dict = json.load(f)

dataPath = Path(conf_dict['dataPath'])

# fixme:
#    Note that the function is imported from
#    model_specific_helpers which means that you have to provide
#    your version of this fuction which will most likely return different
#    variables
npp, rh, croot, csoil, cveg, clitter = get_example_site_vars(dataPath)

# will become a test later
epa_0 = EstimatedParameters(
    beta_leaf   =0.15,                                      #0
    beta_wood   =0.6,                                       #1
    Ca          =12,                                          #2
    clay        =0.5,                                         #3
    leachedwtr30  = 50,                                       #4
    silt_clay   =0.5,                                          #5
    k_leaf      =   (1/(2*365) + 1/(0.3*365))/2,               #6
    k_wood      =(1/(60*365) + 1/365)/2,                         #7
    k_root      =(1/(10*365) + 1/(365*0.8))/2,                    #8
    C_leaf_0      =((cveg[0]-croot[0])*0.1 + (cveg[0]-croot[0])*0.9)/2,      #9
    C_abvstrlit_0_frac    = 0.3,                                  #10
    C_abvmetlit_0_frac    = 0.3,                                   #11
    C_blwstrlit_0_frac    = 0.3,                                       #12
    C_surfacemic_0_frac   = 0.3,                                   #13
    C_soilmic_0_frac      = 0.3,                               #14
    C_slow_0_frac         = 0.3                                      #15
)

cpa = UnEstimatedParameters(
    C_Veg_0 = cveg[0],
    C_root_0 = croot[0],
    C_litter_0 = clitter[0],
    npp = npp,
    number_of_months = 120,
    C_soil_0 = csoil[0],
    rh_0=rh[0]
)

# -

StateVariables._fields

construct_V0(cpa,epa_0)


# +
    
def make_daily_iterator_sym(
        mvs,
        V_init,
        mpa,
        func_dict
    ):
        
        #from bgc_md2.models.cable_yuanyuan.source import mvs 
        B_func, u_func = make_B_u_funcs(mvs,mpa,func_dict)  
        
        def f(it,V):
            X = V[0:9]
            co2 = V[9]
            b = u_func(it,X)
            B = B_func(it,X)
            X_new = X + b + B@X

            # we also compute the respired co2 in every (daily) timestep
            # and use this part of the solution later to sum up the monthly amount
            co2_new = -np.sum(B @ X) # fixme add computer for respirattion
            
            V_new = np.concatenate((X_new,np.array([co2_new]).reshape(1,1)), axis=0)
            return V_new
    
        return TimeStepIterator2(
                initial_values=V_init,
                f=f#,
                #max_it=max(day_indices)+1
        )


# -

it_sym = make_daily_iterator_sym(
            mvs,
            V_init=construct_V0(
                self.cpa,
                self.epa0
            ),
            mpa=mpa,
            func_dict={
                #'NPP':npp_func
            },
        )
 #       n=5
  #      res= np.zeros((n,10))
   #     res_sym = copy(res)
    #    for i in range(n):
     #       res[i,:]=it.__next__().reshape(10,)
      #      res_sym[i,:]=it_sym.__next__().reshape(10,)
      #  self.assertTrue(
       #     np.allclose(
        #        res,
         #       res_sym
          #  )
        #)


par_dict={
 C_soil_microbe : 3.4,
 Ca = 12,
 clay=0.5,
 k_C_abvmetlit,
 k_C_abvstrlit,
 k_C_belowmetlit,
 k_C_belowstrlit,
 k_C_leaf,
 k_C_root,
 k_C_slowsom,
 k_C_soil_microbe,
 k_C_surface_microbe,
 k_C_wood,
 leached,
 leachedwtr30,
 silt_clay}
}

srm=mvs.get_SmoothReservoirModel()

from sympy import Symbol,var
var("a b c")
h=a**2+b**2+c**2
h

par_dict={a:2,b:3,c:5}

h.subs(par_dict)


