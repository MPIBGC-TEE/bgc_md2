# +
# To do:
# - check corner cases for empty Influx dict in the plot
# - check the corner case for no flux from a node in the matplotlib_plot
# +
# adjust the output to full width
#from IPython.display import HTML
#display(HTML("<style>.container { width:100% !important; }</style>"))

# make changes to imported files immidiately available 
# avoiding the need to reload (in most cases)
# #%load_ext autoreload
# #%autoreload 2

# +
from pathlib import Path
from sympy import var, Symbol, Function, Tuple, lambdify, exp, Abs, Matrix, diff
#from ComputabilityGraphs.CMTVS import CMTVS
#from ComputabilityGraphs.helpers import module_computers
##from bgc_md2.helper import module_computers
from bgc_md2.resolve.mvars import (
#    InFluxesBySymbol,
#    OutFluxesBySymbol,
#    InternalFluxesBySymbol,
    TimeSymbol,
#    StateVariableTuple,
)
#from importlib import import_module
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
# the substrates
# in our example only C (OC_1 =OC_C)
i_s=[1,2]

# assume different electron acceptors e.g. oxygen O
js=[1]

# assume different functional groups of micro organisms
# They each have given stochiometric ratios between the elements
# [C    H   N    O    P    S]
# [a    b   c    d    e    f] with charge z
# In our examples from the paper we have  
# [a=1, b=1.8, c=0.2,  d=0.5, e=0 f=0, z=0
ks=[1]

# for every pair of Substrate and Electron acceptor we have one reaction that involves different chemical species
# e.g. for 
# 1*C  +2*O-  >1* CO_2 we would have 3 species 
# sp1=C,  sp2=O,  sp3=CO_2  with coefficient 
# y1=-1, y2=-2, y3=1 
# We group them into coefficant for reactants (negative) and products (positive) 
# since this decides if they are related to influxes or outfluxes of the pools 
sym_dict={
    **{
        f"OC_{i}": f"mass of substrate {i}"
        for i in i_s
    }
    ,
    **{
        f"EA_{j}": f"mass of electron acceptor {j}"
        for j in js
    }
    ,
    **{
        f"B_{k}": f"mass of functional group {k}"
        for k in ks
    }
    ,
    **{
        f"HCO3": f"mass of in HCO3"
    }
    ,
    **{
        f"mu_{i}_{j}_{k}":
        f"growth rate of functional group B{k} due to reaction combining  OC_{i} and EA_{j}"
        for i in i_s
        for j in js
        for k in ks
    }
    ,
    **{
        f"y_OC_{i}_{j}":
        f"stocheometric vector component for element for reaction combining OC_{i} and EA_{j}"
        for i in i_s
        for j in js
    }
    ,
    **{
        f"y_EA_{i}_{j}":
        f"stocheometric vector for element for reaction combining {i} and {j}"
        for i in i_s
        for j in js
    }
    ,
    **{
        f"y_HCO3_{i}_{j}":
        f"stocheometric vector for element rjction combining {i} and {j}"
        for i in i_s
        for j in js
    }

}
# Make symbols from  the strings that we can later use in expressions  
# B_1, B_2,...
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)
# Note:
# The use of exec to execute python code strings seems a bit funny 
# The normal way woud be: 
# B_1 = Symbol("B_1")
# B_2 = Symbol("B_2")
# ...
# we do it this way, because we want the dictionary for documentation later...

## We will also use some symbolic functions ("symbols" with an argument) 
#sym_func_dict={
#}
#for k in sym_func_dict.keys():
#    code=k+" = Function('{0}')".format(k)
#    exec(code)
## Note:
## Again we want the dictionary for later

t=TimeSymbol("t")
state_var_tuple=Tuple(
    *[ 
        Symbol(f"OC_{i}")
        for i in i_s
    ]
    +
    [
        Symbol(f"EA_{j}")
        for j in js
    ] 
    +
    [ 
        Symbol("HCO3")
    ]
    +
    [
        Symbol(f"B_{k}")
        for k in ks
    ]
)    
#rhs_sym=Tuple((Symbol(
#state_var_tuple=Tuple(x,y)
#rhs_sym=Tuple(x,p*y)
#y0=(1,1)
rhs_sym=Tuple(
    *(
        [
            sum(
                [
                    (
                        Symbol(f"y_OC_{i}_{j}")
                        *Symbol(f"mu_{i}_{j}_{k}")
                        * Symbol(f"B_{k}")
                    )
                    for j in js
                    for k in ks
                ]
            )
            for i in i_s
        ]
        +
        [
            sum(
                [
                    (
                        Symbol(f"y_EA_{i}_{j}")
                        * Symbol(f"mu_{i}_{j}_{k}")
                        * Symbol(f"B_{k}")
                    )
                    for i in i_s
                    for k in ks
                ]
            )
            for j in js
        ]
        +
        [
            sum(
                [
                    (
                        Symbol(f"y_HCO3_{i}_{j}")
                        *Symbol(f"mu_{i}_{j}_{k}")
                        * Symbol(f"B_{k}")
                    )
                    for j in js
                    for i in i_s
                    for k in ks
                ]
            )
        ]
        +
        [
            sum(
                [
                    (
                        Symbol(f"mu_{i}_{j}_{k}")
                        * Symbol(f"B_{k}")
                    )
                    for j in js
                    for i in i_s
                ]
            )
            for k in ks
        ]
    )
)

start_value_dict={
    **{
        Symbol(f"OC_{i}"): 1.0
        for i in i_s
    }
    ,
    **{
        Symbol(f"EA_{j}"): 2.0
        for j in js
    }
    ,
    Symbol("HCO3"): 4
    ,
    **{
        Symbol(f"B_{k}"): 3.0
        for k in ks
    }
}
y0=state_var_tuple.subs(start_value_dict)
par_dict = {
    **{
        Symbol(f"mu_{i}_{j}_{k}"): 0.1 
        for i in i_s
        for j in js
        for k in ks
    }
    ,
    **{
        Symbol(f"y_OC_{i}_{j}"): 0.1
        for i in i_s
        for j in js
    }
    ,
    **{
        Symbol(f"y_EA_{i}_{j}"): 0.1
        for i in i_s
        for j in js
    }
    ,
    **{
        Symbol(f"y_HCO3_{i}_{j}"): 0.1
        for i in i_s
        for j in js
    }
}        

solve_ivp(
    fun=lambdify((t,state_var_tuple),rhs_sym.subs(par_dict),"numpy"),
    t_span=(0,5),
    y0=y0,
)
# until now we treated the mu's as parameters 
# now we use eq 6 and 7 to replace them
subs_dict_2 = {
    **{
        Symbol(f"mu_{i}_{j}_{k}"): 
        (
            Symbol(f"mu_max_{i}_{j}_{k}")
            *
            exp(
                -Abs(Symbol(f"y_OC_{i}_{j}"))
                /
                (
                    Symbol(f"V_h")
                    *Symbol(f"OC_{i}")
                )
            )
            *
            exp(
                -Abs(Symbol(f"y_EA_{i}_{j}"))
                /
                (
                    Symbol(f"V_h")
                    *Symbol(f"EA_{j}")
                )
            )
        )
        for i in i_s
        for j in js
        for k in ks
    }
}
rhs_sym_2=rhs_sym.subs(subs_dict_2)
# instead of the mus we now have mu_max_{i}_{j}_{k} and Vh
par_dict_2 = {
    Symbol("V_h"):5 # fake value
    ,
    **{
        Symbol(f"mu_max_{i}_{j}_{k}"): 0.1 # fake value
        for i in i_s
        for j in js
        for k in ks
    }
    ,
    **{
        Symbol(f"y_OC_{i}_{j}"): -0.1 # fake value, negative because it is a reactant
        for i in i_s
        for j in js
    }
    ,
    **{
        Symbol(f"y_EA_{i}_{j}"): -0.1 # fake value, negative because it is a reactant
        for i in i_s
        for j in js
    }
    ,
    **{
        Symbol(f"y_HCO3_{i}_{j}"): -0.1 ## fake value
        for i in i_s
        for j in js
    }
}        
sol2=solve_ivp(
    fun=lambdify((t,state_var_tuple),rhs_sym_2.subs(par_dict_2),"numpy"),
    t_span=(0,50),
    y0=y0,
    dense_output=True
).sol
# -

plot_times=np.linspace(sol2.t_min,sol2.t_max,100)
fig=plt.figure()
ax=fig.add_subplot()
values=sol2(plot_times)
for i in range(values.shape[0]):
    ax.plot(
        plot_times,
        values[i,:],
        label=str(state_var_tuple[i])
    ) 
ax.legend()
fig.savefig(f"{Path(__file__).stem}.pdf")




