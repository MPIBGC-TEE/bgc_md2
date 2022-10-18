# # Problems with $\xi$
# ## 1. $\xi$ is not unambigously defined.
# ### 1. Example for a one pool system

from IPython.display import display_pretty,display
from sympy import var,Function,sin,cos,diag,simplify
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
from bgc_md2.display_helpers import mass_balance_equation

# +
#define all the symbols you are going to use in the equations
var("""
    C_leaf 
    I_leaf
""")
#I_leaf=Function("I_leaf") # We can also represent functions symbolically
#I_wood=Function("I_wood") 
t=TimeSymbol("t") # the symbol used for time since it has a special role

# formulate the model
r_leaf=(sin(t)+2)
mvs = CMTVS(
    {
        StateVariableTuple( # the pool names in your preferred order
            (
                C_leaf,
            )
        ), 
        t, 
        InFluxesBySymbol({
            C_leaf: I_leaf,#(t), 
        }),
        OutFluxesBySymbol({
            C_leaf: r_leaf * C_leaf
        }),
        InternalFluxesBySymbol(),
    },
    bgc_md2_computers()
)
mass_balance_equation(mvs)
#mvs.get_CompartmentalMatrix()
#dh.
# -


# The time dependent rate $r_{leaf}(t)$ can now be separated into $r_leaf(t)=k_1 \xi_1(t)=k_2 \xi_2(t)$ where the $\xi(t)_i \; $ are  arbitrary multiples of $1+\frac{sin(t)}{2}$. 

# +
xi_1=2*sin(t)+4
K_1=simplify(r_leaf/xi_1)

xi_2=1/2*sin(t)+1
K_2=simplify(r_leaf/xi_2)

display(xi_1,K_1)
display(xi_2,K_2)
# -

# The different $K_i$ will lead to different "base line residence times" although the model as a whole has not changed at all.
# Now consider that this model is compared to another model. If there is no additional condition to define $\xi$ then the following consequences ensue: 
# 1. A comparison by size of  $\xi_{Model_a}(t)$ versus $\xi_{Model_b}(t)$ does not make sense. Since we can scale either participant by an arbitrary factor( resulting in a different $K$  
# 2. Any method trying to attribute the uncertainty to a different variable

# +
#define all the symbols you are going to use in the equations
var("""
    C_leaf 
    C_wood
    k_leaf2wood
    k_wood2leaf
    I_leaf
    I_wood
""")
#I_leaf=Function("I_leaf") # We can also represent functions symbolically
#I_wood=Function("I_wood") 
t=TimeSymbol("t") # the symbol used for time since it has a special role

# formulate the model
mvs = CMTVS(
    {
        StateVariableTuple( # the pool names in your preferred order
            (
                C_leaf,
                C_wood
            )
        ), 
        t, 
        InFluxesBySymbol({
            C_leaf: I_leaf,#(t), 
            C_wood: I_wood#(t)
        }),
        OutFluxesBySymbol({
            C_leaf: 20*(sin(t)+1)* C_leaf,
            C_wood: 1.0/10*(cos(t)+1) * C_wood
        }),
        InternalFluxesBySymbol({
            #(C_leaf, C_wood): k_leaf2wood* C_leaf, 
            #(C_wood, C_leaf): k_wood2leaf * C_wood
        }),
    },
    bgc_md2_computers()

)
# -


#start to query the model description..
M=mvs.get_CompartmentalMatrix()
M
#M.inverse_LU()

#lets define some xi
Xi_1=diag(*[sin(t)+1,5*(cos(t)+1)])
Xi_2=diag(*[20*(sin(t)+1),(cos(t)+1)])
display(Xi_1,Xi_2)

smr=mvs.get_SmoothReservoirModel()
_,A,K,_,_=smr.xi_T_N_u_representation(factor_out_xi=False)
display(A,K)

#display(A*K,M)
K_1=simplify(Xi_1.inverse_LU()*-K)
K_2=simplify(Xi_2.inverse_LU()*-K)
display(K_1,K_2)
display(A*Xi_1*K_1,A*Xi_2*K_2,-M)
#from bgc_md2.helper import compartmental_graph
#compartmental_graph(mvs)


# Now we compute the "baseline transittime"
u=sum(mvs.get_InputTuple())
beta=mvs.get_InputTuple()/u
CT_1=(A*K_1).inverse_LU()*beta
CT_2=(A*K_2).inverse_LU()*beta
display(CT_1,CT_2)

par_dict={
    k_leaf2wood: 0.1,
    k_wood2leaf: 0.2,
    I_leaf: 3,
    I_wood:4
}

display(CT_1.subs(par_dict),CT_2.subs(par_dict))


