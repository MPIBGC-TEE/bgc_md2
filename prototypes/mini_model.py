from sympy import var,Function
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)

#define all the symbols you are going to use in the equations
var("""
    C_leaf 
    C_wood
    k_leaf2wood
    k_wood2leaf
    k_leaf2out
    k_wood2out
""")
I_leaf=Function("I_leaf") # We can also represent functions symbolically
I_wood=Function("I_wood") 
t=TimeSymbol("t") # the symbol used for time since it has a special role


# +
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
            C_leaf: I_leaf(t), 
            C_wood: I_wood(t)
        }),
        OutFluxesBySymbol({
            C_leaf: k_leaf2out * C_leaf,
            C_wood: k_wood2out * C_wood
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
M.inverse_LU()

mvs.get_InputTuple()

mvs.get_StateVariableTuple()

from bgc_md2.helper import compartmental_graph
compartmental_graph(mvs)


from bgc_md2.display_helpers import mass_balance_equation
mass_balance_equation(mvs)

# for comparison the century model as found in our database
from bgc_md2.models.Parton1987SoilSciSocAmJ.source_by_name import mvs as mvs_century


mvs_century.get_CompartmentalMatrix()

mvs_century.get_InputTuple()

compartmental_graph(mvs_century)

BI=mvs_century.get_BibInfo()
BI.sym_dict


