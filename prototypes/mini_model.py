from sympy import var,Function,Symbol
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)

# +
# obtain dictionary of symbols for equations from txt file
params_in = open(r"mini_model_sym_dict.txt", 'r')
sym_dict = {}
for line in params_in:
    k, v = line.strip().split('=')
    sym_dict[k.strip()] = v.strip()
    
params_in.close()

# -

#define all the symbols you are going to use in the equations
sym_dict={
    "AGLIVC": "Above ground live grass C pool ", 
    "BGLIVC": "Below ground live grass root C pool ",
    "STDEDC": "Standing dead C pool ", 
    "STRUCC_1": "Surface structural C pool",
    "STRUCC_2": "Below ground structural C pool",
    "SOM1C_1": "Surface microbe C pool",
    "SOM1C_2": "Active organic C pool",
    "SOM2C": "Slow organic C pool",
    "SOM3C": "Passive soil organic matter C pool",
    "METABC_1": "Surface metabolic C ",
    "METABC_2": "Below ground metabolic C ",
    "CROOTC": "C in forest system coarse root component",
    "FBRCHC": "C in forest system fine branch component",
    "FROOTC": "C in forest system fine root component",
    "RLEAVC": "C in forest system leaf component",
    "RLWODC": "C in forest system large wood component",
    "WOOD1C": "C in dead fine branch component of forest system",
    "WOOD2C": "C in dead large wood component of forest system",
    "WOOD3C": "C in dead coarse roots component of forest system",
    "cprodc_A": "fraction of Gpp that goes to the above ground part",
    "cprodc_B": "fraction of Gpp that goes to the below ground part",
}

# +
#define all the symbols you are going to use in the equations

for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)
    
func_dict = {
    "Gpp_grass": "Total C production for grass in g/m^2*month FIXME"
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)    
    
#I_wood=Function("I_wood") 
t=TimeSymbol("t") # the symbol used for time since it has a special role
e = Symbol("e")   # for exponential functions

cprodc_B = 1 - cprodc_A
k_STRUCC_1 = dec1_1 * defac_0 * (e**(-3*strlig_1))

# formulate the model
mvs = CMTVS(
    {
        StateVariableTuple( # the pool names in your preferred order
            (
                AGLIVC,
                BGLIVC,
                STDEDC,
                STRUCC_1,
                STRUCC_2,
                SOM1C_1,
                SOM1C_2,
                SOM2C,
                SOM3C,
                METABC_1,
                METABC_2,
                CROOTC,
                FBRCHC,
                FROOTC,
                RLEAVC,
                RLWODC,
                WOOD1C,
                WOOD2C,
                WOOD3C,
            )
        ), 
        t, 
        InFluxesBySymbol({
            AGLIVC: Gpp_grass(t)*cprodc_A, 
            BGLIVC: Gpp_grass(t)*cprodc_B
        }),
        OutFluxesBySymbol({
            STDEDC: fallrt * STDEDC,
            STRUCC_1: (k_STRUCC_1) * (0.3 * strlig_1 + (0.55 - 0.55 * strlig_1)) * STRUCC_1, 
            STRUCC_2: (dec1_2 * defac_0 * (e**(-3*strlig_2))) * STRUCC_2  # anerb=1 see litdec.F
        }),
        InternalFluxesBySymbol({
            (STRUCC_1, SOM1C_1): ((k_STRUCC_1 - strlig_1 * k_STRUCC_1) - (k_STRUCC_1 - strlig_1 * k_STRUCC_1) * 0.55) * STRUCC_1,
            (STRUCC_1, SOM2C): k_STRUCC_1 * (strlig_1 - strlig_1 * 0.3) * STRUCC_1  
        }),
    },
    bgc_md2_computers()

)
# -


#start to query the model description..
M=mvs.get_CompartmentalMatrix()
#M.inverse_LU()

mvs.get_InputTuple()

mvs.get_StateVariableTuple()

from bgc_md2.helper import compartmental_graph
compartmental_graph(mvs)


from bgc_md2.display_helpers import mass_balance_equation
mass_balance_equation(mvs)

# for comparison the century model as found in our database
from bgc_md2.models.Parton1987SoilSciSocAmJ.source_by_name import mvs as mvs_century


mvs.computable_mvar_types

mvs_century.get_InputTuple()

compartmental_graph(mvs_century)

mass_balance_equation(mvs_century)

BI=mvs_century.get_BibInfo()
BI.sym_dict

x=Symbol("x")

# +

s=x**2-x

s
# -

type(s)
