from sympy import Symbol, Function
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import module_computers
from bgc_md2.models.BibInfo import BibInfo
from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.resolve.computers as bgc_c

# Make a small dictionary for the variables we will use
sym_dict = {
    'C_wood1': '',
    'C_wood2': '',
    'C_wood3': '',
    'C_wood4': '',
    'C_leaf': '',
    'C_root': '',
    'C_fruit': '',
    'C_litter1': '',
    'C_litter2': '',
    'C_litter3': '',
    'C_litter4': '',
    'C_litter5': '',
    'C_litter6': '',
    'C_som1': '',
    'C_som2': '',
    'C_som3': '',
    'C_som4': '',
    'r_C_wood1_2_C_wood3': '',
    'r_C_wood1_2_C_litter1': '',
    'r_C_wood2_2_C_wood4': '',
    'r_C_wood2_2_C_litter2': '',
    'r_C_wood3_2_C_litter1': '',
    'r_C_wood4_2_C_litter2': '',
    'r_C_leaf_2_C_litter3': '',
    'r_C_leaf_2_C_litter5': '',
    'r_C_root_2_C_litter4': '',
    'r_C_root_2_C_litter6': '',
    'r_C_fruit_2_C_litter3': '',
    'r_C_fruit_2_C_litter5': '',
    'r_C_litter1_2_C_som1': '',
    'r_C_litter1_2_C_som2': '',
    'r_C_litter2_2_C_som2': '',
    'r_C_litter2_2_C_som3': '',
    'r_C_litter3_2_C_som1': '',
    'r_C_litter3_2_C_som3': '',
    'r_C_litter4_2_C_som1': '',
    'r_C_litter4_2_C_som2': '',
    'r_C_litter5_2_C_som1': '',
    'r_C_litter6_2_C_som2': '',
    'r_C_som1_2_C_som3': '',
    'r_C_som2_2_C_som3': '',
    'r_C_som2_2_C_som4': '',
    'r_C_som3_2_C_som2': '',
    'r_C_som3_2_C_som4': '',
    'r_C_som4_2_C_som2': '',
    'r_C_som4_2_C_som3': '',
    'r_C_litter1_rh': '',
    'r_C_litter2_rh': '',
    'r_C_litter3_rh': '',
    'r_C_litter4_rh': '',
    'r_C_litter5_rh': '',
    'r_C_litter6_rh': '',
    'r_C_som1_rh': '',
    'r_C_som2_rh': '',
    'r_C_som3_rh': '',
    'r_C_som4_rh': '',
    'Ts': '',
    'mrso': '',
    't': '',
    'T_0': '',
    'E': '',
    'KM': '',
    'beta_wood1': '',
    'beta_wood2': '',
    'beta_leaf': '',
    'beta_root': '',
}
for k in sym_dict.keys():
    code = k + " = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument)
func_dict = {
    'xi': 'a scalar function of temperature and moisture and thereby ultimately of time',
    'NPP': '',
}
for k in func_dict.keys():
    code = k + " = Function('{0}')".format(k)
    exec(code)

t = TimeSymbol("t")
beta_fruit = 1.0 - (beta_leaf + beta_wood1 + beta_wood2 + beta_root)
mvs = CMTVS(
    {
        t,
        StateVariableTuple((
            C_wood1,
            C_wood2,
            C_wood3,
            C_wood4,
            C_leaf,
            C_root,
            C_fruit,
            C_litter1,
            C_litter2,
            C_litter3,
            C_litter4,
            C_litter5,
            C_litter6,
            C_som1,
            C_som2,
            C_som3,
            C_som4,
        )),
        InFluxesBySymbol(
            {
                C_leaf: NPP(t) * beta_leaf,
                C_root: NPP(t) * beta_root,
                C_wood1: NPP(t) * beta_wood1,
                C_wood2: NPP(t) * beta_wood2,
                C_fruit: NPP(t) * beta_fruit,
            }
        ),
        OutFluxesBySymbol(
            {
                C_litter1: r_C_litter1_rh * C_litter1 * xi(t),
                C_litter2: r_C_litter2_rh * C_litter2 * xi(t),
                C_litter3: r_C_litter3_rh * C_litter3 * xi(t),
                C_litter4: r_C_litter4_rh * C_litter4 * xi(t),
                C_litter5: r_C_litter5_rh * C_litter5 * xi(t),
                C_litter6: r_C_litter6_rh * C_litter6 * xi(t),
                C_som1: r_C_som1_rh * C_som1 * xi(t),
                C_som2: r_C_som2_rh * C_som2 * xi(t),
                C_som3: r_C_som3_rh * C_som3 * xi(t),
                C_som4: r_C_som4_rh * C_som4 * xi(t),
            }
        ),
        InternalFluxesBySymbol(
            {
                (C_wood1, C_wood3): r_C_wood1_2_C_wood3 * C_wood1,
                (C_wood1, C_litter1): r_C_wood1_2_C_litter1 * C_wood1,
                (C_wood2, C_wood4): r_C_wood2_2_C_wood4 * C_wood2,
                (C_wood2, C_litter2): r_C_wood2_2_C_litter2 * C_wood2,
                (C_wood3, C_litter1): r_C_wood3_2_C_litter1 * C_wood3,
                (C_wood4, C_litter2): r_C_wood4_2_C_litter2 * C_wood4,
                (C_leaf, C_litter3): r_C_leaf_2_C_litter3 * C_leaf,
                (C_leaf, C_litter5): r_C_leaf_2_C_litter5 * C_leaf,
                (C_root, C_litter4): r_C_root_2_C_litter4 * C_root,
                (C_root, C_litter6): r_C_root_2_C_litter6 * C_root,
                (C_fruit, C_litter3): r_C_fruit_2_C_litter3 * C_fruit,
                (C_fruit, C_litter5): r_C_fruit_2_C_litter5 * C_fruit,

                (C_litter1, C_som1): r_C_litter1_2_C_som1 * C_litter1 * xi(t),
                (C_litter1, C_som2): r_C_litter1_2_C_som2 * C_litter1 * xi(t),
                (C_litter2, C_som2): r_C_litter2_2_C_som2 * C_litter2 * xi(t),
                (C_litter2, C_som3): r_C_litter2_2_C_som3 * C_litter2 * xi(t),
                (C_litter3, C_som1): r_C_litter3_2_C_som1 * C_litter3 * xi(t),
                (C_litter3, C_som3): r_C_litter3_2_C_som3 * C_litter3 * xi(t),
                (C_litter4, C_som1): r_C_litter4_2_C_som1 * C_litter4 * xi(t),
                (C_litter4, C_som2): r_C_litter4_2_C_som2 * C_litter4 * xi(t),
                (C_litter5, C_som1): r_C_litter5_2_C_som1 * C_litter5 * xi(t),
                (C_litter6, C_som2): r_C_litter6_2_C_som2 * C_litter6 * xi(t),

                (C_som1, C_som3): r_C_som1_2_C_som3 * C_som1 * xi(t),
                (C_som2, C_som3): r_C_som2_2_C_som3 * C_som2 * xi(t),
                (C_som2, C_som4): r_C_som2_2_C_som4 * C_som2 * xi(t),
                (C_som3, C_som2): r_C_som3_2_C_som2 * C_som3 * xi(t),
                (C_som3, C_som4): r_C_som3_2_C_som4 * C_som3 * xi(t),
                (C_som4, C_som2): r_C_som4_2_C_som2 * C_som4 * xi(t),
                (C_som4, C_som3): r_C_som4_2_C_som3 * C_som4 * xi(t),
            }
        ),
        BibInfo(  # Bibliographical Information
            name="OCN",
            longName="",
            version="1",
            entryAuthor="Song Wang",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        ),
    },

    computers=module_computers(bgc_c)
)
