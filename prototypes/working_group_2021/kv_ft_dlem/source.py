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

# +
# Make a small dictionary for the variables we will use
sym_dict={
    'C_leaf': '',
    'C_wood': '',
    'C_root': '',
    'C_aom1': '',
    'C_aom2': '',
    'C_smb1': '',
    'C_smb2': '',
    'C_smr': '',
    'C_nom': '',
    'C_dom': '',
    'C_psom': '',
    'r_C_leaf_2_C_aom1': '',
    'r_C_leaf_2_C_aom2': '',
    'r_C_wood_2_C_aom1': '',
    'r_C_wood_2_C_aom2': '',
    'r_C_root_2_C_aom1': '',
    'r_C_root_2_C_aom2': '',
    'r_C_aom1_2_C_smb1': '',
    'r_C_aom1_2_C_smb2': '',
    'r_C_aom1_2_C_dom': '',
    'r_C_aom1_2_C_nom': '',
    'r_C_aom2_2_C_smb1': '',
    'r_C_aom2_2_C_smb2': '',
    'r_C_aom2_2_C_dom': '',
    'r_C_smb2_2_C_smr': '',
    'r_C_smr_2_C_smb1': '',
    'r_C_dom_2_C_smb1': '',
    'r_C_dom_2_C_nom': '',
    'r_C_nom_2_C_dom': '',
    'r_C_nom_2_C_psom': '',
    'r_C_smb1_2_C_psom': '',
    'r_C_smb1_2_C_nom': '',
    'r_C_nom_2_C_smb1': '',
    'r_C_psom_2_C_smb1': '',
    'r_C_aom1_rh': '',
    'r_C_aom2_rh': '',
    'r_C_smb1_rh': '',
    'r_C_smb2_rh': '',
    'r_C_smr_rh': '',
    'r_C_dom_rh': '',
    'r_C_nom_rh': '',
    'r_C_psom_rh': '',
    'tsl': '',
    'mrso': '',
    't': '',
    'Theta_sat': '',
    'Theta_fc': '',
    'beta_leaf': '',
    'beta_wood': '',
}
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    #"I_vl": "Influx into vegetation leaf pool",
    #"k_vl_o": "out flux rate of leaf pool",
    'xi': 'a scalar function of temperature and moisture and thereby ultimately of time',
    'NPP': '',
   # 'tsl': '',
   # 'mrso': '',
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")

mvs = CMTVS(
    {
        t,
        StateVariableTuple((
            #vl, 
            #vw, 
            C_leaf,
            C_wood,
            C_root,
            C_aom1,
            C_aom2,
            C_smb1,
            C_smb2,
            C_smr,
            C_nom,
            C_dom,
            C_psom,
        )),
        InFluxesBySymbol(
            {
                #vl: I_vl, vw: I_vw
                C_leaf: NPP(t) * beta_leaf, 
                C_root: NPP(t) * (1.0-beta_leaf-beta_wood), 
                C_wood: NPP(t) * beta_wood
            }
        ),
        OutFluxesBySymbol(
            {
                #vl: k_vl_o * vl, vw: k_vw_o * vw
                C_aom1: r_C_aom1_rh*C_aom1*xi(t),
                C_aom2: r_C_aom2_rh*C_aom2*xi(t),
                C_smb1: r_C_smb1_rh*C_smb1*xi(t),
                C_smb2: r_C_smb2_rh*C_smb2*xi(t),
                C_smr: r_C_smr_rh*C_smr*xi(t),
                C_nom: r_C_nom_rh*C_nom*xi(t),
                C_dom: r_C_dom_rh*C_dom*xi(t),
                C_psom: r_C_psom_rh*C_psom*xi(t),
            }
        ),
        InternalFluxesBySymbol(
            {
                #(vl, vw): k_vl_2_vw * vl, (vw, vl): k_vw_2_vl * vw
                (C_leaf, C_aom1): r_C_leaf_2_C_aom1*C_leaf,
                (C_leaf, C_aom2): r_C_leaf_2_C_aom2*C_leaf,
                (C_wood, C_aom1): r_C_wood_2_C_aom1*C_wood,
                (C_wood, C_aom2): r_C_wood_2_C_aom2*C_wood,
                (C_root, C_aom1): r_C_root_2_C_aom1*C_root,
                (C_root, C_aom2): r_C_root_2_C_aom2*C_root,
                (C_aom1, C_smb1): r_C_aom1_2_C_smb1 * C_aom1*xi(t),
                (C_aom1, C_smb2): r_C_aom1_2_C_smb2 * C_aom1*xi(t),
                (C_aom1, C_dom): r_C_aom1_2_C_dom * C_aom1*xi(t),
                (C_aom1, C_nom): r_C_aom1_2_C_nom * C_aom1*xi(t),
                (C_aom2, C_smb1): r_C_aom2_2_C_smb1 * C_aom2*xi(t),
                (C_aom2, C_smb2): r_C_aom2_2_C_smb2 * C_aom2*xi(t),
                (C_aom2, C_dom): r_C_aom2_2_C_dom * C_aom2*xi(t),
                (C_smb2, C_smr): r_C_smb2_2_C_smr * C_smb2*xi(t),                
                (C_smr, C_smb1): r_C_smr_2_C_smb1 * C_smr*xi(t), 
                (C_dom, C_smb1): r_C_dom_2_C_smb1 * C_dom*xi(t), 
                (C_dom, C_nom): r_C_dom_2_C_nom * C_dom*xi(t),                
                (C_nom, C_dom): r_C_nom_2_C_dom * C_nom*xi(t),  
                (C_nom, C_psom): r_C_nom_2_C_psom * C_nom*xi(t),   
                (C_smb1, C_psom): r_C_smb1_2_C_psom * C_smb1*xi(t),
                (C_smb1, C_nom): r_C_smb1_2_C_nom * C_smb1*xi(t),
                (C_nom, C_smb1): r_C_nom_2_C_smb1 * C_nom*xi(t), 
                (C_psom, C_smb1): r_C_psom_2_C_smb1 * C_psom*xi(t),
            }
        ),
        BibInfo(# Bibliographical Information
            name="Visit",
            longName="",
            version="1",
            entryAuthor="Kostiantyn Viatkin",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        ),
    },


    computers=module_computers(bgc_c)
)
