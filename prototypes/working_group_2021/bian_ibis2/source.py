from sympy import  Symbol, Function
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
sym_dict={    
    # Carbon Pools
    "C_leaf": "vegetation leaf pool",
    "C_wood": "",
    "C_root": "",
    "C_mll": "Metabolic Leaf Litter pool",
    "C_mwl": "Metabolic Wood Litter pool",
    "C_mrl": "Metabolic Root Litter pool",
    "C_sll": "Structural Leaf Litter Pool",
    "C_swl": "Structural Wood Litter Pool",
    "C_srl": "Structural Root Litter Pool",
    "C_lll": "Lignin Leaf Litter Pool",
    "C_lwl": "Lignin Wood Litter Pool",
    "C_lrl": "Lignin Root Litter Pool",
    "C_mic": "microbial",
    "C_prot": "protected som",
    "C_nonprot": "non-protected som",
    "C_pass": "passive som", 
    
    "r_C_mll_rh": "",
    "r_C_mwl_rh": "",
    "r_C_mrl_rh": "",
    "r_C_sll_rh": "",
    "r_C_swl_rh": "",
    "r_C_srl_rh": "",
    "r_C_lll_rh": "",
    "r_C_lwl_rh": "",
    "r_C_lrl_rh": "",
    "r_C_mic_rh": "",
    "r_C_prot_rh": "",
    "r_C_nonprot_rh": "",
    "r_C_pass_rh": "",
    
    "r_C_leaf_2_C_mll": "f4_1",
    "r_C_wood_2_C_mwl": "f5_2",
    "r_C_root_2_C_mrl": "f6_3",
    "r_C_leaf_2_C_sll": "f7_1",
    "r_C_wood_2_C_swl": "f8_2",
    "r_C_root_2_C_srl": "f9_3",
    "r_C_leaf_2_C_lll": "f10_1",
    "r_C_wood_2_C_lwl": "f11_2",
    "r_C_root_2_C_lrl": "f12_3",
    
    "r_C_mll_2_C_mic": "f13_4",
    "r_C_mwl_2_C_mic": "f13_5",
    "r_C_mrl_2_C_mic": "f13_6",
    "r_C_sll_2_C_mic": "f13_7",
    "r_C_swl_2_C_mic": "f13_8",
    "r_C_srl_2_C_mic": "f13_9",
    "r_C_prot_2_C_mic": "f13_14",
    "r_C_nonprot_2_C_mic": "f13_15",
    "r_C_pass_2_C_mic": "f13_16",
    
    "r_C_lll_2_C_prot": "f14_10",
    "r_C_lwl_2_C_prot": "f14_11",
    "r_C_lrl_2_C_prot": "f14_12",
    "r_C_mic_2_C_prot": "f14_13",
    
    "r_C_lll_2_C_nonprot": "f15_10",
    "r_C_lwl_2_C_nonprot": "f15_11",
    "r_C_lrl_2_C_nonprot": "f15_12",
    "r_C_mic_2_C_nonprot": "f15_13",
 
    "r_C_prot_2_C_pass": "f16_14",
    "r_C_nonprot_2_C_pass": "f16_15",
    
    "tas": "",
    "mrso": "",
    "t": "",
    "beta_leaf": "",
    "beta_wood": "",
    "beta_root": "",
}
# Make symbols from  the strings that we can later use in expressions  
# vl, vw,...
for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument) 
func_dict={
    # "I_vl": "Influx into vegetation leaf pool",
    # "r_vl_o": "out flux rate of leaf pool",
    "xi": "Environmental scalar as a function of time",
    "NPP": "",
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")
beta_root = 1.0 - (beta_leaf+beta_wood)

mvs = CMTVS(
    {
        t,
        StateVariableTuple((
            C_leaf,
            C_wood,
            C_root,
            C_mll,
            C_mwl,
            C_mrl,
            C_sll,
            C_swl,
            C_srl,
            C_lll,
            C_lwl,
            C_lrl,
            C_mic,
            C_prot,
            C_nonprot,
            C_pass,  
        )),
        InFluxesBySymbol(
            {
            # Recieving Pool: Input * Allocation
            C_leaf: NPP(t) * beta_leaf,
            C_wood: NPP(t) * beta_wood,
            C_root: NPP(t) * beta_root
            }
        ),
        OutFluxesBySymbol({
            # define fluxes leaving the system
            # fluxes leaving the system: FluxRate * DonorPool * Environment scalers
            
            C_mll: r_C_mll_rh * C_mll * xi(t),
            C_mwl: r_C_mwl_rh * C_mwl * xi(t),
            C_mrl: r_C_mrl_rh * C_mrl * xi(t),
            C_sll: r_C_sll_rh * C_sll * xi(t),
            C_swl: r_C_swl_rh * C_swl * xi(t),
            C_srl: r_C_srl_rh * C_srl * xi(t),
            C_lll: r_C_lll_rh * C_lll * xi(t),
            C_lwl: r_C_lwl_rh * C_lwl * xi(t),
            C_lrl: r_C_lrl_rh * C_lrl * xi(t),
            C_mic: r_C_mic_rh * C_mic * xi(t),
            C_prot: r_C_prot_rh * C_prot * xi(t),
            C_nonprot: r_C_nonprot_rh * C_nonprot * xi(t),
            C_pass: r_C_pass_rh * C_pass * xi(t) 
        }),
        InternalFluxesBySymbol({
            # define tranfer matrix
            #(Donor Pool, recieving pool): FluxRate * DonorPool
            (C_leaf, C_mll): r_C_leaf_2_C_mll * C_leaf,
            (C_leaf, C_sll): r_C_leaf_2_C_sll * C_leaf,
            (C_leaf, C_lll): r_C_leaf_2_C_lll * C_leaf,
            
            (C_wood, C_mwl): r_C_wood_2_C_mwl * C_wood,
            (C_wood, C_swl): r_C_wood_2_C_swl * C_wood,
            (C_wood, C_lwl): r_C_wood_2_C_lwl * C_wood,
            
            (C_root, C_mrl): r_C_root_2_C_mrl * C_root,
            (C_root, C_srl): r_C_root_2_C_srl * C_root,
            (C_root, C_lrl): r_C_root_2_C_lrl * C_root,

            (C_mll, C_mic): r_C_mll_2_C_mic * C_mll * xi(t),
            (C_mwl, C_mic): r_C_mwl_2_C_mic * C_mwl * xi(t),
            (C_mrl, C_mic): r_C_mrl_2_C_mic * C_mrl * xi(t),
            (C_sll, C_mic): r_C_sll_2_C_mic * C_sll * xi(t),
            (C_swl, C_mic): r_C_swl_2_C_mic * C_swl * xi(t),
            (C_srl, C_mic): r_C_srl_2_C_mic * C_srl * xi(t),
            (C_prot, C_mic): r_C_prot_2_C_mic * C_prot * xi(t),
            (C_nonprot, C_mic): r_C_nonprot_2_C_mic * C_nonprot * xi(t),
            (C_pass, C_mic): r_C_pass_2_C_mic * C_pass * xi(t),

            (C_lll, C_prot): r_C_lll_2_C_prot * C_lll * xi(t),
            (C_lwl, C_prot): r_C_lwl_2_C_prot * C_lwl * xi(t),
            (C_lrl, C_prot): r_C_lrl_2_C_prot * C_lrl * xi(t),
            (C_mic, C_prot): r_C_mic_2_C_prot * C_mic * xi(t),
            
            (C_lll, C_nonprot): r_C_lll_2_C_nonprot * C_lll * xi(t),
            (C_lwl, C_nonprot): r_C_lwl_2_C_nonprot * C_lwl * xi(t),
            (C_lrl, C_nonprot): r_C_lrl_2_C_nonprot * C_lrl * xi(t),
            (C_mic, C_nonprot): r_C_mic_2_C_nonprot * C_mic * xi(t),

            (C_prot, C_pass): r_C_prot_2_C_pass * C_prot * xi(t),
            (C_nonprot, C_pass): r_C_nonprot_2_C_pass * C_nonprot * xi(t),
        }),
        BibInfo(# Bibliographical Inoformation
            name="IBIS",
            longName="",
            version="1.0",
            entryAuthor="Chenyu Bian",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict   
        ),
    },

    computers=module_computers(bgc_c)
)
