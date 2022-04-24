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
sym_dict={
    'C_NWT': 'Non-Woody Tree Parts',
    'C_AGWT': 'Aboveground-Woody Tree Parts',
    'C_TR': 'Tree Roots',
    'C_GVF': 'Ground Vegetation Foliage',
    'C_GVR': 'Ground Vegetation Roots',
    'C_AGML': 'AG: Metabolic Litter',
    'C_AGSL': 'AG: Structural Litter',
    'C_AGMS': 'AG: Microbial Soil',
    'C_YHMS': 'AG: Young Humus Soil',
    'C_BGDL': 'BG: Decomposable Litter',
    'C_BGRL': 'BG: Resistant Litter',
    'C_BGMS': 'BG: Microbial Soil',
    'C_SHMS': 'BG: Stabilized Humus Soil',

    'r_C_AGML_rh': '',
    'r_C_AGSL_rh': '',
    'r_C_AGMS_rh': '',
    'r_C_YHMS_rh': '',
    'r_C_BGDL_rh': '',
    'r_C_BGRL_rh': '',
    'r_C_BGMS_rh': '',
    'r_C_SHMS_rh': '',
    
    'r_C_NWT_2_C_AGML': '',
    'r_C_NWT_2_C_AGSL': '',
    'r_C_AGWT_2_C_AGSL': '',
    'r_C_TR_2_C_BGDL': '',
    'r_C_TR_2_C_BGRL': '',
    'r_C_GVF_2_C_AGML': '',
    'r_C_GVF_2_C_AGSL': '',
    'r_C_GVR_2_C_BGDL': '',
    'r_C_GVR_2_C_BGRL': '',
    
    'r_C_AGML_2_C_AGMS': '',
    'r_C_AGSL_2_C_AGMS': '',
    'r_C_AGSL_2_C_YHMS': '',
    'r_C_AGMS_2_C_YHMS': '',
    'r_C_YHMS_2_C_AGMS': '',
    'r_C_YHMS_2_C_SHMS': '',
    
    'r_C_BGDL_2_C_SHMS': '',
    'r_C_BGRL_2_C_SHMS': '',
    'r_C_BGRL_2_C_BGMS': '',
    'r_C_SHMS_2_C_BGMS': '',
    'r_C_BGMS_2_C_SHMS': '',
    
    'tas': '',
    'mrso': '',
    't': '',
#    'T_0': '',
#    'E': '',
#    'KM': '',
    'fco': '',
    'fml': '',
    'fd': '',
    'fwt': '',
    'fgv': '',
    'beta_NWT': '',
    'beta_AGWT': '',
    'beta_TR': '',
    'beta_GVF': '',
    'beta_GVR': '',
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
}
for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")
beta_GVR = 1.0- (beta_NWT+beta_AGWT+beta_TR+beta_GVF)
mvs = CMTVS(
    {
        t,
        StateVariableTuple((
            C_NWT,
	        C_AGWT,
	        C_TR,
	        C_GVF,
	        C_GVR,
	        C_AGML,
	        C_AGSL,
	        C_AGMS,
	        C_YHMS,
	        C_BGDL,
	        C_BGRL,
	        C_BGMS,
	        C_SHMS,            
        )),
        InFluxesBySymbol(
            {
                C_NWT: NPP(t) * beta_NWT, 
                C_AGWT: NPP(t) * beta_AGWT, 
                C_TR: NPP(t) * beta_TR,
                C_GVF: NPP(t) * beta_GVF,
                C_GVR: NPP(t) * beta_GVR,
            }
        ),
        OutFluxesBySymbol(
            {
                C_AGML: r_C_AGML_rh*C_AGML*xi(t),
                C_AGSL: r_C_AGSL_rh*C_AGSL*xi(t),
                C_AGMS: r_C_AGMS_rh*C_AGMS*xi(t),
                C_YHMS: r_C_YHMS_rh*C_YHMS*xi(t),
                C_BGDL: r_C_BGDL_rh*C_BGDL*xi(t),
                C_BGRL: r_C_BGRL_rh*C_BGRL*xi(t),
                C_BGMS: r_C_BGMS_rh*C_BGMS*xi(t),
                C_SHMS: r_C_SHMS_rh*C_SHMS*xi(t),
            }
        ),
        InternalFluxesBySymbol(
            {
                (C_NWT,C_AGML): r_C_NWT_2_C_AGML*C_NWT,
                (C_NWT,C_AGSL): r_C_NWT_2_C_AGSL*C_NWT,
                (C_AGWT,C_AGSL): r_C_AGWT_2_C_AGSL*C_AGWT,
                (C_TR,C_BGDL): r_C_TR_2_C_BGDL*C_TR,
                (C_TR,C_BGRL): r_C_TR_2_C_BGRL*C_TR,
                (C_GVF,C_AGML): r_C_GVF_2_C_AGML*C_GVF,
                (C_GVF,C_AGSL): r_C_GVF_2_C_AGSL*C_GVF,
                (C_GVR,C_BGDL): r_C_GVR_2_C_BGDL*C_GVR,
                (C_GVR,C_BGRL): r_C_GVR_2_C_BGRL*C_GVR,
                (C_AGML,C_AGMS): r_C_AGML_2_C_AGMS*C_AGML*xi(t),
                (C_AGSL,C_AGMS): r_C_AGSL_2_C_AGMS*C_AGSL*xi(t),
                (C_AGSL,C_YHMS): r_C_AGSL_2_C_YHMS*C_AGSL*xi(t),
                (C_AGMS,C_YHMS): r_C_AGMS_2_C_YHMS*C_AGMS*xi(t),
                (C_YHMS,C_AGMS): r_C_YHMS_2_C_AGMS*C_YHMS*xi(t),
                (C_YHMS,C_SHMS): r_C_YHMS_2_C_AGMS*C_YHMS*xi(t),
                (C_BGDL,C_SHMS): r_C_BGDL_2_C_SHMS*C_BGDL*xi(t),
                (C_BGRL,C_SHMS): r_C_BGRL_2_C_SHMS*C_BGRL*xi(t),
                (C_BGRL,C_BGMS): r_C_BGRL_2_C_BGMS*C_BGRL*xi(t),
                (C_SHMS,C_BGMS): r_C_SHMS_2_C_BGMS*C_SHMS*xi(t),
                (C_BGMS,C_SHMS): r_C_BGMS_2_C_SHMS*C_BGMS*xi(t),
            }
        ),
        BibInfo(# Bibliographical Information
            name="ISAM",
            longName="",
            version="1",
            entryAuthor="Cuijuan Liao",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        ),
    },


    computers=module_computers(bgc_c)
)

