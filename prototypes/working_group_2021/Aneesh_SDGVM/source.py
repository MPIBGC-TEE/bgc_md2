import sys
sys.path.insert(0,'..') # necessary to import general_helpers
from sympy import Symbol, Function 
from bgc_md2.models.BibInfo import BibInfo
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
        'beta_leaf': 'NPP partitioning to leaf', 
        'beta_wood': 'NPP partitioning to wood',
        'beta_root': 'NPP partitioning to root',
         #'k_C_leaf': 'Turn over rate of leaf pool',
         #'k_C_wood':'Turn over rate of wood pool',
         #'k_C_root': 'Turn over rate of root pool',
         #'k_C_abvmetlit':'Turn over rate of aboveground metabolic litter pool',
         #'k_C_abvstrlit' :'Turn over rate of aboveground structural litter pool',
         #'k_C_belowstrlit' :'Turn over rate of belowground structural litter pool',
         #'k_C_belowmetlit':'Turn over rate of belowground metabolic litter pool',
         #'k_C_surface_microbe' :'Turn over rate of surface microbe pool',
         #'k_C_soil_microbe' :'Turn over rate of soil microbial pool',
         #'k_C_slowsom' :'Turn over rate of slow soil pool',
         #'k_C_passsom':'Turn over rate of passive soil pool',
         #'ls_aboveground': 'Fraction of structural C that is lignin in aboveground',
         #'ls_belowground': 'Fraction of structural C that is lignin in belowground',
         #'f_leaf2abvstrlit': 'Transfer coefficient of C from leaf to aboveground structural litter',
         #'f_leaf2abvmetlit':'Transfer coefficient of C from leaf to aboveground metabolic litter',
         #'f_wood2abvstrlit':'Transfer coefficient of C from wood to aboveground structural litter',
         #'f_wood2abvmetlit':'Transfer coefficient of C from wood to aboveground metabolic litter',
         #'f_root2belowstrlit':'Transfer coefficient of C from root to belowground structural litter',
         #'f_root2belowmetlit':'Transfer coefficient of C from root to belowground metabolic litter',
         #'f_abvstrlit2surface_microbe':'Transfer coefficient of C from aboveground structural litter to surface microbe',
         #'f_abvstrlit2slowsom':'Transfer coefficient of C from aboveground structural litter to slow soil',
         #'f_abvmetlit2surface_microbe':'Transfer coefficient of C from aboveground metabolic litter to surface microbe',
         #'f_belowstrlit2soil_microbe':'Transfer coefficient of C from belowground structural litter to soil microbe',
         #'f_belowstrlit2slowsom':'Transfer coefficient of C from belowground structural litter to slow soil',
         #'f_belowmetlit2soil_microbe':'Transfer coefficient of C from belowground metabolic litter to soil microbe',
         #'f_slowsom2soil_microbe':'Transfer coefficient of C from soil soil to soil microbe',
         #'f_passsom2soil_microbe':'Transfer coefficient of C from passive soil to soil microbe',
         #'f_surface_microbe2slowsom':'Transfer coefficient of C from surface microbe to slow soil',
         #'f_soil_microbe2passsom':'Transfer coefficient of C from soil microbe to passive soil',
         #'f_soil_microbe2slowsom':'Transfer coefficient of C from soil microbe to slow soil',
         #'f_slowsom2passsom':'Transfer coefficient of C from slow soil to passive soil',
        'Ca': 'Need to find (unknown)',
        'clay': 'Clay content (fraction)',
        'silt_clay': 'Silt plus clay content (fraction)',
        'leachedwtr30': 'H20 leached below 30cm (cm/m)',
        'r_C_leaf2abvstrlit':'',
        'r_C_leaf2abvmetlit':'',
        'r_C_wood2abvstrlit':'', 
        'r_C_wood2abvmetlit' :'',
        'r_C_root2belowstrlit':'',  
        'r_C_root2belowmetlit':'',
        'r_C_abvstrlit2surface_microbe':'',
        'r_C_abvmetlit2surface_microbe':'', 
        'r_C_abvstrlit2slowsom':'',
        'r_C_belowstrlit2soil_microbe':'',  
        'r_C_belowmetlit2soil_microbe':'',
        'r_C_belowstrlit2slowsom':'' ,
        'r_C_surface_microbe2slowsom':'',
        'r_C_soil_microbe2slowsom' :'',''
        'r_C_slowsom2soil_microbe':'',
        'r_C_soil_microbe2passsom':'',
        'r_C_slowsom2passsom':'',
        'r_C_passsom2soil_microbe':'',
        'r_C_leached': '',
        'r_C_abvstrlit_rh':'',
        'r_C_abvmetlit_rh':'',
        'r_C_belowstrlit_rh':'',
        'r_C_belowmetlit_rh':'',
        #'r_C_surface_microbe_rh':'',
        'r_C_slowsom_rh':'',
        'r_C_passsom_rh':'',
        'r_C_soil_microbe_rh' :''
}

for k in sym_dict.keys():
    code=k+" = Symbol('{0}')".format(k)
    exec(code)

# some we will also use some symbols for functions (which appear with an argument)
func_dict={
    'xi': 'a scalar function of temperature and moisture and thereby ultimately of time',
    'NPP': '',
}

for k in func_dict.keys():
    code=k+" = Function('{0}')".format(k)
    exec(code)

t=TimeSymbol("t")
beta_root=1-(beta_leaf+beta_wood)
mvs = CMTVS(
    {
        StateVariableTuple((C_leaf, 
                            C_wood, 
                            C_root, 
                            C_abvstrlit,
                            C_abvmetlit, 
                            C_belowstrlit, 
                            C_belowmetlit, 
                            C_surface_microbe, 
                            C_soil_microbe, 
                            C_slowsom, 
                            C_passsom)),
        t,
        InFluxesBySymbol(
            {
                C_leaf: NPP(t)* beta_leaf, 
                C_wood: NPP(t)* beta_wood, 
                C_root:NPP(t)*  beta_root
            }),
        OutFluxesBySymbol(
            {
                C_abvstrlit: r_C_abvstrlit_rh * C_abvstrlit,
                C_abvmetlit: r_C_abvmetlit_rh * C_abvmetlit, 
                C_belowstrlit: r_C_belowstrlit_rh * C_belowstrlit, 
                C_belowmetlit: r_C_belowmetlit_rh * C_belowmetlit, 
                C_surface_microbe: r_C_leached* C_surface_microbe, 
                C_soil_microbe: r_C_soil_microbe_rh * C_soil_microbe, 
                C_slowsom: r_C_slowsom_rh * C_slowsom, 
                C_passsom: r_C_passsom_rh * C_passsom
            }
        ),
        InternalFluxesBySymbol(
            {
                                (C_leaf, C_abvstrlit): r_C_leaf2abvstrlit* C_leaf, 
                                (C_leaf, C_abvmetlit): r_C_leaf2abvmetlit* C_leaf, 
                                (C_wood, C_abvstrlit): r_C_wood2abvstrlit* C_wood, 
                                (C_wood, C_abvmetlit): r_C_wood2abvmetlit* C_wood, 
                                (C_root, C_belowstrlit): r_C_root2belowstrlit* C_root, 
                                (C_root, C_belowmetlit): r_C_root2belowmetlit * C_root, 
                                (C_abvstrlit , C_surface_microbe ): r_C_abvstrlit2surface_microbe* C_abvstrlit, 
                                (C_abvstrlit , C_slowsom ): r_C_abvstrlit2slowsom*C_abvstrlit,
                                (C_abvmetlit, C_surface_microbe ): r_C_abvmetlit2surface_microbe* C_abvmetlit, 
                                (C_belowstrlit, C_soil_microbe): r_C_belowstrlit2soil_microbe*C_belowstrlit, 
                                (C_belowmetlit , C_soil_microbe): r_C_belowmetlit2soil_microbe*C_belowmetlit, 
                                (C_belowstrlit, C_slowsom): r_C_belowstrlit2slowsom* C_belowstrlit, 
                                (C_surface_microbe , C_slowsom): r_C_surface_microbe2slowsom*C_surface_microbe, 
                                (C_soil_microbe, C_slowsom): r_C_soil_microbe2slowsom*C_soil_microbe, 
                                (C_slowsom, C_soil_microbe): r_C_slowsom2soil_microbe*C_slowsom, 
                                (C_soil_microbe, C_passsom): r_C_soil_microbe2passsom*C_soil_microbe, 
                                (C_slowsom , C_passsom): r_C_slowsom2passsom*C_slowsom,
                               (C_passsom, C_soil_microbe): r_C_passsom2soil_microbe * C_passsom
            }
        ),
        BibInfo(# Bibliographical Information
            name="SDGVM",
            longName="",
            version="1",
            entryAuthor="Aneesh",
            entryAuthorOrcid="",
            entryCreationDate="",
            doi="",
            sym_dict=sym_dict,
            func_dict=func_dict
        ),
    },

    computers=module_computers(bgc_c)
)
