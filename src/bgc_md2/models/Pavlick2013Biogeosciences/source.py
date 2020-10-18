from sympy import var, ImmutableMatrix, exp, Max 
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
)
from ..BibInfo import BibInfo 
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet

sym_dict = {
        'C_A': 'Carbon in stored assimilates' #"gC*m^{-2}" 
        ,'C_L': 'Carbon in leaves' #"gC*m^{-2}" 
        ,'C_R': 'Carbon in fine roots' #"gC*m^{-2}" 
        ,'C_WL': 'Carbon in aboveground wood (branches and stems)' #"gC*m^{-2}" 
        ,'C_WR': 'Carbon in belowground wood (coarse roots)' #"gC*m^{-2}" 
        ,'C_S': 'Carbon in seeds (reproductive tisses)' #"gC*m^{-2}" 
        ,'SLA': 'Specific Leaf Area ' #"m^2*gC^{-1}" exprs: #Function of leaf turnover rate. See page 23
        ,'LAI': 'Leaf Area Index' #"m^2*m^{-2}" 
        ,'k': 'Light extinction coefficient'
        ,'f_VEG': 'Fractional vegetative cover'
        ,'GPP': 'Gross Primary Production' #"gC*m^{-2}*d^{-1}" # exprs: #Function of LAI, see equations in page 24
        ,'RES_a': 'Autotrophic respiration' #"gC*m^{-2}*d^{-1}" 
        ,'NPP': 'Net Primary Production' #"gC*m^{-2}*d^{-1}" 
        ,'t_1': 'Growth response time to moisture conditions'
        ,'t_2': 'Growth response time to temperature conditions'
        ,'t_3': 'Critical temperature for growth'
        ,'t_4': 'Germination fraction'
        ,'t_5': 'Allocation to reproduction'
        ,'t_6': 'Allocation to aboveground growth'
        ,'t_7': 'Allocation to belowground growth'
        ,'t_8': 'Allocation to storage'
        ,'t_9': 'Relative allocation to aboveground structure'
        ,'t_10': 'Relative allocation to belowground structure'
        ,'f_GERM': ''
        ,'gamma_GERM': 'Germination fraction'
        ,'k_GERM': ''
        ,'p': ''
        ,'GERM': '"Germination of carbon from C_S to C_A. Occurs when germination conditions are favourable (f_GERM = 1) and C_S > 0"'
        ,'f_SEED': ''
        ,'f_GROW': 'Growing conditions are controlled by environmental conditions, specifically, soil wetness and near-surface air temperature' # exprs: #See A3 in page 22
        ,'A_S': 'Allocation fraction to seeds'
        ,'A_L': 'Allocation fraction to leaves'
        ,'A_R': 'Allocation fraction to fine roots'
        ,'A_WL': 'Allocation fraction to aboveground wood'
        ,'A_WR': 'Allocation fraction to belowground wood'
        ,'tau_S': 'Seeds turnover rate' #days 
        ,'tau_L': 'Stem turnover rate' #days 
        ,'tau_R': 'Fine roots turnover rate' #days 
        ,'tau_WL': 'Aboveground wood turnover rate' #days 
        ,'tau_WR': 'Belowground wood turnover rate' #days 
        ,'C_RES_S': 'Growth respiration coefficient' #"gC*gC^-1"
        ,'C_RES_L': 'Growth respiration coefficient' #"gC*gC^-1"
        ,'C_RES_R': 'Growth respiration coefficient' #"gC*gC^-1"
        ,'C_RES_WL': 'Growth respiration coefficient' #"gC*gC^-1"
        ,'C_RES_WR': 'Growth respiration coefficient' #"gC*gC^-1"
}

for name in sym_dict.keys():
    var(name)

LAI = C_L*SLA
k =  0.5
f_VEG = 1 - exp(-k*LAI)
NPP = GPP - RES_a
gamma_GERM = 10**(4*(t_4**(-4))) #"days^{-1}"
GERM = f_GERM*gamma_GERM*(C_S/Max(p,k_GERM)) 
A_S = f_SEED*(t_5/(t_5+t_6+t_7+t_8))
A_L = f_GROW*(1-t_9)*(t_6/(t_5+t_6+t_7+t_8))
A_R = f_GROW*(1-t_10)*(t_7/(t_5+t_6+t_7+t_8))
A_WL = f_GROW*f_VEG*t_9*(t_6/(t_5+t_6+t_7+t_8))
A_WR = f_GROW*f_VEG*t_10*(t_7/(t_5+t_6+t_7+t_8)) #sum of A_S:A_WR < 1

x = StateVariableTuple((C_A, C_S, C_L, C_R, C_WL, C_WR))
u = NPP
b = (1,0,0,0,0,0)
Input = InputTuple(u*ImmutableMatrix(b))
A = CompartmentalMatrix(
[[-((A_S*(1-C_RES_S))+(A_L*(1-C_RES_L))+(A_R*(1-C_RES_R))+(A_WL*(1-C_RES_WL))+(A_WR*(1-C_RES_WR))),f_GERM*gamma_GERM*(1/Max(p,k_GERM)),0,0,0,0]
                               ,[(A_S*(1-C_RES_S)),-(f_GERM*gamma_GERM*(1/Max(p,k_GERM)))-(1/tau_S),0,0,0,0]
                               ,[(A_L*(1-C_RES_L)),0,-(1/tau_L),0,0,0]
                               ,[(A_R*(1-C_RES_R)),0,0,-(1/tau_R),0,0]
                               ,[(A_WL*(1-C_RES_WL)),0,0,0,-(1/tau_WL),0]
                               ,[(A_WR*(1-C_RES_WR)),0,0,0,0,-1/tau_WR]]
)
t = TimeSymbol("t") # "day"

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="JeDi-DGVM",
        longName="The Jena Diversity-Dynamic Global Vegetation Model", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="",
        doi="10.5194/bg-10-4137-2013",
        sym_dict=sym_dict
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_A, C_S, C_L, C_R, C_WL, C_WR)),
})
