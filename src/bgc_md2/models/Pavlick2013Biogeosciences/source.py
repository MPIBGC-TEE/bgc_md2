from sympy import var, ImmutableMatrix, exp, Max, Piecewise 
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
    NumericParameterization,
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
        ,'C_LIT': 'Fine litter carbon'
        ,'C_CWD': 'Woddy litter carbon'
        ,'C_SOIL': 'Soil carbon'
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
        ,'t_11': 'Turnover time of structural pools'
        ,'t_12': 'Turnover time of leaf and fine root pools'
        ,'t_13': 'Senescence response time to productivity conditions'
        ,'t_14': 'Relative senescence aboveground'
        ,'t_15': 'Plant nitrogen status'
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
        ,'tau_A': 'Stored assimilates turnover rate' #days 
        ,'tau_SEN': 'Time constant'
        ,'tau_NPP': 'Time constant'
        ,'f_SEN': 'Senescence'
        ,'f_NPP': 'Time-averaged NPP'
        ,'tau_L0': 'Base turnover time for leaf and fine root pools' #days 
        ,'tau_L': 'Leaves turnover rate' #days 
        ,'tau_S': 'Seeds turnover rate' #days 
        ,'tau_R': 'Fine roots turnover rate' #days 
        ,'tau_WL': 'Aboveground wood turnover rate' #days 
        ,'tau_WR': 'Belowground wood turnover rate' #days 
        ,'C_RES_S': 'Growth respiration coefficient' #"gC*gC^-1"
        ,'C_RES_L': 'Growth respiration coefficient' #"gC*gC^-1"
        ,'C_RES_R': 'Growth respiration coefficient' #"gC*gC^-1"
        ,'C_RES_WL': 'Growth respiration coefficient' #"gC*gC^-1"
        ,'C_RES_WR': 'Growth respiration coefficient' #"gC*gC^-1"
        ,'T': 'Temperature'
        ,'Q_10h': 'Sensitivity of heterotrophic respiration to air temperature'
        ,'c_lit_atm': 'Fraction of fine litter decomposition to atmosphere'
        ,'c_cwd_atm': 'Fraction of woody litter decomposition to atmosphere'
        ,'tau_LIT': 'Turnover time of fine litter at 20 °C'
        ,'tau_CWD': 'Turnover time of woody litter at 20 °C'
        ,'tau_SOIL': 'Turnover time of soil carbon at 20 °C'
        ,'LIT_L': 'Litter flux'
        ,'LIT_R': 'Litter flux'
        ,'LIT_A': 'Litter flux'
        ,'LIT_S': 'Litter flux'
        ,'LIT_WL': 'Litter flux'
        ,'LIT_WR': 'Litter flux'
        ,'DEC_LIT': 'Decomposition flux'
        ,'DEC_CWD': 'Decomposition flux'
        ,'DEC_SOIL': 'Decomposition flux'
        ,'t': 'Time symbol'
        ,'delta_t': 'time step' #fixme?
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
tau_WL = tau_WR
tau_WR = 365*((79*t_11)+1)
tau_L0 = (365/12)*10**(2*t_12)
tau_NPP = 10**((5*t_13)-2)
f_NPP = (NPP + tau_NPP * f_NPP * (t - delta_t))/(1+tau_NPP)
f_SEN = Piecewise((0,(f_NPP >= 0) | (NPP >= 0)),(1, (f_NPP < 0) & (NPP < 0)))
tau_L = ((1/tau_L0)+((f_SEN*t_14)/tau_SEN))**-1
tau_R = ((1/tau_L0)+((f_SEN*(1-t_14))/tau_SEN))**-1
LIT_L = C_L/tau_L # See equation A24 (page 25). This model includes plant biodiversity, so each grid cell has different traits (k). The LIT_ equations actually are community aggregated fluxes from the turnover of the various vegetation tissue pools according to the relative abundance of each trait
LIT_R = C_R/tau_R 
LIT_A = C_A/tau_A 
LIT_WL = C_WL/tau_WL 
LIT_WR = C_WR/tau_WR 
LIT_S = C_S/tau_S
DEC_LIT = (Q_10h**((T-20)/10)) * C_LIT / tau_LIT
DEC_CWD = (Q_10h**((T-20)/10)) * C_CWD / tau_CWD
DEC_SOIL = (Q_10h**((T-20)/10)) * C_SOIL / tau_SOIL

x = StateVariableTuple((C_A, C_S, C_L, C_R, C_WL, C_WR, C_LIT, C_CWD, C_SOIL))
u = NPP
b = (1,0,0,0,0,0)
Input = InputTuple((u*ImmutableMatrix(b),0,0,0))
A = CompartmentalMatrix(
[[-((A_S*(1-C_RES_S))+(A_L*(1-C_RES_L))+(A_R*(1-C_RES_R))+(A_WL*(1-C_RES_WL))+(A_WR*(1-C_RES_WR))),f_GERM*gamma_GERM*(1/Max(p,k_GERM)),0,0,0,0,0,0,0]
,[(A_S*(1-C_RES_S)),-(f_GERM*gamma_GERM*(1/Max(p,k_GERM)))-(1/tau_S),0,0,0,0,0,0,0]
,[(A_L*(1-C_RES_L)),0,-(1/tau_L),0,0,0,0,0,0]
,[(A_R*(1-C_RES_R)),0,0,-(1/tau_R),0,0,0,0,0]
,[(A_WL*(1-C_RES_WL)),0,0,0,-(1/tau_WL),0,0,0,0]
,[(A_WR*(1-C_RES_WR)),0,0,0,0,-1/tau_WR,0,0,0]
,[LIT_A/C_A, LIT_S/C_S, LIT_L/C_L, LIT_R/C_R,           0,           0,-DEC_LIT/C_LIT,             0, 0]
,[        0,         0,         0,         0, LIT_WL/C_WL, LIT_WR/C_WR,             0,-DEC_CWD/C_CWD, 0]
,[        0,         0,         0,         0,           0,           0,((1-c_lit_atm)*DEC_LIT)/C_LIT,((1-c_cwd_atm)*DEC_CWD)/C_CWD,-DEC_SOIL/C_SOIL]
])
t = TimeSymbol("t") # "day"

np1 = NumericParameterization(
    par_dict={
    c_lit_atm: 0.77
    ,c_cwd_atm: 0.2
    ,tau_LIT: 2.05 # "yr"
    ,tau_CWD: 60 # "yr"
    ,tau_SOIL: 100 # "yr"
    ,Q_10h: 1.4
},
    func_dict = frozendict({})
)

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
