from sympy import var, ImmutableMatrix, Piecewise, exp, Min, Max
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
        'C_L': 'Amount of carbon for the leaf' # "kgC*m^{-2}" 
        ,'C_S': 'Amount of carbon for the stem' # "kgC*m^{-2}" 
        ,'C_R': 'Amount of carbon for the root' # "kgC*m^{-2}" 
        ,'R_gL': 'Growth respiration flux for the leaves'
        ,'R_mL': 'Maintenance respiration flux for the leaves'
        ,'R_gS': 'Growth respiration flux for the stem'
        ,'R_mS': 'Maintenance respiration flux for the stem'
        ,'R_gR': 'Growth respiration flux for the root'
        ,'R_mR': 'Maintenance respiration flux for the root'
        ,'R_hD': 'Heterotrophic respiration from litter (debris)'
        ,'R_hH': 'Heterotrophic respiration from soil carbon (humus)'
        ,'G': 'Carbon gain via photosynthesis (Gross Primary Productivity, GPP)'
        ,'N': 'Net primary Productivity (NPP)'
        ,'LAI': 'Leaf Area Index'
        ,'k_n': 'PFT-dependent light extinction coefficient'
        ,'L': 'Light availability (scalar index between 0 and 1)'
        ,'theta_i': 'Volumetric soil moisture content'
        ,'theta_field': 'Field capacity'
        ,'theta_wilt': 'Wilting point'
        ,'W_i': 'Availability of water in soil layer i. Weighted by the fraction of roots present in each soil layer'
        ,'W': 'Averaged soil water availability index'
        ,'epsilon_L': 'PFT-dependent parameter for leaf'
        ,'epsilon_S': 'PFT-dependent parameter for stem'
        ,'epsilon_R': 'PFT-dependent parameter for root'
        ,'omega': 'PFT-dependent parameter'
        ,'a_S': 'Stem allocation fraction'
        ,'a_R': 'Root allocation fration'
        ,'a_L': 'Leaf allocation fraction'
        ,'A_S': 'Amount of carbon allocated to the stem'
        ,'A_R': 'Amount of carbon allocated to the root '
        ,'T_air': 'Temperature of the air' # "°C"
        ,'T_cold': 'Cold temperature threshold for a PFT below which leaf loss begins to occur' # "°C"
        ,'b_T': '"Parameter that describes sensitivity of leaf loss to temp. below the T$_{cold}$"'
        ,'beta_T': 'Temperature measure (varies between 0 and 1)'
        ,'gamma_N': 'Loss rate (normal turnover)' # day^{-1} 
        ,'gamma_W': 'Loss rate under drought stress' # day^{-1} 
        ,'gamma_Tmax': 'Maximum loss rate of specified PFT'
        ,'gamma_T': 'Loss rate under cold stress' # day^{-1} 
        ,'gamma_S': 'Stem turnover rate'  # year^{-1} 
        ,'gamma_R': 'Root turnover rate' # year^{-1} 
        ,'D_L': 'Litter loss from the leaves'
        ,'D_S': 'Litter loss from the stem'
        ,'D_R': 'Litter loss from the root'
}

for name in sym_dict.keys():
    var(name)

N=G-(R_gL+R_gS+R_gR)-(R_mL+R_mS+R_mR)
L=(exp(-k_n*LAI))
#L=(Piecewise(((exp(-k_n*LAI)),FOR TREES AND CROPS),((Max(0,1-(LAI/4.5))),FOR GRASSES)))
W_i=(Max(0,(Min(1,((theta_i-theta_wilt)/(theta_field-theta_wilt))))))
epsilon_R=1-epsilon_L-epsilon_S
a_S=((epsilon_S+(omega*(1-L)))/(1+(omega*(2-L-W))))
a_R=((epsilon_R+(omega*(1-W)))/(1+(omega*(2-L-W))))
#a_R=((epsilon_R+(omega*(1-W)))/(1+(omega*(1-L-W)))) # For grasses
a_L=1-a_S-a_R 
#a_L=((epsilon_L+(omega*L))/(1+(omega*(1-L-W)))) # For grasses
A_S=(Piecewise((a_S*G,N<0),(a_S*N+R_gS+R_mS,N>0)))
A_R=(Piecewise((a_R*G,N<0),(a_R*N+R_gR+R_mR,N>=0)))
#beta_T=(Piecewise((1,T_air>= T_cold),(Piecewise((((T_air-T_cold-5)/5),T_air>(T_cold-5)),(0,T_air<=(T_cold-5))),T_cold>T_air))) # Had to comment it out because of problems of piecewise in matrix
gamma_T=gamma_Tmax*(1-beta_T)**b_T
D_L=(gamma_N+gamma_W+gamma_T)*C_L
D_S=gamma_S*C_S
D_R=gamma_R*C_R

x = StateVariableTuple((C_L, C_S, C_R))
b = 1
u = (C_L, C_S, C_R)
Input = InputTuple(ImmutableMatrix(u)*1) #Fixme: does input always have to be a Tuple? what happens when, in this case, "f_v = u + A*x"? 
A = CompartmentalMatrix([[-(gamma_N+gamma_W+gamma_T),0,0],
                               [0, -gamma_S-R_gS-R_mS, 0],
                               [0, 0, -gamma_R-R_gR-R_mR]]
)
# 
t = TimeSymbol("t")


#    parameter_sets:
#        - "Original dataset of the publication":
#            values: {k_n: 0.5, omega: 0.8, epsilon_L: 0.35, epsilon_S: 0.1, epsilon_R: 0.55}
#            desc: Eastern US and Germany, cold broadleaf deciduous
#            doi: 10.1111/j.1365-2486.2004.00890.x

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="CTEM",
        longName="Canadian Terrestrial Ecosystem Model", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="21/1/2016",
        doi="10.1111/j.1365-2486.2004.00890.x",
        sym_dict=sym_dict
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(b),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(u),
    VegetationCarbonStateVariableTuple((C_L, C_S, C_R)),
})
