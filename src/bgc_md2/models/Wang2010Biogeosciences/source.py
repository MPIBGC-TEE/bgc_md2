#import numpy as np
from sympy import var, symbols, Symbol, ImmutableMatrix, Min#, Rational
#from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
#    NumericParameterization,
#    NumericStartValueDict,
#    NumericSimulationTimes,
   )
from ..BibInfo import BibInfo 
from bgc_md2.helper import MVarSet

sym_dict={
        'C_leaf': 'Plant (carbon) pool Leaf' #unit: "gC*m^{-2}" 
        ,'C_root': 'Plant (carbon) pool Root' #unit: "gC*m^{-2}" 
        ,'C_wood': 'Plant (carbon) pool Wood' #unit: "gC*m^{-2}" 
        ,'C_j1': 'Amount of carbon in metabolic litter' #unit: "gC*m^{-2}" 
        ,'C_j2': 'Amount of carbon in structural litter' #unit: "gC*m^{-2}" 
        ,'C_j3': 'Amount of carbon in coarse woody debris (cwd) litter' #unit: "gC*m^{-2}" 
        ,'C_k1': 'Amount of carbon in microbial soil' #unit: "gC*m^{-2}" 
        ,'C_k2': 'Amount of carbon in slow soil' #unit: "gC*m^{-2}" 
        ,'C_k3': 'Amount of carbon in passive soil' #unit: "gC*m^{-2}" 
        ,'Delta_t': 'Time step of model integration' #FixMe: is it the same as TimeSymbol in this case?
        ,'N_min': 'Amount of mineral N in soil' #unit: "gN*m^{-2}"
        ,'P_lab': 'Amount of labile P in soil' #unit: "gP*m^{-2}"
        ,'F_nupmin': 'Minimum amount of N uptake required to sustain a given NPP'
        ,'F_pupmin': 'Minimum amount of P uptake required to sustain a given NPP'
        ,'x_nup': 'Nitrogen uptake limitation on NPP'
        ,'x_pup': 'Phosphorus uptake limitation on NPP'
        ,'x_npup': 'Nutrient uptake limiting factor'
        ,'n_leaf': 'N:C ratio of leaf biomass'
        ,'p_leaf': 'P:C ratio of leaf biomass'
        ,'k_n': 'Empirical constant'
        ,'k_p': 'Empirical constant'
        ,'x_nleaf':''
        ,'x_pleaf':''
        ,'x_npleaf': 'Nutrient concentration limiting factor'
        ,'F_cmax': 'Nutrient unlimited NPP' #unit: "gC*m^{-2}*d^{-1}"
        ,'F_c': 'Net Primary Productivity (flux)'
        ,'a_leaf': 'Fraction of NPP allocated to plant pool Leaf'
        ,'a_root': 'Fraction of NPP allocated to plant pool Root'
        ,'a_wood': 'Fraction of NPP allocated to plant pool Wood'
        ,'mu_leaf': 'Turnover rate of plant pool Leaf' #unit: "year^{-1}" # In table with parameter values. In Appendix B the unit is "d^{-1}" 
        ,'mu_root': 'Turnover rate of plant pool Root' #unit: "year^{-1}" # In table with parameter values. In Appendix B the unit is "d^{-1}" 
        ,'mu_wood': 'Turnover rate of plant pool Wood' #unit: "year^{-1}" # In table with parameter values. In Appendix B the unit is "d^{-1}" 
        ,'m_n': 'N limitation on litter decomposition, varies from 0 to 1' 
        ,'mu_j1': 'Metabolic litter turnover rate' #"d^{-1}"
        ,'mu_j2': 'Structural litter turnover rate' #"d^{-1}"
        ,'mu_j3': 'cwd litter turnover rate' #"d^{-1}"
        ,'mu_k1': 'Microbial soil turnover rate' #"d^{-1}"
        ,'mu_k2': 'Slow soil turnover rate' #"d^{-1}"
        ,'mu_k3': 'Passive soil turnover rate' #"d^{-1}"
        ,'b1l': 'Fraction of metabolic litter fall coming from leaves' 
        ,'b1r': 'Fraction of metabolic litter fall coming from root' 
        ,'b2l': 'Fraction of structural litter fall coming from leaves' 
        ,'b2r': 'Fraction of structural litter fall coming from root' 
        ,'b3w': 'Fraction of cwd litter fall coming from wood' 
        ,'c11': 'Fraction of microbial soil coming from metabolic litter' 
        ,'c12': 'Fraction of microbial soil coming from structural litter' 
        ,'c21': 'Fraction of slow soil fall coming from metabolic litter' 
        ,'c22': 'Fraction of slow soil coming from structural litter' 
        ,'c23': 'Fraction of slow soil coming from cwd litter' 
        ,'c33': 'Fraction of passive soil coming from cwd litter' 
        ,'d12': 'Fraction of microbial soil coming from slow soil' 
        ,'d13': 'Fraction of microbial soil coming from passive soil' 
        ,'d21': 'Fraction of slow soil coming from microbial soil' 
        ,'d23': 'Fraction of slow soil coming from passive soil' 
        ,'d31': 'Fraction of passive soil coming from microbial soil' 
        ,'d32': 'Fraction of passive soil coming from slow soil' 
}
for name in sym_dict.keys():
    var(name)
#a_leaf + a_root + a_wood = 1
#b1l + b2l = 1
#b1r + b2r = 1
x_nup = Min(1,(N_min/(F_nupmin*Delta_t)))
x_pup = Min(1,(P_lab/(F_pupmin*Delta_t)))
x_npup = Min(x_nup,x_pup)
x_nleaf = (n_leaf/(n_leaf+k_n))
x_pleaf = (p_leaf/(p_leaf+k_p))
x_npleaf = Min(x_nleaf,x_pleaf)
F_c = x_npleaf*x_npup*F_cmax #unit: "gC*m^{-2}*d^{-1}" 
u = F_c
x = StateVariableTuple((C_leaf, C_root, C_wood, C_j1, C_j2, C_j3, C_k1, C_k2, C_k3))
b = (a_leaf, a_root, a_wood)
Input = InputTuple((u * ImmutableMatrix(b),0,0,0,0,0,0))
A = CompartmentalMatrix([
 [-mu_leaf*(b1l+b2l),         0        ,      0     ,          0         ,          0         ,           0        ,        0       ,        0       ,        0       ]
,[         0        ,-mu_root*(b1r+b2r),      0     ,          0         ,          0         ,           0        ,        0       ,        0       ,        0       ]
,[         0        ,         0        ,-mu_wood*b3w,          0         ,          0         ,           0        ,        0       ,        0       ,        0       ]
,[    mu_leaf*b1l   ,    mu_root*b1r   ,      0     ,-mu_j1*m_n*(c11+c21),          0         ,           0        ,        0       ,        0       ,        0       ]
,[    mu_leaf*b2l   ,    mu_root*b2r   ,      0     ,          0         ,-mu_j2*m_n*(c12+c22),           0        ,        0       ,        0       ,        0       ]
,[         0        ,         0        , mu_wood*b3w,          0         ,          0         ,-mu_j3*m_n*(c23+c33),        0       ,        0       ,        0       ]
,[         0        ,         0        ,      0     ,    mu_j1*m_n*c11   ,    mu_j2*m_n*c12   ,           0        ,-mu_k1*(d21+d31),    mu_k2*d12   ,    mu_k3*d13   ]
,[         0        ,         0        ,      0     ,    mu_j1*m_n*c21   ,    mu_j2*m_n*c22   ,    mu_j3*m_n*c23   ,    mu_k1*d21   ,-mu_k2*(d12+d32),    mu_k3*d23   ]
,[         0        ,         0        ,      0     ,          0         ,          0         ,    mu_j3*m_n*c33   ,    mu_k1*d31   ,    mu_k2*d32   ,-mu_k3*(d13+d23)]
])
### When not explicit in equations B1, B2, B3 (Appendix B), followed the fluxes from Fig 1. However, the soil fluxes are not clear.
t = TimeSymbol("t") #'day' # or 'year'? incongruent turnover units

# Commented out the following lines because original publication only has 3 parameter values
# In original publication: 
## See table 1 for parameter values; a_(leaf,wood,root) and 1/mu_(leaf,wood,root) are based on CASA
## 1/mu_(leaf,wood,root) = Mean residence time of plant tissue
## 1/n_max,leaf  and  1/p_max,leaf : maximal leaf N:C and P:C ratios = 1.2*mean (min is obtained by 0.8*mean)  
#    # mean estimates were obtained from Glopnet datasets for each biome (Wright et al., 2004)
#np1 = NumericParameterization(
#    par_dict={#"Evergreen needle leaf forest":
#    Delta_t: 1, 
#    k_n: 0.01, #"gN*(gC)^{-1}"
#    k_p: 0.0006, #"gP*(gC)^{-1}"
#    a_leaf: 0.42, 
#    a_root: 0.25, 
#    a_wood: 0.33, 
#    mu_leaf: 'Rational(1,2)', 
#    mu_root: 'Rational(1,18)', 
#    mu_wood: 'Rational(1,70)', 
#    n_leaf: 'Rational(1,42)',#"gN/gC"  
#    p_leaf: 'Rational(1,408)'#"gP/gC" 
#},
#    par_dict={#"Evergreen broadleaf forest":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.25, a_root: 0.65, a_wood: 0.1, mu_leaf: 'Rational(1,1.5)', mu_root: 'Rational(1,10)', mu_wood: 'Rational(1,60)', n_leaf: 'Rational(1,21)', p_leaf: 'Rational(1,400)'
#},
#    par_dict={#"Deciduous needle leaf forest":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.4, a_root: 0.3, a_wood: 0.3, mu_leaf: 'Rational(1,0.8)', mu_root: 'Rational(1,10)', mu_wood: 'Rational(1,80)', n_leaf: 'Rational(1,50)', p_leaf: 'Rational(1,405)'
#},
#    par_dict={#"Deciduous broadleaf forest":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.3, a_root: 0.5, a_wood: 0.2, mu_leaf: 'Rational(1,8)', mu_root: 'Rational(1,10)', mu_wood: 'Rational(1,40)', n_leaf: 'Rational(1,21)', p_leaf: 'Rational(1,333)'
#},
#    par_dict={#"Mixed forest":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.35, a_root: 0.25, a_wood: 0.4, mu_leaf: 'Rational(1,1.2)', mu_root: 'Rational(1,10)', mu_wood: 'Rational(1,50)', n_leaf: 'Rational(1,28)', p_leaf: 'Rational(1,278)'
#},
#    par_dict={#"Shrub land (open and close shrubland)":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.4, a_root: 0.45, a_wood: 0.15, mu_leaf: 'Rational(1,1.2)', mu_root: 'Rational(1,5)', mu_wood: 'Rational(1,40)', n_leaf: 'Rational(1,33)', p_leaf: 'Rational(1,293)'
#},
#    par_dict={#"Woddy savannah":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.3, a_root: 0.6, a_wood: 0.1, mu_leaf: 'Rational(1,1.5)', mu_root: 'Rational(1,5)', mu_wood: 'Rational(1,40)', n_leaf: 'Rational(1,21)', p_leaf: 'Rational(1,354)'
#},
#    par_dict={#"Savannah":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.2, a_root: 0.7, a_wood: 0.1, mu_leaf: 'Rational(1,1.5)', mu_root: 'Rational(1,3)', mu_wood: 'Rational(1,40)', n_leaf: 'Rational(1,21)', p_leaf: 'Rational(1,492)'
#},
#    par_dict={#"Grassland":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.3, a_root: 0.7, a_wood: 0, mu_leaf: 1, mu_root: 'Rational(1,3)', mu_wood: 1, n_leaf: 'Rational(1,42)', p_leaf: 'Rational(1,833)'
#},
#    par_dict={#"Crop land (cropland mosaic was aggregated into this term)":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.3, a_root: 0.7, a_wood: 0, mu_leaf: 1, mu_root: 'Rational(1,0.9)', mu_wood: 1, n_leaf: 'Rational(1,21)', p_leaf: 'Rational(1,333)'
#},
#    par_dict={#"Barren or sparse vegetation":
#    Delta_t: 1, k_n: 0.01, k_p: 0.0006, a_leaf: 0.2, a_root: 0.6, a_wood: 0.2, mu_leaf: 1, mu_root: 'Rational(1,4)', mu_wood: 'Rational(1,5)', n_leaf: 'Rational(1,17)', p_leaf: 'Rational(1,167)'
#},
#    func_dict=frozendict({})
#    # state_var_units=gram/meter**2,
#    # time_unit=day
#)
#nsv1 = NumericStartValueDict({
#    C_leaf: , #"g*m^{-2}"
#    C_root: , #"g*m^{-2}"
#    C_wood: #"g*m^{-2}"
#})
#
#ntimes = NumericSimulationTimes(np.arange(, , ))

mvs=MVarSet({
    BibInfo(# Bibliographical Information
        name="CABLE",
        longName="CSIRO Atmosphere Biosphere Land Exchange", 
        version="1",
#        basedOn = "CASA'" # Fung et al. (2005) -> Need to look it up
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="14/3/2016",
        doi="10.5194/bg-7-2261-2010",
        sym_dict=sym_dict
        
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_leaf, C_root, C_wood)),
#    np1
})

