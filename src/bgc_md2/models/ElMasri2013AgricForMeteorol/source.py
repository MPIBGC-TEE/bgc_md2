from sympy import var, ImmutableMatrix, Piecewise, Max, exp, sympify
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
        'C_leaf': 'Amount of carbon for the leaf' #"kgC*m^{-2}" 
        ,'C_stem': 'Amount of carbon for the stem' #"kgC*m^{-2}" 
        ,'C_roots': 'Amount of carbon for the root' #"kgC*m^{-2}" 
        ,'GPP': 'Carbon gain via photosynthesis (Gross Primary Productivity, GPP)' #"KgC*m^{−2}*yr^{−1}"
        # since NPP became an expression it is not part of the symdict ### Q(Vero) But what if we want to see the description?
        ,'NPP': 'Net primary Productivity (NPP)' #"KgC*m^{−2}*yr^{−1}"  #Previous equation was exactly like the Arora... because the paper cited it as a reference in this part, would be interesting to compare both models once we have the "keys" (Vero) 
        ,'R_G': 'Growth respiration'
        ,'R_mtotal': 'Total maintenance respiration'
        ,'R_auto': 'Autotrophic respiration'
        ,'k_leaf': 'Base respiration rate for leaves'
        ,'cn_leaf': 'Carbon:Nitrogen ratio for the leaves'
        ,'k_stem': 'Base respiration rate for stem'
        ,'cn_stem': 'Carbon:Nitrogen ratio for the stem'
        ,'k_roots': 'Base respiration rate for roots'
        ,'cn_roots': 'Carbon:Nitrogen ratio for the roots'
        ,'gt_ls': 'Function of Q_10 and temperature for leaves and stems' #See table A1 for equations
        ,'gt_r': 'Function of Q_10 and temperature for roots' #See table A1 for equations
        ,'teta': 'Leaf phenological status defined as the current fraction of maximum LAI'
        ,'LAI': 'Leaf Area Index'
        ,'Q_10_ls': ''
        ,'Q_10_r': ''
        ,'T_veg': 'Vegetation temperature ºC'
        ,'T_soil': 'Soil temperature ºC'
        ,'T_r': 'Root temperature ºC'
        ,'T_air': 'Air temperature ºC'
        ,'T_cold': 'Temperature threshold for leaf loss because of cold stress ºC'
        ,'R_leaf': 'Leaf respiration'
        ,'R_stem': 'Stem respiration'
        ,'R_roots': 'Roots respiration'
        ,'Allo_fact_stem': 'Allocation factor for stem'
        ,'Allo_fact_roots': 'Allocation factor for roots'
        ,'Allo_fact_leaf': 'Allocation factor for leaves (herbaceous)'
        ,'Allo_fact_roots': 'Allocation factor for roots (herbaceous)'
        ,'Allo_fact_leaf': 'Allocation factor for leaves'
        ,'epsilon_leaf': 'Parameter controlling allocation to leaves'
        ,'epsilon_stem': 'Parameter controlling allocation to stem'
        ,'epsilon_root': 'Parameter controlling allocation to roots'
        ,'omega': 'Allocation parameter'
        ,'ws': 'Water stress factor'
        ,'k_n': 'Light extinction coefficient'
        ,'L': 'Light availability factor'
        ,'L_grass': 'Light availability factor for grasses'
        #,'a_L': 'Parameter introduced by the author of this entry in order to summarize equations on the paper.'
        #,'a_S': 'Parameter introduced by the author of this entry in order to summarize equations on the paper.'  
        #,'a_R': 'Parameter introduced by the author of this entry in order to summarize equations on the paper.'  
        ,'r_wmax': 'Maximum drought leaf loss rate'
        ,'r_tmax': 'Maximum cold leaf loss rate'
        ,'b_w': 'Parameter for drought leaf loss'
        ,'b_t': 'Parameter for cold leaf loss'
        ,'Y_leaf': 'Leaf life span' #"year"
        ,'Y_stem': 'Stem turnover rate' #"year"
        ,'Y_roots': 'Fine root turnover rate' #"year"
}

for name in sym_dict.keys():
    var(name)

x = StateVariableTuple((C_leaf, C_stem, C_roots))
Q_10_ls = 3.22 - (0.046*T_veg)
Q_10_r = 3.22 - (0.046*T_soil)
gt_ls = Q_10_ls**((T_veg-20)/10)
gt_r = Q_10_r**((T_soil-20)/10)
R_leaf = k_leaf * (C_leaf/cn_leaf) * teta * gt_ls
R_stem = k_stem * (C_stem/cn_stem) * gt_ls
R_roots = k_roots * (C_roots/cn_roots) * teta * gt_r
R_mtotal = R_leaf + R_stem + R_roots
R_G = 0.25 * Max(0,GPP-R_mtotal) 
R_auto = R_mtotal + R_G
NPP = GPP-R_mtotal-R_G
L = exp(-k_n*LAI)
L_grass = Max(0,(1-LAI)/4.5)
Allo_fact_stem = (epsilon_stem + (omega*(1-L)))/(1+(omega*(2-L-ws))) 
Allo_fact_roots = (epsilon_root + (omega*(1-ws)))/(1+(omega*(2-L-ws))) 
Allo_fact_leaf = 1 - Allo_fact_stem - Allo_fact_roots
Allo_fact_roots_herb = (epsilon_root + (omega*(1-ws)))/(1+(omega*(1+L-ws)))
Allo_fact_leaf_herb = (epsilon_leaf + (omega*L))/(1+(omega*(1+L-ws)))
#beta_t = Piecewise((1,T_air > T_cold),((T_air -(T_cold - 5))/5,T_cold > T_air > (T_cold - 5)))
beta_t = Piecewise((1,T_air > T_cold),(0,T_air <= (T_cold-5)))
r_n = 1/(Y_leaf*365)
r_w = r_wmax*(1-ws)**b_w
r_t = r_tmax*(1-beta_t)**b_t
a_L = Piecewise((((GPP*Allo_fact_leaf)-R_leaf),NPP<0),((NPP*Allo_fact_leaf),NPP>=0))
a_S = Piecewise((((GPP*Allo_fact_stem)-R_stem),NPP<0),((NPP*Allo_fact_stem),NPP>=0))
a_R = Piecewise((((GPP*Allo_fact_roots)-R_roots),NPP<0),((NPP*Allo_fact_roots),NPP>=0))
# a_L = Piecewise((((GPP*Allo_fact_leaf)-R_leaf),NPP<0),((NPP*Allo_fact_leaf),NPP>=0))
# a_S = Piecewise((((GPP*Allo_fact_stem)-R_stem),NPP<0),((NPP*Allo_fact_stem),NPP>=0))
# a_R = Piecewise((((GPP*Allo_fact_roots)-R_roots),NPP<0),((NPP*Allo_fact_roots),NPP>=0))
# 
# u = VegetationCarbonInputPartitioningTuple((a_L,a_S,a_R))
# 
# b = VegetationCarbonInputScalar(sympify(1))
# b looks very funny to me. Should it not be something more along the lines
# of sum(Iv) where Iv is the input to the vegetation pools (VegetationCarbonInputTuple)
# my take would actually be:
#u = VegetationCarbonInputScalar(
#    Piecewise((GPP,NPP<0),(NPP,NPP>=0))
#)
#b = VegetationCarbonInputPartitioningTuple(
#    (
#        Piecewise((Allo_fact_leaf-R_leaf/GPP,NPP<0),(Allo_fact_leaf,NPP>=0)),
#        Piecewise((Allo_fact_stem-R_stem/GPP,NPP<0),(Allo_fact_stem,NPP>=0)), 
#        Piecewise((Allo_fact_roots-R_roots/GPP,NPP<0),(Allo_fact_roots,NPP>=0))
#    )
#)
# framework computes
#Input = InputTuple(b*u)
Input = InputTuple(
    (
        a_L,
        a_S,
        a_R
    )
)

vcsvt = VegetationCarbonStateVariableTuple((C_leaf, C_stem, C_roots))
# from Input and vcsvt we could compute
#u = a_L + a_S + a_R
#b= VegetationCarbonInputPartitioningTuple(
#    (
#        a_L/u,
#        a_S/u,
#        a_R/u
#    )
#)
#
#
#b=VegetationCarbonInputScalar(
#    Piecewise((GPP,NPP<0),(NPP,NPP>=0))
#)
#u = VegetationCarbonInputPartitioningTuple(
#    (
#        Piecewise((Allo_fact_leaf-R_leaf/GPP,NPP<0),(Allo_fact_leaf,NPP>=0)),
#        Piecewise((Allo_fact_stem-R_stem/GPP,NPP<0),(Allo_fact_stem,NPP>=0)), 
#        Piecewise((Allo_fact_roots-R_roots/GPP,NPP<0),(Allo_fact_roots,NPP>=0))
#    )
#)
# 
#Input = InputTuple(b*u)
#"f_v = u + A*x"
A = CompartmentalMatrix(
[[-(r_n+r_w+r_t),0,0]
,[0, -(1/(Y_stem*365)), 0]
,[0, 0, -(1/(Y_roots*365))]]
)
t = TimeSymbol("t")

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="ISAM",
        longName="Integrated Science Assessment Model", 
        version="",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="3/5/2018",
        doi="10.1111/j.1365-2486.2004.00890.x",
        sym_dict=sym_dict
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    # vegetation carbon partitioning.
    #b,
    #u,
    vcsvt
})
