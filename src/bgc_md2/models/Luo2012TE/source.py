from sympy import var, symbols, Symbol, ImmutableMatrix, diag, Min
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
from bgc_md2.helper import MVarSet

sym_dict={
        'C_f': 'Carbon in foliage',
        'C_r': 'Carbon in roots',
        'C_w': 'Carbon in woody tissue',
        'GPP': 'Photosynthetic rate (Carbon input) at time t ',#unit: "gC*day^{-1}"
        'T':   'Temperature',
        'Q_10': 'Temperature quotient that describes a change in decomposition rate for evey 10°C difference in temperature',
        'W': 'Volumetric soil moisture',
        'f_W': 'Function of W',
        'f_T': 'Function of T',
        'epsilon_t': 'Environmental scalar #unit: km^2',
        'eta_f': 'Fixed partitioning ratio (fraction) of available carbon allocated to foliage',
        'eta_r': 'Fixed partitioning ratio (fraction) of available carbon allocated to roots',
        'eta_w': 'Fixed partitioning ratio (fraction) of available carbon allocated to wood',
        'gamma_f': 'Foliage turnover rate',
        'gamma_r': 'Roots turnover rate',
        'gamma_w': 'Wood turnover rate'
}
for name in sym_dict.keys():
    var(name)

f_W = Min((0.5*W),1)
f_T = Q_10**((T-10)/10)
epsilon_t = f_W*f_T
x = StateVariableTuple((C_f, C_r, C_w ))
u = GPP*epsilon_t
b = (eta_f, eta_w, eta_r)

Input = InputTuple(u * ImmutableMatrix(b))
A = CompartmentalMatrix(
    diag(-gamma_f, -gamma_w, -gamma_r)
)
t = TimeSymbol("t") #unit: "day"

np1 = NumericParameterization(
    par_dict={
    },
    func_dict=frozendict({})
)

mvs=MVarSet({
    BibInfo(# Bibliographical Information
        name="",
        longName="", 
        version="",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="24/3/2016",
#        doi="",
#        further_references=BibInfo(doi=""),
#        modApproach: process based,
#        partitioningScheme: fixed,
#        claimedDynamicPart: "no",
#        spaceScale: global, 
#        #    unit: "1°",
#        timeResolution: monthly,
#        #    unit: month^{-1},
#        bibtex: "@article{Luo2012TE,
#                 address = {Berkeley},
#                 author = {Yiqi Luo and Ensheng Weng and Yuanhe Yang},
#                 booktitle = {Encyclopedia of Theoretical Ecology},
#                 editor = {Alan Hastings and Louis Gross},
#                 pages = {219-229},
#                 publisher = {University of California Press},
#                 title = {Ecosystem Ecology},
#                 year = {2012}
#                }",
#        
#        abstract: "Ecosystem ecology is a subdiscipline of ecology that focuses on exchange of energy and materials between organisms and the environment. The materials that are commonly studied in ecosystem ecology include water, carbon, nitrogen, phosphorus, and other elements that organisms use as nutrients. The source of energy for most ecosystems is solar radiation. In this entry, material cy-cling and energy exchange are generally described before the carbon cycle is used as an example to illustrate our quantitative and theoretical understanding of ecosystem ecology.",
        sym_dict=sym_dict
        
    ),
    A,  # the overall compartmental matrix
    Input,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_f, C_w, C_r)),
    np1
})
