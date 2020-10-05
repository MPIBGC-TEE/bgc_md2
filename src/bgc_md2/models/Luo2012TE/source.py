import numpy as np
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
    NumericStartValueDict,
    NumericSimulationTimes,
   )
from ..BibInfo import BibInfo 
from bgc_md2.helper import MVarSet

sym_dict={
        'C_f': 'Carbon in foliage',
        'C_r': 'Carbon in roots',
        'C_w': 'Carbon in woody tissue',
        'GPP': 'Photosynthetic rate (Carbon input) at time t ',
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
t = TimeSymbol("t")

# "Original parameters of the publication. Parameter value of GPP corresponds to an annual average"
# T, Q_10, and W are variables that should be looked at in a data set. What I have here are invented values
np1 = NumericParameterization(
    par_dict={
    Q_10: 1, 
    W: 4.2, 
    T: 25, 
    GPP: 3370, #"gC*day^{-1}"
    eta_f: 0.14,
    eta_r: 0.26,
    eta_w: 0.14, 
    gamma_f: 0.00258, 
    gamma_w: 0.0000586, 
    gamma_r: 0.00239},
    func_dict=frozendict({})
    # state_var_units=gram/kilometer**2,
    # time_unit=day
)
nsv1 = NumericStartValueDict({
    C_f: 250, 
    C_w: 4145, 
    C_r: 192
})

ntimes = NumericSimulationTimes(np.arange(0, 150, 2.5))

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
