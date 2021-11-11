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
from ComputabilityGraphs.CMTVS import CMTVS
from bgc_md2.helper import bgc_md2_computers

sym_dict={
        'x_1': 'Carbon in foliage',
        'x_2': 'Carbon in woody tissue',
        'x_3': 'Carbon in roots',
        'x_4': 'Carbon in metabolic litter',
        'x_5': 'Carbon in structural litter',
        'x_6': 'Carbon in fast soil organic matter (SOM)',
        'x_7': 'Carbon in slow soil organic matter (SOM)',
        'x_8': 'Carbon in passive soil organic matter (SOM)',
        'U': 'Photosynthetic rate (Carbon input) at time t ',
        'T':   'Temperature',
        'Q_10': 'Temperature quotient that describes a change in decomposition rate for evey 10°C difference in temperature',
        'W': 'Volumetric soil moisture',
        'f_W': 'Function of W',
        'f_T': 'Function of T',
        'b_1': 'Fixed partitioning ratio (fraction) of available carbon allocated to foliage',
        'b_2': 'Fixed partitioning ratio (fraction) of available carbon allocated to wood',
        'b_3': 'Fixed partitioning ratio (fraction) of available carbon allocated to roots',
        'c_1': 'Foliage exit rate',
        'c_2': 'Wood exit rate',
        'c_3': 'Roots exit rate',
        'c_4': 'exit rate',
        'c_5': 'exit rate',
        'c_6': 'exit rate',
        'c_7': 'exit rate',
        'c_8': 'exit rate',
        'f_41': 'transfer coefficients from pool j to pool i',
        'f_51': 'transfer coefficients from pool j to pool i',
        'f_52': 'transfer coefficients from pool j to pool i',
        'f_43': 'transfer coefficients from pool j to pool i',
        'f_53': 'transfer coefficients from pool j to pool i',
        'f_64': 'transfer coefficients from pool j to pool i',
        'f_65': 'transfer coefficients from pool j to pool i',
        'f_75': 'transfer coefficients from pool j to pool i',
        'f_76': 'transfer coefficients from pool j to pool i',
        'f_86': 'transfer coefficients from pool j to pool i',
        'f_67': 'transfer coefficients from pool j to pool i',
        'f_87': 'transfer coefficients from pool j to pool i',
        'f_68': 'transfer coefficients from pool j to pool i',
        'A': 'Carbon transfer coefficients between plant, litter, and soil pools',
        'B': 'Vector of partitioning coefficients of the photosynthetically fixed carbon (b_1, b_2, b_3)',
        'C': 'Diagonal matrix with exit rates of carbon from the eight carbon pools (c_)',
        'xi': 'Environmental scalar representing effects of temperature and moisture on the carbon transfer among pools' # xi_(t)
}
for name in sym_dict.keys():
    var(name)

f_W = Min((0.5*W),1)
f_T = Q_10**((T-10)/10)
xi = f_W*f_T
x = StateVariableTuple((x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8))
B = (b_1, b_2, b_3,0,0,0,0,0) 
Input = InputTuple(U * ImmutableMatrix(B))
A = ImmutableMatrix([
 [-1,0,0,0,0,0,0,0]
,[0,-1,0,0,0,0,0,0]
,[0,0,-1,0,0,0,0,0]
,[f_41,0,f_43,-1,0,0,0,0]
,[f_51,f_52,f_53,0,-1,0,0,0]
,[0,0,0,f_64,f_65,-1,f_67,f_68]
,[0,0,0,0,f_75,f_76,-1,0]
,[0,0,0,0,0,f_86,f_87,-1]
]) 
C = ImmutableMatrix(diag(c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8))
CM = CompartmentalMatrix(xi*A*C)
t = TimeSymbol("t")

# "Original parameters of the publication. Parameter value of GPP corresponds to an annual average"
# T, Q_10, and W are variables that should be looked at in a data set. What I have here are invented values
np1 = NumericParameterization(
    par_dict={
    Q_10: 1 
    ,W: 4.2 
    ,T: 25 
    ,U: 3370 #"gC*day^{-1}" # Average estimated value
    ,b_1: 0.14
    ,b_2: 0.14 
    ,b_3: 0.26
    ,c_1: 0.00258 
    ,c_2: 0.0000586 
    ,c_3: 0.00239
    ,c_4: 0.0109
    ,c_5: 0.00095
    ,c_6: 0.0105
    ,c_7: 0.0000995
    ,c_8: 0.0000115
    ,f_41: 0.9 
    ,f_51: 0.1
    ,f_52: 1
    ,f_43: 0.2
    ,f_53: 0.8
    ,f_64: 0.45
    ,f_65: 0.275
    ,f_75: 0.275
    ,f_76: 0.296
    ,f_86: 0.004
    ,f_67: 0.42
    ,f_87: 0.01
    ,f_68: 0.45
},
    func_dict=frozendict({})
    # state_var_units=gram/kilometer**2,
    # time_unit=day
)
nsv1 = NumericStartValueDict({
    x_1: 250 
    ,x_2: 4145
    ,x_3: 192
    ,x_4: 93
    ,x_5: 545
    ,x_6: 146
    ,x_7: 1585
    ,x_8: 300
})

ntimes = NumericSimulationTimes(np.arange(0, 150, 2.5))

mvs=CMTVS(
    {
        BibInfo(# Bibliographical Information
            name="",
            longName="", 
            version="",
            entryAuthor="Verónika Ceballos-Núñez",
            entryAuthorOrcid="0000-0002-0046-1160",
            entryCreationDate="24/3/2016",
            doi="",
    #        further_3eferences=BibInfo(doi=""),
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
        CM,  # the overall compartmental matrix
        Input,  # the overall input
        t,  # time for the complete system
        x,  # state vector of the complete system
        VegetationCarbonInputScalar(U),
        # vegetation carbon partitioning.
        VegetationCarbonInputPartitioningTuple((b_1, b_2, b_3)),
        VegetationCarbonStateVariableTuple((x_1, x_2, x_3)),
        np1,
        nsv1
    },
    bgc_md2_computers()
)
