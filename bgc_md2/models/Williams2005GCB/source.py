import numpy as np
from sympy import var, symbols, Symbol, ImmutableMatrix, exp, Rational
from sympy.physics.units.systems import SI
from sympy.physics.units import (
    Quantity,
    length,
    mass,
    time,
    temperature,
    kilogram,
    meter,
    day,
    kelvin,
)
from sympy.physics.units.systems.si import dimsys_SI

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
    NumericParameterizedSmoothReservoirModel,
    NumericSimulationTimes,
    QuantityStartValueDict,
    QuantitySimulationTimes,
    QuantityParameterization,
    QuantityModelRun,
    QuantityParameterizedSmoothReservoirModel,
)
from bgc_md2.resolve.computers import (
    smooth_reservoir_model_from_input_tuple_and_matrix,
    quantity_parameterization_1,
    numeric_model_run_1,
    numeric_start_value_array_1,
    quantity_model_run_1,
    quantity_start_value_array_1,
)
from ..BibInfo import BibInfo 



t = TimeSymbol("t")

# It is better to store the purely symbolical description (without units or
# even dimensions)
# since it is more versatile. For purposes of comparison we might want to reinterprete
# something that was given originally as mass/aread as mass or vice versa
# The responsibility of the user is to provide consisten parameterisations
# for which the comments about dimensions are helpful and must be correct
# but will not be checked automatically by the system.
# 
# We can retain make the original description available to the framework as 
# part of the bibliographical information in a dictionary which could also be used to define
# the symbols (to avoid duplication) as demonstrated here.

sym_dict={
    "C_f":      "Foliar C mass",
    "C_lab":    "Labile C mass",
    "C_w":	    "Wood C mass",
    "C_r":	    "Fine root C mass",
    "p_10":     "Parameter in exponential term of temperature dependent rate parameter",
    "mint":	    "Dayly minimum temperature",
    "maxt":	    "Dayly maximum temperature",
    "multtl":	"Turnover of labile C (0 = off)		 1 = On)			",
    "multtf":	"Turnover of foliage C (0 = off)		 1 = On)			",
    "LAI":	    "LeafAreaIndex",
    "p_3":	    "Fraction of NPP partitioned to foliage",
    "p_4":	    "Fraction of NPP partitioned to roots",
    "p_5":	    "Turnover rate of foliage",
    "p_6":	    "Turnover rate of wood",
    "p_7":	    "Turnover rate of roots",
    "p_14":	    "Fraction of leaf loss transferred to litter",
    "p_15":	    "Turnover rate of labile carbon",
    "p_16":	    "Fraction of labile transfers respired ",
    "NPP":	    "Net Primary Production per area",
}
# For the code we only use the first column
for name in sym_dict.keys():
    var(name)

# The following parameters would be needed to compute NPP from GPP  since they are not given
#    ("p_2"		, "1"	                            ,"Fraction of GPP respired "), # no parameter values in publication
#    ("GPP"		, "mass*length**(-2)*time**(-1)"	,"Gross Primary Production per area"),
#    ("G"		, "mass*length**(-2)*time**(-1)"	                            ,),

# Temperature sensitive rate parameterd
T_rate = 0.5 * exp(p_10 * ((maxt + mint) / 2))

# state vector
x = StateVariableTuple((C_f, C_lab, C_w, C_r))

# NetPrimaryProduction can not be computted in the absence of G and p_2 with
# the following formula NPP = (1-p_2)*G # CAUTION: there are 2 NPPs, need to
# correct according to ACM model so NPP has been made a parameter instead of an
# expression

# scalar function of photosynthetic inputs
u = NPP

# tuple of partitioning coefficients of photosynthetically fixed carbon
b = ((p_3 * multtl), 0, (1 - p_4), p_4)

# f_v = u*b+A*x

# matrix of cycling rates
A = CompartmentalMatrix(
    [
        [
            (
                (-1)
                * (
                    ((1 - p_14) * p_5 * (1 - p_16) * multtf * T_rate)
                    + ((1 - p_14) * p_5 * p_16 * multtf * T_rate)
                    + (p_5 * p_14 * multtf)
                )
            ),
            (p_15 * (1 - p_16) * multtl * T_rate),
            0,
            0,
        ],
        [
            ((1 - p_14) * p_5 * (1 - p_16) * multtf * T_rate),
            (
                (-p_15 * (1 - p_16) * multtl * T_rate)
                - (p_15 * p_16 * multtl * T_rate)
            ),
            0,
            0,
        ],
        [0, 0, -p_6, 0],
        [0, 0, 0, -p_7],
    ]
)
Input = InputTuple(u * ImmutableMatrix(b))

# Here create the instance explicitly by calling a function from the computers
# sub module that is also available to the framework to compute srm 
# automatically.
# It is not necessarry to create the model instance explicitly, since the framework recognizes the ingredients. 
# But since the models serve as examples for 
srm = smooth_reservoir_model_from_input_tuple_and_matrix(Input, A, t, x)

np1 = NumericParameterization(
    par_dict={
        NPP: Rational(409, 365),  # *gram*/(meter**2*day)
        # Note:
        # The (negative ) parameter value for mint (see below) suggests that
        # maxt and mint are given in the celsius scale.
        # The sympy unit system does not destinguish between absolute
        # temperature and celsius scale.
        # Conceptually 5 deg Celsius describe a temperature DIFFERENCE of 5 Kelvin
        # to the triple point of water.
        # Differences are always mesured in Kelvin
        # So although mint and maxt are given in Kelvin they are not understood
        # as absolute temperatures but as differences to the triple point.
        mint: -4,  # Kelvin (deg Celsius)
        maxt: 5,  # Kelvin (de Celsius)
        multtf: 0,
        multtl: 1,
        # p_2		    :0.47,
        p_3: 0.31,
        p_4: 0.43,
        p_5: 0.0027,  # 1/day
        p_6: 0.00000206,  # 1/day
        p_7: 0.00248,  # 1/day
        #
        # Althouhg p10 is considered to have unit 1/Kelvin
        # inspection of the T_rate expression connects it
        # to Temperature DIFFERENCE from  the triple point
        p_10: 0.0693,  # /kelvin
        p_14: 0.45,
        p_15: 0.001,  # 1/day
        p_16: 0.251,
    },
    func_dict=frozendict({})
    # state_var_units=kilogram/meter**2,
    # time_unit=day
)
nsv1 = NumericStartValueDict({
    C_f: 58,
    C_lab: 60,
    C_w: 770,
    C_r: 102
})

## We create the parameterized model explicitly , although the framework would find the ingredients
nsrm = NumericParameterizedSmoothReservoirModel(srm,np1)



ntimes = NumericSimulationTimes(np.arange(0, 1096, 1))

## We create the model run explicitly again, although the framework would find the ingredients
nsmr= numeric_model_run_1(
    nsrm, 
    numeric_start_value_array_1(nsv1,x),
    ntimes
)



qp1 = quantity_parameterization_1(
    np1,
    state_var_units=(
        kilogram/meter**2,
        kilogram/meter**2,
        kilogram/meter**2,
        kilogram/meter**2
    ),
    time_unit=day
)
qsv1 = QuantityStartValueDict({
    C_f: 58*kilogram/meter**2,
    C_lab: 60*kilogram/meter**2,
    C_w: 770*kilogram/meter**2,
    C_r: 102*kilogram/meter**2
})
qsrm = QuantityParameterizedSmoothReservoirModel(srm,qp1)
qtimes = QuantitySimulationTimes(np.arange(0, 1096, 1)*day)
## We create the model run explicitly again, although the framework would find the ingredients
qsmr= numeric_model_run_1(
    qsrm, 
    quantity_start_value_array_1(qsv1,x),
    qtimes
)

# Open questions regarding the translation
# - The variable G seems to be identical with GPP but
specialVars = {
    BibInfo(# Bibliographical Information
        name="DALEC",
        longName="Data Assimilation Linked Ecosystem model", 
        version="1",
        entryAuthor="Verónika Ceballos-Núñez",
        entryAuthorOrcid="0000-0002-0046-1160",
        entryCreationDate="16/9/2016",
        doi="10.1111/j.1365-2486.2004.00891.x",
        further_references=BibInfo(doi="10.1016/j.agrformet.2009.05.002"),
        #  Also from PDF in Reflex experiment
        sym_dict=sym_dict
        
    ),
    #
    # the following variables constitute the compartmental system:
    #
    A,  # the overall compartmental matrix
    Input,  # the overall imput
    t,  # time for the complete system
    x,  # state vector of the the complete system
    #
    # the following variables constitute a numerical model run :
    #
    np1,
    nsv1,
    ntimes,
    # Again the model run could be created explicitly (in several ways)
    #
    # the following variables constitute a quantiy model run :
    #
    qp1,
    qsv1,
    qtimes,
    # Again the model run could be created explicitly (in several ways)
    #
    # the following variables give a more detailed view with respect to
    # carbon and vegetation variables
    # This imformation can be used to extract the part
    # of a model that is concerned with carbon and vegetation
    # in the case of this model all of the state variables
    # are vegetation and carbon related but for an ecosystem model
    # for different elements there could be different subsystems
    # e.g. consisting of  Nitrogen Soil state variables
    #
    # vegetation carbon
    VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((C_f, C_lab, C_w, C_r))
}
