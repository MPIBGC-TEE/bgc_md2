from sympy import symbols ,Symbol,ImmutableMatrix,exp
from sympy.physics.units.systems import SI
from sympy.physics.units import Quantity, length, mass,time,temperature
from sympy.physics.units.systems.si import dimsys_SI

from bgc_md2.resolve.mvars import (
        CompartmentalMatrix
        ,InputTuple
        ,TimeSymbol
        ,StateVariableTuple
        ,VegetationCarbonInputScalar
        ,VegetationCarbonInputPartitioningTuple
        ,VegetationCarbonStateVariableTuple
        )


#class DescribedQuantity(Quantity):
#    def __new__(cls,name,dimension,description,**kwargs):
#        # build a normal sympy Quantity (works like Symbol)
#        obj=super().__new__(cls,name=name)
#
#
#        # add the additional text field
#        SI.set_quantity_dimension(obj,dimension)
#
#        # add the additional text field
#        obj.description = description
#        return obj

def describedQuantity(name, dimension, description):
    obj=Quantity(name=name)
    SI.set_quantity_dimension(obj,dimension)
    obj.description = description
    return obj


# define Quantities (as you would symbols, but with an additional description)
t=TimeSymbol('t')

C_f     =  describedQuantity("C_f"      , mass/length**2                ,"Foliar C mass")
C_lab	=  describedQuantity("C_lab"	, mass/length**2                ,"Labile C mass")
C_w		=  describedQuantity("C_w"		, mass/length**2                ,"Wood C mass")
C_r		=  describedQuantity("C_r"		, mass/length**2                ,"Fine root C mass")
p_10	=  describedQuantity("p_10"	    , 1/temperature                 ,"Parameter in exponential term of temperature dependent rate parameter ")
mint	=  describedQuantity("mint"	    , temperature                   ,"Dayly minimum temperature")
maxt	=  describedQuantity("maxt"	    , temperature                   ,"Dayly maximum temperature")
multtl	=  describedQuantity("multtl"	, 1/time                        ,"Turnover of labile C (0 = off, 1 = On)")
multtf	=  describedQuantity("multtf"	, 1/time                        ,"Turnover of foliage C (0 = off, 1 = On)")
LAI		=  describedQuantity("LAI"		, 1                             ,"LeafAreaIndex")
p_3		=  describedQuantity("p_3"		, 1	                            ,"Fraction of NPP partitioned to foliage")
p_4		=  describedQuantity("p_4"		, 1	                            ,"Fraction of NPP partitioned to roots")
p_5		=  describedQuantity("p_5"		, 1/time                        ,"Turnover rate of foliage")
p_6		=  describedQuantity("p_6"		, 1/time           	            ,"Turnover rate of wood")
p_7		=  describedQuantity("p_7"		, 1/time           	            ,"Turnover rate of roots")
p_14	=  describedQuantity("p_14"	    , 1	                            ,"Fraction of leaf loss transferred to litter")
p_15	=  describedQuantity("p_15"	    , 1/time                        ,"Turnover rate of labile carbon")
p_16	=  describedQuantity("p_16"	    , 1	                            ,"Fraction of labile transfers respired ")
NPP		=  describedQuantity("NPP"		, mass*length**(-2)*time**(-1)  ,"Net Primary Production per area")

#C_f     =  Quantity("C_f"       , mass/length**2                )   #"Foliar C mass")
#C_lab	=  Quantity("C_lab"	    , mass/length**2                )   #"Labile C mass")
#C_w		=  Quantity("C_w"		, mass/length**2                )   #"Wood C mass")
#C_r		=  Quantity("C_r"		, mass/length**2                )   #"Fine root C mass")
#p_10	=  Quantity("p_10"	    , 1/temperature                 )   #"Parameter in exponential term of temperature dependent rate parameter ")
#mint	=  Quantity("mint"	    , temperature                   )   #"Dayly minimum temperature")
#maxt	=  Quantity("maxt"	    , temperature                   )   #"Dayly maximum temperature")
#multtl	=  Quantity("multtl"	, 1/time                        )   #"Turnover of labile C (0 = off, 1 = On)")
#multtf	=  Quantity("multtf"	, 1/time                        )   #"Turnover of foliage C (0 = off, 1 = On)")
#LAI		=  Quantity("LAI"		, 1                             )   #"LeafAreaIndex")
#p_3		=  Quantity("p_3"		, 1	                            )   #"Fraction of NPP partitioned to foliage")
#p_4		=  Quantity("p_4"		, 1	                            )   #"Fraction of NPP partitioned to roots")
#p_5		=  Quantity("p_5"		, 1/time                        )   #"Turnover rate of foliage")
#p_6		=  Quantity("p_6"		, 1/time           	            )   #"Turnover rate of wood")
#p_7		=  Quantity("p_7"		, 1/time           	            )   #"Turnover rate of roots")
#p_14	=  Quantity("p_14"	    , 1	                            )   #"Fraction of leaf loss transferred to litter")
#p_15	=  Quantity("p_15"	    , 1/time                        )   #"Turnover rate of labile carbon")
#p_16	=  Quantity("p_16"	    , 1	                            )   #"Fraction of labile transfers respired ")
#NPP		=  Quantity("NPP"		, mass*length**(-2)*time**(-1)  )   #"Net Primary Production per area")


# The following parameters would be needed to compute NPP from GPP  since they are not given
#    ("p_2"		, "1"	                            ,"Fraction of GPP respired "), # no parameter values in publication
#    ("GPP"		, "mass*length**(-2)*time**(-1)"	,"Gross Primary Production per area"),
#    ("G"		, "mass*length**(-2)*time**(-1)"	                            ,),

# Temperature sensitive rate parameter
T_rate = 0.5*exp(p_10*0.5*(maxt+mint))

# state vector
x = StateVariableTuple((C_f, C_lab, C_w, C_r))

# NetPrimaryProduction can not be computted in the absence of G and p_2 with
# the following formula NPP = (1-p_2)*G # CAUTION: there are 2 NPPs, need to
# correct according to ACM model so NPP has been made a parameter instead of an
# expression

#scalar function of photosynthetic inputs
u = NPP

# tuple of partitioning coefficients of photosynthetically fixed carbon
b = ((p_3*multtl), 0, (1-p_4), p_4) 

#f_v = u*b+A*x

# matrix of cycling rates
A=CompartmentalMatrix(
        [
            [((-1)*(((1-p_14)*p_5*(1-p_16)*multtf*T_rate)+((1-p_14)*p_5*p_16*multtf*T_rate)+(p_5*p_14*multtf))),    (p_15*(1-p_16)*multtl*T_rate),                                  0,      0],
            [((1-p_14)*p_5*(1-p_16)*multtf*T_rate),                                                                 ((-p_15*(1-p_16)*multtl*T_rate)-(p_15*p_16*multtl*T_rate)),     0,      0],
            [0,                                                                                                     0,                                                           -p_6,      0],
            [0,                                                                                                     0,                                                              0,   -p_7]
        ]
)
I = InputTuple(u*ImmutableMatrix(b))

#model_run_data:
#    parameter_sets:
#        - "param_vals":
#            desc: "NPP value was given per year. p~14~, p~15~ and p~16~ values not given in publication"
#            values: {NPP: 'Rational(409,365)', mint: -4, maxt: 5, multtf: 0, multtl: 1, p_2: 0.47, p_3: 0.31, p_4: 0.43, p_5: 0.0027, p_6: 0.00000206, p_7: 0.00248, p_10: 0.0693, p_14: 0.45, p_15: 0.001, p_16: 0.25}
#            doi: 10.1111/j.1365-2486.2004.00891.x
#    initial_values:
#        - "init_vals":
#            desc: C_lab values not provided in publication
#            values: {C_f: 58, C_lab: 60, C_w: 770, C_r: 102}
#            doi: 10.1111/j.1365-2486.2004.00891.x
#    run_times:
#        - RT1:
#            start: 0
#            end: 1096 # See drivers data set
#            step_size: 1
#
#    possible_combinations:
#        - ["param_vals", "init_vals", RT1]
#

# Open questions regarding the translation
# - The variable G seems to be identical with GPP but
specialVars={
    A, # the overall compartmental matrix
    I, # the overall imput
    t, # time for the complete system
    x, # state vector of the the complete system
    # 
    ## the following variables give a more detailed view with respect to 
    ## carbon and vegetation variables
    ## This imformation can be used to extract the part
    ## of a model that is concerned with carbon and vegetation
    ## in the case of this model all of the state variables 
    ## are vegetation and carbon related but for an ecosystem model
    ## for different elements there could be different subsystems 
    ## e.g. consisting of  Nitrogen Soil state variables 
    VegetationCarbonInputScalar(u),   # vegetation carbon 
    VegetationCarbonInputPartitioningTuple(b), # vegetation carbon partitioning.
    VegetationCarbonStateVariableTuple((C_f, C_lab, C_w, C_r))
}

