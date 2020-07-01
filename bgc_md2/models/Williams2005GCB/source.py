from sympy import symbols ,Symbol,ImmutableMatrix
from sympy.physics.units.systems import SI
from sympy.physics.units import Quantity, length, mass
from sympy.physics.units.systems.si import dimsys_SI

from bgc_md2.resolve.mvars import (
        CompartmentalMatrix
        ,TimeSymbol
        ,StateVariableTuple
        )


class DescribedQuantity(Quantity):
    def __new__(cls,*args,description,**kwargs):
        # build a normal sympy Quantity (works like Symbol)
        obj=super().__new__(cls,*args,**kwargs)


        # add the additional text field
        obj.description = description
        return obj


# define Quantities (as you would symbols, but with an additional description)
C_f=DescribedQuantity("C_f",description="Foliar C mass")
C_lab=DescribedQuantity("C_lab",description="Labile C mass")
C_w=DescribedQuantity("C_w",description="Wood C mass")
C_r=DescribedQuantity("C_r",description="Fine root C mass")


# It is possible but not required to set the dimension If provided it can be
# used later to derive the dimensions of expressions containing the variable or
# to check the dimensional correctness of the these expressions.
for q in [C_f]:#,C_lab,C_w,C_r]:
    SI.set_quantity_dimension(q,mass/length**2)

#    - phenology_parameters:
#        - t:
#            desc: time 
#            unit: "day"
#            type: variable
#            key: time_symbol
#        - p_10:
#            range: "0.05 - 0.2"
#            desc: Parameter in exponential term of temperature dependent rate parameter 
#            key:
#            type: parameter
#        - mint:
#            desc: Dayly minimum temperature
#            type: variable
#            key: Temperature
#        - maxt:
#            desc: Dayly maximum temperature
#            type: variable
#            key: Temperature
#        - T_rate:
#            desc: Temperature sensitive rate parameter
#            type: variable
#            exprs: "T_rate = 0.5*exp(p_10*0.5*(maxt+mint))"
#        - multtl:
#            desc: "Turnover of labile C (0 = off, 1 = On)"
#        - multtf:
#            desc: "Turnover of foliage C (0 = off, 1 = On)"
#    - respiration_parameters:
#        - p_2:
#            range: "0.2 - 0.7"
#            desc: Fraction of GPP respired 
#            key:
#            type: parameter
#        - p_16:
#            range: "0.01 - 0.5"
#            desc: Fraction of labile transfers respired 
#            unit: "day^{-1}"
#            key:
#            type: parameter
#    - photosynthetic_parameters:
##        - LAI:
##            desc: 
##            key:
##            type:
##            unit:
##            exprs:
##        - GPP:
##            desc: 
##            key:
##            type:
##            unit:
##            exprs:
##        - G:
##            desc: 
##            key:
##            type:
##            unit:
##            exprs:
#        - NPP:
#            desc: 
##            exprs: "NPP = (1-p_2)*G" # CAUTION: there are 2 NPPs, need to correct according to ACM model
#            unit: "gC*m^{-2}*day^{-1}" 
#            type: variable
#            key: NPP
#    - partitioning_coefficients:
#        - p_3:
#            range: "0.01 - 0.5"
#            desc: Fraction of NPP partitioned to foliage
#            key: "part_foliage"
#            type: parameter
#        - p_4:
#            range: "0.01 - 0.5" 
#            desc: Fraction of NPP partitioned to roots
#            key: "part_roots"
#            type: parameter
#    - cycling_rates:
#        - p_5:
#            range: "1*10**(-4) - 0.1"
#            desc: Turnover rate of foliage
#            unit: "day^{-1}"
#            key: "cyc_foliage"
#            type: parameter
#        - p_6:
#            range: "1*10**(-6) - 0.01"
#            desc: Turnover rate of wood
#            unit: "day^{-1}"
#            key: "cyc_wood"
#            type: parameter
#        - p_7:
#            range: "1*10**(-4) - 0.01"
#            desc: Turnover rate of roots
#            unit: "day^{-1}"
#            key: "cyc_roots"
#            type: parameter
#        - p_14:
#            range: "0.2 - 0.7"
#            desc: Fraction of leaf loss transferred to litter 
#            key: "cyc_foliage" #?
#            type: parameter
#        - p_15:
#            range: "1*10**(-4) - 0.1"
#            desc: Turnover rate of labile carbon
#            unit: "day^{-1}"
#            key: "cyc_labile" # Need to include it in list or search for cyc_NSC?
#            type: parameter
#    - components:
#        - x:
#            desc: vector of states of vegetation
#            key: state_vector
#            exprs: "x = Matrix(4,1,[C_f, C_lab, C_w, C_r])"
#        - u:
#            desc: scalar function of photosynthetic inputs
#            key: scalar_func_phot
#            exprs: "u = NPP"
#        - b:
#            desc: vector of partitioning coefficients of photosynthetically fixed carbon
#            key: part_coeff
#            exprs: "b = Matrix(4,1,[(p_3*multtl),0,(1-p_4), p_4])"
#        - A:
#            desc: matrix of cycling rates
#            key: cyc_matrix
#            exprs: "A = Matrix([[((-1)*(((1-p_14)*p_5*(1-p_16)*multtf*T_rate)+((1-p_14)*p_5*p_16*multtf*T_rate)+(p_5*p_14*multtf))), (p_15*(1-p_16)*multtl*T_rate), 0, 0],
#                                [((1-p_14)*p_5*(1-p_16)*multtf*T_rate), ((-p_15*(1-p_16)*multtl*T_rate)-(p_15*p_16*multtl*T_rate)), 0, 0],
#                                [0, 0, -p_6, 0],
#                                [0, 0, 0, -p_7]])"
#        - f_v:
#            desc: the right hand side of the ode
#            key: state_vector_derivative
#            exprs: "f_v = u*b+A*x"
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
#A = ImmutableMatrix([[((-1)*(((1-p_14)*p_5*(1-p_16)*multtf*T_rate)+((1-p_14)*p_5*p_16*multtf*T_rate)+(p_5*p_14*multtf))), (p_15*(1-p_16)*multtl*T_rate), 0, 0],
#                                [((1-p_14)*p_5*(1-p_16)*multtf*T_rate), ((-p_15*(1-p_16)*multtl*T_rate)-(p_15*p_16*multtl*T_rate)), 0, 0],
#                                [0, 0, -p_6, 0],
#                                [0, 0, 0, -p_7]])
specialVars={
    CompartmentalMatrix([[C_f]]),
    TimeSymbol('t'),
    StateVariableTuple((C_f,))#, C_lab, C_w, C_r))
    #,'SmoothReservoirModel':srm
}
