from bgc_md2.resolve.mvars import NumericStartValueDict
from sympy import Symbol, var
from functools import lru_cache
from CompartmentalSystems.smooth_model_run import SmoothModelRun
from CompartmentalSystems.TimeStepIterator import TimeStepIterator2
import CompartmentalSystems.helpers_reservoir as hr
import netCDF4 as nc
import numpy as np
from pathlib import Path
from collections import namedtuple
from functools import reduce
from bgc_md2.models.cable_yuanyuan.source import mvs 

def construct_matrix_func_sym(pa):
    # we create a parameterdict for the fixed values
    # and extend it by the parameters provided 
    symbol_names = mvs.get_BibInfo().sym_dict.keys()   
    for name in symbol_names:
        var(name)
    parDict = {
        clay: 0.2028,
        silt: 0.2808,
        lig_wood: 0.4,
        f_wood2CWD: 1,
        f_metlit2mic: 0.45,
    #    NPP: npp_in
    }
    model_params = {Symbol(k): v for k,v in pa._asdict().items()}
    parDict.update(model_params)
    B_func = hr.numerical_array_func(
            state_vector = mvs.get_StateVariableTuple(),
            time_symbol=mvs.get_TimeSymbol(),
            expr=mvs.get_CompartmentalMatrix(),
            parameter_dict=parDict,
            func_dict={}
    )
    # in the general nonautonomous nonlinear B_func is a function of t,x
    # although for this example it does not depend on either t, nor x.
    return B_func 
