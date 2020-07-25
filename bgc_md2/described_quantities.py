import numpy as np
from sympy import (
    Symbol,
    # symbols,
    Function,
    prod,
    sin,
    cos,
    pi,
    lambdify,
    simplify,
    factor
)
from sympy.physics.units import Quantity
from sympy.physics.units.systems import SI
from sympy.physics.units import (
    day,
    year,
    kilogram,
    gram
)
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun

def describedQuantity(name, dimension, description):
    obj=Quantity(name=name)
    SI.set_quantity_dimension(obj,dimension)
    obj.description = description
    #obj = Symbol(name)
    return obj

def to_number(q, targetUnit):
    q_s = simplify(q)
    return 0 if q_s == 0 else float(simplify((q/targetUnit)))

