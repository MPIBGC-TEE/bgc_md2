from sympy.physics.units import Quantity
from sympy.physics.units.systems import SI

def describedQuantity(name, dimension, description):
    obj=Quantity(name=name)
    SI.set_quantity_dimension(obj,dimension)
    obj.description = description
    #obj = Symbol(name)
    return obj

