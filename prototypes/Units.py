from sympy.physics.units import Quantity, length, time, day, second, meter, kilometer
from sympy.physics.units.systems import SI
from sympy.physics.units.systems.si import dimsys_SI

# only fixed qunatities (numbers times unit like in a parameter set)

a = Quantity("a")
SI.set_quantity_dimension(a, length)

t = Quantity("t")
SI.set_quantity_dimension(t, time)

res = a ** 2 / t
# we can now determine the physical dimension of res
print(SI.get_dimensional_expr(res))


# In parameter dicts we can describe values along with units
a_val = Quantity("a_val")
SI.set_quantity_dimension(a_val, length)
SI.set_quantity_scale_factor(a_val, 5 * meter)

t_val = Quantity("t_val")
SI.set_quantity_dimension(t_val, time)
SI.set_quantity_scale_factor(t_val, 6 * second)

parameter_dict = {a: a_val, t: t_val}

res_val = res.subs(parameter_dict)
# we can now determine the physical dimension of res_val
# and check it against the expected
print(SI.get_dimensional_expr(res_val))
# dimsys_SI.equivalent_dims(res,res_val)
