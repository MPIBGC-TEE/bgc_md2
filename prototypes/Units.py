import matplotlib.pyplot as plt
import numpy as np
from sympy import Symbol, Function, prod, sin, pi, lambdify, simplify
from sympy.physics.units import (
        convert_to,
        Quantity,
        length,
        time,
        second,
        minute,
        hour,
        day,
        kelvin,
        meter,
        kilometer,
        kilogram)
from sympy.physics.units.systems import SI
from sympy.physics.units.systems.si import dimsys_SI
# only fixed qunatities (numbers times unit like in a parameter set)

from bgc_md2.described_quantities import describedQuantity
from bgc_md2.resolve.mvars import ParameterDict

# ####################################################################################################
#
# Exampe 1
#
a = Quantity("a")
SI.set_quantity_dimension(a, length)

t = Quantity("t")
SI.set_quantity_dimension(t, time)

res = a**2/t
# we can now determine the physical dimension of res
print(SI.get_dimensional_expr(res))


# In parameter dicts we can describe values along with units
a_val = Quantity("a_val")
SI.set_quantity_dimension(a_val, length)
SI.set_quantity_scale_factor(a_val, 5*meter)

parameter_dict = {a: 5*meter, t: 4*second}

res_val = res.subs(parameter_dict)
# we can now determine the physical dimension of res_val
# and check it against the expected
print(SI.get_dimensional_expr(res_val))
# dimsys_SI.equivalent_dims(res, res_val)

# ####################################################################################################
#
# Exampe 2
#
# plot two iterables with units

ts = [t*hour for t in range(1, 10)]
ys = [meter/t for t in ts]


def scalar_unit(expr):
    # extract unit if not 0
    if simplify(expr)==0:
        raise Exception("""
            ant determine the unit of 0, please provide an expression that does
            not evaluate to zero
            """)
    return simplify(prod((1.0*expr.n()).args[1:]))


def iterable_unit(exprs):
    # try to guess unit from an iterable if not zero
    simp_exprs = map(simplify, exprs)
    non_zero_exprs = [e for e in filter(
            lambda exp: 0 != exp,
            simp_exprs
    )]

    if len([e for e in non_zero_exprs]) < 1:
        raise Exception("""
            All expressions in the iterable evaluate to 0.
            I cannot determine the unit of 0, please provide at least one expression that does
            not evaluate to zero
            """)

    units = tuple(set(map(scalar_unit, non_zero_exprs)))
    if len(units) > 1:
        raise Exception("The components have different units:{units}".format(units=units))

    return tuple(units)[0]



def to_number(q, targetUnit):
    q_s = simplify(q)
    return 0 if q_s == 0 else simplify((q/targetUnit).n())


def plot_with_units(ax, xt, yt):
    xs, x_unit = xt
    ys, y_unit = yt
    x_nums = [to_number(x, x_unit) for x in xs]
    y_nums = [to_number(y, y_unit) for y in ys]
    ax.set_xlabel(str(x_unit))
    ax.set_ylabel(str(y_unit))
    ax.plot(x_nums, y_nums)


def auto_plot_with_units(ax, xs, ys):
    x_unit = iterable_unit(xs)
    y_unit = iterable_unit(ys)
    plot_with_units(ax, (xs, x_unit), (ys, y_unit))


fig = plt.figure()
# ax = fig.add_axes((0, 0, 1, 1))
ax = fig.add_subplot(1, 1, 1)
# The units can be prescribed by the user (hopefully correctly)
plot_with_units(ax, (ts, day), (ys, kilometer/second))
fig.savefig("example2.pdf")

fig = plt.figure()
# or determined automatically
ax = fig.add_subplot(1, 1, 1)
auto_plot_with_units(ax, ts, ys)
fig.savefig("example2_auto.pdf")


# ####################################################################################################
#
# Exampe 3
#
# apply a function defined with units

def v_unit(t):
    # Note:
    # At the moment t has to be a scalar
    cond_1 = simplify(t) == 0
    cond_2 = dimsys_SI.equivalent_dims(
        SI.get_dimensional_expr(t),
        time
    )
    assert(cond_1 | cond_2)
    omega = 2*pi/second
    # phi = 0
    phi = pi/8
    V_0 = 20*meter/second
    V_range = 5*meter/second
    return V_0+V_range*sin(omega*t+phi)


ts = [second*t for t in np.linspace(0, float(2*pi), 100)]
ysf = [v_unit(t) for t in ts]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
auto_plot_with_units(ax, ts, ysf)
fig.savefig("example3_auto.pdf")


# ####################################################################################################
#
# Exampe 4
#
# substitute and lambify


# We can transform Expressions to numerical functions
# This also works for Expressions of Quanteties and functions that
# contain units


v = Function('v')
m = Quantity('m')
# create an expression containing a function
E = m/2*v(t)**2


# print(v_unit(3*day))

# substitute paramters
E_parUnit = E.subs({m: 5*kilogram})

# lambify the expression to a function
tup = (t,)
E_funcUnit = lambdify(tup, E_parUnit, {'v': v_unit})

ysE = [E_funcUnit(t) for t in ts]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
auto_plot_with_units(ax, ts, ysE)
fig.savefig("example4_auto.pdf")

