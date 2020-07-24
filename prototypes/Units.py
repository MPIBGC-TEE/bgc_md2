import matplotlib.pyplot as plt
import numpy as np
import time as tm
import unittest
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
from sympy.physics.units import (
    convert_to,
    Quantity,
    length,
    mass,
    time,
    second,
    minute,
    hour,
    day,
    year,
    kelvin,
    meter,
    kilometer,
    gram,
    giga,
    kilogram
)
from sympy.physics.units.systems import SI
from sympy.physics.units.systems.si import dimsys_SI
# only fixed qunatities (numbers times unit like in a parameter set)

# from inspect import signature

from bgc_md2.described_quantities import describedQuantity
from bgc_md2.resolve.mvars import ParameterDict

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun
def scalar_unit(expr):
    # extract unit if not 0
    if simplify(expr) == 0:
        raise Exception("""
            ant determine the unit of 0, please provide an expression that does
            not evaluate to zero
            """)
    return simplify(prod((1.0*expr.n()).args[1:]))


def scalar_value(expr):
    return expr.n().args[0]


def iterable_unit(exprs):
    # try to guess unit from an iterable if not zero
    simp_exprs = map(simplify, exprs)
    non_zero_exprs = [e for e in filter(
            lambda exp: 0 != exp,
            simp_exprs
    )]

    if len([e for e in non_zero_exprs]) < 1:
        raise Exception("""
            All expressions in the iterable evaluate to 0.  I cannot determine
            the unit of 0, please provide at least one expression that does not
            evaluate to zeroif __name__ == '__main__':
    unittest.main()
            """)

    units = tuple(set(map(scalar_unit, non_zero_exprs)))
    if len(units) > 1:
        raise Exception(
            "The components have different units:{units}".format(units=units)
        )

    return tuple(units)[0]


def to_number(q, targetUnit):
    q_s = simplify(q)
    return 0 if q_s == 0 else float(simplify((q/targetUnit)))


def plot_with_units(ax, xt, yt):
    print("##############################")
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

class TestUnitExamples(unittest.TestCase):
    # ####################################################################################################
    #
    # Exampe 1
    #
    def test_dimensions_of_expressions(self):
        a = Quantity("a")
        SI.set_quantity_dimension(a, length)

        t = Quantity("t")
        SI.set_quantity_dimension(t, time)

        res = a**2/t
        # we can now determine the physical dimension of res
        res_dim = SI.get_dimensional_expr(res)
        self.assertTrue(
            dimsys_SI.equivalent_dims(res_dim, length**2/time)
        )

        # In parameter dicts we can describe values along with units
        a_val = Quantity("a_val")
        SI.set_quantity_dimension(a_val, length)
        SI.set_quantity_scale_factor(a_val, 5*meter)

        parameter_dict = {a: 5*meter, t: 4*second}

        res_val = res.subs(parameter_dict)
        # we can now determine the physical dimension of res_val
        # and check it against the expected
        print(SI.get_dimensional_expr(res_val))
        self.assertTrue(
            dimsys_SI.equivalent_dims(
                SI.get_dimensional_expr(res),
                SI.get_dimensional_expr(res_val)
            )
        )

    # ####################################################################################################
    #
    # Exampe 2
    #
    # plot two iterables with units
    def test_scalar_unit_and_value(self):
        t = 4*hour
        self.assertEqual(scalar_unit(t), hour)
        self.assertEqual(scalar_value(t), 4)

        # zero has no unit
        with self.assertRaises(Exception):
            scalar_unit(0)

    def test_iterable_unit_and_value(self):
        self.assertEqual(iterable_unit([0,4*hour,5*hour]), hour)

        # zeros without values have no unit
        with self.assertRaises(Exception):
            iterable_unit([0, 0])

        # different units
        with self.assertRaises(Exception):
            iterable_unit([0*hour, 0*second])

    def test_to_number(self):
        t = 4*hour
        self.assertEqual(to_number(t, hour), 4)
        self.assertEqual(to_number(t, minute), 240)


    def test_plot(self):
        ts = [t*hour for t in range(1, 10)]
        ys = [meter/t for t in ts]

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


    # ###################################################################################################
    #
    # Exampe 4
    #
    # functions of quanteties substitute and lambify
    def test_numerification(self):
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

        # Note
        
        ts = [second*t for t in np.linspace(0, float(2*pi), 100)]
        ysf = [v_unit(t) for t in ts]

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        auto_plot_with_units(ax, ts, ysf)
        fig.savefig("example3_auto.pdf")
        t = Quantity("t")
        SI.set_quantity_dimension(t, time)

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

        #ysE = [E_funcUnit(t) for t in ts]
        #
        #fig = plt.figure()
        #ax = fig.add_subplot(1, 1, 1)
        #auto_plot_with_units(ax, ts, ysE)
        #fig.savefig("example4_auto.pdf")

        # ####################################################################################################
        #
        # Exampe 5
        # Here we assume that parameters, arguments and return values of functions
        # have units attached.
        # This is natural from sympys perspective and allows computations involving
        # quanteties ( numbers with units ) to be completely transparent.
        # E.g. a function or expression receiving a
        # length argument  will always compute the correct result since the argument
        # carries its unit with it. Also the result can be expressed in any unit
        # compatible with its dimension.
        #
        # However numerical computations in scipy like solving an ode require numbers
        # or functions that consume and return numbers (as opposed to quanteties)
        # To perform those computations we temporarily have to remove units from the
        # arguments before the computation and attach units to the results.
        # We can automatically convert 'numerify' a function working on quanteties
        # (f_quant) to  a function on number (f_num) if we choose the units for
        # arguments and return values.
        # The resulting purely numerical function f_num represents
        # f_quant faithfully for those units and for those units only.
        # This has several consequences:
        # - Along with f_num the untis of the arguments and the return value have to
        #   be remembered externally since they are no intrinsic part of f_num.
        # - The numerical representations of different functions and paramters
        #   should be consistent. E.g in f_quant(g_quant(x_qunat)) x_qant,f_quant
        #   and g_quant should be "numerified" making sure that x_num represents
        #   the original quantity with respect to the unit that f_num expects and
        #   that f_num returns its result w.r.t. the unit g_num expects...
        #   This is possible with the help of the unit system.
        # - Unfortunately the "numerification" of functions is computationally
        #   expensive as the numeric function f_num everytime it is called it will
        #   attach units call f_quant and remove units from the result.


        def numerify(f_quant, *units):
            def f_num(*num_args):
                target_unit = units[-1]
                u_args = tuple(num_arg*units[ind] for ind,num_arg in enumerate(num_args))
                res_quant = f_quant(*u_args)
                #res_quant_target_unit = convert_to(res_quant, target_unit)
                #res_num = factor(res_quant_target_unit, target_unit)/target_unit
                res_num = simplify(res_quant/target_unit)
                return float(res_num)

            return f_num


        C_0  = describedQuantity("C_0", mass, "")
        C_1  = describedQuantity("C_1", mass, "")
        t    = describedQuantity("t", time, "")
        k_01 = describedQuantity("k_01", 1/time, "")
        k_10 = describedQuantity("k_10", 1/time, "")
        k_0o = describedQuantity("k_0o", 1/time, "")

        k_1o = Function('k_1o')

        state_variables = [C_0, C_1]  # order is important
        inputs = {
            0: sin(t)+2,  # input to pool 0
            1: cos(t)+2   # input to pool 1
            }
        outputs = {
            0: k_0o*C_0**3,  # output from pool 0
            1: k_1o(t)*C_1**3   # output from pool 0
            }
        internal_fluxes = {
            (0, 1): k_01*C_0*C_1**2,  # flux from pool 0  to pool 1
            (1, 0): k_10*C_0*C_1  # flux from pool 1 to pool 0
            }
        time_symbol = t
        srm = SmoothReservoirModel(
            state_variables,
            time_symbol,
            inputs,
            outputs,
            internal_fluxes
        )

        par_dict_quant = {
            k_01: 1/100*1/day,
            k_10: 1/100*1/day,
            k_0o: 1/2*1/day,
        }


        def k_1o_func_quant(t):
            omega = 2*pi/day
            phi = pi/8
            V_0 = 20 * kilogram/day
            V_range = 5*kilogram/day
            u_res = V_0+V_range*sin(omega*t+phi)
            return u_res


        times_quant = np.linspace(0, 20, 16)*year
        start_values_quant = np.array([1*gram, 2*kilogram])

        # Note that the time units of the parameters the function and the time array
        # are different.  Also note that the components of the startvalue tuple are
        # given with respect to two different units (gigagram and kilogram)  and the
        # k_1o_func_quant uses kilogram/second as the unit for the return value.
        #
        # The different units are handled correctly by sympy, becuase the same quantety
        # can be represented with respect to differnt units (1*day=24*hour)
        #
        # If we convert ot purely numerical values and functions, we have to consider
        # the consistency of the whole ensemble The easiest way to achieve this is to
        # normalize all times to an arbitrary but common unit of time (e.g. year), and
        # all masses to an arbitrary unit of mass (e.g. kilogram)

        # create a numerical function expecting its argument to be a number measuring
        # years and returning a number to be interpreted as kilogram/year
        k_1o_func_num = numerify(k_1o_func_quant,year,kilogram/year)

        # create a par_dict of numbers (since all parameter in the dictionary are rates
        # we can deal with the whole dict at once
        par_dict_num = {k:to_number(v,1/year) for k,v in par_dict_quant.items()}

        # extract the times in years
        times_num = np.array([to_number(t,year) for t in times_quant])

        # crete numerical start values
        start_values_num = np.array(
            [to_number(v, kilogram) for v in start_values_quant]
        )
        n_smr = SmoothModelRun(
            srm,
            par_dict_num,
            start_values_num,
            times_num,
            func_set={k_1o: k_1o_func_num})
        before = tm.time()
        sol_num = n_smr.solve() # extremely slow
        after = tm.time()
        print(before-after)



        def k_1o_func_num_manual(t):
            omega = 2*pi
            phi = pi/8
            V_0 = 20 
            V_range = 5
            u_res = V_0+V_range*sin(omega*t+phi)
            return u_res

        n_smr = SmoothModelRun(
            srm,
            par_dict_num,
            start_values_num,
            times_num,
            func_set={k_1o: k_1o_func_num_manual})
        before = tm.time()
        sol_num = n_smr.solve()
        after = tm.time()
        print(before-after)

        #k_1o_func_quant_from_num = quantify(k_1o_func_num_manual,second,kilogram/second)



#suite=unittest.defaultTestLoader(TestUnitExamples)
#suite.run(unittest.TestResult())
tr=unittest.TextTestRunner()
#tr.run(TestUnitExamples('test_dimensions_of_expressions'))
#tr.run(TestUnitExamples('test_scalar_unit_and_value'))
#tr.run(TestUnitExamples('test_to_number'))
#tr.run(TestUnitExamples('test_plot'))
tr.run(TestUnitExamples('test_numerification'))

if __name__ == '__main__':
    unittest.main()
