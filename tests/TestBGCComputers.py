import unittest
from sympy import (
    Symbol,
    symbols,
    var,
    sympify,
    Function,
    prod,
    sin,
    cos,
    pi,
    lambdify,
    simplify,
)

import numpy as np

from bgc_md2.resolve.mvars import (
    InFluxesBySymbol,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
    VegetationCarbonInputTuple,
    NumericParameterization,
    NumericStartValueDict,
    NumericStartValueArray,
    NumericSimulationTimes,
    NumericParameterizedSmoothReservoirModel,
    QuantityParameterization,
    QuantityParameterizedSmoothReservoirModel,
    QuantitySimulationTimes,
    QuantityModelRun,
)
from bgc_md2.resolve.computers import (
        quantity_model_run_1,
        vegetation_carbon_input_tuple_1,
        vegetation_carbon_input_tuple_2
)
from sympy.physics.units import days, day, kilogram, gram

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun


class TestBGCComputers(unittest.TestCase):
    def setUp(self):
        C_0, C_1, t, k_01, k_10, k_0o = (
            Symbol(s) for s in ("C_0", "C_1", "t", "k_01", "k_10", "k_0o")
        )

        k_1o = Function("k_1o")

        state_variables = [C_0, C_1]  # order is important
        # input to pool 0  # input to pool 1
        inputs = {0: sin(t) + 2, 1: cos(t) + 2}
        outputs = {
            0: k_0o * C_0 ** 3,  # output from pool 0
            1: k_1o(t) * C_1 ** 3,  # output from pool 0
        }
        internal_fluxes = {
            (0, 1): k_01 * C_0 * C_1 ** 2,  # flux from pool 0  to pool 1
            (1, 0): k_10 * C_0 * C_1,  # flux from pool 1 to pool 0
        }
        time_symbol = t
        self.srm = SmoothReservoirModel(
            state_variables, time_symbol, inputs, outputs, internal_fluxes
        )

    def test_quantity_model_run_1(self):
        C_0, C_1, t, k_01, k_10, k_0o = (
            Symbol(s) for s in ("C_0", "C_1", "t", "k_01", "k_10", "k_0o")
        )
        k_1o = Function("k_1o")
        
        par_dict = {
            k_01: 1 / 100,  # 1/year
            k_10: 1 / 100,  # 1/year
            k_0o: 1 / 2,  # 1/year
        }

        def k_1o_func(t):
            omega = 2 * pi  # 1/year
            phi = pi / 8
            V_0 = 20  # kilogram/year
            V_range = 5  # kilogram/year
            u_res = V_0 + V_range * sin(omega * t + phi)
            return u_res

        para = QuantityParameterization(
                par_dict=par_dict,
                func_dict={k_1o: k_1o_func},
                state_var_units=(kilogram,kilogram),
                time_unit=days
        )
        qpm = QuantityParameterizedSmoothReservoirModel(srm=self.srm, parameterization=para)
        # The parameterdict and the functions, possibly even the matrix/flux expressions
        # have implicit unit assumption, which the user is required to maintain consistently.
        # To make modelruns comparable it is important to remember the units for which this
        # consistency can be guaranteed.

        # A model run can then adapt the units of times and masses since it
        times_quant = QuantitySimulationTimes(np.linspace(0, 20, 16) * day)
        start_values_quant = [1 * kilogram, 2 * gram]  # kg

        qmr = quantity_model_run_1(qpm, start_values_quant, times_quant)
        #qmr = QuantityModelRun(qpm, start_values_quant, times_quant)
        print(qmr.solve())

class TestVegetationCarbonComputers(unittest.TestCase):
    def test_vegetation_carbon_input_tuple_1(self):
        for name in ("a_L","a_S","a_R"):
            var(name)
        u = VegetationCarbonInputScalar(sympify(1))
        b = VegetationCarbonInputPartitioningTuple((a_L,a_S,a_R))
        res_1 = u*b
        print(type(res_1))

        #res_2 = VegetationCarbonInputTuple(res_1)
        #print(type(res_2))

        #res_3 = vegetation_carbon_input_tuple_1(u,b)
    

    def test_vegetation_carbon_input_tuple_2(self):
        Labile = Symbol("Labile")
        Leaf = Symbol("Leaf")
        Root = Symbol("Root")
        Wood = Symbol("Wood")
        gpp_to_labile = Symbol("gpp_to_labile")
        gpp_to_leaf = Symbol("gpp_to_leaf")
        gpp_to_root = Symbol("gpp_to_root")
        gpp_to_wood = Symbol("gpp_to_wood")
        d = {
                𝙻𝚊𝚋𝚒𝚕𝚎: gpp_to_labile,
                Leaf: gpp_to_leaf,
                𝚁𝚘𝚘𝚝: gpp_to_root,
                𝚆𝚘𝚘𝚍: gpp_to_wood
        }
        u = InFluxesBySymbol(d)
        vcsv = VegetationCarbonStateVariableTuple((Labile, Leaf, Root, Wood))
        res = vegetation_carbon_input_tuple_2(u, vcsv)
        self.assertEqual(type(res),VegetationCarbonInputTuple)
