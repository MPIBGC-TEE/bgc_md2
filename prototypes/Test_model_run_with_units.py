import unittest

import matplotlib.pyplot as plt
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
from sympy.physics.units import (
    day,
    year,
    kilogram
)

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun

class ModelRunWithUnits(SmoothModelRun):
    def __init__(
        self,
        srm,
        par_dict,
        start_values,
        times,
        func_dict,
        state_var_unit,
        time_unit,
    ):
        super().__init__(
            srm,
            par_dict,
            start_values,
            times,
            func_dict
        )
        self.state_var_unit = state_var_unit
        self.time_unit = time_unit

    def solve(self):
        sol_num = super().solve()
        return sol_num*self.state_var_unit


class TestModelRunWithUnits(unittest.TestCase):

    def test_solve(self):
        C_0,  C_1, t, k_01, k_10, k_0o = (
                Symbol(s) for s in ("C_0", "C_1", "t", "k_01", "k_10", "k_0o")
        )

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

        # the solution of a modelrun that received quatities as parameters
        # should be a quantity and have a unit
        # The quantities of the results could be infered from parameters and functions as quanteties but this is very
        # expensive computationally. 
        # To minimize computational cost we require the user to provide all expressions, parameters, functions, start values and times 
        # with respect to consistent units for time and statevariables.
        # From these tow units all computable quanteties can be inferred
        par_dict = {
            k_01: 1/100,  # 1/year
            k_10: 1/100,  # 1/year
            k_0o: 1/2     # 1/year  
        }
        times = np.linspace(0, 20, 16)  # year
        start_values= np.array([1, 2])  # kg

        def k_1o_func(t):
            omega = 2*pi    # 1/year
            phi = pi/8
            V_0 = 20        # kilogram/year
            V_range = 5     # kilogram/year
            u_res = V_0+V_range*sin(omega*t+phi)
            return u_res

        smrwu = ModelRunWithUnits(
            srm,
            par_dict,
            start_values,
            times,
            {k_1o: k_1o_func},
            kilogram,
            year
        )
        print(smrwu.solve())
