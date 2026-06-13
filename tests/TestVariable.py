import unittest

import numpy as np

from bgc_md2.Variable import Variable, FluxVariable


class TestVariable(unittest.TestCase):
    def test_init(self):
        unit = "gC/m^2"
        var = Variable(data=np.arange(10), unit=unit)
        self.assertEqual(var.unit, "g/m^2")

        unit = "gC14/m^2"
        var = Variable(data=np.arange(10), unit=unit)
        self.assertEqual(var.unit, "g/m^2")

    def test_data_mult(self):
        var_masked_data = np.ma.masked_array(
            data=np.arange(9).reshape(3, 3),
            mask=[[0, 0, 0], [1, 0, 0], [0, 0, 0]],
            fill_value=-1,
        )
        var = Variable(data=var_masked_data, unit="gC/m^3")

        dz_masked_data = np.ma.masked_array(
            data=np.cumsum(np.arange(3)), mask=[0, 0, 1], fill_value=-1
        )
        dz = Variable(data=dz_masked_data, unit="m")

        res = var.data_mult(dz, given_axes=0)
        res2 = np.ma.masked_array(
            data=[[0, 0, 0], [3, 4, 5], [18, 21, 25]],
            mask=[[0, 0, 0], [1, 0, 0], [1, 1, 1]],
            fill_value=-1,
        )
        self.assertTrue(np.all(res.data == res2))
        self.assertEqual(res.unit, "m-2.kg")
        self.assertTrue(np.all(res.data.filled() == res2.filled()))
        res2.mask[2, 2] = 0
        self.assertFalse(np.all(res.data.filled() == res2.filled()))

        res = var.data_mult(dz, given_axes=1)
        res2 = np.ma.masked_array(
            data=[[0, 1, 6], [0, 4, 15], [0, 7, 24]],
            mask=[[0, 0, 1], [1, 0, 1], [0, 0, 1]],
            fill_value=-1,
        )
        self.assertTrue(np.all(res.data == res2))
        self.assertEqual(res.unit, "m-2.kg")
        self.assertTrue(np.all(res.data.filled() == res2.filled()))
        res2.mask[2, 2] = 0
        self.assertFalse(np.all(res.data.filled() == res2.filled()))

        var_data = np.ma.masked_array(
            data=np.arange(12).reshape(3, 2, 2), mask=False, fill_value=-1
        )
        var = Variable(data=var_data, unit="g/m^2")
        area_data = np.ma.masked_array(
            data=np.arange(4).reshape(2, 2), mask=[[0, 0], [1, 0]], fill_value=-1
        )
        area = Variable(data=area_data, unit="km^2")
        res = var.data_mult(area, given_axes=(1, 2))
        res2 = np.ma.masked_array(
            data=np.array([[[0, 1], [4, 9]], [[0, 5], [12, 21]], [[0, 9], [20, 33]]],),
            mask=[[[0, 0], [1, 0]], [[0, 0], [1, 0]], [[0, 0], [1, 0]]],
            fill_value=-1,
        )
        self.assertTrue(np.all(res.data.filled() == res2.filled()))
        self.assertEqual(res.unit, "t")
        res2.set_fill_value(-2)
        self.assertFalse(np.all(res.data.filled() == res2.filled()))

    def test_aggregateInTime(self):
        ## one-dimensioanl data (time-like)
        data = np.arange(10)
        unit = "kg"
        var = Variable(data=data, unit=unit)

        nstep = 1
        var_agg = var.aggregateInTime(nstep)
        res2 = np.arange(10)
        self.assertTrue(np.all(var_agg.data == res2))

        nstep = 2
        var_agg = var.aggregateInTime(nstep)
        res2 = np.array([0, 2, 4, 6, 8, 9])
        self.assertTrue(np.all(var_agg.data == res2))

        nstep = 3
        var_agg = var.aggregateInTime(nstep)
        res2 = np.array([0, 3, 6, 9])
        self.assertTrue(np.all(var_agg.data == res2))

        ## multi-dimensional
        data = np.arange(20).reshape((10, 2))
        unit = "kg"
        var = Variable(data=data, unit=unit)

        nstep = 3
        var_agg = var.aggregateInTime(nstep)
        res2 = [[0, 1], [6, 7], [12, 13], [18, 19]]
        self.assertTrue(np.all(var_agg.data == res2))

        ## test mask
        masked_data = np.ma.masked_array(
            data=np.arange(10), mask=[0, 0, 1, 1, 0, 1, 1, 0, 0, 0], fill_value=-1
        )
        var = Variable(data=masked_data, unit=unit)

        nstep = 1
        var_agg = var.aggregateInTime(nstep)
        res2 = np.ma.masked_array(
            data=np.arange(10), mask=[0, 0, 1, 1, 0, 1, 1, 0, 0, 0], fill_value=-1
        )
        self.assertTrue(np.all(var_agg.data.filled() == res2.filled()))

        nstep = 2
        var_agg = var.aggregateInTime(nstep)
        res2 = np.array([0, -1, 4, -1, 8, 9])
        self.assertTrue(np.all(var_agg.data.filled() == res2))

    def test_convert(self):
        masked_data = np.ma.masked_array(data=np.arange(2), mask=[1, 0], fill_value=-1)
        unit = "d"
        var = Variable(data=masked_data, unit=unit)
        tar_unit = "s"
        var.convert(tar_unit)
        self.assertTrue(np.all(var.data.filled() == [-1, 86400]))
        self.assertEqual(var.unit, "s")

        tar_unit = "m"
        with self.assertRaises(ValueError):
            var.convert(tar_unit)

    def test_add(self):
        masked_data = np.ma.masked_array(
            data=np.arange(3), mask=[0, 1, 0], fill_value=-1
        )
        unit = "g"
        var1 = Variable(data=masked_data, unit=unit)

        masked_data2 = np.ma.masked_array(data=[0, 2, 4], mask=[1, 0, 0], fill_value=-2)
        var2 = Variable(data=masked_data2, unit=unit)
        res = var1 + var2
        self.assertTrue(np.all(res.data.filled() == [-1, -1, 6]))
        self.assertEqual(var1.unit, res.unit)

        var3 = Variable(data=masked_data, unit="kg")
        with self.assertRaises(AssertionError):
            var1 + var3


class TestFluxVariable(unittest.TestCase):
    def test_aggregateInTime(self):
        ## one-dimensional
        masked_data = np.ma.masked_array(
            data=np.arange(10), mask=[0, 0, 1, 1, 0, 1, 1, 0, 0, 0], fill_value=-1
        )
        unit = "kg"
        fv = FluxVariable(data=masked_data, unit=unit)

        nstep = 1
        fv_agg = fv.aggregateInTime(nstep)
        res2 = [0, 1, -1, -1, 4, -1, -1, 7, 8, 9]
        self.assertTrue(np.all(fv_agg.data.filled() == res2))

        nstep = 2
        fv_agg = fv.aggregateInTime(nstep)
        res2 = np.array([1, -1, -1, -1, 17])
        self.assertTrue(np.all(fv_agg.data.filled() == res2))

        nstep = 3
        fv_agg = fv.aggregateInTime(nstep)
        res2 = np.array([-1, -1, -1, 9])
        self.assertTrue(np.all(fv_agg.data.filled() == res2))

        ## multi-dimensional
        masked_data = np.ma.masked_array(
            data=np.arange(20).reshape((10, 2)),
            mask=[
                [0, 0],
                [0, 1],
                [1, 0],
                [0, 0],
                [0, 0],
                [0, 0],
                [1, 0],
                [0, 0],
                [0, 0],
                [0, 1],
            ],
            fill_value=-1,
        )
        unit = "kg"
        fv = FluxVariable(data=masked_data, unit=unit)

        nstep = 3
        fv_agg = fv.aggregateInTime(nstep)
        res2 = [[-1, -1], [24, 27], [-1, 45], [18, -1]]
        self.assertTrue(np.all(fv_agg.data.filled() == res2))

    def test_absolute_and_relative_error(self):
        v_right = Variable(data=np.arange(10).reshape(2, 5), unit="m")
        v_wrong = Variable(data=np.arange(10).reshape(2, 5) + 1, unit="m")
        abs_err = v_wrong.absolute_error(v_right)
        self.assertTrue(np.all(abs_err.data == np.ones(10).reshape(2, 5)))
        self.assertEqual(abs_err.unit, "m")
        rel_err = v_wrong.relative_error(v_right)
        res = 100 / np.ma.arange(10).reshape(2, 5)
        self.assertTrue(np.all(rel_err.data == res.data))
        self.assertEqual(rel_err.unit, "%")

        ## convertible units
        v_right = Variable(data=np.arange(10).reshape(2, 5) * 1000, unit="m")
        v_wrong = Variable(data=np.arange(10).reshape(2, 5) + 1, unit="km")
        abs_err = v_wrong.absolute_error(v_right)
        self.assertTrue(np.all(abs_err.data == np.ones(10).reshape(2, 5) * 1000))
        self.assertEqual(abs_err.unit, "m")
        rel_err = v_wrong.relative_error(v_right)
        res = 100 / np.ma.arange(10).reshape(2, 5)
        self.assertTrue(np.all(rel_err.data == res.data))
        self.assertEqual(rel_err.unit, "%")

        ## unconvertible units
        v_right = Variable(data=np.arange(10).reshape(2, 5) * 1000, unit="m")
        v_wrong = Variable(data=np.arange(10).reshape(2, 5) + 1, unit="km/s")
        with self.assertRaises(ValueError):
            abs_err = v_wrong.absolute_error(v_right)
        with self.assertRaises(ValueError):
            rel_err = v_wrong.relative_error(v_right)

        ## incompatible shapes
        v_right = Variable(data=np.arange(10) * 1000, unit="m")
        v_wrong = Variable(data=np.arange(10).reshape(2, 5) + 1, unit="km/s")
        with self.assertRaises(AssertionError):
            abs_err = v_wrong.absolute_error(v_right)
        with self.assertRaises(AssertionError):
            rel_err = v_wrong.relative_error(v_right)


################################################################################


if __name__ == "__main__":
    unittest.main()
