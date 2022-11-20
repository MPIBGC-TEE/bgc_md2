import numpy as np
import netCDF4 as nc
from unittest import TestCase, skip
from numpy.core.fromnumeric import shape
from pathlib import Path
from time import time, sleep
from testinfrastructure.InDirTest import InDirTest
import general_helpers as gh
import matplotlib.pyplot as plt


class Test_general_helpers(InDirTest):

    # def test_make_fluxrates_from_kf(,xi_d):

    @skip
    # the test is now obsolete since we 
    # decided to work with a continuous function
    # and a constant factor of 30 between days and months.

    def test_month_2_day_index(self):
        self.assertEqual(month_2_day_index([0]), [0])
        self.assertEqual(month_2_day_index([1]), [31])
        self.assertEqual(month_2_day_index([2]), [59])
        self.assertEqual(month_2_day_index([3]), [90])
        self.assertEqual(month_2_day_index([1, 3]), [31, 90])

    def test_day_2_month_index(self):
        # note that days are counted from zero so day 29 is January 30.
        self.assertEqual(gh.day_2_month_index(0), 0)
        self.assertEqual(gh.day_2_month_index(30), 1)
        self.assertEqual(gh.day_2_month_index(31), 1)
        self.assertEqual(gh.day_2_month_index(60), 2)

    # def test_month_2_day_index_vm(self):
    #    self.assertEqual(
    #            gh.month_2_day_index_vm([0]),
    #            [0]
    #    )
    #    self.assertEqual(
    #            gh.month_2_day_index_vm([1]),
    #            [31]
    #    )
    #    self.assertEqual(
    #            gh.month_2_day_index_vm([2]),
    #            [59]
    #    )
    #    self.assertEqual(
    #            gh.month_2_day_index_vm([3]),
    #            [90]
    #    )
    #    self.assertEqual(
    #            gh.month_2_day_index_vm([1,3]),
    #            [31,90]
    #    )

    # def test_day_2_month_index_vm(self):
    #    # note that days are counted from zero so day 30 is January 31.
    #    self.assertEqual(gh.day_2_month_index_vm( 0), 0)
    #    self.assertEqual(gh.day_2_month_index_vm(30), 0)
    #    self.assertEqual(gh.day_2_month_index_vm(31), 1)
    #    self.assertEqual(gh.day_2_month_index_vm(60), 2)

    def test_pixel_area_on_unit_sphere(self):
        # we test that the pixel areas for a whole
        # grid sum up to 1
        # lat_len=181
        # lon_len=360
        lat_len = 19
        lon_len = 360
        lats = np.ma.masked_array(np.linspace(-90, 90, lat_len))
        lons = np.ma.masked_array(np.linspace(-179.5, 179.5, lon_len))
        delta_lat = (lats.max() - lats.min()) / (len(lats) - 1)
        delta_lon = (lons.max() - lons.min()) / (len(lons) - 1)

        puaf_sym = gh.make_pixel_area_on_unit_spehre(delta_lat, delta_lon, sym=True)
        puaf = gh.make_pixel_area_on_unit_spehre(delta_lat, delta_lon)

        # assert identical values
        lw = np.array(
            [puaf(lats[lat_ind]) for lat_ind in range(len(lats))], dtype=np.float64
        )
        lws = np.array(
            [puaf_sym(lats[lat_ind]) for lat_ind in range(len(lats))], dtype=np.float64
        )
        # check for equivalence of our manually solved integral with sympy
        # (found typos this way....)
        self.assertTrue(np.allclose(lw, lws))

        # check for symetrie (the weighs for pos and negative latitudes should
        # be the same
        pres = np.array([puaf(lat) for lat in lats], dtype=np.float64)
        nres = np.array([puaf(-lat) for lat in lats], dtype=np.float64)
        self.assertTrue(np.allclose(pres, nres))

        # having established the equality of the two ways to compute the
        # weights we check that they add up to 1
        self.assertTrue(
            np.allclose(
                np.sum(
                    [
                        np.sum([puaf(lats[lat_ind]) for lon_ind in range(len(lons))])
                        for lat_ind in range(len(lats))
                    ],
                    dtype=np.float64,
                ),
                4 * np.pi,  # surface area of the unit sphere
            )
        )

    def test_global_average(self):
        # we create data similar to cmip6 and trendy. These are masked arrays:
        # the lat lon combinations that have no landpoints are marked with
        # False in the mask which is a boolean array of the same shape as
        # the actual value array) the
        #
        # The following commented lines show how to get the
        # lats and longs are as in jon_yibs
        # lat_name = 'latitude'
        # lon_name = 'longitude'
        # lats=nds.variables[lat_name].__array__()
        # longs=nds.variables[lon_name].__array__()
        # target_var_names=set(nds.variables.keys()).difference([lat_name,lon_name,'time'])
        # arr = target_var_names[0]
        lat_len_1 = 3
        lat_len_2 = 3
        lat_len_3 = 3
        lat_len = sum([lat_len_1, lat_len_2, lat_len_3])
        lon_len_1 = 120
        lon_len_2 = 120
        lon_len_3 = 120
        lon_len = sum([lon_len_1, lon_len_2, lon_len_3])
        time_len = 3
        lats = np.ma.masked_array(np.linspace(-90, 90, lat_len))
        lons = np.ma.masked_array(np.linspace(-179.5, 179.5, lon_len))
        arr = np.ma.ones(shape=(time_len, lat_len, lon_len))

        res = gh.global_mean(lats, lons, arr)
        self.assertEqual(res.shape, (time_len,))

        # we assert that the average of a constant field is the constant value
        self.assertTrue(np.allclose(res, np.ones((time_len,))))

        # now we test the same properties with a masked array
        # as e.g the CMIP6 files that Kostia uses
        # dataPath = Path(conf_dict['dataPath'])
        # ds=nc.Dataset('cLeaf_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc')
        # pwd
        # dataPath
        # ds=nc.Dataset(dataPath.joinpath('cLeaf_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'))
        # cLeaf = ds.variables['cLeaf'].__array__()
        # in this array many fields are masked

        # note that in python boolean values are equivalent to integers
        # 1 = True , 0 = False
        mask = np.stack(
            [
                np.ones(shape=(time_len, lat_len_1, lon_len)),
                np.zeros(shape=(time_len, lat_len_2, lon_len)),
                np.ones(shape=(time_len, lat_len_3, lon_len)),
            ],
            axis=1,
        )
        arr = np.ma.array(np.ones(shape=(time_len, lat_len, lon_len)), mask=mask)
        res = gh.global_mean(lats, lons, arr)
        self.assertTrue(np.allclose(res, np.ones((time_len,))))
        # now we block out some longitudes
        mask = np.stack(
            [
                np.ones(shape=(time_len, lat_len, lon_len_1)),
                np.zeros(shape=(time_len, lat_len, lon_len_2)),
                np.ones(shape=(time_len, lat_len, lon_len_3)),
            ],
            axis=2,
        )
        arr = np.ma.array(np.ones(shape=(time_len, lat_len, lon_len)), mask=mask)
        # we assert that the average of a constant field is the constant value
        self.assertTrue(
            np.allclose(gh.global_mean(lats, lons, arr), np.ones((time_len,)))
        )

    def test_get_nan_pixels(self):
        n_t = 2
        n_lats = 3
        n_lons = 4
        ref = np.zeros((n_t, n_lats, n_lons))
        ref[0, 2, 3] = np.nan
        ref[1, 1, 3] = np.nan
        mask = np.zeros((n_t, n_lats, n_lons))
        ma_arr = np.ma.array(ref, mask=mask)
        # diskless would require different filenames in different
        # tests making them depend on each other
        ds = nc.Dataset("example.nc", "w")  # ,diskless=True,persist=False)
        time = ds.createDimension("time", size=n_t)
        lat = ds.createDimension("lat", size=n_lats)
        lon = ds.createDimension("lon", size=n_lons)
        test = ds.createVariable("test", np.float64, ["time", "lat", "lon"])
        test[:, :, :] = ma_arr

        self.assertEqual(((1, 3), (2, 3)), gh.get_nan_pixels(test))
        ds.close()

    def test_get_nan_pixel_mask_3D(self):
        n_t = 2
        n_lats = 3
        n_lons = 4
        arg = np.zeros((n_t, n_lats, n_lons))
        arg[0, 2, 3] = np.nan
        arg[1, 1, 3] = np.nan
        mask = np.zeros((n_t, n_lats, n_lons), dtype=np.bool_)
        mask[:, 0, 0] = True
        ma_arr = np.ma.array(arg, mask=mask)
        # diskless would require different filenames in different
        # tests making them depend on each other
        ds = nc.Dataset("example.nc", "w")  # ,diskless=True,persist=False)
        time = ds.createDimension("time", size=n_t)
        lat = ds.createDimension("lat", size=n_lats)
        lon = ds.createDimension("lon", size=n_lons)
        test_var = ds.createVariable("test_var", np.float64, ["time", "lat", "lon"])
        test_var[:, :, :] = ma_arr

        ref_mask = np.zeros((n_lats, n_lons), dtype=np.bool_)  # 2 dimensional
        ref_mask[0, 0] = True
        ref_mask[2, 3] = True
        ref_mask[1, 3] = True

        res = gh.get_nan_pixel_mask(test_var)
        self.assertTrue((ref_mask == res).all())
        ds.close()

        ma_arr = np.ma.array(arg, mask=False)
        # diskless would require different filenames in different
        # tests making them depend on each other
        ds = nc.Dataset("example1.nc", "w")  # ,diskless=True,persist=False)
        time = ds.createDimension("time", size=n_t)
        lat = ds.createDimension("lat", size=n_lats)
        lon = ds.createDimension("lon", size=n_lons)
        test_var = ds.createVariable("test_var", np.float64, ["time", "lat", "lon"])
        test_var[:, :, :] = ma_arr
        ref_mask = np.zeros((n_lats, n_lons), dtype=np.bool_)  # 2 dimensional
        ref_mask[2, 3] = True
        ref_mask[1, 3] = True
        res = gh.get_nan_pixel_mask(test_var)
        self.assertTrue((ref_mask == res).all())
        ds.close()

    def test_get_nan_pixel_mask_4D(self):
        # 4-D-example
        n_t = 2
        n_d = 3
        n_lats = 4
        n_lons = 5
        arg = np.zeros((n_t, n_d, n_lats, n_lons))
        arg[0, 0, 2, 3] = np.nan
        arg[1, 2, 1, 3] = np.nan
        mask = np.zeros((n_t, n_d, n_lats, n_lons), dtype=np.bool_)
        mask[:, 0, 0, 0] = True
        ma_arr = np.ma.array(arg, mask=mask)
        # diskless would require different filenames in different
        # tests making them depend on each other
        ds = nc.Dataset("example.nc", "w")  # ,diskless=True,persist=False)
        time = ds.createDimension("time", size=n_t)
        depth = ds.createDimension("depth", size=n_d)
        lat = ds.createDimension("lat", size=n_lats)
        lon = ds.createDimension("lon", size=n_lons)
        test_var = ds.createVariable(
            "test_var", np.float64, ["time", "depth", "lat", "lon"]
        )
        test_var[:, :, :, :] = ma_arr

        ref_mask = np.zeros((n_lats, n_lons), dtype=np.bool_)  # 2 dimensional
        # NAN contribution
        ref_mask[2, 3] = True
        ref_mask[1, 3] = True
        # mask contribution
        ref_mask[0, 0] = True

        res = gh.get_nan_pixel_mask(test_var)
        self.assertTrue(np.all(ref_mask == res))
        ds.close()

    def test_covariance_allocation_sum(self):
        ts = np.linspace(0, 2 * np.pi, 10000)
        x1s = np.sin(ts)
        x2s = np.ones_like(ts)
        ys = 1 * x1s + x2s

        res = gh.product_attribution(ys, x1s, x2s)
        print(np.sum(res))

    def test_covariance_allocation_product(self):
        ts = np.linspace(0, 2 * np.pi, 10000)
        x1s = np.sin(ts) + 1.1  # must be positive for application of log
        x2s = np.ones_like(ts)
        ys = x1s * x2s
        res = gh.product_attribution(ys, x1s, x2s)
        self.assertEqual(res, (1, 0))

        # divide equally
        x1s = np.sin(ts) + 1.1  # must be positive for application of log
        x2s = x1s
        ys = x1s * x2s
        res = gh.product_attribution(ys, x1s, x2s)
        self.assertTrue(np.allclose(res, (0.5, 0.5)))

    def test_covariance_allocation_product_temporal_variation(self):
        # now we create a test case with a  possibly scaled \fac and
        # check its influence on the attribution
        # first if xis==base_rates
        ts = np.linspace(0, 2 * np.pi, 10000)
        xis = np.sin(ts) + 1.1  # must be positive for application of log
        base_rates = xis
        rts = xis * base_rates
        res = gh.product_attribution(rts, xis, base_rates)

        fac = 5
        xis_ = 1 / fac * xis
        base_rates_ = fac * base_rates
        res = gh.product_attribution(rts, xis_, base_rates_)
        self.assertTrue(np.allclose(np.sum(res), 1))

        xis_ = 1 / fac * np.sin(ts) + 1.1  # must be positive for application of log
        res_ = gh.product_attribution(rts, xis, base_rates)
        self.assertTrue(np.allclose(res, res_))

        # now we create a test case with a  possibly scaled \fac and
        # check its influence on the attribution
        # first if xis==base_rates
        xis = np.sin(ts) + 1.1  # must be positive for application of log
        base_rates = ts / ts[-1] + 1
        rts = xis * base_rates
        res = gh.product_attribution(rts, xis, base_rates)
        print(res)
        self.assertTrue(np.allclose(np.sum(res), 1))

        fac = 5
        xis_ = 1 / fac * xis
        base_rates_ = fac * base_rates
        rts_ = xis_ * base_rates_
        res_ = gh.product_attribution(rts_, xis_, base_rates_)
        print(res)
        self.assertTrue(np.allclose(res, res_))

        xis_ = 1 / fac * np.sin(ts) + 1.1  # must be positive for application of log
        res_ = gh.product_attribution(rts, xis, base_rates)
        self.assertTrue(np.allclose(res, res_))

    def test_covariance_allocation_product_inter_model_variation(self):
        # now we create a test case with two models where we 
        # provide as scaled xi and baseline version of the second model  
        # that is otherwise identical
        # and check the influence on the attribution
        # first if x1s == x2s
        xi1 = 1
        xi2 = 2
        xis = np.array([xi1, xi2])  
        br1 = 1
        br2 = 2
        base_rates = np.array([br1, br2])
        rts = xis * base_rates
        res = gh.product_attribution(rts, xis, base_rates)

        # now we assume that we change  xi2 and br2 to xi2_ and br_
        fac = 5
        xi2_ = 1 / fac * xi2
        xis_ = np.array([xi1, xi2_])
        br2_ = fac * br2
        base_rates_ = np.array([br1, br2_])
        res_ = gh.product_attribution(rts, xis_, base_rates_)
        self.assertTrue(np.allclose(np.sum(res), 1))
        print(res)
        print(res_)
        #self.assertTrue(np.allclose(res, res_))

    def test_globalmean_var(self):
        # we create data similar to cmip6 and trendy. These are masked arrays:
        # the lat lon combinations that have no landpoints are marked with
        # False in the mask which is a boolean array of the same shape as
        # the actual value array) the
        #
        # The following commented lines show how to get the
        # lats and longs are as in jon_yibs
        # lat_name = 'latitude'
        # lon_name = 'longitude'
        # lats=nds.variables[lat_name].__array__()
        # longs=nds.variables[lon_name].__array__()
        # target_var_names=set(nds.variables.keys()).difference([lat_name,lon_name,'time'])
        # arr = target_var_names[0]
        lat_len_1 = 3
        lat_len_2 = 3
        lat_len_3 = 3
        lat_len = sum([lat_len_1, lat_len_2, lat_len_3])
        lon_len_1 = 120
        lon_len_2 = 120
        lon_len_3 = 120
        lon_len = sum([lon_len_1, lon_len_2, lon_len_3])
        time_len = 3
        lats = np.ma.masked_array(np.linspace(-90, 90, lat_len))
        lons = np.ma.masked_array(np.linspace(-179.5, 179.5, lon_len))
        arr = np.ma.ones(shape=(time_len, lat_len, lon_len))
        vn = "test_var"
        trunk = "{}_1".format(vn)
        ds = nc.Dataset(Path("{}.nc".format(trunk)), "w", diskless=True, persist=False)
        time = ds.createDimension("time", size=time_len)
        lat = ds.createDimension("lat", size=lat_len)
        lon = ds.createDimension("lon", size=lon_len)
        test_var = ds.createVariable(vn, np.float64, ["time", "lat", "lon"])
        test_var[:, :, :] = arr

        res_arr = gh.global_mean(lats, lons, arr)
        res_var = gh.global_mean_var(
            lats, lons, np.zeros((lat_len, lon_len), dtype=np.bool_), arr
        )
        self.assertTrue((res_arr == res_var).all())
        cache_path = Path("{}_gm.nc".format(trunk))

        gh.write_global_mean_cache(cache_path, res_var, vn)
        # def get_cached_global_mean(gm_path, vn):
        #    return nc.Dataset(str(gm_path)).variables[vn].__array__()

        res_cache = gh.get_cached_global_mean(cache_path, vn)
        # from IPython import embed;embed()
        self.assertTrue((res_cache == res_var).all())

        ds.close()

        # now we test the same properties with a masked array
        # as e.g the CMIP6 files that Kostia uses
        # dataPath = Path(conf_dict['dataPath'])
        # ds=nc.Dataset('cLeaf_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc')
        # pwd
        # dataPath
        # ds=nc.Dataset(dataPath.joinpath('cLeaf_Lmon_MIROC-ES2L_1pctCO2-bgc_r1i1p1f2_gn_185001-199912.nc'))
        # cLeaf = ds.variables['cLeaf'].__array__()
        # in this array many fields are masked

        # note that in python boolean values are equivalent to integers
        # 1 = True , 0 = False
        mask = np.stack(
            [
                np.ones(shape=(time_len, lat_len_1, lon_len)),
                np.zeros(shape=(time_len, lat_len_2, lon_len)),
                np.ones(shape=(time_len, lat_len_3, lon_len)),
            ],
            axis=1,
        )
        arr = np.ma.array(np.ones(shape=(time_len, lat_len, lon_len)), mask=mask)
        vn = "test_var"
        trunk = "{}_2".format(vn)
        ds = nc.Dataset(Path("{}.nc".format(trunk)), "w", diskless=True, persist=False)
        time = ds.createDimension("time", size=time_len)
        lat = ds.createDimension("lat", size=lat_len)
        lon = ds.createDimension("lon", size=lon_len)
        test_var = ds.createVariable(vn, np.float64, ["time", "lat", "lon"])
        test_var[:, :, :] = arr
        res_arr = gh.global_mean(lats, lons, arr)
        res_var = gh.global_mean_var(lats, lons, mask.sum(axis=0), arr)
        self.assertTrue((res_arr == res_var).all())
        cache_path = Path("{}_gm.nc".format(trunk))

        gh.write_global_mean_cache(cache_path, res_var, vn)
        res_cache = gh.get_cached_global_mean(cache_path, vn)
        self.assertTrue((res_cache == res_var).all())

        ds.close()

        # now we block out some longitudes
        mask = np.stack(
            [
                np.ones(shape=(time_len, lat_len, lon_len_1)),
                np.zeros(shape=(time_len, lat_len, lon_len_2)),
                np.ones(shape=(time_len, lat_len, lon_len_3)),
            ],
            axis=2,
        )
        arr = np.ma.array(np.ones(shape=(time_len, lat_len, lon_len)), mask=mask)
        vn = "test_var"
        trunk = "{}_3".format(vn)
        ds = nc.Dataset(Path("{}.nc".format(trunk)), "w", diskless=True, persist=False)
        time = ds.createDimension("time", size=time_len)
        lat = ds.createDimension("lat", size=lat_len)
        lon = ds.createDimension("lon", size=lon_len)
        test_var = ds.createVariable(vn, np.float64, ["time", "lat", "lon"])
        test_var[:, :, :] = arr
        res_arr = gh.global_mean(lats, lons, arr)
        res_var = gh.global_mean_var(lats, lons, mask.sum(axis=0), arr)
        self.assertTrue((res_arr == res_var).all())
        cache_path = Path("{}_gm.nc".format(trunk))

        gh.write_global_mean_cache(cache_path, res_var, vn)
        res_cache = gh.get_cached_global_mean(cache_path, vn)
        self.assertTrue((res_cache == res_var).all())

        ds.close()

    # def test_common_mask(self):
    #    def i2c_maker(lat_0,lon_0,step_lat,step_lon):
    #        def i2c(i_lat,i_lon):
    #            return(
    #                lat_0+step_lat*i_lat,
    #                lon_0+step_lon*i_lon
    #            )
    #        return i2c
    #    step_lat1,step_lon1=(1,1)
    #    i2c_1=i2c_maker(
    #        step_lat1/2,
    #        step_lon1/2,
    #        step_lat1,
    #        step_lon1,
    #    )
    #
    #
    #    step_lat2,step_lon2=(.7,1)
    #    i2c_2=i2c_maker(
    #        .5+step_lat2/2,
    #        step_lon2/2,
    #        step_lat2,
    #        step_lon2,
    #    )

    #    m1=np.array([0,1,0,0],dtype=np.bool_).reshape(4,1)
    #    m2=np.array([1,1,0,0,0],dtype=np.bool_).reshape(5,1)
    #
    #    m_res=gh.combine_masks([m1,m2],[i2c_1,i2c_2])
    #
    #    m_ref=[
    #            (1.0, 2.0, 0.0, 1.0),
    #            (0.5, 1.2, 0.0, 1.0),
    #            (1.2, 1.8999999999999997, 0.0, 1.0)
    #    ]
    #    print(m_ref==m_res)
    #    self.assertTrue(m_res==m_ref)

    #

    #    common_mask=m_res
    #    m_ref=np.array([1,1,0,0]).reshape(4,1)
    #    res= gh.project(
    #        m1.shape,
    #        i2c_1,
    #        common_mask
    #    )
    #    self.assertTrue(
    #            (
    #                res==m_ref
    #            ).all()
    #    )

    #
    # def test_open_interval_intersect(self):
    #    self.assertTrue(gh.open_interval_intersect((0,1),(0.5,1.5)))
    #    self.assertTrue(gh.open_interval_intersect((0.5,1.5),(0,1)))
    #    self.assertTrue(gh.open_interval_intersect((0,1),(0,1)))
    #    self.assertTrue(gh.open_interval_intersect((0,1),(-1,2)))
    #    self.assertFalse(gh.open_interval_intersect((0,1),(1,2)))
    #    self.assertFalse(gh.open_interval_intersect((1,2),(0,1)))

    # def test_pixel_intersect(self):
    #    self.assertTrue(
    #        gh.pixel_intersect(
    #            gh.boundaries(0,1,0,1),
    #            gh.boundaries(0,1,0,1)
    #        )
    #    )
    #    self.assertTrue(
    #        gh.pixel_intersect(
    #            gh.boundaries(0,1,0,1),
    #            gh.boundaries(.5,1.5,0,1)
    #        )
    #    )
    #    self.assertTrue(
    #        gh.pixel_intersect(
    #            gh.boundaries(0,1,0,1),
    #            gh.boundaries(.5,.7,.5,.7)
    #        )
    #    )
    #    self.assertFalse(
    #        gh.pixel_intersect(
    #            gh.boundaries(0,1,0,1),
    #            gh.boundaries(1,1.5,0,1)
    #        )
    #    )
    #    self.assertFalse(
    #        gh.pixel_intersect(
    #            gh.boundaries(0,1,0,1),
    #            gh.boundaries(1,1.5,1,1.5)
    #        )
    #    )
    def test_transform_maker(self):
        step_lat = 10
        step_lon = 10
        lat_0 = 5
        lon_0 = 5
        tr = gh.transform_maker(
            lat_0,
            lon_0,
            step_lat,
            step_lon,
        )
        n_lat = int(180 / step_lat)
        n_lon = int(360 / step_lon)
        # check that i2lat rejects arguments >n_lat-1
        with self.assertRaises(IndexError):
            lat = tr.i2lat(n_lat)

        # check that i2lon rejects arguments >n_lon-1
        with self.assertRaises(IndexError):
            lon = tr.i2lon(n_lon)

        ## test the inverse property (under too large indices)
        # for i in range(0,n_lat):
        #    lat=tr.i2lat(i)
        #    ii=tr.lat2i(lat)
        #    #print("i={i},lat={lat},ii={ii}".format(i=i,lat=lat,ii=ii))
        #    self.assertEqual(ii,i)

        # for i in range(0,n_lon):
        #    lon=tr.i2lon(i)
        #    ii=tr.lon2i(lon)
        #    #print("i={i},lon={lon},ii={ii}".format(i=i,lon=lon,ii=ii))
        #    self.assertEqual(tr.lon2i(tr.i2lon(i)),i)

    def test_project_2_self_b(self):
        # step_lat = 90
        # lat_0 = -90+45
        step_lat = -90
        lat_0 = 90 - 45
        step_lon = 90
        lon_0 = -180 + 45
        itr = gh.transform_maker(
            lat_0,
            lon_0,
            step_lat,
            step_lon,
        )
        # here we use the identical transformation
        ctr = gh.identicalTransformers()
        sym_tr = gh.SymTransformers(itr=itr, ctr=ctr)
        cm_1 = gh.CoordMask(
            np.array([[0, 1, 0, 0], [0, 0, 0, 0]], dtype=np.bool_).reshape(2, 4), sym_tr
        )

        # self projection
        res = gh.project_2(cm_1, cm_1)
        f = plt.figure()
        ax = f.add_subplot(3, 1, 1)
        cm_1.plot_dots(ax)
        ax = f.add_subplot(3, 1, 2)
        res.plot_dots(ax)
        f.savefig("cm_1.pdf")
        self.assertTrue((res.index_mask == cm_1.index_mask).all())

    def test_project_2_self_c(self):
        # step_lat = 60
        # lat_0 = -90+30

        step_lat = -60
        lat_0 = 90 - 30

        step_lon = 90
        lon_0 = -180 + 45
        itr = gh.transform_maker(
            lat_0,
            lon_0,
            step_lat,
            step_lon,
        )
        # here we use the identical transformation
        ctr = gh.identicalTransformers()
        sym_tr = gh.SymTransformers(itr=itr, ctr=ctr)
        n_lat = int(180 / step_lat)
        n_lon = int(360 / step_lon)
        cm_1 = gh.CoordMask(
            np.array(
                [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]], dtype=np.bool_
            ).reshape(n_lat, n_lon),
            sym_tr,
        )

        # self projection
        res = gh.project_2(cm_1, cm_1)
        f = plt.figure()
        ax = f.add_subplot(2, 1, 1)
        cm_1.plot_dots(ax)
        ax = f.add_subplot(2, 1, 2)
        res.plot_dots(ax)
        f.savefig("cm_1.pdf")
        self.assertTrue((res.index_mask == cm_1.index_mask).all())

    def test_project_2_higher_res_target(self):
        # here we use the identical transformation
        ctr = gh.identicalTransformers()
        step_lat = 60
        step_lon = 90
        # lat_0 = 0
        lat_0 = -90 + 30
        lon_0 = -180 + 45
        itr_1 = gh.transform_maker(
            lat_0,
            lon_0,
            step_lat,
            step_lon,
        )
        sym_tr_1 = gh.SymTransformers(itr=itr_1, ctr=ctr)
        n_lat = int(180 / step_lat)
        n_lon = int(360 / step_lon)
        cm_1 = gh.CoordMask(
            np.array(
                [
                    [0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                ],
                dtype=np.bool_,
            ).reshape(n_lat, n_lon),
            sym_tr_1,
        )

        step_lat = 45
        step_lon = 60
        # lat_0 = 0
        lat_0 = -90 + step_lat / 2.0
        lon_0 = -180 + step_lon / 2.0
        itr_2 = gh.transform_maker(
            lat_0,
            lon_0,
            step_lat,
            step_lon,
        )
        # here we use the identical transformation
        sym_tr_2 = gh.SymTransformers(itr=itr_2, ctr=ctr)
        n_lat = int(180 / step_lat)
        n_lon = int(360 / step_lon)
        cm_2 = gh.CoordMask(
            np.array(
                [
                    [0, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ],
                dtype=np.bool_,
            ).reshape(n_lat, n_lon),
            sym_tr_2,
        )
        res = gh.project_2(target=cm_2, source=cm_1)
        f = plt.figure()
        ax = f.add_subplot(3, 1, 1)
        cm_1.plot_dots(ax)
        ax = f.add_subplot(3, 1, 2)
        cm_2.plot_dots(ax)
        ax = f.add_subplot(3, 1, 3)
        res.plot_dots(ax)
        f.savefig("cms.pdf")
        # raise "the plot shows that something goes definitely wrong"
        # m_ref=np.array([1,1,0,0]).reshape(4,1)
        # self.assertTrue((res==m_ref).all())

    def test_project_2_lower_res_target(self):
        # here we use the identical transformation
        ctr = gh.identicalTransformers()
        step_lat = 45
        step_lon = 60
        # lat_0 = 0
        lat_0 = -90 + step_lat / 2.0
        lon_0 = -180 + step_lon / 2.0
        itr_1 = gh.transform_maker(
            lat_0,
            lon_0,
            step_lat,
            step_lon,
        )
        # here we use the identical transformation
        sym_tr_1 = gh.SymTransformers(itr=itr_1, ctr=ctr)
        n_lat = int(180 / step_lat)
        n_lon = int(360 / step_lon)
        cm_1 = gh.CoordMask(
            np.array(
                [
                    [1, 0, 0, 0, 0, 0],
                    [1, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ],
                dtype=np.bool_,
            ).reshape(n_lat, n_lon),
            sym_tr_1,
        )
        step_lat = 60
        step_lon = 90
        # lat_0 = 0
        lat_0 = -90 + 30
        lon_0 = -180 + 45
        itr_2 = gh.transform_maker(
            lat_0,
            lon_0,
            step_lat,
            step_lon,
        )
        sym_tr_2 = gh.SymTransformers(itr=itr_2, ctr=ctr)
        n_lat = int(180 / step_lat)
        n_lon = int(360 / step_lon)
        cm_2 = gh.CoordMask(
            np.array(
                [
                    [0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                ],
                dtype=np.bool_,
            ).reshape(n_lat, n_lon),
            sym_tr_2,
        )

        res = gh.project_2(target=cm_2, source=cm_1)
        f = plt.figure()
        ax = f.add_subplot(3, 1, 1)
        cm_1.plot_dots(ax)
        ax = f.add_subplot(3, 1, 2)
        cm_2.plot_dots(ax)
        ax = f.add_subplot(3, 1, 3)
        res.plot_dots(ax)
        f.savefig("cms.pdf")
        # raise "the plot shows that something goes definitely wrong"
        # m_ref=np.array([1,1,0,0]).reshape(4,1)
        # self.assertTrue((res==m_ref).all())

        ## shifted and different size
        #
        # i2lat_2, i2lon_2, lat2i_2,lon2i_2= gh.transform_maker(
        #    .5+step_lat2/2,
        #    step_lon2/2,
        # )

        # cm_2=gh.CoordMask(
        #    np.array([1,1,0,0,0],dtype=np.bool_).reshape(5,1),
        #    45+(180/5)/2,
        #    (360/1)/2
        # )
        #

        # res=gh.project_2(
        #    target=cm_1,
        #    source=cm_2
        # )
        # print(res)
        # m_ref=np.array([1,1,0,0]).reshape(4,1)

        ##self.assertTrue(
        ##        (
        ##            res==m_ref
        ##        ).all()
        ##)

    def test_InfiniteIterator(self):
        I = (np.array([1, 1]).reshape(2, 1),)
        k = 0.5

        def f(i, X):
            return X + (I - k * X)

        itr = gh.InfiniteIterator(x0=np.array([1, 1]).reshape(2, 1), func=f)

        # make sure that [] has no side effects

        # first we get a single value
        # the results will be tuples of length 1
        # of 2,1 arrays
        result_1 = itr[0]
        result_2 = itr[0]
        self.assertTrue(np.all((np.all(result_1[0] == result_2[0]))))

        # tuples of arrays
        results_1 = np.stack(itr[0:10])
        results_2 = np.stack(itr[0:10])
        self.assertTrue(np.all(results_1 == results_2))

    def test_TraceTupleIterator(self):
        def B_func(it, X):
            return -0.5 * np.eye(2)

        def I_func(it, X):
            return 2 * np.ones(shape=(2, 1))

        X_0 = np.array([2, 1]).reshape(2, 1)
        traced_functions = {
            # in our case X has
            "test": lambda i, x1, x2: x1
            ** 2
        }
        trace_tuple_instance = gh.make_trace_tuple_func(traced_functions)
        V_init = trace_tuple_instance(
            X_0,
            # in Yiqi's nomenclature: dx/dt=I-Bx
            # instead of           : dx/dt=I+Bx
            # as assumed by B_u_func
            -B_func(0, X_0),
            I_func(0, X_0),
            0,
        )

        # define the function with V_{i+1}=f(i,V_i)
        def f(it: int, V: gh.TraceTuple) -> gh.TraceTuple:
            X = V.X
            I = I_func(it, X)
            B = B_func(it, X)
            X_new = X + I + B @ X
            return trace_tuple_instance(X_new, -B, I, it)

        itr = gh.TraceTupleIterator(x0=V_init, func=f)
        # assert side effect free [ ] application
        results_1 = itr[0:10]
        results_2 = itr[0:10]
        self.assertTrue(
            np.all(
                tuple(
                    np.all(
                        results_1.__getattribute__(name)
                        == results_2.__getattribute__(name)
                    )
                    for name in results_1._fields
                )
            )
        )
        # test correct averaging
        parts = gh.partitions(0, 10, 2)
        res = itr.averaged_values(parts)
        ref = results_1.averages(parts)
        self.assertTrue(res == ref)

    def test_read_or_create(self):
        cachePath = Path("cache.nc")
        var_name = "test"

        def caw(path):
            sleep(2)
            n = 5
            ds = nc.Dataset(path, "w", persist=True)
            lat = ds.createDimension("lat", size=n)
            lon = ds.createDimension("lon", size=n)
            test = ds.createVariable(var_name, np.float64, ["lat", "lon"])
            var = np.diag(np.arange(0, n))
            test[:, :] = var
            return var

        def r(path):
            ds = nc.Dataset(path)
            return ds.variables[var_name][:, :].data

        path = Path("cache", "test.nc")
        before = time()
        res_1 = gh.read_or_create(path=path, create_and_write=caw, read=r)
        after_1 = time()
        res_2 = gh.read_or_create(path=path, create_and_write=caw, read=r)
        after_2 = time()
        self.assertTrue((res_1 == res_2).all())
        self.assertTrue(after_1 - before > after_2 - after_1)
