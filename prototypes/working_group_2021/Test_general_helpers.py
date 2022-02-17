import numpy as np
from unittest import TestCase, skip
from testinfrastructure.InDirTest import InDirTest
from general_helpers import (
        day_2_month_index, 
        month_2_day_index,
        months_by_day_arr,
        TimeStepIterator2,
        respiration_from_compartmental_matrix,
        global_mean,
        make_pixel_area_on_unit_spehre
)
class Test_general_helpers(InDirTest):
    def test_month_2_day_index(self):
        self.assertEqual(
                month_2_day_index([0]),
                [0]
        ) 
        self.assertEqual(
                month_2_day_index([1]),
                [31]
        ) 
        self.assertEqual(
                month_2_day_index([2]),
                [59]
        ) 
        self.assertEqual(
                month_2_day_index([3]),
                [90]
        ) 
        self.assertEqual(
                month_2_day_index([1,3]),
                [31,90]
        ) 
    
    def test_day_2_month_index(self):
        # note that days are counted from zero so day 30 is January 31.
        self.assertEqual(day_2_month_index( 0), 0) 
        self.assertEqual(day_2_month_index(30), 0) 
        self.assertEqual(day_2_month_index(31), 1) 
        self.assertEqual(day_2_month_index(60), 2) 


    def test_pixel_area_on_unit_sphere(self):
        # we test that the pixel areas for a whole
        # grid sum up to 1
        #lat_len=181
        #lon_len=360
        lat_len=19
        lon_len=360
        lats=np.ma.masked_array(np.linspace(-90,90,lat_len))
        lons=np.ma.masked_array(np.linspace(-179.5,179.5,lon_len))
        delta_lat=(lats.max()- lats.min())/(len(lats)-1)
        delta_lon=(lons.max() -lons.min())/(len(lons)-1)

        puaf_sym= make_pixel_area_on_unit_spehre(delta_lat, delta_lon,sym=True)
        puaf= make_pixel_area_on_unit_spehre(delta_lat, delta_lon)

        # assert identical values
        lw=np.array(
            [
                puaf(lats[lat_ind])
                for lat_ind in range(len(lats))
            ],
            dtype=np.float64
        )
        lws=np.array(
            [
                puaf_sym(lats[lat_ind])
                for lat_ind in range(len(lats))
            ],
            dtype=np.float64
        )
        # check for equivalence of our manually solved integral with sympy
        # (found typos this way....)
        self.assertTrue(np.allclose(lw,lws))

        # check for symetrie (the weighs for pos and negative latitudes should
        # be the same
        pres = np.array([puaf(lat)  for lat in lats ],dtype=np.float64)
        nres = np.array([puaf(-lat)  for lat in lats ],dtype=np.float64)
        self.assertTrue(np.allclose(pres,nres))
       
        
        # having established the equality of the two ways to compute the
        # weights we check that they add up to 1
        self.assertTrue(
            np.allclose(
                np.sum(
                    [
                        np.sum(
                            [   
                                puaf(lats[lat_ind])
                                for lon_ind in range(len(lons))    
                            ]
                        )
                        for lat_ind in range(len(lats))    
                    ],
                    dtype=np.float64
                ),
                4*np.pi #surface area of the unit sphere
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
        lat_len_1=3
        lat_len_2=3
        lat_len_3=3
        lat_len=sum(
            [
                lat_len_1,
                lat_len_2,
                lat_len_3
            ]
        )
        lon_len_1=120
        lon_len_2=120
        lon_len_3=120
        lon_len=sum(
            [
                lon_len_1,
                lon_len_2,
                lon_len_3
            ]
        )
        time_len=3
        lats=np.ma.masked_array(np.linspace(-90,90,lat_len))
        lons=np.ma.masked_array(np.linspace(-179.5,179.5,lon_len))
        arr=np.ma.ones(shape=(time_len,lat_len,lon_len))
        
        res = global_mean(lats,lons,arr)
        self.assertEqual(
            res.shape,
            (time_len,)
        )

        # we assert that the average of a constant field is the constant value
        self.assertTrue(
            np.allclose(
                res,
                np.ones((time_len,))
            )
        )

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
        mask=np.stack(
            [
                np.ones(shape=(time_len,lat_len_1,lon_len)),
                np.zeros(shape=(time_len,lat_len_2,lon_len)),
                np.ones(shape=(time_len,lat_len_3,lon_len))
            ],
            axis=1
        )
        arr=np.ma.array(
            np.ones(shape=(time_len,lat_len,lon_len)),
            mask=mask
        )
        res = global_mean(lats,lons,arr)
        self.assertTrue(
            np.allclose(
                res,
                np.ones((time_len,))
            )
        )
        # now we block out some longitudes
        mask=np.stack(
            [
                np.ones(shape=(time_len,lat_len,lon_len_1)),
                np.zeros(shape=(time_len,lat_len,lon_len_2)),
                np.ones(shape=(time_len,lat_len,lon_len_3))
            ],
            axis=2
        )
        arr=np.ma.array(
            np.ones(shape=(time_len,lat_len,lon_len)),
            mask=mask
        )
        # we assert that the average of a constant field is the constant value
        self.assertTrue(
            np.allclose(
                global_mean(lats,lons,arr),
                np.ones((time_len,))
            )
        )

     
