import unittest
import numpy as np

from bgc_md2.ModelDataObject import ModelDataObjectException,\
                                    ModelDataObject,\
                                    readVariable,\
                                    StockDensityVariable2StockVariable,\
                                    getStockVariable_from_Density,\
                                    FluxRateDensityVariable2FluxRateVariable,\
                                    FluxRateVariable2FluxVariable,\
                                    getFluxVariable_from_DensityRate,\
                                    getFluxVariable_from_Rate
from bgc_md2.Variable import Variable, StockVariable, FluxVariable
from example_MDOs import MDO_3pools_3layers_2x2,\
                         MDO_3pools_nolayers_nogrid,\
                         MDO_3pools_nolayers_2x2,\
                         MDO_3pools_3layers_nogrid,\
                         MDO_3pools_3layers_2x2_discretizable


class TestModelDataObject(unittest.TestCase):

    def setUp(self):
        self.mdo = MDO_3pools_3layers_2x2()
        self.lat_arr = np.arange(len(self.mdo.dataset['lat'][:]))
        self.lon_arr = np.arange(len(self.mdo.dataset['lon'][:]))

        ## incomplete models
        mdos_inc = [
            MDO_3pools_nolayers_nogrid(),
            MDO_3pools_nolayers_2x2(),
            MDO_3pools_3layers_nogrid()
        ]

        mdos_inc[0].nr_layers = None
        mdos_inc[0].lat_arr   = None
        mdos_inc[0].lon_arr   = None
        
        mdos_inc[1].nr_layers = None
        mdos_inc[1].lat_arr   = np.arange(2)
        mdos_inc[1].lon_arr   = np.arange(2)
     
        mdos_inc[2].nr_layers = 3
        mdos_inc[2].lat_arr   = None
        mdos_inc[2].lon_arr   = None

        self.mdos_inc = mdos_inc

        ## discretizable MDOs
        self.mdo_discretizable = MDO_3pools_3layers_2x2_discretizable()


    def tearDown(self):
        self.mdo.dataset.close()
        for mdo in self.mdos_inc:
            mdo.dataset.close()
        self.mdo_discretizable.dataset.close()


    def test_readVariable(self):
        ds = self.mdo.dataset

        var = readVariable(
            dataset       = ds,
            variable_name = 'CWDC',
            min_time      = 0,
            max_time      = 3,
            nr_layers     = 2,
            lat_arr       = self.lat_arr,
            lon_arr       = [1],
            data_shift    = 1
        )
        
        res2 = [[[[13.],
                  [15.]],

                 [[17.],
                  [19.]]],


                [[[25.],
                  [27.]],

                 [[29.],
                  [31.]]]]
      
        self.assertTrue(isinstance(var, Variable)) 
        self.assertTrue(np.all(var.data == res2))
        self.assertEqual(var.unit, 'g/m^3')
 
        sv = readVariable(
            ReturnClass   = StockVariable,
            dataset       = ds,
            variable_name = 'CWDC',
            min_time      = 0,
            max_time      = 3,
            nr_layers     = 2,
            lat_arr       = self.lat_arr,
            lon_arr       = [1],
            data_shift    = 1
        )
        
        self.assertTrue(isinstance(sv, StockVariable)) 

        ##### no layers, no grid #####

        mdo = self.mdos_inc[0]
        ds = mdo.dataset
        ms = mdo.model_structure

        sv = readVariable(
            ReturnClass   = StockVariable,
            dataset       = ds,
            variable_name = 'CWDC',
            nr_layers     = ms.get_nr_layers('CWD'),
            min_time      = 0,
            max_time      = 3,
            data_shift    = 1
        )
        self.assertTrue(isinstance(sv, StockVariable))
        self.assertEqual(sv.unit, 'g/m^2')
        self.assertTrue(np.all(sv.data.shape == np.array([2,1,1,1])))

        ##### no layers, 2x2 #####

        mdo = self.mdos_inc[1]
        ds = mdo.dataset
        ms = mdo.model_structure

        sv = readVariable(
            ReturnClass   = StockVariable,
            dataset       = ds,
            variable_name = 'CWDC',
            nr_layers     = ms.get_nr_layers('CWD'),
            lat_arr       = mdo.lat_arr,
            lon_arr       = mdo.lon_arr,
            min_time      = 0,
            max_time      = 3,
            data_shift    = 1
        )
        self.assertTrue(isinstance(sv, StockVariable))
        self.assertEqual(sv.unit, 'g/m^2')
        self.assertTrue(np.all(sv.data.shape == np.array([2,1,2,2])))

        ##### 3 layers, no grid #####

        mdo = self.mdos_inc[2]
        ds = mdo.dataset
        ms = mdo.model_structure

        sv = readVariable(
            ReturnClass   = StockVariable,
            dataset       = ds,
            variable_name = 'CWDC',
            nr_layers     = ms.get_nr_layers('CWD'),
            min_time      = 0,
            max_time      = 3,
            data_shift    = 1
        )
        self.assertTrue(isinstance(sv, StockVariable))
        self.assertEqual(sv.unit, 'g/m^3')
        self.assertTrue(np.all(sv.data.shape == np.array([2,3,1,1])))


    def test_StockDensityVariable2StockVariable(self):
        ds = self.mdo.dataset
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')

        sdv = readVariable(
            ReturnClass   = StockVariable,
            dataset       = ds,
            variable_name = 'CWDC',
            min_time      = 0,
            max_time      = 100,
            nr_layers     = nr_layers,
            lat_arr       = self.lat_arr,
            lon_arr       = self.lon_arr,
            data_shift    = 0
        )
        dz = Variable(
            data = ds['dz'][:],
            unit = ds['dz'].units
        )

        sv = StockDensityVariable2StockVariable(sdv, dz)
        self.assertTrue(isinstance(sv, StockVariable))
        self.assertEqual(sv.unit, 'm-2.kg')
        res = sv.data[4,...] * 1000
        res2 = [[[ 48.,  49.],
                 [ 50.,  51.]],

                [[104., 106.],
                 [108., 110.]],

                [[168., 171.],
                 [174., 177.]]]
        self.assertTrue(np.allclose(res, res2))


    def test_getStockVariable_from_Density(self):
        mdo = self.mdo
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')

        sv = getStockVariable_from_Density(
            mdo = mdo,
            variable_name = 'CWDC',
            min_time      = 0,
            max_time      = mdo.max_time,
            nr_layers     = nr_layers,
            lat_arr       = self.lat_arr,
            lon_arr       = self.lon_arr,
            dz            = mdo.get_dz('CWD'),
            data_shift    = 0
        )

        self.assertTrue(isinstance(sv, StockVariable))
        self.assertEqual(sv.unit, mdo.stock_unit)

        res = sv.data[2,...]
        res2 = [[[ 48.,   49.],
                 [ 50.,   51.]],

                [[104.,  106.],
                 [108.,  110.]],

                [[168., 171.],
                 [174., 177.]]]
        self.assertTrue(np.allclose(res, res2))


    def test_FluxRateDensityVariable2FluxRateVariable(self):
        ds = self.mdo.dataset
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')

        frdv = readVariable(
            ReturnClass   = FluxVariable,
            dataset       = ds,
            variable_name = 'fire_CWD',
            min_time      = 0,
            max_time      = 100,
            nr_layers     = nr_layers,
            lat_arr       = self.lat_arr,
            lon_arr       = self.lon_arr,
            data_shift    = 1
        )
        dz = Variable(
            data = ds['dz'][:],
            unit = ds['dz'].units
        )

        frv = FluxRateDensityVariable2FluxRateVariable(frdv, dz)

        self.assertTrue(isinstance(frv, FluxVariable))
        self.assertEqual(frv.unit, 'm-2.kg.s-1')
        res = frv.data[3,...] * 1000
        res2 = np.array(
               [[[ 48.,  49.],
                 [ 50.,  51.]],

                [[104., 106.],
                 [108., 110.]],

                [[168., 171.],
                 [174., 177.]]]
        )
        self.assertTrue(np.allclose(res, res2*1e-03))


    def test_FluxRateVariable2FluxVariable(self):
        ds = self.mdo.dataset
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')
        
        frdv = readVariable(
            ReturnClass   = FluxVariable,
            dataset       = ds,
            variable_name = 'fire_CWD',
            min_time      = 0,
            max_time      = 100,
            nr_layers     = nr_layers,
            lat_arr       = self.lat_arr,
            lon_arr       = self.lon_arr,
            data_shift    = 1
        )
        dz = Variable(
            data = ds['dz'][:],
            unit = ds['dz'].units
        )
        frv = FluxRateDensityVariable2FluxRateVariable(frdv, dz)

        time = Variable(
            data = ds['time'][:],
            unit = ds['time'].units
        )
        fv = FluxRateVariable2FluxVariable(frv, time)
    
        self.assertTrue(isinstance(fv, FluxVariable))
        self.assertEqual(fv.unit, 'm-2.kg')
        res = fv.data[3,...]  
        res2 = 86.4 * np.array(
               [[[ 48.,  49.],
                 [ 50.,  51.]],

                [[104., 106.],
                 [108., 110.]],

                [[168., 171.],
                 [174., 177.]]]
        )
        self.assertTrue(np.allclose(res, res2*1e-03))


    def test_getFluxVariable_from_DensityRate(self):
        mdo = self.mdo
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')

        fv = getFluxVariable_from_DensityRate(
            mdo = mdo,
            variable_name = 'fire_CWD',
            nr_layers     = nr_layers,
            lat_arr       = self.lat_arr,
            lon_arr       = self.lon_arr,
            dz            = mdo.get_dz('CWD'),
            min_time      = 0,
            max_time      = mdo.max_time,
            data_shift    = 1
        )

        self.assertTrue(isinstance(fv, FluxVariable))
        self.assertEqual(fv.unit, mdo.stock_unit)

        res = fv.data[1,...]
        res2 = 86.4 * np.array(
               [[[ 84.,  86.],
                 [ 88.,  90.]],

                [[184., 188.],
                 [192., 196.]],

                [[300., 306.],
                 [312., 318.]]]
        )
        self.assertTrue(np.allclose(res, res2))


    def test_getFluxVariable_from_Rate(self):
        mdo = self.mdo
        nr_layers = self.mdo.model_structure.get_nr_layers('Litter')

        fv = getFluxVariable_from_Rate(
            mdo = mdo,
            variable_name = 'LITR_flux_down_tb',
            min_time      = 0,
            max_time      = mdo.max_time,
            nr_layers     = nr_layers,
            lat_arr       = self.lat_arr,
            lon_arr       = self.lon_arr,
            data_shift    = 1
        )

        self.assertTrue(isinstance(fv, FluxVariable))
        self.assertEqual(fv.unit, mdo.stock_unit)

        res = fv.data[1,...]
        res2 = 86.4 * np.array(
               [[[ 84.,  86.],
                 [ 88.,  90.]],

                [[ 92.,  94.],
                 [ 96.,  98.]],

                [[100., 102.],
                 [104,  106.]]]
        )
        self.assertTrue(np.allclose(res, res2))


    ##### Class methods #####

    
    ## still to check calendar conversion
    def test_init(self):
        mdo = self.mdo
        self.assertEqual(mdo.stock_unit, 'g/m^2')
        self.assertTrue(np.all(mdo.time.data == 365+5+np.arange(8)))
        self.assertTrue(
            np.all(mdo.time_agg.data == 365+5+np.array([0,2,4,6,7]))
        )
 
        ## test time conversion
        mdo_time = ModelDataObject(
            model_structure = mdo.model_structure, 
            dataset         = mdo.dataset,
            nstep           = mdo.nstep,
            stock_unit      = mdo.stock_unit,
            cftime_unit_src = 'days since 1850-01-01 00:00:00',
            cftime_unit_tar = 'hours since 1850-01-01 00:00:00',
            calendar_src    = 'noleap',
            udtime_unit     = 'hr',
            calendar        = 'noleap',
            max_time        = mdo.max_time
        )
 
        xs = mdo.load_stocks(
            func       = getStockVariable_from_Density,
            lat_arr    = self.lat_arr,
            lon_arr    = self.lon_arr,
            data_shift = 0
        )
 
        self.assertTrue(isinstance(xs, StockVariable))
        self.assertEqual(xs.unit, mdo.stock_unit)
        self.assertTrue(np.all(xs.data.shape == (5,2,2,9)))
 
        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_stocks = xs.data[...,CWD_pools]
        res2 = np.array(
            [[[[  0.,   8.,  24.],
               [  1.,  10.,  27.]],
 
              [[  2.,  12.,  30.],
               [  3.,  14.,  33.]]],
            
            
             [[[ 24.,  56.,  96.],
               [ 25.,  58.,  99.]],
            
              [[ 26.,  60., 102.],
               [ 27.,  62., 105.]]],
            
            
             [[[ 48., 104., 168.],
               [ 49., 106., 171.]],
            
              [[ 50., 108., 174.],
               [ 51., 110., 177.]]],
            
            
             [[[ 72., 152., 240.],
               [ 73., 154., 243.]],
            
              [[ 74., 156., 246.],
               [ 75., 158., 249.]]],
            
            
             [[[ 84., 176., 276.],
               [ 85., 178., 279.]],
            
              [[ 86., 180., 282.],
               [ 87., 182., 285.]]]]
        )
        self.assertTrue(np.allclose(CWD_stocks, res2))
        
 
        us = mdo.load_external_input_fluxes(
            func       = getFluxVariable_from_DensityRate,
            lat_arr    = self.lat_arr,
            lon_arr    = self.lon_arr,
            data_shift = 1
        )
 
        self.assertTrue(isinstance(us, FluxVariable))
        self.assertEqual(us.unit, mdo.stock_unit)
        self.assertTrue(np.all(us.data.shape == (4,2,2,9)))
 
        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_inputs = us.data[...,CWD_pools]
        res2 = np.array(
            [[[[  9336.384,  22819.968,  40450.752],
               [  9854.784,  23856.768,  42005.952]],
            
              [[ 10373.184,  24893.568,  43561.152],
               [ 10891.584,  25930.368,  45116.352]]],
            
            
             [[[ 21777.984,  47703.168,  77775.552],
               [ 22296.384,  48739.968,  79330.752]],
            
              [[ 22814.784,  49776.768,  80885.952],
               [ 23333.184,  50813.568,  82441.152]]],
            
            
             [[[ 34219.584,  72586.368, 115100.352],
               [ 34737.984,  73623.168, 116655.552]],
            
              [[ 35256.384,  74659.968, 118210.752],
               [ 35774.784,  75696.768, 119765.952]]],
            
            
             [[[ 21775.392,  45624.384,  71546.976],
               [ 22034.592,  46142.784,  72324.576]],
            
              [[ 22293.792,  46661.184,  73102.176],
               [ 22552.992,  47179.584,  73879.776]]]]
        )
        self.assertTrue(np.allclose(CWD_inputs, res2))


    def test_get_dz(self):
        mdo = self.mdo
        ds = mdo.dataset
        dz = mdo.get_dz('CWD')

        self.assertEqual(dz.name, 'dz')
        self.assertEqual(dz.unit, 'm')
        res2 = ds['dz']    
        self.assertTrue(np.all(dz.data == res2[:]))

        ## no layers, no grid
        mdo = self.mdos_inc[0]
        dz = mdo.get_dz('CWD')
        self.assertTrue(np.all(dz.data == np.array([1])))
        self.assertEqual(dz.unit, '1')

        ## test manually given dz variable
        mdo = self.mdo
        dz = Variable(
            name = 'dz',
            data = np.cumsum(np.arange(10)),
            unit = 'km'
        )
        mdo_dz = ModelDataObject(
            model_structure = mdo.model_structure, 
            dataset         = mdo.dataset,
            dz_var_names    = {'dz': dz},
            nstep           = mdo.nstep,
            stock_unit      = mdo.stock_unit,
            cftime_unit_src = 'days since 1850-01-01 00:00:00',
            calendar_src    = 'noleap',
            max_time        = mdo.max_time
        )
        
        dz = mdo_dz.get_dz('CWD')
        self.assertEqual(dz.name, 'dz')
        self.assertEqual(dz.unit, 'km')
        res2 = np.cumsum(np.arange(10))
        self.assertTrue(np.all(dz.data == res2[:]))


    def test_load_stocks(self):
        mdo = self.mdo

        xs = mdo.load_stocks(
            func       = getStockVariable_from_Density,
            lat_arr    = self.lat_arr,
            lon_arr    = self.lon_arr,
            data_shift = 0
        )

        self.assertTrue(isinstance(xs, StockVariable))
        self.assertEqual(xs.unit, mdo.stock_unit)
        self.assertTrue(np.all(xs.data.shape == (5,2,2,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_stocks = xs.data[...,CWD_pools]
        res2 = np.array(
            [[[[  0.,   8.,  24.],
               [  1.,  10.,  27.]],

              [[  2.,  12.,  30.],
               [  3.,  14.,  33.]]],
            
            
             [[[ 24.,  56.,  96.],
               [ 25.,  58.,  99.]],
            
              [[ 26.,  60., 102.],
               [ 27.,  62., 105.]]],
            
            
             [[[ 48., 104., 168.],
               [ 49., 106., 171.]],
            
              [[ 50., 108., 174.],
               [ 51., 110., 177.]]],
            
            
             [[[ 72., 152., 240.],
               [ 73., 154., 243.]],
            
              [[ 74., 156., 246.],
               [ 75., 158., 249.]]],
            
            
             [[[ 84., 176., 276.],
               [ 85., 178., 279.]],
            
              [[ 86., 180., 282.],
               [ 87., 182., 285.]]]]
        )
        self.assertTrue(np.allclose(CWD_stocks, res2))

        pool_nr = mdo.model_structure.get_pool_nr('Soil', 2)
        res = xs.data[...,pool_nr]
        res2 = np.array(
            [[[ 24.6,  27.6],
              [ 30.6,  33.6]],

             [[ 96.6,  99.6],
              [102.6, 105.6]],
            
             [[168.6, 171.6],
              [174.6, 177.6]],
            
             [[240.6, 243.6],
              [246.6, 249.6]],
            
             [[276.6, 279.6],
              [282.6, 285.6]]]
        )
        self.assertTrue(np.allclose(res, res2))

        ## incomplete MDOs
        data_shapes = [(5,1,1,3), (5,2,2,3), (5,1,1,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            xs = mdo.load_stocks(
                func       = getStockVariable_from_Density,
                lat_arr    = mdo.lat_arr,
                lon_arr    = mdo.lon_arr,
                data_shift = 0
            )
    
            self.assertTrue(isinstance(xs, StockVariable))
            self.assertEqual(xs.unit, mdo.stock_unit)
            self.assertTrue(np.all(xs.data.shape == data_shapes[nr]))


    def test_load_external_input_fluxes(self):
        mdo = self.mdo

        us = mdo.load_external_input_fluxes(
            func       = getFluxVariable_from_DensityRate,
            lat_arr    = self.lat_arr,
            lon_arr    = self.lon_arr,
            data_shift = 1
        )

        self.assertTrue(isinstance(us, FluxVariable))
        self.assertEqual(us.unit, mdo.stock_unit)
        self.assertTrue(np.all(us.data.shape == (4,2,2,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_inputs = us.data[...,CWD_pools]
        res2 = np.array(
            [[[[  9336.384,  22819.968,  40450.752],
               [  9854.784,  23856.768,  42005.952]],
            
              [[ 10373.184,  24893.568,  43561.152],
               [ 10891.584,  25930.368,  45116.352]]],
            
            
             [[[ 21777.984,  47703.168,  77775.552],
               [ 22296.384,  48739.968,  79330.752]],
            
              [[ 22814.784,  49776.768,  80885.952],
               [ 23333.184,  50813.568,  82441.152]]],
            
            
             [[[ 34219.584,  72586.368, 115100.352],
               [ 34737.984,  73623.168, 116655.552]],
            
              [[ 35256.384,  74659.968, 118210.752],
               [ 35774.784,  75696.768, 119765.952]]],
            
            
             [[[ 21775.392,  45624.384,  71546.976],
               [ 22034.592,  46142.784,  72324.576]],
            
              [[ 22293.792,  46661.184,  73102.176],
               [ 22552.992,  47179.584,  73879.776]]]]
        )
        self.assertTrue(np.allclose(CWD_inputs, res2))

        ## incomplete MDOs
        data_shapes = [(4,1,1,3), (4,2,2,3), (4,1,1,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            us = mdo.load_external_input_fluxes(
                func       = getFluxVariable_from_DensityRate,
                lat_arr    = mdo.lat_arr,
                lon_arr    = mdo.lon_arr,
                data_shift = 1
            )
    
            self.assertTrue(isinstance(us, FluxVariable))
            self.assertEqual(us.unit, mdo.stock_unit)
            self.assertTrue(np.all(us.data.shape == data_shapes[nr]))


    def test_load_external_output_fluxes(self):
        mdo = self.mdo

        rs = mdo.load_external_output_fluxes(
            func       = getFluxVariable_from_DensityRate,
            lat_arr    = self.lat_arr,
            lon_arr    = self.lon_arr,
            data_shift = 1
        )

        self.assertTrue(isinstance(rs, FluxVariable))
        self.assertEqual(rs.unit, mdo.stock_unit)
        self.assertTrue(np.all(rs.data.shape == (4,2,2,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_outputs = rs.data[...,CWD_pools]
        res2 = np.array(
            [[[[ 6222.528, 15209.856, 26961.984],
               [ 6568.128, 15901.056, 27998.784]],
            
              [[ 6913.728, 16592.256, 29035.584],
               [ 7259.328, 17283.456, 30072.384]]],
            
            
             [[[14516.928, 31798.656, 51845.184],
               [14862.528, 32489.856, 52881.984]],
            
              [[15208.128, 33181.056, 53918.784],
               [15553.728, 33872.256, 54955.584]]],
            
            
             [[[22811.328, 48387.456, 76728.384],
               [23156.928, 49078.656, 77765.184]],
            
              [[23502.528, 49769.856, 78801.984],
               [23848.128, 50461.056, 79838.784]]],
            
            
             [[[14516.064, 30414.528, 47695.392],
               [14688.864, 30760.128, 48213.792]],
            
              [[14861.664, 31105.728, 48732.192],
               [15034.464, 31451.328, 49250.592]]]]
        )
        self.assertTrue(np.allclose(CWD_outputs, res2))

        ## incomplete MDOs
        data_shapes = [(4,1,1,3), (4,2,2,3), (4,1,1,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            rs = mdo.load_external_output_fluxes(
                func       = getFluxVariable_from_DensityRate,
                lat_arr    = mdo.lat_arr,
                lon_arr    = mdo.lon_arr,
                data_shift = 1
            )
    
            self.assertTrue(isinstance(rs, FluxVariable))
            self.assertEqual(rs.unit, mdo.stock_unit)
            self.assertTrue(np.all(rs.data.shape == data_shapes[nr]))


    def test_load_horizontal_fluxes(self):
        mdo = self.mdo

        HFs = mdo.load_horizontal_fluxes(
            func       = getFluxVariable_from_DensityRate,
            lat_arr    = self.lat_arr,
            lon_arr    = self.lon_arr,
            data_shift = 1
        )
    
        self.assertTrue(isinstance(HFs, FluxVariable))
        self.assertEqual(HFs.unit, mdo.stock_unit)
        self.assertTrue(np.all(HFs.data.shape == (4,2,2,9,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        LITR_pools = mdo.model_structure.get_pool_nrs('Litter')

        CWD_TO_LITR = HFs.data[...,LITR_pools,CWD_pools]
        res2 = np.array(
            [[[[ 3110.4,  7603.2, 13478.4],
               [ 3283.2,  7948.8, 13996.8]],
            
              [[ 3456. ,  8294.4, 14515.2],
               [ 3628.8,  8640. , 15033.6]]],
            
            
             [[[ 7257.6, 15897.6, 25920. ],
               [ 7430.4, 16243.2, 26438.4]],
            
              [[ 7603.2, 16588.8, 26956.8],
               [ 7776. , 16934.4, 27475.2]]],
            
            
             [[[11404.8, 24192. , 38361.6],
               [11577.6, 24537.6, 38880. ]],
            
              [[11750.4, 24883.2, 39398.4],
               [11923.2, 25228.8, 39916.8]]],
            
            
             [[[ 7257.6, 15206.4, 23846.4],
               [ 7344. , 15379.2, 24105.6]],
            
              [[ 7430.4, 15552. , 24364.8],
               [ 7516.8, 15724.8, 24624. ]]]]
        )
        self.assertTrue(np.allclose(CWD_TO_LITR, res2))

        LITR_TO_CWD = HFs.data[...,CWD_pools,LITR_pools]
        self.assertFalse(np.any(LITR_TO_CWD))

        ## incomplete MDOs
        data_shapes = [(4,1,1,3,3), (4,2,2,3,3), (4,1,1,9,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            Fs = mdo.load_horizontal_fluxes(
                func       = getFluxVariable_from_DensityRate,
                lat_arr    = mdo.lat_arr,
                lon_arr    = mdo.lon_arr,
                data_shift = 1
            )
    
            self.assertTrue(isinstance(Fs, FluxVariable))
            self.assertEqual(Fs.unit, mdo.stock_unit)
            self.assertTrue(np.all(Fs.data.shape == data_shapes[nr]))


    def test_load_vertical_fluxes(self):
        mdo = self.mdo

        VFs, runoffs_up, runoffs_down = mdo.load_vertical_fluxes(
            func       = getFluxVariable_from_Rate,
            lat_arr    = self.lat_arr,
            lon_arr    = self.lon_arr,
            data_shift = 1
        )
    
        self.assertTrue(isinstance(VFs, FluxVariable))
        self.assertEqual(VFs.unit, mdo.stock_unit)
        self.assertTrue(np.all(VFs.data.shape == (4,2,2,9,9)))
        
        self.assertTrue(isinstance(runoffs_up, FluxVariable))
        self.assertEqual(runoffs_up.unit, mdo.stock_unit)
        self.assertTrue(np.all(runoffs_up.data.shape == (4,2,2,9)))

        self.assertTrue(isinstance(runoffs_down, FluxVariable))
        self.assertEqual(runoffs_down.unit, mdo.stock_unit)
        self.assertTrue(np.all(runoffs_down.data.shape == (4,2,2,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        LITR_pools = mdo.model_structure.get_pool_nrs('Litter')
        SOIL_pools = mdo.model_structure.get_pool_nrs('Soil')

        VFs_CWD = VFs.data[...,CWD_pools,CWD_pools]
        self.assertFalse(np.any(VFs_CWD))
        runoffs_up_CWD = runoffs_up.data[...,CWD_pools]
        self.assertFalse(np.any(runoffs_up_CWD))
        runoffs_down_CWD = runoffs_down.data[...,CWD_pools]
        self.assertFalse(np.any(runoffs_down_CWD))
        
        VFs_LITR = VFs.data[...,LITR_pools,:][...,LITR_pools]
        res2 = np.array(
            [[[[    0.  ,  1140.48 ,    0.  ],
                [11404.8,      0.  ,  1209.6 ],
                [    0. ,  12096.  ,     0.  ]],
            
               [[    0. ,   1157.76,     0.  ],
                [11577.6,      0.  ,  1226.88],
                [    0. ,  12268.8 ,     0.  ]]],
            
            
              [[[    0. ,   1175.04,     0.  ],
                [11750.4,      0.  ,  1244.16],
                [    0. ,  12441.6 ,     0.  ]],
            
               [[    0. ,   1192.32,     0.  ],
                [11923.2,      0.  ,  1261.44],
                [    0. ,  12614.4 ,     0.  ]]]]
        )
        self.assertTrue(np.allclose(VFs_LITR[2,...], res2))

        runoffs_up_SOIL = runoffs_up.data[...,SOIL_pools]
        res2 = np.array(
            [[[[ 311.04,    0.,      0.  ],
               [ 328.32,    0.,      0.  ]],
            
              [[ 345.6 ,    0.,      0.  ],
               [ 362.88,    0.,      0.  ]]],
            
            
             [[[ 725.76,    0.,      0.  ],
               [ 743.04,    0.,      0.  ]],
            
              [[ 760.32,    0.,      0.  ],
               [ 777.6 ,    0.,      0.  ]]],
            
            
             [[[1140.48,    0.,      0.  ],
               [1157.76,    0.,      0.  ]],
            
              [[1175.04,    0.,      0.  ],
               [1192.32,    0.,      0.  ]]],
            
            
             [[[ 725.76,    0.,      0.  ],
               [ 734.4 ,    0.,      0.  ]],
            
              [[ 743.04,    0.,      0.  ],
               [ 751.68,    0.,      0.  ]]]]
        )
        self.assertTrue(np.allclose(runoffs_up_SOIL, res2))

        runoffs_down_LITR = runoffs_down.data[...,LITR_pools]
        res2 = np.array(
            [[[[    0.,      0.,   4492.8],
               [    0.,      0.,   4665.6]],
            
              [[    0.,      0.,   4838.4],
               [    0.,      0.,   5011.2]]],
            
            
             [[[    0.,      0.,   8640. ],
               [    0.,      0.,   8812.8]],
            
              [[    0.,      0.,   8985.6],
               [    0.,      0.,   9158.4]]],
            
            
             [[[    0.,      0.,  12787.2],
               [    0.,      0.,  12960. ]],
            
              [[    0.,      0.,  13132.8],
               [    0.,      0.,  13305.6]]],
            
            
             [[[    0.,      0.,   7948.8],
               [    0.,      0.,   8035.2]],
            
              [[    0.,      0.,   8121.6],
               [    0.,      0.,   8208. ]]]]
        )
        self.assertTrue(np.allclose(runoffs_down_LITR, res2))

        VFs_SOIL = VFs.data[...,SOIL_pools,:][...,SOIL_pools]
        res2 = np.array(
             [[[[    0. ,   1209.6 ,     0.  ],
                [12096. ,      0.  ,  1278.72],
                [    0. ,  12787.2 ,     0.  ]],
            
               [[    0. ,   1226.88,     0.  ],
                [12268.8,      0.  ,  1296.  ],
                [    0. ,  12960.  ,     0.  ]]],
            
            
              [[[    0. ,   1244.16,     0.  ],
                [12441.6,      0.  ,  1313.28],
                [    0. ,  13132.8 ,     0.  ]],
            
               [[    0. ,   1261.44,     0.  ],
                [12614.4,      0.  ,  1330.56],
                [    0. ,  13305.6 ,     0.  ]]]]
        )
        self.assertTrue(np.allclose(VFs_SOIL[2,...], res2))

        runoffs_up_LITR = runoffs_up.data[...,LITR_pools]
        self.assertFalse(np.any(runoffs_up_LITR))
        runoffs_down_SOIL = runoffs_down.data[...,SOIL_pools]
        self.assertFalse(np.any(runoffs_down_SOIL))


        ## incomplete MDOs
        data_shapes    = [(4,1,1,3,3), (4,2,2,3,3), (4,1,1,9,9)]
        runoffs_shapes = [(4,1,1,3),   (4,2,2,3)  , (4,1,1,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            VFs, runoffs_up, runoffs_down = mdo.load_vertical_fluxes(
                func       = getFluxVariable_from_Rate,
                lat_arr    = mdo.lat_arr,
                lon_arr    = mdo.lon_arr,
                data_shift = 1
            )
    
            self.assertTrue(isinstance(VFs, FluxVariable))
            self.assertEqual(VFs.unit, mdo.stock_unit)
            self.assertTrue(np.all(VFs.data.shape == data_shapes[nr]))
            self.assertTrue(np.all(runoffs_up.data.shape == runoffs_shapes[nr]))
            self.assertTrue(
                np.all(runoffs_down.data.shape == runoffs_shapes[nr])
            )


    def test_load_xs_us_Fs_rs(self):
        mdo = self.mdo_discretizable
        xs, us, Fs, rs = mdo.load_xs_us_Fs_rs(
            lat_index = 0,
            lon_index = 0
        )
        self.assertTrue(isinstance(xs, StockVariable))
        self.assertEqual(xs.unit, 'g/m^2')
        self.assertEqual(xs.data.shape, (5,1,1,9))
        res2 = np.ones((3,3))
        dz = (np.arange(3)+1).reshape((1,3))
        res2 = (res2*dz).reshape((9,))
        self.assertTrue(np.allclose(xs.data[0,0,0,:], res2*1e5))
 
        self.assertTrue(isinstance(us, FluxVariable))
        self.assertEqual(us.unit, 'g/m^2')
        self.assertEqual(us.data.shape, (4,1,1,9))
        res2 = np.ones((3,3))*1e-02 * 86400
        dz = (np.arange(3)+1).reshape((1,3))
        res2 = (res2*dz).reshape((9,))
        res2[-3:] = 0
        self.assertTrue(np.allclose(us.data[2,0,0,:], res2*2))
        
        self.assertTrue(isinstance(Fs, FluxVariable))
        self.assertEqual(Fs.unit, 'g/m^2')
        self.assertEqual(Fs.data.shape, (4,1,1,9,9))
        res2 = np.array([0,1,0]) *1e-03 * 86400 * 2 * 2
        self.assertTrue(np.allclose(Fs.data[2,0,0,4,:3], res2))

        self.assertTrue(isinstance(rs, FluxVariable))
        self.assertEqual(rs.unit, 'g/m^2')
        self.assertEqual(rs.data.shape, (4,1,1,9))
        res2 = np.ones((3,3))*1e-03 *2 * 86400
        dz = (np.arange(3)+1).reshape((1,3))
        res2 = (res2*dz).reshape((9,))
        self.assertTrue(np.allclose(rs.data[3,0,0,:], res2))
   
 
    def test_create_discrete_model_run(self):
        mdo = self.mdo_discretizable
        
        with self.assertRaises(ModelDataObjectException) as context:
            dmr = mdo.create_discrete_model_run()
        msg = 'Data structure not understood'
        self.assertEqual(str(context.exception), msg)
        
        dmr = mdo.create_discrete_model_run(
            lat_index = 0,
            lon_index = 0
        )
        self.assertEqual(dmr.nr_pools, 9)
        soln = dmr.solve()
        self.assertTrue(np.allclose(dmr.times, [0,2,4,6,7]))

        ## check missing data
        mdo = self.mdos_inc[2]
        dmr = mdo.create_discrete_model_run()
        self.assertTrue(dmr is None)


    def test_create_model_run(self):
        mdo = self.mdo_discretizable
        
        with self.assertRaises(ModelDataObjectException) as context:
            mr_pwc = mdo.create_model_run()
        msg = 'Data structure not understood'
        self.assertEqual(str(context.exception), msg)
        
        mr_pwc = mdo.create_model_run(
            lat_index = 0,
            lon_index = 0
        )
        self.assertEqual(mr_pwc.model.nr_pools, 9)
        soln = mr_pwc.solve()
        self.assertTrue(np.allclose(mr_pwc.data_times, [0,2,4,6,7]))

        ## check missing data
        mdo = self.mdos_inc[2]
        mr_pwc = mdo.create_discrete_model_run()
        self.assertTrue(mr_pwc is None)



################################################################################


if __name__ == '__main__':
    unittest.main()

