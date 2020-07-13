import unittest
import numpy as np
import xarray as xr
from bgc_md2.models.CARDAMOM.CARDAMOMlib import load_mdo

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
from example_MDOs import MDO_3pools_3layers,\
                         MDO_3pools_nolayers_nogrid,\
                         MDO_3pools_3layers_nogrid,\
                         MDO_3pools_3layers_discretizable

import os.path
THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestModelDataObject(unittest.TestCase):

    def setUp(self):
        self.mdo = MDO_3pools_3layers()
        ## incomplete models
        mdos_inc = [
            MDO_3pools_nolayers_nogrid(),
            MDO_3pools_3layers_nogrid()
        ]

        mdos_inc[0].nr_layers = None
        mdos_inc[1].nr_layers = 3

        self.mdos_inc = mdos_inc

        ## discretizable MDOs
        self.mdo_discretizable = MDO_3pools_3layers_discretizable()


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
            nr_layers     = 2,
            data_shift    = 1
        )
        res2 = np.array(
            [[3.0 , 4.0],
             [6.0 , 7.0],
             [9.0 , 10.0],
             [12.0, 13.0],
             [15.0, 16.0],
             [18.0, 19.0],
             [21.0, 22.0],
             [24.0, 25.0],
             [27.0, 28.0]]
        ) 
        self.assertTrue(isinstance(var, Variable)) 
        self.assertTrue(np.allclose(var.data,res2))
        self.assertEqual(var.unit, 'g/m^3')
 
        sv = readVariable(
            ReturnClass   = StockVariable,
            dataset       = ds,
            variable_name = 'CWDC',
            nr_layers     = 2,
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
            data_shift    = 1
        )
        self.assertTrue(isinstance(sv, StockVariable))
        self.assertEqual(sv.unit, 'g/m^2')
        self.assertTrue(np.all(sv.data.shape == np.array([9,1])))

        ##### 3 layers, no grid #####

        mdo = self.mdos_inc[1]
        ds = mdo.dataset
        ms = mdo.model_structure

        sv = readVariable(
            ReturnClass   = StockVariable,
            dataset       = ds,
            variable_name = 'CWDC',
            nr_layers     = ms.get_nr_layers('CWD'),
            data_shift    = 1
        )
        self.assertTrue(isinstance(sv, StockVariable))
        self.assertEqual(sv.unit, 'g/m^3')
        self.assertTrue(np.all(sv.data.shape == np.array([9,3])))


    def test_StockDensityVariable2StockVariable(self):
        ds = self.mdo.dataset
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')

        sdv = readVariable(
            ReturnClass   = StockVariable,
            dataset       = ds,
            variable_name = 'CWDC',
            nr_layers     = nr_layers,
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
        res2 = np.array([12, 26, 42])
        self.assertTrue(np.allclose(res, res2))


    def test_getStockVariable_from_Density(self):
        mdo = self.mdo
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')

        sv = getStockVariable_from_Density(
            mdo = mdo,
            variable_name = 'CWDC',
            nr_layers     = nr_layers,
            dz            = mdo.get_dz('CWD'),
            data_shift    = 0
        )

        self.assertTrue(isinstance(sv, StockVariable))
        self.assertEqual(sv.unit, mdo.stock_unit)

        res = sv.data[2,...]
        res2 = np.array([12, 13*2, 14*3])
        self.assertTrue(np.allclose(res, res2))


    def test_FluxRateDensityVariable2FluxRateVariable(self):
        ds = self.mdo.dataset
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')

        frdv = readVariable(
            ReturnClass   = FluxVariable,
            dataset       = ds,
            variable_name = 'fire_CWD',
            nr_layers     = nr_layers,
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
        res2 = np.array([12, 26, 42])
        self.assertTrue(np.allclose(res, res2*1e-03))


    def test_FluxRateVariable2FluxVariable(self):
        ds = self.mdo.dataset
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')
        
        frdv = readVariable(
            ReturnClass   = FluxVariable,
            dataset       = ds,
            variable_name = 'fire_CWD',
            nr_layers     = nr_layers,
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
        res2 = 86.4 * np.array([12, 26, 42])
        self.assertTrue(np.allclose(res, res2*1e-03))


    def test_getFluxVariable_from_DensityRate(self):
        mdo = self.mdo
        nr_layers = self.mdo.model_structure.get_nr_layers('CWD')

        fv = getFluxVariable_from_DensityRate(
            mdo = mdo,
            variable_name = 'fire_CWD',
            nr_layers     = nr_layers,
            dz            = mdo.get_dz('CWD'),
            data_shift    = 1
        )

        self.assertTrue(isinstance(fv, FluxVariable))
        self.assertEqual(fv.unit, mdo.stock_unit)

        res = fv.data[1,...]
        res2 = 86.4 * np.array([21, 23*2, 25*3])
        self.assertTrue(np.allclose(res, res2))


    def test_getFluxVariable_from_Rate(self):
        mdo = self.mdo
        nr_layers = self.mdo.model_structure.get_nr_layers('Litter')

        fv = getFluxVariable_from_Rate(
            mdo = mdo,
            variable_name = 'LITR_flux_down_tb',
            nr_layers     = nr_layers,
            data_shift    = 1
        )

        self.assertTrue(isinstance(fv, FluxVariable))
        self.assertEqual(fv.unit, mdo.stock_unit)

        res = fv.data[1,...]
        res2 = 86.4 * np.array([9+12, 10+13, 11+14])
        self.assertTrue(np.allclose(res, res2))


    ##### Class methods #####

    
    ## still to check calendar conversion
    def test_init(self):
        mdo = self.mdo
        self.assertEqual(mdo.stock_unit, 'g/m^2')
        self.assertTrue(np.all(mdo.time.data == np.arange(10)))
        self.assertTrue(
            np.all(mdo.time_agg.data == np.array([0,2,4,6,8,9]))
        )
 
        xs = mdo.load_stocks(
            func       = getStockVariable_from_Density,
            data_shift = 0
        )
 
        self.assertTrue(isinstance(xs, StockVariable))
        self.assertEqual(xs.unit, mdo.stock_unit)
        self.assertTrue(np.all(xs.data.shape == (6,9)))
 
        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_stocks = xs.data[...,CWD_pools]
        res2 = np.array(
            [[ 0.0,  2.0,  6.0],
             [ 6.0, 14.0, 24.0],
             [12.0, 26.0, 42.0],
             [18.0, 38.0, 60.0],
             [24.0, 50.0, 78.0],
             [27.0, 56.0, 87.0]]
        )
        self.assertTrue(np.allclose(CWD_stocks, res2))
        
 
        us = mdo.load_external_input_fluxes(
            func       = getFluxVariable_from_DensityRate,
            data_shift = 1
        )
 
        self.assertTrue(isinstance(us, FluxVariable))
        self.assertEqual(us.unit, mdo.stock_unit)
        self.assertTrue(np.all(us.data.shape == (5,9)))
 
        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_inputs = us.data[...,CWD_pools]
        res2 = np.array(
            [[ 2337.984,  5712.768, 10124.352],
             [ 5448.384, 11933.568, 19455.552],
             [ 8558.784, 18154.368, 28786.752],
             [11669.184, 24375.168, 38117.952],
             [ 7000.992, 14520.384, 22558.176]]
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
        time = Variable(
            data = np.arange(10),
            unit = 'd'
        )
        mdo_dz = ModelDataObject(
            model_structure = mdo.model_structure, 
            dataset         = mdo.dataset,
            dz_var_names    = {'dz': dz},
            nstep           = mdo.nstep,
            stock_unit      = mdo.stock_unit,
            time            = time
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
            data_shift = 0
        )

        self.assertTrue(isinstance(xs, StockVariable))
        self.assertEqual(xs.unit, mdo.stock_unit)
        self.assertTrue(np.all(xs.data.shape == (6,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_stocks = xs.data[:,CWD_pools]
        res2 = np.array(
            [[0.0 ,  2.0,  6.0],
             [6.0 , 14.0, 24.0],
             [12.0, 26.0, 42.0],
             [18.0, 38.0, 60.0],
             [24.0, 50.0, 78.0],
             [27.0, 56.0, 87.0]]
        )
        self.assertTrue(np.allclose(CWD_stocks, res2))

        pool_nr = mdo.model_structure.get_pool_nr('Soil', 2)
        res = xs.data[:,pool_nr]
        res2 = 3* np.array(
            [ 2.2,
              8.2,
             14.2,
             20.2,
             26.2,
             29.2]  
        )
        self.assertTrue(np.allclose(res, res2))

        ## incomplete MDOs
        data_shapes = [(6,3), (6,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            xs = mdo.load_stocks(
                func       = getStockVariable_from_Density,
                data_shift = 0
            )
    
            self.assertTrue(isinstance(xs, StockVariable))
            self.assertEqual(xs.unit, mdo.stock_unit)
            self.assertTrue(np.all(xs.data.shape == data_shapes[nr]))


    def test_load_external_input_fluxes(self):
        mdo = self.mdo

        us = mdo.load_external_input_fluxes(
            func       = getFluxVariable_from_DensityRate,
            data_shift = 1
        )

        self.assertTrue(isinstance(us, FluxVariable))
        self.assertEqual(us.unit, mdo.stock_unit)
        self.assertTrue(np.all(us.data.shape == (5,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_inputs = us.data[...,CWD_pools]
        res2 = np.array(
            [[ 2337.984,  5712.768, 10124.352],
             [ 5448.384, 11933.568, 19455.552],
             [ 8558.784, 18154.368, 28786.752],
             [11669.184, 24375.168, 38117.952],
             [ 7000.992, 14520.384, 22558.176]]
        )
        self.assertTrue(np.allclose(CWD_inputs, res2))

        ## incomplete MDOs
        data_shapes = [(5,3), (5,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            us = mdo.load_external_input_fluxes(
                func       = getFluxVariable_from_DensityRate,
                data_shift = 1
            )
    
            self.assertTrue(isinstance(us, FluxVariable))
            self.assertEqual(us.unit, mdo.stock_unit)
            self.assertTrue(np.all(us.data.shape == data_shapes[nr]))


    def test_load_external_output_fluxes(self):
        mdo = self.mdo

        rs = mdo.load_external_output_fluxes(
            func       = getFluxVariable_from_DensityRate,
            data_shift = 1
        )

        self.assertTrue(isinstance(rs, FluxVariable))
        self.assertEqual(rs.unit, mdo.stock_unit)
        self.assertTrue(np.all(rs.data.shape == (5,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        CWD_outputs = rs.data[...,CWD_pools]
        res2 = np.array(
            [[1556.928,  3805.056,  6744.384],
             [3630.528,  7952.256, 12965.184],
             [5704.128, 12099.456, 19185.984],
             [7777.728, 16246.656, 25406.784],
             [4666.464,  9678.528, 15036.192]]
        )
        self.assertTrue(np.allclose(CWD_outputs, res2))

        ## incomplete MDOs
        data_shapes = [(5,3), (5,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            rs = mdo.load_external_output_fluxes(
                func       = getFluxVariable_from_DensityRate,
                data_shift = 1
            )
    
            self.assertTrue(isinstance(rs, FluxVariable))
            self.assertEqual(rs.unit, mdo.stock_unit)
            self.assertTrue(np.all(rs.data.shape == data_shapes[nr]))


    def test_load_horizontal_fluxes(self):
        mdo = self.mdo

        HFs = mdo.load_horizontal_fluxes(
            func       = getFluxVariable_from_DensityRate,
            data_shift = 1
        )
    
        self.assertTrue(isinstance(HFs, FluxVariable))
        self.assertEqual(HFs.unit, mdo.stock_unit)
        self.assertTrue(np.all(HFs.data.shape == (5,9,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        LITR_pools = mdo.model_structure.get_pool_nrs('Litter')

        CWD_TO_LITR = HFs.data[:,LITR_pools,:][:,:,CWD_pools]
        res2 = np.array(
            [[[777.6,    0.0,    0.0],
              [  0.0, 1900.8,    0.0],
              [  0.0,    0.0, 3369.6]],
            
             [[1814.4,    0.0,    0.0],
              [   0.0, 3974.4,    0.0],
              [   0.0,    0.0, 6480.0]],
            
             [[2851.2,    0.0,    0.0],
              [   0.0, 6048.0,    0.0],
              [   0.0,    0.0, 9590.4]],
            
             [[3888.0,    0.0,     0.0],
              [   0.0, 8121.6,     0.0],
              [   0.0,    0.0, 12700.8]],
            
             [[2332.8,    0.0,    0.0],
              [   0.0, 4838.4,    0.0],
              [   0.0,    0.0, 7516.8]]]
        )
        self.assertTrue(np.allclose(CWD_TO_LITR, res2))

        LITR_TO_CWD = HFs.data[:,CWD_pools,:][:,:,LITR_pools]
        self.assertFalse(np.any(LITR_TO_CWD))

        ## incomplete MDOs
        data_shapes = [(5,3,3), (5,9,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            Fs = mdo.load_horizontal_fluxes(
                func       = getFluxVariable_from_DensityRate,
                data_shift = 1
            )
    
            self.assertTrue(isinstance(Fs, FluxVariable))
            self.assertEqual(Fs.unit, mdo.stock_unit)
            self.assertTrue(np.all(Fs.data.shape == data_shapes[nr]))


    def test_load_vertical_fluxes(self):
        mdo = self.mdo

        VFs, runoffs_up, runoffs_down = mdo.load_vertical_fluxes(
            func       = getFluxVariable_from_Rate,
            data_shift = 1
        )
    
        self.assertTrue(isinstance(VFs, FluxVariable))
        self.assertEqual(VFs.unit, mdo.stock_unit)
        self.assertTrue(np.all(VFs.data.shape == (5,9,9)))
        
        self.assertTrue(isinstance(runoffs_up, FluxVariable))
        self.assertEqual(runoffs_up.unit, mdo.stock_unit)
        self.assertTrue(np.all(runoffs_up.data.shape == (5,9)))

        self.assertTrue(isinstance(runoffs_down, FluxVariable))
        self.assertEqual(runoffs_down.unit, mdo.stock_unit)
        self.assertTrue(np.all(runoffs_down.data.shape == (5,9)))

        CWD_pools = mdo.model_structure.get_pool_nrs('CWD')
        LITR_pools = mdo.model_structure.get_pool_nrs('Litter')
        SOIL_pools = mdo.model_structure.get_pool_nrs('Soil')

        VFs_CWD = VFs.data[CWD_pools,CWD_pools]
        self.assertFalse(np.any(VFs_CWD))
        runoffs_up_CWD = runoffs_up.data[...,CWD_pools]
        self.assertFalse(np.any(runoffs_up_CWD))
        runoffs_down_CWD = runoffs_down.data[...,CWD_pools]
        self.assertFalse(np.any(runoffs_down_CWD))

        VFs_LITR = VFs.data[:,LITR_pools,:][...,LITR_pools]
        res2 = np.array(
            [[   0.0,  285.12,   0.0],
             [2851.2,    0.0 , 302.4],
             [   0.0, 3024.0 ,   0.0]]
        )
        self.assertTrue(np.allclose(VFs_LITR[2,...], res2))

        runoffs_up_SOIL = runoffs_up.data[:,SOIL_pools]
        res2 = np.array(
           [[ 77.76, 0.0, 0.0],
            [181.44, 0.0, 0.0],
            [285.12, 0.0, 0.0],
            [388.80, 0.0, 0.0],
            [233.28, 0.0, 0.0]]
        )
        self.assertTrue(np.allclose(runoffs_up_SOIL, res2))

        runoffs_down_LITR = runoffs_down.data[:,LITR_pools]
        res2 = np.array(
            [[0.0, 0.0, 1123.2],
             [0.0, 0.0, 2160.0],
             [0.0, 0.0, 3196.8],
             [0.0, 0.0, 4233.6],
             [0.0, 0.0, 2505.6]]
        )
        self.assertTrue(np.allclose(runoffs_down_LITR, res2))

        VFs_SOIL = VFs.data[:,SOIL_pools,:][...,SOIL_pools]
        res2 = np.array(
            [[   0.0,  302.4,   0.0 ],
             [3024.0,    0.0, 319.68],
             [   0.0, 3196.8,   0.0]]
        )
        self.assertTrue(np.allclose(VFs_SOIL[2,...], res2))

        runoffs_up_LITR = runoffs_up.data[...,LITR_pools]
        self.assertFalse(np.any(runoffs_up_LITR))
        runoffs_down_SOIL = runoffs_down.data[...,SOIL_pools]
        self.assertFalse(np.any(runoffs_down_SOIL))


        ## incomplete MDOs
        data_shapes    = [(5,3,3), (5,9,9)]
        runoffs_shapes = [(5,3),   (5,9)]
        for nr, mdo in enumerate(self.mdos_inc):
            VFs, runoffs_up, runoffs_down = mdo.load_vertical_fluxes(
                func       = getFluxVariable_from_Rate,
                data_shift = 1
            )
    
            self.assertTrue(isinstance(VFs, FluxVariable))
            self.assertEqual(VFs.unit, mdo.stock_unit)
            self.assertTrue(np.all(VFs.data.shape == data_shapes[nr]))
            self.assertTrue(np.all(runoffs_up.data.shape == runoffs_shapes[nr]))
            self.assertTrue(
                np.all(runoffs_down.data.shape == runoffs_shapes[nr])
            )


    def test_load_xs_Us_Fs_Rs(self):
        mdo = self.mdo_discretizable
        xs, Us, Fs, Rs = mdo.load_xs_Us_Fs_Rs()

        self.assertTrue(isinstance(xs, StockVariable))
        self.assertEqual(xs.unit, 'g/m^2')
        self.assertEqual(xs.data.shape, (6,9))
        res2 = np.ones((3,3))
        dz = (np.arange(3)+1).reshape((1,3))
        res2 = (res2*dz).reshape((9,))
        self.assertTrue(np.allclose(xs.data[0,:], res2*1e5))
 
        self.assertTrue(isinstance(Us, FluxVariable))
        self.assertEqual(Us.unit, 'g/m^2')
        self.assertEqual(Us.data.shape, (5,9))
        res2 = np.ones((3,3))*1e-02 * 86400
        dz = (np.arange(3)+1).reshape((1,3))
        res2 = (res2*dz).reshape((9,))
        res2[-3:] = 0
        self.assertTrue(np.allclose(Us.data[2,:], res2*2))
        
        self.assertTrue(isinstance(Fs, FluxVariable))
        self.assertEqual(Fs.unit, 'g/m^2')
        self.assertEqual(Fs.data.shape, (5,9,9))
        res2 = np.array([0,1,0]) *1e-03 * 86400 * 2 * 2
        self.assertTrue(np.allclose(Fs.data[2,4,:3], res2))

        self.assertTrue(isinstance(Rs, FluxVariable))
        self.assertEqual(Rs.unit, 'g/m^2')
        self.assertEqual(Rs.data.shape, (5,9))
        res2 = np.ones((3,3))*1e-03 *2 * 86400
        dz = (np.arange(3)+1).reshape((1,3))
        res2 = (res2*dz).reshape((9,)) * 2 # nstep=2
        self.assertTrue(np.allclose(Rs.data[3,:], res2))
   
 
    def test_create_discrete_model_run(self):
        mdo = self.mdo_discretizable
        
        dmr = mdo.create_discrete_model_run()
        self.assertEqual(dmr.nr_pools, 9)
        soln = dmr.solve()
        self.assertTrue(np.allclose(dmr.times, [0,2,4,6,8,9]))

        ## check missing data
        mdo = self.mdos_inc[1]
        dmr = mdo.create_discrete_model_run()
        self.assertTrue(dmr is None)

#    @unittest.skip #mm just to test the github workflow which had other issues
    def test_create_model_run(self):
        mdo = self.mdo_discretizable
        
        mr_pwc = mdo.create_model_run()
        self.assertEqual(mr_pwc.model.nr_pools, 9)
        soln = mr_pwc.solve()
        self.assertTrue(np.allclose(mr_pwc.times, [0,2,4,6,8,9]))

        ## check missing data
        mdo = self.mdos_inc[1]
        mr_pwc = mdo.create_discrete_model_run()
        self.assertTrue(mr_pwc is None)

    @unittest.skip
    def test_get_stock(self):
        dataset = xr.open_dataset('~/Desktop/CARDAMOM/cardamom_for_holger.nc')
        ds = dataset.isel(ens=0, lat=0, lon=0)
        mdo = load_mdo(ds)
        xs, Us, Fs, Rs = mdo.load_xs_Us_Fs_Rs()

        mr = mdo.create_model_run()
#         print(soln[:, 0]-self.get_stock(pwc_mr_fd, 'Labile'))
#        print(xs.data[:, :])
#        print(mdo.get_stock(mr, 'CWD'))
        self.assertTrue(
            np.all(xs.data[:, 0] == mdo.get_stock(mr, 'CWD').data)
        ) 

    def test_load_datafile(self):
        my_data_path = os.path.join(THIS_DIR, 'datafile.txt')
        with open(my_data_path, 'r') as datafile:
            pass


################################################################################


if __name__ == '__main__':
    unittest.main()

