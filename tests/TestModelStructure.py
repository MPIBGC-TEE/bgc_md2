import unittest

import numpy as np
from sympy import symbols

from bgc_md2.ModelStructure import ModelStructure, ModelStructureException


class TestModelStructure(unittest.TestCase):

    def setUp(self):
        nr_layers = 10
        dz_var_name = 'dz'
        
        pool_structure = [
            {
                'pool_name': 'Vegetation 1',
                    'stock_var': 'Veg1C',
                    'nr_layers': 1,
                    'dz_var': None
            },
            {
                'pool_name': 'CWD',        
                    'stock_var': 'CWDC_vr',
                    'nr_layers': nr_layers, 
                    'dz_var'   : dz_var_name
            },
            {
                'pool_name': 'Litter 1',   
                    'stock_var': 'LITR1C_vr',
                    'nr_layers': nr_layers, 
                    'dz_var'   : dz_var_name
            },
            {
                'pool_name': 'Litter 2',   
                    'stock_var': 'LITR2C_vr',
                    'nr_layers': nr_layers,
                    'dz_var'   : dz_var_name
            },
            {
                'pool_name': 'Litter 3',   
                    'stock_var': 'LITR3C_vr',
                    'nr_layers': nr_layers,
                    'dz_var': dz_var_name
            },
            {
                'pool_name': 'Soil 1',     
                    'stock_var': 'SOIL1C_vr',
                    'nr_layers': nr_layers,
                    'dz_var'   : dz_var_name
            },
            {
                'pool_name': 'Soil 2',     
                    'stock_var': 'SOIL2C_vr',
                    'nr_layers': nr_layers,
                    'dz_var'   : dz_var_name
            },
            {
                'pool_name': 'Soil 3',     
                    'stock_var': 'SOIL3C_vr',
                    'nr_layers': nr_layers,
                    'dz_var'   : dz_var_name
            }
        ]
        
        external_input_structure = {
            'CWD':      ['fire_mortality_c_to_cwdc_col',
                         'gap_mortality_c_to_cwdc_col',
                         'harvest_c_to_cwdc_col'],
            'Litter 1': ['gap_mortality_c_to_litr_met_c_col',
                         'harvest_c_to_litr_met_c_col',
                         'm_c_to_litr_met_fire_col',
                         'phenology_c_to_litr_met_c_col'],
            'Litter 2': ['gap_mortality_c_to_litr_cel_c_col',
                         'harvest_c_to_litr_cel_c_col',
                         'm_c_to_litr_cel_fire_col',
                         'phenology_c_to_litr_cel_c_col'],
            'Litter 3': ['gap_mortality_c_to_litr_lig_c_col',
                         'harvest_c_to_litr_lig_c_col',
                         'm_c_to_litr_lig_fire_col',
                         'phenology_c_to_litr_lig_c_col']
        }
        
        horizontal_structure = {
            ('CWD',      'Litter 2'): 'CWDC_TO_LITR2C_vr',
            ('CWD',      'Litter 3'): 'CWDC_TO_LITR3C_vr',
            ('Litter 1', 'Soil 1'):   'LITR1C_TO_SOIL1C_vr',
            ('Litter 2', 'Soil 1'):   'LITR2C_TO_SOIL1C_vr',
            ('Litter 3', 'Soil 2'):   'LITR3C_TO_SOIL2C_vr',
            ('Soil 1',   'Soil 2'):   'SOIL1C_TO_SOIL2C_vr',
            ('Soil 1',   'Soil 3'):   'SOIL1C_TO_SOIL3C_vr',
            ('Soil 2',   'Soil 1'):   'SOIL2C_TO_SOIL1C_vr',
            ('Soil 2',   'Soil 3'):   'SOIL2C_TO_SOIL3C_vr',
            ('Soil 3',   'Soil 1'):   'SOIL3C_TO_SOIL1C_vr'
        }
        
        vertical_structure = {
            'CWD': 
                {'to_below'  : 'CWD_diffus_down',
                 'from_below': 'CWD_diffus_up',
                 'to_above'  : None,
                 'from_above': None},
             'Litter 1':
                {'to_below'  : 'LITR1_diffus_down',
                 'from_below': 'LITR1_diffus_up',
                 'to_above'  : None,
                 'from_above': None},
             'Litter 2':
                {'to_below'  : 'LITR2_diffus_down',
                 'from_below': 'LITR2_diffus_up',
                 'to_above'  : None,
                 'from_above': None},
             'Litter 3':
                {'to_below'  : 'LITR3_diffus_down',
                 'from_below': 'LITR3_diffus_up',
                 'to_above'  : None,
                 'from_above': None},
             'Soil 1':
                {'to_below'  : 'SOIL1_diffus_down',
                 'from_below': 'SOIL1_diffus_up',
                 'to_above'  : None,
                 'from_above': None},
             'Soil 2':
                {'to_below'  : 'SOIL2_diffus_down',
                 'from_below': 'SOIL2_diffus_up',
                 'to_above'  : None,
                 'from_above': None},
             'Soil 3':
                {'to_below'  : 'SOIL3_diffus_down',
                 'from_below': 'SOIL3_diffus_up',
                 'to_above'  : None,
                 'from_above': None},
        } 
        
        external_output_structure = {
            'CWD'     : ['CWD_HR_L2_vr', 'CWD_HR_L3_vr'],
            'Litter 1': ['LITR1_HR_vr'],
            'Litter 2': ['LITR2_HR_vr'],
            'Litter 3': ['LITR3_HR_vr'],
            'Soil 1':   ['SOIL1_HR_S2_vr', 'SOIL1_HR_S3_vr'],
            'Soil 2':   ['SOIL2_HR_S1_vr', 'SOIL2_HR_S3_vr'],
            'Soil 3':   ['SOIL3_HR_vr']
        }

        
        ms = ModelStructure(
            pool_structure = pool_structure,
            external_input_structure = external_input_structure,
            horizontal_structure = horizontal_structure,
            vertical_structure = vertical_structure,
            external_output_structure = external_output_structure
        )
        ms.nr_layers = nr_layers

        self.ms = ms


    def test_init(self):
        ## test invalid vertical structure
        invalid_vertical_structure = {
         'CWD': 
            {'to_below'  : 'CWD_diffus_down',
             'from_below': None,
             'to_above'  : 'CWD_diffus_up',
             'from_above': None}
        }

        ms = self.ms
        with self.assertRaises(ModelStructureException) as context:
            model_structure = ModelStructure(
                pool_structure            = ms.pool_structure,        
                external_input_structure  = ms.external_input_structure,
                horizontal_structure      = ms.horizontal_structure,
                vertical_structure        = invalid_vertical_structure,
                external_output_structure = ms.external_output_structure
            )


    def test_get_functions(self):
        ms = self.ms

        ## get_nr_layers
        res = ms.get_nr_layers('Soil 2')
        self.assertEqual(res, ms.nr_layers)

        res = ms.get_nr_layers('Vegetation 1')
        self.assertEqual(res, 1)

        with self.assertRaises(KeyError) as context:
            res = ms.get_nr_layers('Inexistent Pool')
        
        ## get_pool_nr
        res = ms.get_pool_nr('Soil 3', 9)
        self.assertEqual(res, 70)

        with self.assertRaises(KeyError):
            res = ms.get_pool_nr('Vegetation 1', 1)

        with self.assertRaises(KeyError) as context:
            res = ms.get_pool_nr('Inexistent Pool', 3)
        
        ## get_pool_name_and_layer_nr
        res = ms.get_pool_name_and_layer_nr(0)
        res2 = {'pool_name': 'Vegetation 1', 'layer_nr': 0}
        self.assertEqual(res, res2)

        res = ms.get_pool_name_and_layer_nr(70)
        res2 = {'pool_name': 'Soil 3', 'layer_nr': 9}
        self.assertEqual(res, res2)

        with self.assertRaises(KeyError):
            res = ms.get_pool_name_and_layer_nr(71)

        ## get_pool_nrs
        res = ms.get_pool_nrs('Vegetation 1')
        self.assertTrue(np.all(res == [0]))

        res = ms.get_pool_nrs('Soil 1')
        self.assertTrue(np.all(res == 41+np.arange(ms.nr_layers)))

        with self.assertRaises(KeyError):
            res = ms.get_pool_nrs('Inexistent pool')

        ## get_pool_nrs_set
        pool_names = ['CWD', 'Soil 1']
        layers = np.arange(5)+2
        res = ms.get_pool_nrs_set(pool_names, layers)
        res2 = [3,4,5,6,7,43,44,45,46,47]
        self.assertTrue(np.all(res == res2))

        pool_names = ['Vegetation 1']
        layers = np.arange(5)+2
        with self.assertRaises(KeyError):
            res = ms.get_pool_nrs_set(pool_names, layers)

        ## get_nr_pools
        res = ms.get_nr_pools()
        self.assertEqual(res, ms.nr_pools)


    def test_get_flux_var_names(self):
        self.assertEqual(
            self.ms.get_flux_var_names()[7],
            'gap_mortality_c_to_litr_cel_c_col'
        )


################################################################################


if __name__ == '__main__':
    unittest.main()


