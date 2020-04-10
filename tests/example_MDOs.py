import numpy as np
from netCDF4 import Dataset

from bgc_md2.ModelDataObject import ModelDataObject
from bgc_md2.ModelStructure import ModelStructure
from bgc_md2.Variable import Variable


def MDO_3pools_3layers():
    ## set up test model strcture
    nr_layers = 3
    dz_var_name = 'dz'

    pool_structure = [
        {
            'pool_name': 'CWD',        
                'stock_var': 'CWDC',
                'nr_layers': nr_layers, 
                'dz_var'   : dz_var_name
        },
        {
            'pool_name': 'Litter',   
                'stock_var': 'LITRC',
                'nr_layers': nr_layers, 
                'dz_var'   : dz_var_name
        },
        {
            'pool_name': 'Soil',   
                'stock_var': 'SOILC',
                'nr_layers': nr_layers,
                'dz_var'   : dz_var_name
        }
    ]

    external_input_structure = {
        'CWD':    ['fire_CWD',
                   'gap_CWD',
                   'harvest_CWD'],
        'Litter': ['gap_LITR',
                   'harvest_LITR',
                   'm_c_LITR',
                   'phenology_LITR'],
    }

    horizontal_structure = {
        ('CWD',    'Litter'): ['CWD_TO_LITR'],
        ('Litter', 'Soil')  : ['LITR_TO_SOIL'],
    }

    vertical_structure = {
         'Litter':
            {'to_below'  : ['LITR_flux_down_tb'],
             'from_below': ['LITR_flux_up_fb'],
             'to_above'  : [],
             'from_above': []},
         'Soil':
            {'to_below'  : [],
             'from_below': [],
             'to_above'  : ['SOIL_flux_up_ta'],
             'from_above': ['SOIL_flux_down_fa']},
    } 

    external_output_structure = {
        'CWD'   : ['CWD_HR1', 'CWD_HR2'],
        'Litter': ['LITR_HR'],
        'Soil' :  ['SOIL_HR1', 'SOIL_HR2'],
    }

    model_structure = ModelStructure(
        pool_structure            = pool_structure,        
        external_input_structure  = external_input_structure,
        horizontal_structure      = horizontal_structure,
        vertical_structure        = vertical_structure,
        external_output_structure = external_output_structure
    )

    ## set up a test dataset
    time = np.arange(10)
    dz = np.arange(nr_layers)+1

    ds = Dataset('test_dataset.nc', 'w', diskless=True, persist=False)
    ds.createDimension('time', len(time))
    ds.createDimension('level', nr_layers)

    time_var = ds.createVariable('time', 'f8', ('time',))
    time_var[:] = time
    time_var.units = 'd'

    dz_var = ds.createVariable('dz', 'f4', ('level'))
    dz_var[:] = dz
    dz_var.units = 'm'


    ## introduce some test variables to test dataset

    ## stocks
    var = ds.createVariable(
        'CWDC',
        'f8',
        ('time','level')
    )
    data = np.arange(len(time)*nr_layers)\
        .reshape((len(time),nr_layers))
    var[...] = data
    var.units = 'gC/m^3'
    var.cell_methods = 'time: instantaneous'

    var = ds.createVariable(
        'LITRC',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.1)\
        .reshape((len(time),nr_layers))
    var[...] = data
    var.units = 'gC/m^3'
    var.cell_methods = 'time: instantaneous'

    var = ds.createVariable(
        'SOILC',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.2)\
        .reshape((len(time),nr_layers))
    var[...] = data.copy()
    var.units = 'gC/m^3'
    var.cell_methods = 'time: instantaneous'


    ## external input fluxes
    var = ds.createVariable(
        'fire_CWD',
        'f8',
        ('time','level')
    )
    data = np.arange(len(time)*nr_layers)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'gap_CWD',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.01)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'harvest_CWD',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.02)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'gap_LITR',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.10)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'harvest_LITR',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.11)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'm_c_LITR',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.12)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'phenology_LITR',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.13)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'


    ## horizontal fluxes
    var = ds.createVariable(
        'CWD_TO_LITR',
        'f8',
        ('time','level')
    )
    data = np.arange(len(time)*nr_layers)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'LITR_TO_SOIL',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.1)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'


    ## vertical fluxes
    var = ds.createVariable(
        'LITR_flux_down_tb',
        'f8',
        ('time','level')
    )
    data = np.arange(len(time)*nr_layers)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^2/s'
    var.cell_methods = 'time: mean'
    
    var = ds.createVariable(
        'LITR_flux_up_fb',
        'f8',
        ('time','level')
    )
    data = np.arange(len(time)*nr_layers)\
        .reshape((len(time),nr_layers))*1e-04
    var[...] = data
    var.units = 'gC/m^2/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'SOIL_flux_down_fa',
        'f8',
        ('time','level')
    )
    data = np.arange(len(time)*nr_layers)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^2/s'
    var.cell_methods = 'time: mean'
    
    var = ds.createVariable(
        'SOIL_flux_up_ta',
        'f8',
        ('time','level')
    )
    data = np.arange(len(time)*nr_layers)\
        .reshape((len(time),nr_layers))*1e-04
    var[...] = data
    var.units = 'gC/m^2/s'
    var.cell_methods = 'time: mean'
    
    
    ## external_output_fluxes
    var = ds.createVariable(
        'CWD_HR1',
        'f8',
        ('time','level')
    )
    data = np.arange(len(time)*nr_layers)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'CWD_HR2',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.01)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'LITR_HR',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.10)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'SOIL_HR1',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.20)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'SOIL_HR2',
        'f8',
        ('time','level')
    )
    data = (np.arange(len(time)*nr_layers)+0.21)\
        .reshape((len(time),nr_layers))*1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'


    ## set up a ModelDataObject
    nstep            = 2
    stock_unit       = 'g/m^2'
#    cftime_unit_src  = 'days since 1851-01-01 00:00:00'
#    calendar_src     = 'noleap'
#    time_shift       = 5
#    max_time         = 8


    time = Variable(
        data = time,
        unit = 'd'
    )
    mdo = ModelDataObject(
        model_structure = model_structure, 
        dataset         = ds,
        nstep           = nstep,
        stock_unit      = stock_unit,
        time            = time
#        cftime_unit_src = cftime_unit_src,
#        calendar_src    = calendar_src,
#        time_shift      = time_shift,
#        max_time        = max_time
    )
    return mdo


################################################################################


def MDO_3pools_nolayers_nogrid():
    ## set up test model strcture

    pool_structure = [
        {
            'pool_name': 'CWD',        
                'stock_var': 'CWDC',
        },
        {
            'pool_name': 'Litter',   
                'stock_var': 'LITRC',
        },
        {
            'pool_name': 'Soil',   
                'stock_var': 'SOILC',
        }
    ]

    external_input_structure = {
        'CWD':    ['fire_CWD'],
    }

    horizontal_structure = {
        ('CWD','Litter'): ['CWD_TO_LITR'],
    }

    external_output_structure = {
        'CWD': ['CWD_HR'],
    }

    model_structure = ModelStructure(
        pool_structure            = pool_structure,        
        external_input_structure  = external_input_structure,
        horizontal_structure      = horizontal_structure,
        external_output_structure = external_output_structure
    )

    ## set up a test dataset
    time = np.arange(10)

    ds = Dataset('test_dataset2.nc', 'w', diskless=True, persist=False)
    ds.createDimension('time', len(time))

    time_var = ds.createVariable('time', 'f8', ('time',))
    time_var[:] = time
    time_var.units = 'd'


    ## introduce some test variables to test dataset

    ## stocks
    var = ds.createVariable(
        'CWDC',
        'f8',
        ('time',)
    )
    data = np.arange(len(time))
    var[...] = data
    var.units = 'gC/m^2'
    var.cell_methods = 'time: instantaneous'

    var = ds.createVariable(
        'LITRC',
        'f8',
        ('time',)
    )
    data = np.arange(len(time))+0.1
    var[...] = data
    var.units = 'gC/m^2'
    var.cell_methods = 'time: instantaneous'

    var = ds.createVariable(
        'SOILC',
        'f8',
        ('time',)
    )
    data = np.arange(len(time))+0.2
    var[...] = data.copy()
    var.units = 'gC/m^2'
    var.cell_methods = 'time: instantaneous'


    ## external input fluxes
    var = ds.createVariable(
        'fire_CWD',
        'f8',
        ('time',)
    )
    data = np.arange(len(time))*1e-03
    var[...] = data
    var.units = 'gC/m^2/s'
    var.cell_methods = 'time: mean'


    ## horizontal fluxes
    var = ds.createVariable(
        'CWD_TO_LITR',
        'f8',
        ('time',)
    )
    data = np.arange(len(time))*1e-03
    var[...] = data
    var.units = 'gC/m^2/s'
    var.cell_methods = 'time: mean'


    ## external_output_fluxes
    var = ds.createVariable(
        'CWD_HR',
        'f8',
        ('time',)
    )
    data = np.arange(len(time))
    var[...] = data
    var.units = 'gC/m^2/s'
    var.cell_methods = 'time: mean'


    ## set up a ModelDataObject
    nstep            = 2
    stock_unit       = 'g/m^2'
#    cftime_unit_src  = 'days since 1850-01-01 00:00:00'
#    calendar_src     = 'noleap'
#    max_time         = 8

    time = Variable(
        data = time,
        unit = 'd'
    )
    mdo = ModelDataObject(
        model_structure = model_structure, 
        dataset         = ds,
        nstep           = nstep,
        stock_unit      = stock_unit,
        time            = time
#        cftime_unit_src = cftime_unit_src,
#        calendar_src    = calendar_src,
#        max_time        = max_time
    )
    return mdo


################################################################################


def MDO_3pools_3layers_nogrid():
    ## set up test model strcture
    
    nr_layers = 3
    dz_var_name = 'dz'

    pool_structure = [
        {
            'pool_name': 'CWD',        
                'stock_var': 'CWDC',
                'nr_layers': nr_layers,
                'dz_var': dz_var_name
        },
        {
            'pool_name': 'Litter',   
                'stock_var': 'LITRC',
                'nr_layers': nr_layers,
                'dz_var': dz_var_name
        },
        {
            'pool_name': 'Soil',   
                'stock_var': 'SOILC',
                'nr_layers': nr_layers,
                'dz_var': dz_var_name
        }
    ]

    external_input_structure = {
        'CWD':    ['fire_CWD'],
    }

    horizontal_structure = {
        ('CWD','Litter'): ['CWD_TO_LITR']
    }

    external_output_structure = {
        'CWD': ['CWD_HR'],
    }

    model_structure = ModelStructure(
        pool_structure            = pool_structure,        
        external_input_structure  = external_input_structure,
        horizontal_structure      = horizontal_structure,
        external_output_structure = external_output_structure
    )

    ## set up a test dataset
    time = np.arange(10)
    dz = np.arange(nr_layers)+1

    ds = Dataset('test_dataset4.nc', 'w', diskless=True, persist=False)
    ds.createDimension('time', len(time))
    ds.createDimension('level', nr_layers)

    time_var = ds.createVariable('time', 'f8', ('time',))
    time_var[:] = time
    time_var.units = 'd'

    dz_var = ds.createVariable('dz', 'f4', ('level'))
    dz_var[:] = dz
    dz_var.units = 'm'


    ## introduce some test variables to test dataset

    ## stocks
    var = ds.createVariable(
        'CWDC',
        'f8',
        ('time', 'level')
    )
    data = np.arange(len(time)*nr_layers)
    var[...] = data.reshape(len(time), nr_layers)
    var.units = 'gC14/m^3'
    var.cell_methods = 'time: instantaneous'

    var = ds.createVariable(
        'LITRC',
        'f8',
        ('time', 'level')
    )
    data = np.arange(len(time)*nr_layers)+0.1
    var[...] = data.reshape(len(time), nr_layers)
    var.units = 'gC14/m^3'
    var.cell_methods = 'time: instantaneous'

    var = ds.createVariable(
        'SOILC',
        'f8',
        ('time', 'level')
    )
    data = np.arange(len(time)*nr_layers)+0.2
    var[...] = data.reshape(len(time), nr_layers)
    var.units = 'gC14/m^3'
    var.cell_methods = 'time: instantaneous'


    ## external input fluxes
    var = ds.createVariable(
        'fire_CWD',
        'f8',
        ('time', 'level')
    )
    masked_data = np.ma.masked_array(
        data = np.arange(len(time)*nr_layers)*1e-03,
        mask = np.zeros(len(time)*nr_layers)
    )
    masked_data.mask[21] = 1
    var[...] = masked_data.reshape(len(time), nr_layers)
    var.units = 'gC14/m^3/s'
    var.cell_methods = 'time: mean'


    ## horizontal fluxes
    var = ds.createVariable(
        'CWD_TO_LITR',
        'f8',
        ('time', 'level')
    )
    data = np.arange(len(time)*nr_layers)*1e-03
    var[...] = data.reshape(len(time), nr_layers)
    var.units = 'gC14/m^3/s'
    var.cell_methods = 'time: mean'


    ## external_output_fluxes
    var = ds.createVariable(
        'CWD_HR',
        'f8',
        ('time', 'level')
    )
    data = np.arange(len(time)*nr_layers)
    var[...] = data.reshape(len(time), nr_layers)
    var.units = 'gC14/m^3/s'
    var.cell_methods = 'time: mean'


    ## set up a ModelDataObject
    nstep            = 2
    stock_unit       = 'g/m^2'
#    cftime_unit_src  = 'days since 1900-01-01 00:00:00'
#    calendar_src     = 'noleap'
#    max_time         = 8

    time = Variable(
        data = time,
        unit = 'd'
    )
    mdo = ModelDataObject(
        model_structure = model_structure, 
        dataset         = ds,
        nstep           = nstep,
        stock_unit      = stock_unit,
        time            = time
#        cftime_unit_src = cftime_unit_src,
#        calendar_src    = calendar_src,
#        max_time        = max_time
    )
    return mdo


##### discretizable models #####


def MDO_3pools_3layers_discretizable():
    ## set up test model strcture
    nr_layers = 3
    dz_var_name = 'dz'

    pool_structure = [
        {
            'pool_name': 'CWD',        
                'stock_var': 'CWDC',
                'nr_layers': nr_layers, 
                'dz_var'   : dz_var_name
        },
        {
            'pool_name': 'Litter',   
                'stock_var': 'LITRC',
                'nr_layers': nr_layers, 
                'dz_var'   : dz_var_name
        },
        {
            'pool_name': 'Soil',   
                'stock_var': 'SOILC',
                'nr_layers': nr_layers,
                'dz_var'   : dz_var_name
        }
    ]

    external_input_structure = {
        'CWD':    ['input_CWD'],
        'Litter': ['input_LITR']
    }

    horizontal_structure = {
        ('CWD',    'Litter'): ['CWD_TO_LITR'],
        ('Litter', 'Soil')  : ['LITR_TO_SOIL'],
    }

    external_output_structure = {
        'CWD'   : ['CWD_HR'],
        'Litter': ['LITR_HR'],
        'Soil'  : ['SOIL_HR'],
    }

    model_structure = ModelStructure(
        pool_structure            = pool_structure,        
        external_input_structure  = external_input_structure,
        horizontal_structure      = horizontal_structure,
#        vertical_structure        = vertical_structure,
        external_output_structure = external_output_structure
    )

    ## set up a test dataset
    time = np.arange(10)
    dz = np.arange(nr_layers)+1

    ds = Dataset(
        'test_dataset_discretizable1.nc', 
        'w', 
        diskless=True, 
        persist=False
    )
    ds.createDimension('time', len(time))
    ds.createDimension('level', nr_layers)

    time_var = ds.createVariable('time', 'f8', ('time',))
    time_var[:] = time
    time_var.units = 'd'

    dz_var = ds.createVariable('dz', 'f4', ('level'))
    dz_var[:] = dz
    dz_var.units = 'm'


    ## introduce some test variables to test dataset

    ## stocks
    var = ds.createVariable(
        'CWDC',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers)) * 1e05
    var[...] = data
    var.units = 'gC/m^3'
    var.cell_methods = 'time: instantaneous'

    var = ds.createVariable(
        'LITRC',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers)) * 1e05
    var[...] = data
    var.units = 'gC/m^3'
    var.cell_methods = 'time: instantaneous'

    var = ds.createVariable(
        'SOILC',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers,)) *1e05
    var[...] = data.copy()
    var.units = 'gC/m^3'
    var.cell_methods = 'time: instantaneous'


    ## external input fluxes
    var = ds.createVariable(
        'input_CWD',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers)) * 1e-02
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'input_LITR',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers)) * 1e-02
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    ## horizontal fluxes
    var = ds.createVariable(
        'CWD_TO_LITR',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers)) * 1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'LITR_TO_SOIL',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers)) * 1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    ## external_output_fluxes
    var = ds.createVariable(
        'CWD_HR',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers)) * 2 * 1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'LITR_HR',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers)) * 2 * 1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    var = ds.createVariable(
        'SOIL_HR',
        'f8',
        ('time','level')
    )
    data = np.ones((len(time),nr_layers)) * 2 * 1e-03
    var[...] = data
    var.units = 'gC/m^3/s'
    var.cell_methods = 'time: mean'

    ## set up a ModelDataObject
    nstep            = 2
    stock_unit       = 'g/m^2'
#    cftime_unit_src  = 'days since 1850-01-01 00:00:00'
#    calendar_src     = 'noleap'
#    max_time         = 8

    time = Variable(
        data = time,
        unit = 'd'
    )
    mdo = ModelDataObject(
        model_structure = model_structure, 
        dataset         = ds,
        nstep           = nstep,
        stock_unit      = stock_unit,
        time            = time
#        cftime_unit_src = cftime_unit_src,
#        calendar_src    = calendar_src,
#        max_time        = max_time
    )
    return mdo

