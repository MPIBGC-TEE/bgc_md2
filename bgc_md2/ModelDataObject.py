import numpy as np
from cf_units import Unit
from cftime import utime
from netCDF4 import Dataset
from sympy import symbols

from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR

from bgc_md2.SmoothModelRunFromData import SmoothModelRunFromData as SMRFD
from bgc_md2.Variable import FixDumbUnits, Variable, StockVariable, FluxVariable


class ModelDataObjectException(Exception):
    def __str__(self): return self.args[0]


def readVariable(**keywords):
    ReturnClass   = keywords.get('ReturnClass', Variable)
    dataset       = keywords['dataset']
    variable_name = keywords['variable_name']
    min_time      = keywords['min_time']
    max_time      = keywords['max_time']
    nr_layers     = keywords['nr_layers']
    lat_arr       = keywords.get('lat_arr', None)
    lon_arr       = keywords.get('lon_arr', None)
    data_shift    = keywords['data_shift']

    var = dataset[variable_name]

    ## check right output of data
    try:
        if ReturnClass == StockVariable:
            if var.cell_methods != 'time: instantaneous':
                raise(ModelDataException('Stock data is not instantaneous'))

        if ReturnClass == FluxVariable:
            if var.cell_methods != 'time: mean':
                raise(ModelDataException('Flux data is not as a mean'))
    except AttributeError:
        s = "'cell_methods' not specified'"
        raise(ModelDataObjectException(s))

    ## read variable depending on dimensions
    ndim = var.ndim
    if (lat_arr is not None) and (lon_arr is not None):
        if ndim == 1 + 1 + 1 + 1:
            data = var[(min_time+data_shift):max_time,:nr_layers,lat_arr,lon_arr]
        elif ndim == 1 + 1 + 1:
            data = var[(min_time+data_shift):max_time,lat_arr,lon_arr]
            data = np.expand_dims(data, axis=1)
        else:
            raise(ModelDataObjectException('Data structure not understood'))
    elif (lat_arr is None) and (lon_arr is None):
        if ndim == 1 + 1:
            data = var[(min_time+data_shift):max_time,:nr_layers]
            data = np.expand_dims(data, axis=2)
            data = np.expand_dims(data, axis=3)
        elif ndim == 1:
            data = var[(min_time+data_shift):max_time]
            data = np.expand_dims(data, axis=1)
            data = np.expand_dims(data, axis=2)
            data = np.expand_dims(data, axis=3)
        else:
            raise(ModelDataObjectException('Data structure not understood'))
    else:
        raise(ModelDataObjectException('Data structure not understood'))

    sdv = ReturnClass(
        data = data,
        unit = var.units
    )

    return sdv


def StockDensityVariable2StockVariable(sdv, dz):
    depth_axis = 1
    sv = sdv.data_mult(dz, depth_axis)

    return sv


def getStockVariable_from_Density(**keywords):
    mdo   = keywords['mdo']
    dz    = keywords['dz']

    dataset    = mdo.dataset
    nstep     = mdo.nstep
    stock_unit = mdo.stock_unit

    sdv = readVariable(        
        ReturnClass = StockVariable,
        dataset = dataset, 
        **keywords
    )
    sv = StockDensityVariable2StockVariable(sdv, dz)

    sv_agg = sv.aggregateInTime(nstep)
    sv_agg.convert(stock_unit)
    return sv_agg


def FluxRateDensityVariable2FluxRateVariable(frdv, dz):
    depth_axis = 1
    frv = frdv.data_mult(dz, depth_axis)

    return frv


def FluxRateVariable2FluxVariable(frv, time):
    time_axis = 0
    dt_data = np.diff(time.data)
    dt = Variable(
        data = dt_data,
        unit = time.unit
    )
    fv = frv.data_mult(dt, time_axis)

    return fv


def getFluxVariable_from_DensityRate(**keywords):
    mdo   = keywords['mdo']
    dz    = keywords['dz']

    dataset    = mdo.dataset
    time       = mdo.time
    nstep      = mdo.nstep
    stock_unit = mdo.stock_unit

    frdv = readVariable(
        ReturnClass = FluxVariable,
        dataset     = dataset,
        **keywords
    )
    frv = FluxRateDensityVariable2FluxRateVariable(frdv, dz)
    fv = FluxRateVariable2FluxVariable(frv, time)

    fv_agg = fv.aggregateInTime(nstep)
    fv_agg.convert(stock_unit)
    return fv_agg


def getFluxVariable_from_Rate(**keywords):
    mdo   = keywords['mdo']
 
    dataset    = mdo.dataset
    time       = mdo.time
    nstep      = mdo.nstep
    stock_unit = mdo.stock_unit

    frv = readVariable(
        ReturnClass = FluxVariable,
        dataset = dataset, 
        **keywords
    )
    fv = FluxRateVariable2FluxVariable(frv, time)

    fv_agg = fv.aggregateInTime(nstep)
    fv_agg.convert(stock_unit)
    return fv_agg


################################################################################


class ModelDataObject(object):
    def __init__(self, **keywords):
        self.model_structure = keywords['model_structure']

        ds_filename     = keywords.get('ds_filename', None)
        self.dataset    = Dataset(ds_filename) if ds_filename is not None\
                            else keywords['dataset']
        stock_unit           = keywords['stock_unit']
        cftime_unit_src      = keywords['cftime_unit_src']
        calendar_src         = keywords['calendar_src']
        time_shift           = keywords.get('time_shift', 0)
        min_time             = keywords.get('min_time', 0)
        max_time             = keywords.get('max_time', None)
        self.nstep           = keywords.get('nstep', 1)
        self.dz_var_names    = keywords.get('dz_var_names', dict())

        udtime_unit     = 'd'
        cftime_unit_tar = 'days since 1850-01-01 00:00:00'
        calendar_tar    = 'noleap'

        self.stock_unit = FixDumbUnits(stock_unit)

        time_data = self.dataset['time'][min_time:max_time]
        time_data += time_shift

        unit_src = utime(cftime_unit_src, calendar_src)
        unit_tar = utime(
            cftime_unit_tar,
            calendar = calendar_tar
        )
        time_data_tar = np.array(
            [unit_tar.date2num(unit_src.num2date(t))\
                 for t in time_data]
        )
        
        self.time = Variable(
            data = time_data_tar,
            unit = udtime_unit
        )
        self.calendar = calendar_tar

        self.min_time = min_time
        self.max_time = max_time

        self.time_agg = self.time.aggregateInTime(self.nstep)
        self.cftime_unit = cftime_unit_tar


    @property
    def plot_times(self):
        return 1850.+self.time_agg.data/365.


    def get_dz(self, pool_name):
        ms = self.model_structure
        dataset = self.dataset

        for item in ms.pool_structure:
            if item['pool_name'] == pool_name:
                dz_var_name = item.get('dz_var', None)

        if dz_var_name is not None:
            dz = self.dz_var_names.get(
                dz_var_name,
                None
            )
            if dz is None:
                dz_var = dataset[dz_var_name]
            
                dz = Variable(
                    name = dz_var_name,
                    data = dz_var[...],
                    unit = dz_var.units
                )
        else:
            dz = Variable(
                name = 'dz_default',
                data = np.array([1]),
                unit = '1'
            )

        return dz


    def load_stocks(self, **keywords):
        func       = keywords['func']
        lat_arr    = keywords['lat_arr']
        lon_arr    = keywords['lon_arr']
 
        ms         = self.model_structure
        min_time   = self.min_time
        max_time   = self.max_time
        time_agg   = self.time_agg
        stock_unit = self.stock_unit

        nlat = 1
        if lat_arr is not None: nlat = len(lat_arr)
        nlon = 1
        if lon_arr is not None: nlon = len(lon_arr)
        xs_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data), nlat, nlon, ms.nr_pools)
            ),
            mask = False
        )
    
        for item in ms.pool_structure:
            pool_name     = item['pool_name']
            variable_name = item['stock_var']
            print(pool_name, variable_name, flush=True)
            nr_layers     = ms.get_nr_layers(pool_name)
            dz = self.get_dz(pool_name)

            sv_pool_agg = func(        
                mdo           = self,
                variable_name = variable_name,
                min_time      = min_time,
                max_time      = max_time,
                nr_layers     = nr_layers,
                dz            = dz,
                **keywords
            )

            for ly in range(nr_layers):
                pool_nr = ms.get_pool_nr(pool_name, ly)
                data = sv_pool_agg.data[:,ly,...]
                xs_data[...,pool_nr] = data

        xs = StockVariable(
            name = 'stocks',
            data = xs_data,
            unit = stock_unit
        )
        return xs

    
    def _load_external_fluxes(self, **keywords):
        func           = keywords['func']
        flux_structure = keywords['flux_structure']
        lat_arr        = keywords['lat_arr']
        lon_arr        = keywords['lon_arr']

        ms         = self.model_structure
        time_agg   = self.time_agg
        stock_unit = self.stock_unit

        nlat = 1
        if lat_arr is not None: nlat = len(lat_arr)
        nlon = 1
        if lon_arr is not None: nlon = len(lon_arr)
        fs_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, nlat, nlon, ms.nr_pools)
            ),
            mask = False
        )
    
        for pool_name, variable_names in flux_structure.items():
            nr_layers = ms.get_nr_layers(pool_name)
            dz = self.get_dz(pool_name)

            fvs_agg = []
            for variable_name in variable_names:
                print(pool_name, variable_name, flush=True)
                fv_agg = func(
                    mdo           = self,
                    variable_name = variable_name,
                    nr_layers     = nr_layers,
                    dz            = dz,
                    **keywords
                )
                fvs_agg.append(fv_agg)        
    
            fv_pool_agg = sum(fvs_agg)
 
            for ly in range(nr_layers):
                pool_nr = ms.get_pool_nr(pool_name, ly)
                data = fv_pool_agg.data[:,ly,...]
                fs_data[...,pool_nr] = data
    
        fs = FluxVariable(
            data = fs_data,
            unit = self.stock_unit
        )
        return fs


    def load_external_input_fluxes(self, **keywords):
        return self._load_external_fluxes(
            flux_structure = self.model_structure.external_input_structure,
            min_time = self.min_time,
            max_time = self.max_time,
            **keywords
        )


    def load_external_output_fluxes(self, **keywords):
        return self._load_external_fluxes(
            flux_structure = self.model_structure.external_output_structure,
            min_time = self.min_time,
            max_time = self.max_time,
            **keywords
        )


    def load_horizontal_fluxes(self, **keywords):
        func       = keywords['func']
        lat_arr    = keywords['lat_arr']
        lon_arr    = keywords['lon_arr']

        ms         = self.model_structure
        min_time   = self.min_time
        max_time   = self.max_time
        time_agg   = self.time_agg
        stock_unit = self.stock_unit

        nlat = 1
        if lat_arr is not None: nlat = len(lat_arr)
        nlon = 1
        if lon_arr is not None: nlon = len(lon_arr)
        HFs_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, nlat, nlon, ms.nr_pools, ms.nr_pools)
            ),
            mask = False
        )
    
        for pools, variable_names in ms.horizontal_structure.items():
            src_pool_name = pools[0]
            tar_pool_name = pools[1]
            print(src_pool_name, tar_pool_name, flush=True)

            src_nr_layers = ms.get_nr_layers(src_pool_name)
            tar_nr_layers = ms.get_nr_layers(tar_pool_name)
            assert(src_nr_layers == tar_nr_layers)
            nr_layers = src_nr_layers

            src_dz = self.get_dz(src_pool_name)
            tar_dz = self.get_dz(tar_pool_name)

            assert(src_dz.name == tar_dz.name)
            dz = src_dz

            fvs_agg = []
            for variable_name in variable_names:
                fv_agg = func(
                    mdo           = self,
                    variable_name = variable_name,
                    min_time      = min_time,
                    max_time      = max_time,
                    nr_layers     = nr_layers,
                    dz            = dz,
                    **keywords
                )
                fvs_agg.append(fv_agg)        

            fv_agg = sum(fvs_agg)

            for ly in range(nr_layers):
                src_pool_nr = ms.get_pool_nr(src_pool_name, ly)
                tar_pool_nr = ms.get_pool_nr(tar_pool_name, ly)
#                print(src_pool_name, tar_pool_name, ly, src_pool_nr, tar_pool_nr, flush=True)
                data = fv_agg.data[:,ly,...]
                HFs_data[...,tar_pool_nr,src_pool_nr] = data

        HFs = FluxVariable(
            name = 'horizontal fluxes',
            data = HFs_data,
            unit = self.stock_unit
        )
        return HFs


    def load_vertical_fluxes(self, **keywords):
        func       = keywords['func']
        lat_arr    = keywords['lat_arr']
        lon_arr    = keywords['lon_arr']

        ms         = self.model_structure
        min_time   = self.min_time
        max_time   = self.max_time
        time_agg   = self.time_agg
        stock_unit = self.stock_unit

        nlat = 1
        if lat_arr is not None: nlat = len(lat_arr)
        nlon = 1
        if lon_arr is not None: nlon = len(lon_arr)
        VFs_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, nlat, nlon, ms.nr_pools, ms.nr_pools)
            ),
            mask = False
        )
        runoffs_up_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, nlat, nlon, ms.nr_pools)
            ),
            mask = False
        )
        runoffs_down_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, nlat, nlon, ms.nr_pools)
            ),
            mask = False
        )
    
        src_ly_shift = {
            'to_below'  :  0,
            'from_below':  1,
            'to_above'  :  0,
            'from_above': -1
        }
        tar_ly_shift = {
            'to_below'  :  1,
            'from_below':  0,
            'to_above'  : -1,
            'from_above':  0
        }

        for pool_name, flux_structure in ms.vertical_structure.items():
            nr_layers = ms.get_nr_layers(pool_name)

            for flux_type in src_ly_shift.keys():
                fvs_agg = []
                variable_names = flux_structure.get(flux_type, None)
                for variable_name in variable_names:
                    print(pool_name, variable_name, flush=True)
                    fv_agg = func(
                        mdo           = self,
                        variable_name = variable_name,
                        min_time      = min_time,
                        max_time      = max_time,
                        nr_layers     = nr_layers,
                        **keywords
                    )
                    fvs_agg.append(fv_agg)

                if fvs_agg != []:
                    fv_agg = sum(fvs_agg)
                    for ly in range(nr_layers):
                        pool_nr = ms.get_pool_nr(pool_name, ly)
    
                        src_ly = ly + src_ly_shift[flux_type]
                        tar_ly = ly + tar_ly_shift[flux_type]
                        if (src_ly in range(nr_layers)) \
                                and (tar_ly in range(nr_layers)):
    
                            src_pool_nr = ms.get_pool_nr(
                                pool_name,
                                src_ly
                            )
                            tar_pool_nr = ms.get_pool_nr(
                                pool_name,
                                tar_ly
                            )
                            data = fv_agg.data[:,ly,...]
                            VFs_data[...,tar_pool_nr,src_pool_nr] = data
                        elif (src_ly in range(nr_layers) \
                                and (tar_ly == -1)):
                            src_pool_nr = ms.get_pool_nr(
                                pool_name,
                                src_ly
                            )
                            data = fv_agg.data[:,ly,...]
                            runoffs_up_data[...,src_pool_nr] = data 
                        elif (src_ly in range(nr_layers) \
                                and (tar_ly == nr_layers)):
                            src_pool_nr = ms.get_pool_nr(
                                pool_name,
                                src_ly
                            )
                            data = fv_agg.data[:,ly,...]
                            runoffs_down_data[...,src_pool_nr] = data 

        VFs = FluxVariable(
            name = 'vertical_fluxes',
            data = VFs_data,
            unit = self.stock_unit
        )
        runoffs_up = FluxVariable(
            name = 'runoffs up',
            data = runoffs_up_data,
            unit = self.stock_unit
        )
        runoffs_down = FluxVariable(
            name = 'runoffs down',
            data = runoffs_down_data,
            unit = self.stock_unit
        )
        return VFs, runoffs_up, runoffs_down


    def load_xs_us_Fs_rs(self, lat_index=None, lon_index = None):
        lat_arr = None
        if lat_index is not None: lat_arr = np.array([lat_index])
        lon_arr = None
        if lon_index is not None: lon_arr = np.array([lon_index])

        xs = self.load_stocks(
            func       = getStockVariable_from_Density,
            lat_arr    = lat_arr,
            lon_arr    = lon_arr,
            data_shift = 0
        )
    
        us = self.load_external_input_fluxes(
            func       = getFluxVariable_from_DensityRate,
            lat_arr    = lat_arr, 
            lon_arr    = lon_arr, 
            data_shift = 1 
        )
    
        HFs = self.load_horizontal_fluxes(
            func       = getFluxVariable_from_DensityRate,
            lat_arr    = lat_arr, 
            lon_arr    = lon_arr, 
            data_shift = 1 
        )
    
        ## we ignore runoffs until we might experience existing ones 
        ## for a model at some point
        ## then we have to decide what to do with them
        VFs, runoffs_up, runoffs_down = self.load_vertical_fluxes(
            func       = getFluxVariable_from_Rate,
            lat_arr    = lat_arr, 
            lon_arr    = lon_arr, 
            data_shift = 1 
        )
#        print(np.where(runoffs_up.data !=0))
#        print(np.where(runoffs_down.data !=0))
#        print(runoffs_down.data[runoffs_down.data!=0])
   
        Fs = HFs + VFs
    
        rs = self.load_external_output_fluxes(
            func       = getFluxVariable_from_DensityRate,
            lat_arr    = lat_arr, 
            lon_arr    = lon_arr, 
            data_shift = 1 
        )

        return xs, us, Fs, rs


    def create_discrete_model_run(self, lat_index=None, lon_index=None):
        out = self.load_xs_us_Fs_rs(lat_index, lon_index)
        xs, us, Fs, rs = out

        start_values = xs.data[0,0,0,:]

        if xs.data.mask.sum() + Fs.data.mask.sum()\
                + us.data.mask.sum() + rs.data.mask.sum() == 0:
            dmr = DMR.reconstruct_from_data(
                self.time_agg.data.filled(),
                start_values.filled(),
                xs.data[:,0,0,...].filled(),
                Fs.data[:,0,0,...].filled(),
                rs.data[:,0,0,...].filled(),
                us.data[:,0,0,...].filled()
            )
        else:
            dmr = None

        return dmr


    def create_model_run(self, lat_index=None, lon_index=None):
        out = self.load_xs_us_Fs_rs(lat_index, lon_index)
        xs, us, Fs, rs = out

        start_values = xs.data[0,0,0,:]

#        print(self.time_agg.data)
#        print(start_values.data)
#        print(xs.data)
#        print(Fs.data)
#        print(rs.data)
#        print(us.data)
#        input()

        times = self.time_agg.data.filled()
#        times = np.arange(len(self.time_agg.data))
        if xs.data.mask.sum() + Fs.data.mask.sum()\
                + us.data.mask.sum() + rs.data.mask.sum() == 0:
            smrfd = SMRFD.reconstruct_from_data(
                symbols('t'),
                times,
                start_values.filled(),
                xs.data[:,0,0,...].filled(),
                Fs.data[:,0,0,...].filled(),
                rs.data[:,0,0,...].filled(),
                us.data[:,0,0,...].filled()
            )
        else:
            smrfd = None

        return smrfd

