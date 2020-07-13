import numpy as np
from netCDF4 import Dataset
from sympy import symbols

# from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
from CompartmentalSystems.discrete_model_run_with_gross_fluxes import\
    DiscreteModelRunWithGrossFluxes as DMRWGF

from CompartmentalSystems.pwc_model_run_fd import PWCModelRunFD as PWCMRFD
from bgc_md2.Variable import (
    FixDumbUnits,
    Variable,
    StockVariable,
    FluxVariable
)


class ModelDataObjectException(Exception):
    def __str__(self): return self.args[0]


def readVariable(**keywords):
    ReturnClass   = keywords.get('ReturnClass', Variable)
    dataset       = keywords['dataset']
    variable_name = keywords['variable_name']
    nr_layers     = keywords['nr_layers']
    data_shift    = keywords['data_shift']

    var = dataset[variable_name]

    ## check right output format of data
    try:
        if ReturnClass == StockVariable:
            if var.cell_methods != 'time: instantaneous':
                raise(ModelDataObjectException(
                    'Stock data is not instantaneous'
                ))

        if ReturnClass == FluxVariable:
            if var.cell_methods != 'time: mean':
                raise(ModelDataObjectException('Flux data is not as a mean'))
    except AttributeError:
        s = "'cell_methods' not specified"
        raise(ModelDataObjectException(s))

    ## read variable depending on dimensions
    ndim = var.ndim
    if ndim == 1 + 1:
        # time and depth
        data = var[data_shift:, :nr_layers]
    elif ndim == 1:
        # only time
        data = var[data_shift:]
        # add artificial depth axis
        data = np.expand_dims(data, axis=1)
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
        self.stock_unit = FixDumbUnits(stock_unit)
        self.nstep           = keywords.get('nstep', 1)
        self.dz_var_names    = keywords.get('dz_var_names', dict())

        self.time = keywords['time']
        self.time_agg = self.time.aggregateInTime(self.nstep)


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
 
        ms         = self.model_structure
        time_agg   = self.time_agg
        stock_unit = self.stock_unit

        xs_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data), ms.nr_pools)
            ),
            mask = False
        )
    
        for item in ms.pool_structure:
            pool_name     = item['pool_name']
            variable_name = item['stock_var']
            #print(pool_name, variable_name, flush=True)
            nr_layers     = ms.get_nr_layers(pool_name)
            dz = self.get_dz(pool_name)

            sv_pool_agg = func(        
                mdo           = self,
                variable_name = variable_name,
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

        ms         = self.model_structure
        time_agg   = self.time_agg
        stock_unit = self.stock_unit

        fs_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, ms.nr_pools)
            ),
            mask = False
        )
    
        for pool_name, variable_names in flux_structure.items():
            nr_layers = ms.get_nr_layers(pool_name)
            dz = self.get_dz(pool_name)

            fvs_agg = []
            for variable_name in variable_names:
                #print(pool_name, variable_name, flush=True)
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
            **keywords
        )


    def load_external_output_fluxes(self, **keywords):
        return self._load_external_fluxes(
            flux_structure = self.model_structure.external_output_structure,
            **keywords
        )


    def load_horizontal_fluxes(self, **keywords):
        func       = keywords['func']

        ms         = self.model_structure
        time_agg   = self.time_agg
        stock_unit = self.stock_unit

        HFs_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, ms.nr_pools, ms.nr_pools)
            ),
            mask = False
        )
    
        for pools, variable_names in ms.horizontal_structure.items():
            src_pool_name = pools[0]
            tar_pool_name = pools[1]
            #print(src_pool_name, tar_pool_name, flush=True)

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

        ms         = self.model_structure
        time_agg   = self.time_agg
        stock_unit = self.stock_unit

        VFs_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, ms.nr_pools, ms.nr_pools)
            ),
            mask = False
        )
        runoffs_up_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, ms.nr_pools)
            ),
            mask = False
        )
        runoffs_down_data = np.ma.masked_array(
            data = np.zeros(
                (len(time_agg.data)-1, ms.nr_pools)
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


    def load_xs_Us_Fs_Rs(self):
        xs = self.load_stocks(
            func       = getStockVariable_from_Density,
            data_shift = 0
        )
    
        Us = self.load_external_input_fluxes(
            func       = getFluxVariable_from_DensityRate,
            data_shift = 1 
        )
    
        HFs = self.load_horizontal_fluxes(
            func       = getFluxVariable_from_DensityRate,
            data_shift = 1 
        )
    
        ## we ignore runoffs until we might experience existing ones 
        ## for a model at some point
        ## then we have to decide what to do with them
        VFs, Runoffs_up, Runoffs_down = self.load_vertical_fluxes(
            func       = getFluxVariable_from_Rate,
            data_shift = 1 
        )
#        print(np.where(runoffs_up.data !=0))
#        print(np.where(runoffs_down.data !=0))
#        print(runoffs_down.data[runoffs_down.data!=0])
   
        Fs = HFs + VFs
    
        Rs = self.load_external_output_fluxes(
            func       = getFluxVariable_from_DensityRate,
            data_shift = 1 
        )

        return xs, Us, Fs, Rs


    def create_discrete_model_run(self, errors=False):
        out = self.load_xs_Us_Fs_Rs()
        xs, Us, Fs, Rs = out
        start_values = xs.data[0,:]

        if xs.data.mask.sum() + Fs.data.mask.sum()\
                + Us.data.mask.sum() + Rs.data.mask.sum() == 0:

            #fixme hm 2020-04-21:
            # which reconstruction version is the right choice?
#            dmr = DMR.reconstruct_from_data(
#                self.time_agg.data.filled(),
#                start_values.filled(),
#                xs.data.filled(),
#                Fs.data.filled(),
#                Rs.data.filled(),
#                Us.data.filled()
#            )

#            dmr = DMR.reconstruct_from_fluxes_and_solution(
            dmr = DMRWGF.reconstruct_from_fluxes_and_solution(
                self.time_agg.data.filled(),
                xs.data.filled(),
                Fs.data.filled(),
                Rs.data.filled(),
                Us.data.filled(),
                Fs.data.filled(),
                Rs.data.filled()
            )
        else:
            dmr = None

        if errors:
            soln_dmr = Variable(
                data = dmr.solve(),
                unit = self.stock_unit
            )
            abs_err = soln_dmr.absolute_error(xs)
            rel_err = soln_dmr.relative_error(xs)

            return dmr, abs_err, rel_err
        else:
            return dmr


    def create_model_run(self, errors=False):
        out = self.load_xs_Us_Fs_Rs()
        xs, Us, Fs, Rs = out

#        print(self.time_agg.data)
#        print(xs.data)
#        print(Fs.data)
#        print(Rs.data)
#        print(Us.data)
#        input()

        times = self.time_agg.data.filled()
#        times = np.arange(len(self.time_agg.data))
        if xs.data.mask.sum() + Fs.data.mask.sum()\
                + Us.data.mask.sum() + Rs.data.mask.sum() == 0:
            pwc_mr_fd = PWCMRFD(
                symbols('t'),
                times,
                xs.data.filled()[0],
#                xs.data.filled(),
                Us.data.filled(),
                Fs.data.filled(),
                Rs.data.filled()
            )
        else:
            pwc_mr_fd = None

        if errors:
            soln = pwc_mr_fd.solve()

            soln_pwc_mr_fd = Variable(
                data = soln,
                unit = self.stock_unit
            )
            abs_err = soln_pwc_mr_fd.absolute_error(xs)
            rel_err = soln_pwc_mr_fd.relative_error(xs)

            print(soln[:, 0]-self.get_stock(pwc_mr_fd, 'Labile'))

            return pwc_mr_fd, abs_err, rel_err
        else:
            return pwc_mr_fd

    def get_stock(self, mr, pool_name, nr_layer=0):
        pool_nr = self.model_structure.get_pool_nr(pool_name, nr_layer)
        soln = mr.solve()
        # make it a Variable
        return Variable(
            data=soln[:, pool_nr],
            unit=self.stock_unit
        )
