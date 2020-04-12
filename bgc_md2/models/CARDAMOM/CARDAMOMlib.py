import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import xarray as xr
from scipy.interpolate import interp1d

from bgc_md2.ModelStructure import ModelStructure
from bgc_md2.ModelDataObject import (
    ModelDataObject,
    getStockVariable_from_Density,
    getFluxVariable_from_DensityRate,
    getFluxVariable_from_Rate
)
from bgc_md2.Variable import Variable

from CompartmentalSystems.discrete_model_run import DMRError

def load_model_structure():
    # labile, leaf, root, wood, litter, and soil

    pool_structure = [
        {'pool_name': 'Labile', 'stock_var': 'Labile'},
        {'pool_name': 'Leaf',   'stock_var': 'Leaf'},
        {'pool_name': 'Root',   'stock_var': 'Root'},
        {'pool_name': 'Wood',   'stock_var': 'Wood'},
        {'pool_name': 'Litter', 'stock_var': 'Litter'},
        {'pool_name': 'Soil',   'stock_var': 'Soil'}
    ]

    external_input_structure = {
        'Labile': ['NPP_to_Labile'],
        'Leaf':   ['NPP_to_Leaf'],
        'Root':   ['NPP_to_Root'],
        'Wood':   ['NPP_to_Wood'],
    }

    horizontal_structure = {
        ('Labile', 'Leaf'):   ['Labile_to_Leaf'],
        ('Leaf',   'Litter'): ['Leaf_to_Litter'],
        ('Wood',   'Soil'):   ['Wood_to_Soil'],
        ('Root',   'Litter'): ['Root_to_Litter'],
        ('Litter', 'Soil'):   ['Litter_to_Soil']
    }

    vertical_structure = {}

    external_output_structure = {
        'Litter': ['Litter_to_RH'],
        'Soil':   ['Soil_to_RH']
    }

    model_structure = ModelStructure(
        pool_structure            = pool_structure,        
        external_input_structure  = external_input_structure,
        horizontal_structure      = horizontal_structure,
        vertical_structure        = vertical_structure,
        external_output_structure = external_output_structure
    )

    return model_structure


def load_mdo(ds):
#    days_per_month = np.array(np.diff(ds.time), dtype='timedelta64[D]').astype(float)

    ms = load_model_structure()            

   
#    ## bring fluxes from gC/m2/day to gC/m2/month
 
    ## all months are supposed to comprise 365.25/12 days
    days_per_month = 365.25/12.0

    time = Variable(
        data = np.arange(len(ds.time)) * days_per_month,
        unit = 'd'
    )

    mdo = ModelDataObject(
        model_structure = ms,
        dataset = ds,
        stock_unit = 'gC/m2',
        time = time
    )

    return mdo


def load_dmr_14C(dmr):
    ## create 14C dmr

    # compute 14C external input
    atm_delta_14C = np.loadtxt(
#        '/home/hmetzler/Desktop/CARDAMOM/C14Atm_NH.csv',
        'C14Atm_NH.csv',
        skiprows  = 1,
        delimiter = ','
    )
    _F_atm_delta_14C = interp1d(
        atm_delta_14C[:,0],
        atm_delta_14C[:,1],
        fill_value = 'extrapolate'
    )
    t_conv = lambda t: 2001 + (15+t)/365.25
    F_atm_delta_14C = lambda t: _F_atm_delta_14C(t_conv(t))

    alpha = 1.18e-12
    us_12C = dmr.external_input_vector
    with np.errstate(divide='ignore'):
        us_14C = alpha * us_12C * (
            1 + 1/1000 * F_atm_delta_14C(dmr.times[:-1]).reshape(-1,1)
        )
    np.nan_to_num(us_14C, posinf=0, copy=False)
    

    # compute 14C start_values
    lamda = 0.0001209681 


#    # V1: assume system at t0 in eq.
#    B0_14C = np.matmul(
#        dmr.Bs[0],
#        np.exp(-lamda*365.25/12)*np.identity(dmr.nr_pools)
#    )
#    start_values_14C = np.linalg.solve(
#        (np.eye(dmr.nr_pools)-B0_14C),
#        us_14C[0]
#    )
#
#    start_values_14C = 30 * alpha * np.ones(6)


#    # V2: problem for pools with no external input
#    ks = np.diag(np.mean(dmr.Bs, axis=0))
#    start_values_14C = np.mean(us_14C, axis=0)/(1-ks*np.exp(-lamda*365.25/12))

    # V1: mean of 14C Bs and mean of 12C us
    B0_14C = np.mean([np.matmul(
        dmr.Bs[k],
        np.exp(-lamda*365.25/12)*np.identity(dmr.nr_pools)
    ) for k in range(len(dmr.Bs))], axis=0)
    start_values_14C = np.linalg.solve(
        (np.eye(dmr.nr_pools)-B0_14C),
        np.mean(us_14C, axis=0)
    )

    dmr_14C = dmr.to_14C_only(
        start_values_14C,
        us_14C
    )   

    return dmr_14C

 
def plot_Delta_14C_in_time(ms, dmr, dmr_14C):
    ## plot Delta 14C profiles
    soln = dmr.solve()
    soln_14C = dmr_14C.solve()

    alpha = 1.18e-12
    t_conv = lambda t: 2001 + (15+t)/365.25
    F_Delta_14C = lambda C14, C12: (C14/C12/alpha - 1) * 1000

    fig, axes = plt.subplots(
        figsize = (14,7),
        nrows   = 2,
        ncols   = 3,
        sharex  = True,
        sharey  = True
    )
   
    for nr, pool_name in enumerate(ms.pool_names):
        ax = axes.flatten()[nr]
        Delta_14C = F_Delta_14C(soln_14C[:,nr], soln[:,nr])
        ax.plot(t_conv(dmr.times), Delta_14C)
        ax.set_title(ms.pool_names[nr])
    
    plt.show()
    plt.close(fig)


def create_dmr_14C_dataset(ms, ds, dmr_14C):
    ## create 14C dataset
    data_vars = {}

    def add_variable(data_vars, var_name_ds, data):   
        if var_name_ds in ds.variables.keys():
            var_ds = ds[var_name_ds] 
            attrs  = var_ds.attrs
            attrs['units'] = 'per mille'
            var = xr.DataArray(
                data   = data,
                coords = var_ds.coords,
                dims   = var_ds.dims,
                attrs  = attrs
            )
            data_vars[var_name_ds] = var


    ## save stocks
    for pool_nr, pool_name  in enumerate(ms.pool_names):
        var_name_ds = pool_name
        data = dmr_14C.solve()[:, pool_nr]
        add_variable(data_vars, var_name_ds, data)


    ## save external input fluxes

    # insert np.nan at time t0
    us_monthly = Variable(
        data = np.concatenate(
            [
                np.nan * np.ones((1,len(ms.pool_names))),
                dmr_14C.external_input_vector
            ],
            axis = 0
        ),
        unit = 'g/(365.25/12 d)'
    )

    # convert to daily flux rates
    us = us_monthly.convert('g/d')
    for pool_nr, pool_name  in enumerate(ms.pool_names):        
        var_name_ds = 'NPP_to_' + pool_name
        data = us.data[:, pool_nr]
        add_variable(data_vars, var_name_ds, data)


    ## save internal fluxes

    # insert np.nan at time t0
    Fs_monthly = Variable(
        data = np.concatenate(
            [
                np.nan * np.ones(
                    (1, len(ms.pool_names), len(ms.pool_names))
                ),
                dmr_14C.internal_flux_matrix
            ],
            axis = 0
        ),
        unit = 'g/(365.25/12 d)'
    )

    # convert to daily flux rates
    Fs = Fs_monthly.convert('g/d')
    for pnr_from, pn_from in enumerate(ms.pool_names):        
        for pnr_to, pn_to in enumerate(ms.pool_names):        
            var_name_ds = pn_from +'_to_' + pn_to
            data = Fs.data[:, pnr_to, pnr_from]
            add_variable(data_vars, var_name_ds, data)


    ## save external output fluxes

    # insert np.nan at time t0
    rs_monthly = Variable(
        data = np.concatenate(
            [
                np.nan * np.ones((1,len(ms.pool_names))),
                dmr_14C.external_output_vector
            ],
            axis = 0
        ),
        unit = 'g/(365.25/12 d)'
    )

    # convert to daily flux rates
    rs = rs_monthly.convert('g/d')
    for pool_nr, pool_name  in enumerate(ms.pool_names):        
        var_name_ds = pool_name + '_to_RH'
        data = rs.data[:, pool_nr]
        add_variable(data_vars, var_name_ds, data)


    ds_dmr_14C = xr.Dataset(
        data_vars = data_vars,
        coords    = ds.coords,
        attrs     = ds.attrs
    )

    return ds_dmr_14C


def load_dmr_14C_dataset(ds):
#    return ds 
    # fake return data struicture on first call (with empty data)
    if ds.ens.values.shape == (0,):
        empty_var = xr.DataArray(
            data   = np.ndarray(dtype=float, shape=(0,0,0)),
            dims   = ('ens', 'lat' , 'lon')
        )
        ds['max_abs_err'] = empty_var
        ds['max_rel_err'] = empty_var
        ds['log'] = xr.DataArray(
            data = np.ndarray(dtype='<U50', shape=(0,0,0)),
            dims = ('ens', 'lat', 'lon')
        )    
        return ds
    

    log = ''
    try:
        mdo = load_mdo(ds)
    
        dmr, abs_err, rel_err = mdo.create_discrete_model_run(errors=True)
        dmr_14C = load_dmr_14C(dmr)
    
    
        ms = mdo.model_structure
    #    plot_Delta_14C_in_time(ms, dmr, dmr_14C)
        ds_dmr_14C = create_dmr_14C_dataset(ms, ds, dmr_14C)
    
        ## add reconstruction error
        var_abs_err = xr.DataArray(
            data  = np.max(abs_err.data),
            attrs = {
                'units':     abs_err.unit,
                'long_name': 'max. abs. error on reconstructed stock sizes'
            }
        )
        ds_dmr_14C['max_abs_err'] = var_abs_err

        var_rel_err = xr.DataArray(
            data  = np.max(rel_err.data),
            attrs = {
                'units':     rel_err.unit,
                'long_name': 'max. rel. error on reconstructed stock sizes'
            }
        )
        ds_dmr_14C['max_rel_err'] = var_rel_err
    
    except DMRError as e:
        log = str(e)
        
        data_vars = {}
        for key, val in ds.data_vars.items():
            if key != 'time':
                data_vars[key] = np.nan * val
        ds_dmr_14C = ds.copy(data=data_vars)
        
        ds_dmr_14C['max_abs_err'] = np.nan
        ds_dmr_14C['max_rel_err'] = np.nan

    ds_dmr_14C['log'] = log
    ds_dmr_14C.close()    

    return ds_dmr_14C


################################################################################


if __name__ == '__main__':
    pass
#    dataset = xr.open_dataset('~/Desktop/CARDAMOM/cardamom_for_holger.nc')
#    ds = dataset.isel(ens=0, lat=0, lon=0)
#    ds_dmr_14C = load_dmr_14C_dataset(ds)
#
#    ds.close()
#    dataset.close()


