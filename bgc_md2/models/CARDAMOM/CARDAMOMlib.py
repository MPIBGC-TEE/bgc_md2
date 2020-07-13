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
from CompartmentalSystems.discrete_model_run_14C import DiscreteModelRun_14C as DMR_14C
from CompartmentalSystems.pwc_model_run_14C import PWCModelRun_14C as PWCMR_14C
from .compute_start_values_14C import compute_start_values_14C


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


def plot_Delta_14C_in_time(ms, soln, soln_14C):
    ## plot Delta 14C profiles
    times = np.arange(soln.shape[0])

    alpha = 1.18e-12
    t_conv = lambda t: 2001 + (15+t)/365.25
    F_Delta_14C = lambda C12, C14: (C14/C12/alpha - 1) * 1000

    fig, axes = plt.subplots(
        figsize = (14,7),
        nrows   = 2,
        ncols   = 3,
        sharex  = True,
        sharey  = True
    )
   
    for nr, pool_name in enumerate(ms.pool_names):
        ax = axes.flatten()[nr]
        Delta_14C = F_Delta_14C(soln[:,nr], soln_14C[:,nr])
        ax.plot(t_conv(times), Delta_14C)
        ax.set_title(ms.pool_names[nr])
    
    plt.show()
    plt.close(fig)


def add_variable(ds, data_vars, var_name_ds, new_var):   
    if var_name_ds in ds.variables.keys():
        # avoid side-effects on ds
        var_ds = ds[var_name_ds].copy(deep=True)

        attrs  = var_ds.attrs
        attrs['units'] = new_var.unit
        var = xr.DataArray(
            data   = new_var.data,
            coords = var_ds.coords,
            dims   = var_ds.dims,
            attrs  = attrs
        )
        data_vars[var_name_ds] = var


def create_Delta_14C_dataset(mdo, ds, mr, mr_14C):
    alpha = 1.18e-12
    t_conv = lambda t: 2001 + (15+t)/365.25
    F_Delta_14C = lambda C12, C14: (C14/C12/alpha - 1) * 1000
    ms = mdo.model_structure

    ## create 14C dataset
    data_vars = {}

    ## save stocks
    for pool_nr, pool_name  in enumerate(ms.pool_names):
        var_name_ds = pool_name
        soln = mr.solve()
        soln_14C = mr_14C.solve()

        new_var = Variable(
            data = F_Delta_14C(soln[:, pool_nr], soln_14C[:, pool_nr]),
            unit = mdo.stock_unit
        )
        add_variable(ds, data_vars, var_name_ds, new_var)


    ## save external input fluxes

    # insert np.nan at time t0
#    raise(Exception('check units and values, use acc variants'))
    Us_monthly = Variable(
        data = np.concatenate(
            [
                np.nan * np.ones((1,len(ms.pool_names))),
                mr.acc_gross_external_input_vector()
            ],
            axis = 0
        ),
        unit = 'g/(365.25/12 d)'
    )
    print(mr, Us_monthly)
    Us_monthly_14C = Variable(
        data = np.concatenate(
            [
                np.nan * np.ones((1,len(ms.pool_names))),
                mr_14C.acc_gross_external_input_vector()
            ],
            axis = 0
        ),
        unit = 'g/(365.25/12 d)'
    )

    # convert to daily flux rates
    us = Us_monthly.convert('g/d')
    us_14C = Us_monthly_14C.convert('g/d')
    for pool_nr, pool_name  in enumerate(ms.pool_names):        
        var_name_ds = 'NPP_to_' + pool_name
        new_var = Variable(
            data = F_Delta_14C(us.data[:, pool_nr], us_14C.data[:, pool_nr]),
            unit = us.unit
        )
        add_variable(ds, data_vars, var_name_ds, new_var)


    ## save internal fluxes

    # insert np.nan at time t0
    Fs_monthly = Variable(
        data = np.concatenate(
            [
                np.nan * np.ones(
                    (1, len(ms.pool_names), len(ms.pool_names))
                ),
                mr.acc_gross_internal_flux_matrix()
            ],
            axis = 0
        ),
        unit = 'g/(365.25/12 d)'
    )
    Fs_monthly_14C = Variable(
        data = np.concatenate(
            [
                np.nan * np.ones(
                    (1, len(ms.pool_names), len(ms.pool_names))
                ),
                mr_14C.acc_gross_internal_flux_matrix()
            ],
            axis = 0
        ),
        unit = 'g/(365.25/12 d)'
    )

    # convert to daily flux rates
    fs = Fs_monthly.convert('g/d')
    fs_14C = Fs_monthly_14C.convert('g/d')
    for pnr_from, pn_from in enumerate(ms.pool_names):        
        for pnr_to, pn_to in enumerate(ms.pool_names):        
            var_name_ds = pn_from +'_to_' + pn_to
            new_var = Variable(
                data = F_Delta_14C(
                    fs.data[:, pnr_to, pnr_from],
                    fs_14C.data[:, pnr_to, pnr_from]
                ),
                unit = fs.unit
            )
            add_variable(ds,data_vars, var_name_ds, new_var)


    ## save external output fluxes

    # insert np.nan at time t0
    Rs_monthly = Variable(
        data = np.concatenate(
            [
                np.nan * np.ones((1,len(ms.pool_names))),
                mr.acc_gross_external_output_vector()
            ],
            axis = 0
        ),
        unit = 'g/(365.25/12 d)'
    )
    Rs_monthly_14C = Variable(
        data = np.concatenate(
            [
                np.nan * np.ones((1,len(ms.pool_names))),
                mr_14C.acc_gross_external_output_vector()
            ],
            axis = 0
        ),
        unit = 'g/(365.25/12 d)'
    )

    # convert to daily flux rates
    rs = Rs_monthly.convert('g/d')
    rs_14C = Rs_monthly_14C.convert('g/d')
    for pool_nr, pool_name  in enumerate(ms.pool_names):        
        var_name_ds = pool_name + '_to_RH'
        new_var = Variable(
            data = F_Delta_14C(rs.data[:, pool_nr], rs_14C.data[:, pool_nr]),
            unit = rs.unit
        )
        add_variable(ds, data_vars, var_name_ds, new_var)


    ds_Delta_14C = xr.Dataset(
        data_vars = data_vars,
        coords    = ds.coords,
        attrs     = ds.attrs
    )

    return ds_Delta_14C


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
    Fa_func = lambda t: alpha * (F_atm_delta_14C(t)/1000+1)
    net_Us_12C = dmr.acc_net_external_input_vector()
    with np.errstate(divide='ignore'):
        net_Us_14C = net_Us_12C * Fa_func(dmr.times[:-1]).reshape(-1,1)
        #acc_net_us_14C = alpha * net_Us_12C * (
        #    1 + 1/1000 * F_atm_delta_14C(dmr.times[:-1]).reshape(-1,1)
        #)
    net_Us_14C = np.nan_to_num(net_Us_14C, posinf=0)
    
    # compute 14C start_values
    start_values_14C = compute_start_values_14C(
        dmr.times,
        dmr.Bs,
        dmr.net_Us,
        Fa_func,
        method = 'D3'
    )
    
    dmr_14C = DMR_14C(
        dmr,
        start_values_14C,
        net_Us_14C
    )   

    return dmr_14C

 
def load_pwc_mr_fd_14C(pwc_mr_fd):
    ## create 14C pwc_mr_fd

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
    Fa_func = lambda t: alpha * (F_atm_delta_14C(t)/1000+1)

    ## compute 14C start_values
    B_func = pwc_mr_fd.B_func()
    Bs_12C = np.array([B_func(t) for t in pwc_mr_fd.times[:-1]]) 
    u_func = pwc_mr_fd.external_input_vector_func()
    us_12C = np.array([u_func(t) for t in pwc_mr_fd.times[:-1]])
    start_values_14C = compute_start_values_14C(
        pwc_mr_fd.times,
        Bs_12C,
        us_12C,
        Fa_func,
        method = 'C3'
    )

    pwc_mr_fd_14C = PWCMR_14C(
        pwc_mr_fd.pwc_mr,
        start_values_14C,
        Fa_func
    )   

    return pwc_mr_fd_14C

 
def load_Delta_14C_dataset(ds, method='continuous'):
    if method not in ['discrete', 'continuous']:
        raise(
            ValueError(
                "method must be either 'discrete' or 'continuous'"
            )
        ) 
    if method == 'discrete':
        raise(
            ValueError(
                """
                discrete Delta 14C computation is impossible because
                we cannot compute the necessary net Delta 14C Us without
                knowledge of the state transition operator
                """
            )
        )

#    return ds 
    # fake return data structure on first call (with empty data)
    if ds.ens.values.shape == (0,):
        empty_var = xr.DataArray(
            data = np.ndarray(dtype=float, shape=(0,0,0)),
            dims = ('ens', 'lat' , 'lon')
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
    
 #       if method == 'discrete':
 #           mr, abs_err, rel_err =\
 #                mdo.create_discrete_model_run(errors=True)
 #           mr_14C = load_dmr_14C(mr)
        if method == 'continuous':
            mr, abs_err, rel_err =\
                 mdo.create_model_run(errors=True)
            var_abs_err = xr.DataArray(
                data  = np.max(abs_err.data),
                attrs = {
                    'units':     abs_err.unit,
                    'long_name': 'max. abs. error on reconstructed stock sizes'
                }
            )
    
            var_rel_err = xr.DataArray(
                data  = np.max(rel_err.data),
                attrs = {
                    'units':     rel_err.unit,
                    'long_name': 'max. rel. error on reconstructed stock sizes'
                }
            )
            print('abs', var_abs_err) 
            print('rel', var_rel_err) 
    
            mr_14C = load_pwc_mr_fd_14C(mr)

        ms = mdo.model_structure
        #plot_Delta_14C_in_time(ms, mr.solve(), mr_14C.solve())
        ds_Delta_14C = create_Delta_14C_dataset(mdo, ds, mr, mr_14C)
    
        ## add reconstruction error
        ds_Delta_14C['max_abs_err'] = var_abs_err
        ds_Delta_14C['max_rel_err'] = var_rel_err
    except DMRError as e:
        log = str(e)
        
        data_vars = {}
        for key, val in ds.data_vars.items():
            if key != 'time':
                data_vars[key] = np.nan * val
        ds_Delta_14C = ds.copy(data=data_vars)
        
        ds_Delta_14C['max_abs_err'] = np.nan
        ds_Delta_14C['max_rel_err'] = np.nan

    ds_Delta_14C['log'] = log
    ds_Delta_14C.close()    

    return ds_Delta_14C


################################################################################


if __name__ == '__main__':
    pass
    dataset = xr.open_dataset('~/Desktop/CARDAMOM/cardamom_for_holger.nc')
    ds = dataset.isel(ens=0, lat=0, lon=0)
#    ds_Delta_14C_dmr = load_Delta_14C_dataset(ds, 'discrete')
    ds_Delta_14C_pwc_mr_fd = load_Delta_14C_dataset(ds, 'continuous')

#    # check for similarity
#    for name, var_dmr in ds_Delta_14C_dmr.data_vars.items():
#        if name not in ['log', 'max_abs_err', 'max_rel_err']:
#            val_dmr = var_dmr.data
#            val_pwc_mr_fd = ds_Delta_14C_pwc_mr_fd[name].data
#            print(name)
#            print(val_dmr)
#            print(val_pwc_mr_fd)
#
#
#            abs_err = np.abs(val_dmr-val_pwc_mr_fd)
#            print(np.nanmax(abs_err))
#            rel_err = abs_err/np.abs(val_dmr)
#            print(np.nanmax(rel_err)*100)
#            rel_err = abs_err/np.abs(val_pwc_mr_fd)
#            print(np.nanmax(rel_err))

    ds_Delta_14C_dmr.close()
#    ds_Delta_14C_pwc_mr_fd.close()
    ds.close()
    dataset.close()


