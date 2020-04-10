import numpy as np
import xarray as xr

from bgc_md2.ModelStructure import ModelStructure
from bgc_md2.ModelDataObject import (
    ModelDataObject,
    getStockVariable_from_Density,
    getFluxVariable_from_DensityRate,
    getFluxVariable_from_Rate
)
from bgc_md2.Variable import Variable


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


def load_dmr(ds):
#    days_per_month = np.array(np.diff(ds.time), dtype='timedelta64[D]').astype(float)

    ms = load_model_structure()            

   
#    ## bring fluxes from gC/m2/day to gC/m2/month
 
    ## all months are supposed to comprise 365.25/12 days
    days_per_month = 365.25/12.0
#    fv_names = ms.get_flux_var_names()
#    for fv_name in fv_names:
#        var = ds[fv_name]
#        var.data *= days_per_month
#        var.attrs['units'] = 'gC/m2'


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

    dmr = mdo.create_discrete_model_run()
    soln_dmr = dmr.solve()
    smrfd = mdo.create_model_run()
    soln_smrfd = smrfd.solve()

    print(np.max(np.abs((soln_dmr-soln_smrfd)/soln_dmr*100)))



################################################################################


if __name__ == '__main__':
    dataset = xr.open_dataset('~/Desktop/CARDAMOM/cardamom_for_holger.nc')

    ds = dataset.isel(ens=0, lat=0, lon=0)
    load_dmr(ds)

    dataset.close()




