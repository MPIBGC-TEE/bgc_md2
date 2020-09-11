import xarray as xr
from .CARDAMOMlib import load_mdo

# dataset = xr.open_dataset('~/Desktop/CARDAMOM/cardamom_for_holger.nc')
#dataset = xr.open_dataset("/home/data/CARDAMOM/cardamom_for_holger.nc")


def create_pwc_model_run_fd(ens, lat, lon):
    #dataset = xr.open_dataset('~/Desktop/CARDAMOM/cardamom_for_holger.nc')
    dataset = xr.open_dataset('/home/data/CARDAMOM/cardamom_for_holger.nc')

    ds = dataset.isel(ens=ens, lat=lat, lon=lon)
    mdo = load_mdo(ds)
    pwc_mr_fd = mdo.create_model_run()

    ds.close()
    dataset.close()
    return pwc_mr_fd


# examples
# create_pwc_model_run_fd(ens=0, lat=0, lon=0)
# create_pwc_model_run_fd(ens=(0, 10, 1))
