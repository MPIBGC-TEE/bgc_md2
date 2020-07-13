import xarray as xr
from CARDAMOMlib import load_mdo

dataset = xr.open_dataset('~/Desktop/CARDAMOM/cardamom_for_holger.nc')

def create_pwc_model_run_fd(ens, lat, lon):
    ds = dataset.isel(ens=ens, lat=lat, lon=lon)
    mdo = load_mdo(ds)
    pwc_mr_fd = mdo.create_model_run()
    return pwc_mr_fd

# example
# create_pwc_model_run_fd(ens=0, lat=0, lon=0)
