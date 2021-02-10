import xarray as xr
from sympy import symbols, Symbol, var

from .CARDAMOMlib import _load_mdo, load_model_structure

from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonStateVariableTuple,
)
from bgc_md2.helper import MVarSet
# dataset = xr.open_dataset('~/Desktop/CARDAMOM/cardamom_for_holger.nc')
#dataset = xr.open_dataset("/home/data/CARDAMOM/cardamom_for_holger.nc")


def create_pwc_model_run_fd(ens, lat, lon):
    #dataset = xr.open_dataset('~/Desktop/CARDAMOM/cardamom_for_holger.nc')
    dataset = xr.open_dataset('/home/data/CARDAMOM/cardamom_for_holger.nc')

    ds = dataset.isel(ens=ens, lat=lat, lon=lon)
    mdo = _load_mdo(ds)
    pwc_mr_fd = mdo.create_model_run()

    ds.close()
    dataset.close()
    return pwc_mr_fd


# examples
# create_pwc_model_run_fd(ens=0, lat=0, lon=0)
# create_pwc_model_run_fd(ens=(0, 10, 1))

t = TimeSymbol('t')
ms = load_model_structure()
x = StateVariableTuple(ms.pool_names)
def make_input_flux(sv_name):
    sv = Symbol(sv_name)
    fluxrate_densities = [ Symbol(name) for name in ms.external_input_structure]
    if sv in fluxrate_densities:
        fl = sum(fluxrate_densities)*sv
    else:
        fl = 0
    return fl
u = InputTuple((make_input_flux(key) for key in ms.pool_names))
mvs = MVarSet({
    #BibInfo(# Bibliographical Information
    #    name="CARDAMOM",
    #    longName="to be ", 
    #    version="2.6",
    #    entryAuthor="",
    #    entryAuthorOrcid="",
    #    entryCreationDate="22/3/2016",
    #    doi="",
    #    #further_references=BibInfo(doi="10.5194/bg-10-2255-2013"),
    #    sym_dict=sym_dict
    #),
    #A,  # the overall compartmental matrix
    u,  # imput tuple
    t,  # time symbol 
    x,  # state vector of the complete system
    #VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    #VegetationCarbonInputPartitioningTuple(b),
    #VegetationCarbonStateVariableTuple((C_il, C_is, C_ir)),
})
