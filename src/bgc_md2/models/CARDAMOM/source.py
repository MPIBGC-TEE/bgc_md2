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
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
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
x = StateVariableTuple(tuple(Symbol(name) for name in ms.pool_names))
gpp = Symbol('gpp')


def make_internal_flux(sv_name_tup):
    source, target = map(Symbol,sv_name_tup)
    fluxrate_densities = [
        Symbol(var_name)
        for var_name in ms.horizontal_structure[sv_name_tup]
    ]
    return sum(fluxrate_densities) * source


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
    #u,  # imput tuple
    InFluxesBySymbol(
        {
            Symbol(name): gpp * Symbol(val[0]) 
            for name, val in ms.external_input_structure.items()
        }
    ),
    OutFluxesBySymbol(
        {
            Symbol(name): Symbol(name) * sum([Symbol(rate) for rate in val])
            for name, val in ms.external_output_structure.items()
        }
    ),
    InternalFluxesBySymbol(
        {
            (Symbol(name_tup[0]),Symbol(name_tup[1])): make_internal_flux(name_tup)
            for name_tup in ms.horizontal_structure.keys()
        }
    ),
    t,  # time symbol 
    x,  # state vector of the complete system
    #VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    #VegetationCarbonInputPartitioningTuple(b),
    #VegetationCarbonStateVariableTuple((C_il, C_is, C_ir)),
})
