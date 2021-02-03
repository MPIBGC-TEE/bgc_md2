import xarray as xr
from sympy import symbols, Symbol, var

from .CARDAMOMlib import _load_mdo, load_model_structure

from ..BibInfo import BibInfo 

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

sym_dict = {
    'Labile':                   'Labile Carbon',
    'Leaf':                     'Leaf carbon pool',
    'Root':                     'Root carbon pool',
    'Wood':                     'Wood carbon pool',
    'Litter':                   'Litter carbon pool',
    'Soil':                     'SOM carbon pool',
    'gpp':                      'Gross primary production',
    'gpp_to_labile':            'Fraction of GPP allocated to labile carbon pool',
    'gpp_to_leaf':              'Fraction of GPP allocated to leaf carbon pool',
    'gpp_to_root':              'Fraction of GPP allocated to root carbon pool',
    'gpp_to_wood':              'Fraction of GPP allocated to wood carbon pool',
    'fire_em_labile':           'external flux rate of labile carbon by fire',
    'fire_em_foliar':           'external flux rate for leaf pool by fire',
    'fire_em_root':             'external flux rate for root pool by fire',
    'fire_em_wood':             'external flux rate for wood pool by fire',
    'fire_em_litter':           'external flux rate for litter pool by fire',
    'fire_em_som':              'external flux rate for soil pool by fire',
    'hetresp_litter':           'heterotrophic respiration rate of litter pool',
    'hetresp_som':              'heterotrophic respiration rate of som pool',
    'labile_to_foliar':         'internal flux rate',
    'fire_labile_to_litter':    'internal flux rate',
    'fire_foliar_to_litter':    'internal flux rate',
    'fire_wood_to_som':         'internal flux rate',
    'wood_to_soilc':            'internal flux rate',
    'fire_root_to_litter':      'internal flux rate',
    'root_to_litter':           'internal flux rate',
    'fire_litter_to_som':       'internal flux rate',
    'litter_to_som':            'internal flux rate',
}

for name in sym_dict.keys():
    var(name)


t = TimeSymbol('t')
ms = load_model_structure()
x = StateVariableTuple(tuple(Symbol(name) for name in ms.pool_names))
gpp = Symbol('gpp')


def make_internal_flux(sv_name_tup):
    source, target = map(Symbol, sv_name_tup)
    fluxrate_densities = [
        Symbol(var_name)
        for var_name in ms.horizontal_structure[sv_name_tup]
    ]
    return sum(fluxrate_densities) * source


mvs = MVarSet({
    BibInfo(# Bibliographical Information
       name="CARDAMOM",
       longName="?",
       version="?",
       entryAuthor="Holger Metzler",
       entryAuthorOrcid="",
       entryCreationDate="02/02/2021",
       doi="",
       #further_references=BibInfo(doi="10.5194/bg-10-2255-2013"),
       sym_dict=sym_dict
    ),
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
    # VegetationCarbonInputScalar(u),
    # vegetation carbon partitioning.
    # VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((Labile, Leaf, Root, Wood)),
})
