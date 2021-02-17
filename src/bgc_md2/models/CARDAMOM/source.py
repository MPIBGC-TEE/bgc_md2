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
    'fire_em_labile':           'external flux of labile carbon by fire',
    'fire_em_foliar':           'external flux for leaf pool by fire',
    'fire_em_root':             'external flux for root pool by fire',
    'fire_em_wood':             'external flux for wood pool by fire',
    'fire_em_litter':           'external flux for litter pool by fire',
    'fire_em_som':              'external flux for soil pool by fire',
    'hetresp_litter':           'heterotrophic respiration rate of litter pool',
    'hetresp_som':              'heterotrophic respiration rate of som pool',
    'labile_to_foliar':         'internal flux ',
    'fire_labile_to_litter':    'internal flux ',
    'fire_foliar_to_litter':    'internal flux ',
    'fire_wood_to_som':         'internal flux ',
    'wood_to_soilc':            'internal flux ',
    'fire_root_to_litter':      'internal flux ',
    'root_to_litter':           'internal flux ',
    'fire_litter_to_som':       'internal flux ',
    'litter_to_som':            'internal flux ',
}

# created from the mo
# to be able to use symbols directly 
for name in sym_dict.keys():
    var(name)


t = TimeSymbol('t')
ms = load_model_structure()
x = StateVariableTuple(tuple(Symbol(name) for name in ms.pool_names))


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
            Symbol(name): sum([Symbol(flux) for flux in val])
            for name, val in ms.external_input_structure.items()
        }
    ),
    # direct description would be
    # InFluxesBySymbol(
    #     {
    #         𝙻𝚊𝚋𝚒𝚕𝚎: 𝚐𝚙𝚙⎯𝚝𝚘⎯𝚕𝚊𝚋𝚒𝚕𝚎,
    #         𝙻𝚎𝚊𝚏: 𝚐𝚙𝚙⎯𝚝𝚘⎯𝚕𝚎𝚊𝚏,
    #         𝚁𝚘𝚘𝚝: 𝚐𝚙𝚙⎯𝚝𝚘_𝚛𝚘𝚘𝚝,
    #         𝚆𝚘𝚘𝚍: 𝚐𝚙𝚙⎯𝚝𝚘⎯𝚠𝚘𝚘𝚍
    #     }
    # ),
    OutFluxesBySymbol(
        {
            Symbol(name): sum([Symbol(flux) for flux in val])
            for name, val in ms.external_output_structure.items()
        }
    ),
    # direct description would be
    # OutFluxesBySymbol(
    #     {
    #         𝙻𝚊𝚋𝚒𝚕𝚎: 𝚏𝚒𝚛𝚎⎯𝚎𝚖⎯𝚕𝚊𝚋𝚒𝚕𝚎,
    #         𝙻𝚎𝚊𝚏: 𝚏𝚒𝚛𝚎⎯𝚎𝚖⎯𝚏𝚘𝚕𝚒𝚊𝚛,
    #         𝚁𝚘𝚘𝚝: 𝚏𝚒𝚛𝚎⎯𝚎𝚖⎯𝚛𝚘𝚘𝚝,
    #         𝚆𝚘𝚘𝚍: 𝚏𝚒𝚛𝚎⎯𝚎𝚖⎯𝚠𝚘𝚘𝚍,
    #         𝙻𝚒𝚝𝚝𝚎𝚛: 𝚏𝚒𝚛𝚎⎯𝚎𝚖⎯𝚕𝚒𝚝𝚝𝚎𝚛 + 𝚑𝚎𝚝𝚛𝚎𝚜𝚙⎯𝚕𝚒𝚝𝚝𝚎𝚛,
    #         𝚂𝚘𝚒𝚕: 𝚏𝚒𝚛𝚎⎯𝚎𝚖⎯𝚜𝚘𝚖 + 𝚑𝚎𝚝𝚛𝚎𝚜𝚙⎯𝚜𝚘𝚖
    #     }
    # ),
    InternalFluxesBySymbol(
        {
            (Symbol(name_tup[0]),Symbol(name_tup[1])): sum([ Symbol(flux) for flux in val])
            for name_tup, val in ms.horizontal_structure.items()
        }
    ),
    # direct description would be
    # InternalFluxesBySymbol(
    #     {
    #         (𝙻𝚊𝚋𝚒𝚕𝚎, 𝙻𝚎𝚊𝚏): 𝚕𝚊𝚋𝚒𝚕𝚎⎯𝚝𝚘⎯𝚏𝚘𝚕𝚒𝚊𝚛,
    #         (𝙻𝚊𝚋𝚒𝚕𝚎, 𝙻𝚒𝚝𝚝𝚎𝚛): 𝚏𝚒𝚛𝚎⎯𝚕𝚊𝚋𝚒𝚕𝚎⎯𝚝𝚘⎯𝚕𝚒𝚝𝚝𝚎𝚛,
    #         (𝙻𝚎𝚊𝚏, 𝙻𝚒𝚝𝚝𝚎𝚛): 𝚏𝚒𝚛𝚎⎯𝚏𝚘𝚕𝚒𝚊𝚛⎯𝚝𝚘⎯𝚕𝚒𝚝𝚝𝚎𝚛 + 𝚕𝚎𝚊𝚏⎯𝚝𝚘⎯𝚕𝚒𝚝𝚝𝚎𝚛,
    #         (𝚆𝚘𝚘𝚍, 𝚂𝚘𝚒𝚕): 𝚏𝚒𝚛𝚎⎯𝚠𝚘𝚘𝚍⎯𝚝𝚘⎯𝚜𝚘𝚖 + 𝚠𝚘𝚘𝚍⎯𝚝𝚘⎯𝚜𝚘𝚒𝚕𝚌,
    #         (𝚁𝚘𝚘𝚝, 𝙻𝚒𝚝𝚝𝚎𝚛): 𝚏𝚒𝚛𝚎⎯𝚛𝚘𝚘𝚝⎯𝚝𝚘⎯𝚕𝚒𝚝𝚝𝚎𝚛 + 𝚛𝚘𝚘𝚝⎯𝚝𝚘⎯𝚕𝚒𝚝𝚝𝚎𝚛,
    #         (𝙻𝚒𝚝𝚝𝚎𝚛, 𝚂𝚘𝚒𝚕): 𝚏𝚒𝚛𝚎⎯𝚕𝚒𝚝𝚝𝚎𝚛⎯𝚝𝚘⎯𝚜𝚘𝚖 + 𝚕𝚒𝚝𝚝𝚎𝚛⎯𝚝𝚘⎯𝚜𝚘𝚖
    #     }
    # ),
    t,  # time symbol
    x,  # state vector of the complete system
    # VegetationCarbonInputScalar(gpp), # ? not sure see ticket
    # vegetation carbon partitioning.
    # VegetationCarbonInputPartitioningTuple(b),
    VegetationCarbonStateVariableTuple((Labile, Leaf, Root, Wood)),
})
