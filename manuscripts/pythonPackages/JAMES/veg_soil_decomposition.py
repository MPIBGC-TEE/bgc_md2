import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
import bgc_md2
from importlib import import_module
import CompartmentalSystems.helpers_reservoir as hr

mf="kv_visit2"
#mf="jon_yib"
model_mod=f'bgc_md2.models.{mf}'
mvs = import_module(f"{model_mod}.source").mvs
part_dict =  {
    frozenset(mvs.get_VegetationCarbonStateVariableTuple()):'green',
    frozenset(mvs.get_SoilCarbonStateVariableTuple()):'brown',
}
gplot=hr.igraph_part_plot(
    mvs.get_StateVariableTuple(),
    mvs.get_InFluxesBySymbol(),#in_fluxes,
    mvs.get_InternalFluxesBySymbol(),#internal_fluxes,
    mvs.get_OutFluxesBySymbol(),#out_fluxes,
    part_dict,
    Path(f"{mf}_veg_soil_decomposition.pdf")
)



