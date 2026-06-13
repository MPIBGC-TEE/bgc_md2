import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import bgc_md2
from importlib import import_module
from pathlib import Path

from ComputabilityGraphs import helpers as h
mf="jon_yib"
model_mod=f'bgc_md2.models.{mf}'
##
mvs = import_module(f"{model_mod}.source").mvs
p=Path("/home/mm/bgc_md2/manuscripts/pythonPackages/JAMES/")
tl, fig, fl = mvs.dep_graph_figure(
    root_type=bgc_md2.resolve.mvars.NumericVegetationCarbonMeanBackwardTransitTimeSolution,
    targetPaths=[
        p.joinpath(*names) 
        for names in [
            ["TypeLegend_2.tex"],
            ["figures","dep_graph_2.pdf"],
            ["ComputerLegend_2.tex"]
        ]
    ]
)

