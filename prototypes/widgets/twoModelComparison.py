# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import ipywidgets as widgets
from ipywidgets import Output, HTML, Button, HBox, VBox

from sympy import latex
import matplotlib.pyplot as plt

from ComputabilityGraphs import CMTVS

from bgc_md2.resolve.mvars import CompartmentalMatrix, InputTuple
import bgc_md2.helper as h
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from ipywidgets import Output,Layout, Button, Box, VBox, Label, HTMLMath
def table(
        records,
        types
    ):
    cws=['20%','30%']
    def line(record):
        srm = record._get_single_value(SmoothReservoirModel)
        with plt.ioff():
            fig = plt.figure()
            rect = 0, 0, 0.8, 1.2  # l, b, w, h
            ax = fig.add_axes(rect)
        
        graph_out = Output(
            layout=Layout(
                height='auto',
                min_width=cws[0],
                description="Compartmental Graph",
            )
        )
        with graph_out:
            ax.clear()
            srm.plot_pools_and_fluxes(ax)
            display(ax.figure)
    
        line = Box(
                children=[ graph_out]+[
                    HTMLMath(
                        value=latex(record._get_single_value(t)),
                        #value=r"Some math and <i>HTML</i>: \(x^2\) and $$\frac{x+1}{x-1}$$",
                        layout=Layout(width=cws[1])
                    ) 
                    for t in types
                ],
                layout= Layout(
                    overflow='scroll hidden',
                    border='3px solid black',
                    width='1500px',
                    height='',
                    flex_flow='row',
                    display='flex'
                )
        )
        return line
    
    return VBox([line(r) for r in records])


# -

display(
    table(
        records=[
            h.CMTVS_from_model_name(name) 
            for name in ["Williams2005GCB"]
        ],
        types = [InputTuple]
    )
)




