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
display(HTML("<style>.container { width:100% !important; }</style>"))


# -

def table(
        records,
        types
    ):
    n=len(types)+1
    cws=['{}%'.format(100/n) for i in range(n)]
    print(cws)
    line_layout= Layout(
        overflow='scroll hidden',
        border='3px solid black',
        #width='1500px',
        width='100%',
        height='',
        flex_flow='row',
        display='flex'
    )
    def firstline():
        names=['graph'] + [ t.__name__ for t in types]
        return Box(
            children=[
                HTML(
                    value="<h3>{}</h3>".format(name),
                    layout=Layout(width=cws[i])
                ) 
                for i,name in enumerate(names)
            ],
            layout = line_layout
        )
    def line(record):
        srm = record._get_single_value(SmoothReservoirModel)
        with plt.ioff():
            fig = plt.figure()
            rect = 0, 0, 0.8, 1.2  # l, b, w, h
            ax = fig.add_axes(rect)
        
        graph_out = Output(
            layout=Layout(
                height='auto',
                width=cws[0],
                description="Compartmental Graph",
            )
        )
        with graph_out:
            ax.clear()
            srm.plot_pools_and_fluxes(ax)
            display(ax.figure)
            
        outs = [ 
            Output(layout=Layout(width=cws[i]))
            for i,t in enumerate(types)
        ]
        for i,t in enumerate(types):
            with outs[i]:
                display(record._get_single_value(t))
                    

        
        line = Box(
                children=[ graph_out]+outs,
                layout = line_layout
        )
        return line
    
    return VBox(
        [firstline()]+
        [line(r) for r in records]
    )

display(
    table(
        records=[
            h.CMTVS_from_model_name(name) 
            for name in ["Williams2005GCB"]
        ],
        types = [InputTuple,CompartmentalMatrix]
    )
)




