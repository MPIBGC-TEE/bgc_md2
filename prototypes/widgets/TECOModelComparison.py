# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.3
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

from bgc_md2.resolve.mvars import CompartmentalMatrix, InputTuple, StateVariableTuple
import bgc_md2.helper as h
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from ipywidgets import Output,Layout, Button, Box, VBox, Label, HTMLMath
display(HTML("<style>.container { width:100% !important; }</style>"))

# +
#h.list_models()

# +
from typing import Callable
from ipywidgets import Layout, Button, Box, VBox, Label, HTMLMath, HTML,HBox
class ExpandBox(VBox):
    
    # instead of content we use contant generating functions
    # Thus the expanded content is only created when the button
    # is pressed and both collapsed and expanded content can 
    # be recreated as often as necessary
    def __init__(
        self,
        button_width_str: str,
        collapsed_content_func:Callable,
        expanded_content_func:Callable
    ):
        box_layout = Layout(
            overflow='scroll hidden',
            border='3px solid black',
            width='100%',
            height='',
            flex_flow='row',
            display='flex'
        )
        super().__init__(layout=box_layout)
        self.button_width_str = button_width_str
        self.collapsed_content_func = collapsed_content_func
        self.expanded_content_func = expanded_content_func
        self.show_collapsed_content()
        
    def triggered_show_expanded(self,button):
        c=Button(
            layout= Layout(height='auto', min_width=self.button_width_str),
            description='collapse',
            button_style='warning'
        )
        c.on_click(self.triggered_show_collapsed)
        self.children = [c]+self.expanded_content_func()
        
    def triggered_show_collapsed(self,button):
        # a wrapper that has the necessary button argument
        self.show_collapsed_content()
        
    def show_collapsed_content(self):
        b=Button(
            layout= Layout(
                height='auto',
                min_width=self.button_width_str,
                max_width=self.button_width_str
            ),
            description='expand',
            button_style='warning'
        )
        b.on_click(self.triggered_show_expanded)
        self.children = [b] + self.collapsed_content_func()

# Now do it programmatically:
def line(
    tup,
    types_compact,
    types_expanded
):
    n=len(types_compact)
    name, record = tup
    button_width=10
    thumbnail_width=10
    remaining_width=100 - (button_width + thumbnail_width)
    def collapsed_c():
        #################### make a thumbnail graph
        srm = record._get_single_value(SmoothReservoirModel)
        with plt.ioff():
            fig = plt.figure(figsize=(0.7,0.7))
            rect = 0, 0, 0.8, 1.2  # l, b, w, h
            ax = fig.add_axes(rect)
        
        graph_out = Output(
            layout=Layout(
                height='auto',
                width="{}%".format(thumbnail_width),
                description="Compartmental Graph",
            )
        )
        with graph_out:
            ax.clear()
            srm.plot_pools_and_fluxes(
                ax,
                thumbnail=True,
                legend=False,
                mutation_scale=15,
                fontsize=16
            )
            display(ax.figure)
            plt.close(fig)
        #################### print the model name
        name_out= HTML("""{}""".format(name))
            
        #################### print the additionally required fields
        # give the additional colums equal width
        cws=['{}%'.format(remaining_width/n) for i in range(n)]
        outs = [ 
            Output(layout=Layout(width=cws[i]))
            for i,t in enumerate(types_compact)
        ]
        for i,t in enumerate(types_compact):
            with outs[i]:
                display(record._get_single_value(t))
                    
        return [ name_out, graph_out ] + outs
    
    
    def expanded_c():
        # here we expand vertically and use the full remaining width
        width_str="{}".format(remaining_width)
        srm = record._get_single_value(SmoothReservoirModel)
        with plt.ioff():
            fig = plt.figure()
            rect = 0, 0, 0.8, 1.2  # l, b, w, h
            ax = fig.add_axes(rect)
        
        graph_out = Output(
            layout=Layout(
                height='auto',
                width=width_str,
                description="Compartmental Graph",
            )
        )
        with graph_out:
            ax.clear()
            srm.plot_pools_and_fluxes(ax)
            display(ax.figure)
            plt.close(fig)
            
        outs = [ 
            Output(layout=Layout(width=width_str))
            for i,t in enumerate(types_expanded)
        ]
        for i,t in enumerate(types_expanded):
            with outs[i]:
                display(record._get_single_value(t))
                    
        return [
            VBox(
                [ graph_out ] + outs
            )
        ]
    
    return  ExpandBox(
        button_width_str="{}%".format(button_width),
        collapsed_content_func=collapsed_c,
        expanded_content_func=expanded_c
    )

def table(tups,types_compact,types_expanded):
    return VBox(
        children=tuple(
            map(
                lambda tups:line(tups,types_compact,types_expanded),
                tups
            )
        )
    )
    
table(
    tups=[
        (name,h.CMTVS_from_model_name(name)) 
        #for name in ["cable_general","TECOmm","TECO","Luo2012TE","TECO_general","cable_general"]
        for name in ["cable_yuanyuan","VISIT_Kostia"]
    ],
    types_compact = [],
    types_expanded = [InputTuple,CompartmentalMatrix]
    #types = [InputTuple,CompartmentalMatrix,StateVariableTuple]
    #types = [InputTuple,InternalFluxesBySymbol]
)

# -


