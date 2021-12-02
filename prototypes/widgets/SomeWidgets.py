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
# %load_ext autoreload
# %autoreload 2
import ipywidgets as widgets
from sympy import Matrix, Symbol
from sympy.printing import latex, mathml

widgets.HTML(
    value="Hello <b>World</b>", placeholder="Some HTML", description="Some HTML",
)
# -


a_f = Symbol("a_f")
widgets.HTMLMath(
    value=r"Some math and <i>HTML</i>: \(x^2\) and $$\frac{x+1}{x-1}$$",
    placeholder="Some HTML",
    description="Some HTML",
)

# +

widgets.HTMLMath(
    value=latex(a_f), placeholder="Some HTML", description="Some HTML",
)

# +

button = widgets.Button(description="Click Me!")
output = widgets.Output()

display(button, output)


def on_button_clicked(b):
    with output:
        print("Button clicked.")


button.on_click(on_button_clicked)
# -
items = [widgets.Label(str(i)) for i in range(8)]
widgets.GridBox(items, layout=widgets.Layout(grid_template_columns="repeat(3, 100px)"))

# +
from ipywidgets import GridspecLayout, AppLayout, Button, Layout, HTMLMath

grid = GridspecLayout(4, 3)

for i in range(4):
    grid[i, 0] = Button(layout=Layout(width="auto", height="auto"))
    grid[i, 1] = HTMLMath(
        value=r"Some math and <i>HTML</i>: \(x^2\) and $$\frac{x+1}{x-1}$$",
        placeholder="Some HTML",
        description="Some HTML"
    )

grid
# -
A = Matrix([[1, 2], [2, Symbol("sigma")]])


out = widgets.Output()
out
from IPython.display import display

for i in range(2):
    display(A)
# with out:
#    display(YouTubeVideo('eWzY2nGfkXk'))

import ipywidgets as widgets

items = [widgets.Label(str(i)) for i in range(4)]
left_box = widgets.VBox([items[0], items[1]])
right_box = widgets.VBox([items[2], items[3]])
widgets.HBox([left_box, right_box])

from IPython.display import Math

Math(latex(a_f))

# +
import matplotlib.pyplot as plt
from sympy import Symbol, Matrix, symbols, diag, zeros, simplify, Function
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

C_0, C_1 = symbols("C_0 C_1")
state_vector = [C_0, C_1]
time_symbol = Symbol("t")
input_fluxes = {0: 3 * C_0}
output_fluxes = {1: 2 * C_0}
internal_fluxes = {(0, 1): 5 * C_0 * C_1, (1, 0): 4 * C_0}
srm = SmoothReservoirModel(
    state_vector, time_symbol, input_fluxes, output_fluxes, internal_fluxes
)

graph_out = widgets.Output()
with graph_out:
    fig = plt.figure()
    rect = 0, 0, 0.8, 1.2  # l, b, w, h

    ax = fig.add_axes(rect)
    ax.clear()
    srm.plot_pools_and_fluxes(ax)


box = widgets.VBox(
    [
        widgets.HTML(
            value="""
                <h1>{name}</h1>
                Overview 
                """.format(
                name="Test"
            )
        ),
        graph_out,
    ]
)

box


# +
from ComputabilityGraphs import CMTVS
from bgc_md2.resolve.mvars import CompartmentalMatrix
import bgc_md2.helper as h
from sympy import latex

import nbformat as nbf
import pkgutil
import matplotlib.pyplot as plt
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

record = h.CMTVS_from_model_name("Williams2005GCB")
target_var = SmoothReservoirModel
srm = record._get_single_value(target_var)

graph_out = widgets.Output()
with graph_out:
    fig = plt.figure()
    rect = 0, 0, 0.8, 1.2  # l, b, w, h

    ax = fig.add_axes(rect)
    ax.clear()
    srm.plot_pools_and_fluxes(ax)


box = widgets.VBox(
    [
        widgets.HTML(
            value="""
                <h1>{name}</h1>
                Overview 
                """.format(
                name="Test"
            )
        ),
        graph_out,
    ]
)

box
# -
from ipywidgets import Layout, Button, Box, VBox, Label, HTMLMath
cw1='60%'
cw2='30%'
items1 = [
    Button(
        layout= Layout(height='auto', min_width=cw1),
        description=str(1),
        button_style='warning'
    ),
    Button(
        layout= Layout(height='auto', min_width=cw2),
        description=str(2),
        button_style='warning'
    ),
]
items2 = [
    HTMLMath(
        value=r"Some math and <i>HTML</i>: \(x^2\) and $$\frac{x+1}{x-1}$$",
        layout= Layout(height='auto', min_width=cw1),
        description=str(1),
    ),
    Button(
        layout= Layout(height='auto', min_width=cw2),
        description=str(2),
        button_style='warning'
    ),
]
box_layout = Layout(overflow='scroll hidden',
                    border='3px solid black',
                    width='500px',
                    height='',
                    flex_flow='row',
                    display='flex')
line1= Box(children=items1, layout=box_layout)
line2= Box(children=items2, layout=box_layout)
VBox([Label('Scroll horizontally:'), line1,line2])





