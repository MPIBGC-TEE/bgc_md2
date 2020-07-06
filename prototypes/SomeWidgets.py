# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import ipywidgets as widgets
from sympy import Matrix, Symbol
from sympy.printing import latex,mathml

widgets.HTML(
    value="Hello <b>World</b>",
    placeholder='Some HTML',
    description='Some HTML',
)
# -


a_f=Symbol('a_f')
widgets.HTMLMath(
    value=r"Some math and <i>HTML</i>: \(x^2\) and $$\frac{x+1}{x-1}$$",
    placeholder='Some HTML',
    description='Some HTML',
)

# +

widgets.HTMLMath(
    value=latex(a_f),
    placeholder='Some HTML',
    description='Some HTML',
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
from ipywidgets import GridspecLayout,AppLayout, Button, Layout,HTMLMath

grid = GridspecLayout(4, 3)

for i in range(4):
    grid[i, 0] = Button(layout=Layout(width='auto', height='auto'))
    grid[i, 1] = HTMLMath(value=r"Some math and <i>HTML</i>: \(x^2\) and $$\frac{x+1}{x-1}$$",
    placeholder='Some HTML',
    description='Some HTML',
)

grid
# -
A=Matrix([[1,2],[2,Symbol("x")]])


out = widgets.Output()
out
from IPython.display import display
for i in range(10):
        display(A)
#with out:
#    display(YouTubeVideo('eWzY2nGfkXk'))

import ipywidgets as widgets
items = [widgets.Label(str(i)) for i in range(4)]
left_box = widgets.VBox([items[0], items[1]])
right_box = widgets.VBox([items[2], items[3]])
widgets.HBox([left_box, right_box])

from IPython.display import Math
Math(latex(a_f))




