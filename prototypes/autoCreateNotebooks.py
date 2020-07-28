# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from IPython.display import HTML
_ = HTML("<style>.container { width:100% !important; }</style>")

# %load_ext autoreload
# %autoreload 2
import ipywidgets as widgets
from IPython.display import display
import bgc_md2.helper as h


# # use a layout widget to build the interactive overview table
#

outerGridBox = widgets.GridspecLayout(3, 1)
model_view_pos = 1

def inspect_model(name):
    outerGridBox[model_view_pos, 0:2] = h.modelVBox(name)

outerGridBox[0, 0:2] = h.ModelListGridBox(inspect_model=inspect_model)

inspect_model(h.list_models()[0])

# make sure that the next cell is not in scroll mode
display(outerGridBox)
# -

# # Some alternatives to create html an markdown programmatically
#

from IPython.display import display, Markdown, Latex,HTML,Math

display(Markdown('*some markdown* $\phi$'))


from bgc_md2.helper import list_models_md,list_models


h.list_models()
"-".join(list_models())

display(Markdown(h.list_models_md()))

display(Markdown('[testVectorFree]("../../../../tmp/test.ipynb")'))
display(HTML("val=4"))


# +
from sympy import symbols,var,latex

var("A_f A_b")


# -

Math(latex(A_f))

display(Math(latex(A_f)))


