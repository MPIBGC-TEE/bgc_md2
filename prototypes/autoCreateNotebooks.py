# ---
# jupyter:
#   jupytext:
#     formats: py:light
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

# %load_ext autoreload
# %autoreload 2
import ipywidgets as widgets
from IPython.display import display
import bgc_md2.helper as h


# # use a layout widget to build the interactive overview table
#

outerGridBox=widgets.GridspecLayout(4,1)

outerGridBox[0,0:2]=h.modelListGridBox()

dd = widgets.Dropdown(
            options=h.list_models(),
            #value='2',
            #description='Inspect (Thsi could also be achieved by a button close to the entry in the list above..):',
            #style={'description_widt':'70%'},
            disabled=False,
)

# +
outerGridBox[1,0:2]=widgets.HBox(
    [
        widgets.Label(
            value=
            'Inspect (The selection of the model to inspect should also be achievable by a button \n close to the entry in the list above..):'
        ),
        dd
    ]
)
modelViewPos=2
#outerGridBox[modelViewPos,0:2]=h.modelGridBox(dd.value)
outerGridBox[modelViewPos,0:2]=h.modelVBox(dd.value)

def updateModelView(x):
    #outerGridBox[modelVielPos,0:2]=h.modelGridBox(dd.value)
    outerGridBox[modelViewPos,0:2]=h.modelVBox(dd.value)
    
dd.observe(updateModelView)
# make sure that the next cell is not in scroll mode
display(outerGridBox)
# -

# # Some alternatives to create html an markdown programmatically
#

from IPython.display import display, Markdown, Latex,HTML,Math

display(Markdown('*some markdown* $\phi$'))


from bgc_md2.helper import list_models_md,list_models,modelListGridBox


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


