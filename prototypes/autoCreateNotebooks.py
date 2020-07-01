# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# %load_ext autoreload
# %autoreload 2
import ipywidgets as widgets
from IPython.display import display
import bgc_md2.helper as h


# # use a layout widget to build the interactive overview table
#

# +
outerGridBox=widgets.GridspecLayout(4,2)
outerGridBox[0,0:2]=h.modelListGridBox()
dd = widgets.Dropdown(
            options=list_models(),
            #value='2',
            #description='Inspect (Thsi could also be achieved by a button close to the entry in the list above..):',
            #style={'description_widt':'70%'},
            disabled=False,
)
outerGridBox[1,0:2]=widgets.HBox(
    [
        widgets.Label(
            value='Inspect (Thsi could also be achieved by a button \n close to the entry in the list above..):'
        ),
        dd
    ]
)
outerGridBox[2,0:2]=h.modelGridBox(dd.value)
def updateModelView(x):
    outerGridBox[2,0:2]=h.modelGridBox(dd.value)
    
dd.observe(updateModelView)
display(outerGridBox)

# -

# # Some alternatives to create html an markdown programmatically
#

from IPython.display import display, Markdown, Latex,HTML

display(Markdown('*some markdown* $\phi$'))


from bgc_md2.helper import list_models_md,list_models,modelTableHtmlWidget,modelTableHtml,modelListGridBox


h.list_models()
"-".join(list_models())

display(Markdown(h.list_models_md()))

display(HTML(h.modelTableHtml()))

display(Markdown('[testVectorFree]("../../../../tmp/test.ipynb")'))

display(h.modelTableHtmlWidget(),output)

# +
outerGridBox=widgets.GridspecLayout(4,2)
outerGridBox[0,0:2]=h.modelListGridBox()
dd = widgets.Dropdown(
            options=list_models(),
            #value='2',
            #description='Inspect (Thsi could also be achieved by a button close to the entry in the list above..):',
            #style={'description_widt':'70%'},
            disabled=False,
)
outerGridBox[1,0:2]=widgets.HBox(
    [
        widgets.Label(
            value='Inspect (Thsi could also be achieved by a button \n close to the entry in the list above..):'
        ),
        dd
    ]
)
outerGridBox[2,0:2]=h.modelGridBox(dd.value)
def updateModelView(x):
    outerGridBox[2,0:2]=h.modelGridBox(dd.value)
    
dd.observe(updateModelView)
display(outerGridBox)

# -

dd.observe


