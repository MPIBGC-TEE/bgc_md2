# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Some alternatives to create html an markdown programmatically
#

from IPython.display import display, Markdown, Latex, HTML, Math

display(Markdown("*some markdown* $\phi$"))


from bgc_md2.helper import list_models_md, list_models


h.list_models()
"-".join(list_models())

display(Markdown(h.list_models_md()))

display(Markdown('[testVectorFree]("../../../../tmp/test.ipynb")'))
display(HTML("val=4"))


# +
from sympy import symbols, var, latex

var("A_f A_b")


# -

Math(latex(A_f))

display(Math(latex(A_f)))
