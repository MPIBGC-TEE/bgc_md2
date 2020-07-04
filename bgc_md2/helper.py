import pkgutil
from . import models
import ipywidgets as widgets
from IPython.display import  display
from pathlib import Path
from sympy import Matrix,Symbol
import nbformat as nbf

from bgc_md2.models.helpers import (
    provided_mvars
    ,computable_mvars
    ,path_dict_to_single_mvar
    ,get_single_mvar_value
)

from bgc_md2.resolve.mvars import (
    TimeSymbol
    ,StateVariableTuple
    ,CompartmentalMatrix
)
def list_models():
    sub_mod_pkgs= [tup[1] for tup in pkgutil.iter_modules(models.__path__) if tup[2]]
    return sub_mod_pkgs

def list_models_md():
    names = list_models()
    links = ["[{name}](/tmp/{name})".format(name=name) for name in names]
    return "".join(links)


def funcmakerInsertLink(grid,names,i):
    def insert_link(b):
        # called for side effect on g
        tmpDirPath= Path("./tmp") #has to be relative for jupyter to open the file
        tmpDirPath.mkdir(exist_ok=True)   
        suffix=".ipynb"
        nbPath=tmpDirPath.joinpath(names[i]+suffix)
        createSingleModelNb(names[i],nbPath)
        grid[i,4]=widgets.HTML(
            value="""
            <a href="{path}" target="_blank">{text}</a>
            """.format(
                    path = nbPath.as_posix(),
                    text = names[i]+suffix
                )
        )

    return insert_link

def createSingleModelNb(modelName,ReportFilePath):
    # this will be specific to the model and depend on the 
    # available information
    nb = nbf.v4.new_notebook()
    text = """\
    # My first automatic Jupyter Notebook
    This is an auto-generated notebook."""
    code = """\
    %pylab inline
    hist(normal(size=2000), bins=50);"""

    nb['cells'] = [nbf.v4.new_markdown_cell(text),
               nbf.v4.new_code_cell(code) ]
    nbf.write(nb,ReportFilePath)

def modelGridBox(modelName):
    nrows=3 
    grid = widgets.GridspecLayout(nrows, 2)
    grid[1,0]=widgets.HTML(value="""
        <h1>{name}</h1>
        This model specific overview page should be depend on the available 
        properties
        """.format(name=modelName))
    return grid

def modelListGridBox():
    names = list_models()
    buttons=list()
    nrows = len(names)
    grid = widgets.GridspecLayout(nrows, 5)
    for i in range(nrows):
        mn=names[i]
        grid[i, 0] = widgets.Text(
                #layout=widgets.Layout(width='auto', height='auto'),
                value = names[i]
        )
        res=get_single_mvar_value(CompartmentalMatrix,mn)
        out = widgets.Output()
        with out:
            display(res)
        grid[i, 2] =out

        b =  widgets.Button(
                        layout=widgets.Layout(width='auto', height='auto'),
                        description="Create notebook \n from template"
                    )
        b.on_click(funcmakerInsertLink(grid,names,i))
        grid[i, 3] = b
    return grid
