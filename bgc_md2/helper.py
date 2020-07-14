import pkgutil
from . import models
import ipywidgets as widgets
from IPython.display import  display,HTML,Markdown,Latex,Math
from pathlib import Path
from sympy import Matrix,Symbol,latex
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



def funcmakerInsertLink(grid,i,j,name):
    def insert_link(b):
        # called for side effect on grid object
        tmpDirPath= Path("./tmp") # has to be relative for jupyter to open the file on click (so the exact location depends on where the notebook server was started)
        tmpDirPath.mkdir(exist_ok=True)   
        suffix=".ipynb"
        nbPath=tmpDirPath.joinpath(name+suffix)
        createSingleModelNb(name,nbPath)
        grid[i,j]=widgets.HTML(
            value="""
            <a href="{path}" target="_blank">{text}</a>
            """.format(
                    path = nbPath.as_posix(),
                    text = name+suffix
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
    # this could probably become a vertical box where line numbers do not have to be given explicitly
    max_row=4 
    max_col=2
    grid = widgets.GridspecLayout(max_row+1, max_col+1)
    grid[0,0:max_col]=widgets.HTML(value="""
        <h1>{name}</h1>
        This model specific overview page should be depend on the available 
        properties
        """.format(name=modelName))
    
    grid[1,0]=widgets.HTML(value=""" CompartmentalMatrix:(add maybe a link to the docs """)
    res=get_single_mvar_value(CompartmentalMatrix,modelName)
    out = widgets.Output()
    with out:
        display(res)
    grid[2,0:max_col-1] =out

    b =  widgets.Button(
                    layout=widgets.Layout(width='auto', height='auto'),
                    description="Create notebook \n from template"
                    )
    b.on_click(funcmakerInsertLink(grid,3,0,modelName))
    grid[0, max_col] = b
    return grid

#################################################################################
def funcmakerInsertLinkInToBox(grid,name):
    def insert_link(b):
        # called for side effect on grid object
        tmpDirPath= Path("./tmp") # has to be relative for jupyter to open the file on click (so the exact location depends on where the notebook server was started)
        tmpDirPath.mkdir(exist_ok=True)   
        suffix=".ipynb"
        nbPath=tmpDirPath.joinpath(name+suffix)
        createSingleModelNb(name,nbPath)
        old_chs=grid.children
        new_chs=( 
            widgets.HTML(
                value="""
                <a href="{path}" target="_blank">{text}</a>
                """.format(
                        path = nbPath.as_posix(),
                        text = name+suffix
                    )
             )
        ,) 
        # we change the children tuple
        grid.children=old_chs + new_chs

    return insert_link

def modelVBox(modelName):
    
    mvs=computable_mvars(modelName)
    res=get_single_mvar_value(CompartmentalMatrix,modelName)
    def output(var):
        res=get_single_mvar_value(var,modelName)
        out = widgets.Output()
        with out:
            display(var.__name__+"=")
            display(Math(latex(res)))
            # The latex could be filtered to display subscripts better
            #display(res)
        return out

    box= widgets.VBox(
        [
            widgets.HTML(value="""
                <h1>{name}</h1>
                This model specific overview page depends on the available 
                properties
                """.format(name=modelName)
            )
            ,
            widgets.HTML(
                "computable_mvars( perhaps as links to the docs or some graph ui ...)"
                +"<ol>\n"
                +"\n".join(["<li> "+var.__name__+" </li>" for var in mvs])
                +"</ol>\n"
            )
        ]
        +
        [ output(var) for var in mvs ]
    )
    b =  widgets.Button(
                    layout=widgets.Layout(width='auto', height='auto'),
                    description="Create notebook \n from template"
                    )
    b.on_click(funcmakerInsertLinkInToBox(box,modelName))
    box.children=box.children+(b,)
    return box 

#################################################################################

def callback_factory_create_notebook(grid, i, name):

    def callback_create_notebook(button):
        # called for side effect on g
        tmp_dir_path = Path('./tmp') # has to be relative for jupyter to open the file
        tmp_dir_path.mkdir(exist_ok=True)
        suffix = '.ipynb'
        nbPath = tmp_dir_path.joinpath(name + suffix)

        createSingleModelNb(name, nbPath)

        grid[i, 9] = widgets.HTML(
            value="""
            <a href="{path}" target="_blank">{text}</a>
            """.format(
                path = nbPath.as_posix(),
                text = name + suffix
            )
        )

    return callback_create_notebook


def modelListGridBox():
    names = list_models()
    grid = widgets.GridspecLayout(len(names), 10)

    for i, name in enumerate(names):
        grid[i, 0] = widgets.Text(
            layout=widgets.Layout(width='auto', height='auto'),
            value=name,
        )

        res = get_single_mvar_value(CompartmentalMatrix, name)
        out = widgets.Output()
        with out:
            display(res)
        grid[i, 1:7] = out

        b = widgets.Button(
            layout=widgets.Layout(width='auto', height='auto'),
            description='Create notebook from template',
        )
        b.on_click(callback_factory_create_notebook(grid, i, name))
        grid[i, 8] = b

    return grid
