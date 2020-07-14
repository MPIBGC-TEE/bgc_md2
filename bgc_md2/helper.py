from . import models
from IPython.display import Math
from IPython.display import display
from bgc_md2.models.helpers import computable_mvars
from bgc_md2.models.helpers import get_single_mvar_value
from bgc_md2.resolve.mvars import CompartmentalMatrix
from pathlib import Path
from sympy import latex
import ipywidgets as widgets
import nbformat as nbf
import pkgutil


def list_models():
    exclude_path = Path('./exclude-models.txt')
    if exclude_path.exists():
        exclude_lines = set(line.strip() for line in open(exclude_path))
        exclude_models = set(name for name in exclude_lines
                             if name and not name.startswith('#'))
    else:
        exclude_models = set()
    sub_mod_pkgs= [tup[1] for tup in pkgutil.iter_modules(models.__path__)
                   if tup[2] and tup[1] not in exclude_models]
    return sub_mod_pkgs


def list_models_md():
    names = list_models()
    links = ["[{name}](/tmp/{name})".format(name=name) for name in names]
    return "".join(links)


def createSingleModelNb(model_name, report_file_path):
    model = Model(model_name)
    nb = nbf.v4.new_notebook()

    text = '# {}'.format(model.name)
    t_mvars = 'Computable mvars:\n' + '\n'.join(
        '1. {}'.format(var) for var in model.mvar_names)
    c_imports = 'import bgc_md2.helper as h'
    c_model = 'model = h.Model({})'.format(repr(model.name))
    c_render = 'for var in model.mvars:\n    model.render(var)'

    nb['cells'] = [
        nbf.v4.new_markdown_cell(text),
        nbf.v4.new_markdown_cell(t_mvars),
        nbf.v4.new_code_cell(c_imports),
        nbf.v4.new_code_cell(c_model),
        nbf.v4.new_code_cell(c_render),
    ]
    nbf.write(nb, report_file_path)


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


class Model:

    def __init__(self, name):
        self.name = name
        self.mvars = computable_mvars(name)
        self.mvar_names = [var.__name__ for var in self.mvars]
        
    def render(self, var, capture=False):
        res = get_single_mvar_value(var, self.name)
        out = widgets.Output()
        with out:
            display(var.__name__ + '=')
            display(Math(latex(res)))
            # The latex could be filtered to display subscripts better
            #display(res)
        if capture:
            return out
        else:
            display(out)


def modelVBox(model_name):
    model = Model(model_name)
    box = widgets.VBox(
        [
            widgets.HTML(value="""
                <h1>{name}</h1>
                This model specific overview page depends on the available 
                properties
                """.format(name=model_name)
            ),
            widgets.HTML(
                "computable_mvars( perhaps as links to the docs or some graph ui ...)"
                +"<ol>\n"
                +"\n".join('<li>{}</li>'.format(var) for var in model.mvar_names)
                +"</ol>\n"
            ),
        ]
        +
        [model.render(var, capture=True) for var in model.mvars]
    )
    b =  widgets.Button(
        layout=widgets.Layout(width='auto', height='auto'),
        description="Create notebook from template"
    )
    b.on_click(funcmakerInsertLinkInToBox(box, model_name))
    box.children += (b,)
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


def modelListGridBox(callback_inspect_model):
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

        button_inspect_model = widgets.Button(
            description='Inspect model',
        )
        button_inspect_model.on_click(
            (lambda name: lambda button: callback_inspect_model(name))(name))

        button_create_notebook = widgets.Button(
            description='Create notebook from template',
        )
        button_create_notebook.on_click(
            callback_factory_create_notebook(grid, i, name))

        grid[i, 8] = widgets.VBox([
            button_inspect_model,
            button_create_notebook,
        ])

    return grid
