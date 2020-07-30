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
import matplotlib.pyplot as plt
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel


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
    c_graph = 'model.graph()'
    c_render = 'for var in model.mvars:\n    model.render(var)'
    
    nb['cells'] = [
        nbf.v4.new_markdown_cell(text),
        nbf.v4.new_markdown_cell(t_mvars),
        nbf.v4.new_code_cell(c_imports),
        nbf.v4.new_code_cell(c_model),
        nbf.v4.new_code_cell(c_graph),
        nbf.v4.new_code_cell(c_render),
    ]
    nbf.write(nb, report_file_path)


def createSingleModelNbFile(model_name):
    tmp_dir = Path('./tmp') # has to be relative for jupyter to open the file on click (so the exact location depends on where the notebook server was started)
    tmp_dir.mkdir(exist_ok=True)
    file_name = model_name + '.ipynb'
    nb_path = tmp_dir.joinpath(file_name)
    createSingleModelNb(model_name, nb_path)
    return file_name, nb_path


class Model:
    # mm @Thomas:
    #
    # 1.
    # in the current implementation the Model has to be reinstanciated
    # when the source.py of a model changes on disk 
    # (For the dot tab completion even the whole class would have to be
    # regenerated.) As a consequence it should be regarded as an ultra
    # light weight shortlived (made on demand) 
    # interface to ipython (which makes working with classes convinient
    # (e.g. by tab completion)
    # see comments on __init__
    #
    # 2.
    # Model as class name is used in the backend packages and 
    # therefore aquired a special connotation.
    # We might have to look for a more specific name

    def __init__(self, name):
        self.name = name
        # mm @Thomas:
        # the next two properties will get stale when the underlying model source.py changes on disk
        # maybe they do not have to be stored ...
        # which would not be a problem it this Model instance is always ready to die
        # and regenerated. 

        # I have no final opinion just a reminder to keep on the lookout
        # for alternative approaches  (see above...)
        # that might be a bit cheaper computationaly by not bundling things into an instance that
        # has to be discarded completely 

        self.mvars = computable_mvars(name) # could be a @property decorated methods
        self.mvar_names = [var.__name__ for var in self.mvars]

    def render(self, var):
        res = get_single_mvar_value(var, self.name)
        display(Math('\\text{' + var.__name__ + '} =' + latex(res)))
            # The latex could be filtered to display subscripts better
            #display(res)

    def graph(self):
        # on demand computation is used
        # I am aware of the possibility of model.computable_mvars
        cmvs = computable_mvars(self.name)
        target_var = SmoothReservoirModel
        if target_var not in cmvs:
            return

        srm = get_single_mvar_value(target_var, self.name)
        fig = plt.figure()
        rect = (0, 0, 0.8, 1.2)  # l, b, w, h
        ax = fig.add_axes(rect)
        ax.clear()
        srm.plot_pools_and_fluxes(ax)
        return ax.figure


#################################################################################

def button_callback(function, *args):

    def callback(button):
        function(*args)

    return callback


class ModelListGridBox(widgets.GridspecLayout):

    def __init__(self, inspection_box):
        self.inspection_box = inspection_box
        self.names = list_models()
        super().__init__(len(self.names), 10)
        self.populate()

    def inspect_model(self, name):
        self.inspection_box.update(name)

    def populate(self):
        for i, name in enumerate(self.names):
            button_inspect_model = widgets.Button(
                description=name,
            )
            button_inspect_model.on_click(
                button_callback(self.inspect_model, name))
            self[i, 0] = button_inspect_model

            res = get_single_mvar_value(CompartmentalMatrix, name)
            out = widgets.Output()
            with out:
                display(res)
            self[i, 1:9] = out


class ModelInspectionBox(widgets.VBox):

    nb_link_box = None

    def create_notebook(self, model_name):
        file_name, nb_path = createSingleModelNbFile(model_name)
        self.nb_link_box.children = (
            widgets.HTML(
                value = """
                <a href="{path}" target="_blank">{text}</a>
                """.format(
                    path = nb_path.as_posix(),
                    text = file_name,
                )
            ),
        )

    def update(self, model_name):
        model = Model(model_name)

        self.children = (
            widgets.HTML(value="""
                <h1>{name}</h1>
                Overview
                """.format(name=model_name)
            ),
            widgets.HTML(
                "computable_mvars( @Thomas perhaps as links to the docs or some graph ui ...)"
                +"<ol>\n"
                +"\n".join('<li>{}</li>'.format(var) for var in model.mvar_names)
                +"</ol>\n"
            ),
        )

        graph = model.graph()
        if graph:
            graph_out = widgets.Output()
            with graph_out:
                display(graph)
            self.children += (graph_out,)
            
        rendered_vars = widgets.Output()
        with rendered_vars:
            for var in model.mvars:
                model.render(var)
        self.children += (rendered_vars,)

        b =  widgets.Button(
            layout=widgets.Layout(width='auto', height='auto'),
            description="Create notebook from template"
        )
        b.on_click(button_callback(self.create_notebook, model_name))
        self.children += (b,)

        self.nb_link_box = widgets.VBox()
        self.children += (self.nb_link_box,)
