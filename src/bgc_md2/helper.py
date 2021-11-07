from IPython.display import Math
from IPython.display import display
from typing import  Tuple, Dict, List, Set, TypeVar
from pathlib import Path
from sympy import latex
from frozendict import frozendict
import ipywidgets as widgets
import nbformat as nbf
import pkgutil
import importlib
import matplotlib.pyplot as plt
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

# from .models.helpers import computable_mvars, get_single_mvar_value
from .resolve.mvars import CompartmentalMatrix, StateVariableTuple
from ComputabilityGraphs.CMTVS import CMTVS
from ComputabilityGraphs.helpers import module_computers
from functools import _lru_cache_wrapper
import inspect
from .resolve import computers as cmod
from . import models

def combine(
        d1: frozendict,
        d2: frozendict
    ) -> frozendict:
    return {
        t[0]: t[1] 
        for t in (
            [(k, v) for k, v in d1.items()]
            +
            [(k, v) for k, v in d2.items()]
        )
    }


def batchSlices(nland, nproc):
    return [
        slice(i*nproc,min((i+1)*nproc,nland))
        for i in range(int(nland/nproc)+1)
    ]

    # fixme mm:

# @Thomas
# If you think that there should be a plumbing class in between (a Model in the Model-View-Controller sense)
# feel free to add it. It is probably not the commented one below anyway...
#

def list_target_models(
    target_classes=frozenset({CompartmentalMatrix, StateVariableTuple}),
    explicit_exclude_models: Set[str] = frozenset()
) -> List[CMTVS]:
    sub_mod_pkgs = list_models(explicit_exclude_models)

    def pred(mn):
        mvs = CMTVS_from_model_name(mn)
        return (target_classes.issubset(mvs.computable_mvar_types()))

    hits = [mod for mod in sub_mod_pkgs if pred(mod)]
    return hits


def list_models(explicit_exclude_models: Set[str] = frozenset()):
    exclude_path = Path("./exclude-models.txt")
    if exclude_path.exists():
        exclude_lines = set(line.strip() for line in open(exclude_path))
        exclude_models_from_file = set(
            name for name in exclude_lines if name and not name.startswith("#")
        )
        exclude_models = exclude_models_from_file.union(
                explicit_exclude_models
        )
    else:
        exclude_models = explicit_exclude_models
    sub_mod_pkgs = [
        tup[1]
        for tup in pkgutil.iter_modules(models.__path__)
        if tup[2] and tup[1] not in exclude_models
    ]
    return sub_mod_pkgs


def list_models_md():
    names = list_models()
    links = ["[{name}](/tmp/{name})".format(name=name) for name in names]
    return "".join(links)


def createSingleModelNb(model_name, report_file_path):
    mvs = CMTVS_from_model_name(model_name)
    nb = nbf.v4.new_notebook()

    text = "# {}".format(model_name)
    # rather let the user do it programmatically
    # t_mvars = "Computable mvars:\n" + "\n".join(
    #     "1. {}".format(var) for var in mvs.computable_mvar_names
    # )
    c_imports = """import bgc_md2.helper as h
import importlib """
    # c_mvs = "mvs = h.CMTVS.from_model_name({})".format(repr(model_name))
    c_mvs = """importlib.invalidate_caches()
mod = importlib.import_module('bgc_md2.models.{}.source')
mvs = mod.mvs""".format(model_name)
    c_graph = "h.compartmental_graph(mvs)"
    c_mvars = "mvs.computable_mvar_types"
    c_render = """for var in mvs.computable_mvar_types():
    res = mvs._get_single_value(var)
    h.latex_render(var,res)"""

    nb["cells"] = [
        nbf.v4.new_markdown_cell(text),
        #nbf.v4.new_markdown_cell(t_mvars),
        nbf.v4.new_code_cell(c_imports),
        nbf.v4.new_code_cell(c_mvs),
        nbf.v4.new_code_cell(c_graph),
        nbf.v4.new_code_cell(c_mvars),
        nbf.v4.new_code_cell(c_render),
    ]
    nbf.write(nb, report_file_path)


def createSingleModelNbFile(model_name):
    tmp_dir = Path(
        "./tmp"
    )  # has to be relative for jupyter to open the file on click (so the exact location depends on where the notebook server was started)
    tmp_dir.mkdir(exist_ok=True)
    file_name = model_name + ".ipynb"
    nb_path = tmp_dir.joinpath(file_name)
    createSingleModelNb(model_name, nb_path)
    return file_name, nb_path

#################################################################################
def funcmakerInsertLinkInToBox(grid, name):
    def insert_link(b):
        # called for side effect on grid object
        tmpDirPath = Path(
            "./tmp"
        )  # has to be relative for jupyter to open the file on click (so the exact location depends on where the notebook server was started)
        tmpDirPath.mkdir(exist_ok=True)
        suffix = ".ipynb"
        nbPath = tmpDirPath.joinpath(name + suffix)
        createSingleModelNb(name, nbPath)
        old_chs = grid.children
        new_chs = (
            widgets.HTML(
                value="""
                <a href="{path}" target="_blank">{text}</a>
                """.format(
                    path=nbPath.as_posix(), text=name + suffix
                )
            ),
        )
        # we change the children tuple
        grid.children = old_chs + new_chs

    return insert_link


#def modelVBox(model_name):
#    mvs = CMTVS.from_model_name(model_name)
#    # on demand computation is used
#    # I am aware of the possibility of mvs.computable_mvars
#    cmvs = computable_mvars(model_name)
#    target_var = SmoothReservoirModel
#    pictlist = []
#    if target_var in cmvs:
#        srm = mvs._get_single_mvar_value(target_var)
#        graph_out = widgets.Output()
#        fig = plt.figure()
#        rect = 0, 0, 0.8, 1.2  # l, b, w, h
#        ax = fig.add_axes(rect)
#        with graph_out:
#            ax.clear()
#            srm.plot_pools_and_fluxes(ax)
#            display(ax.figure)
#        pictlist = [graph_out]
#
#    box = widgets.VBox(
#        [
#            widgets.HTML(
#                value="""
#                <h1>{name}</h1>
#                Overview 
#                """.format(
#                    name=model_name
#                )
#            ),
#            widgets.HTML(
#                "computable_mvars( @Thomas perhaps as links to the docs or some graph ui ...)"
#                + "<ol>\n"
#                + "\n".join("<li>{}</li>".format(var) for var in mvs.computable_mvar_names)
#                + "</ol>\n"
#            ),
#        ]
#        + pictlist
#        + [
#            latex_render(
#                var,
#                mvs._get_single_value(var),
#            )
#            for var in mvs.computable_mvar_types()
#          ]
#    )
#    b = widgets.Button(
#        layout=widgets.Layout(width="auto", height="auto"),
#        description="Create notebook from template",
#    )
#    b.on_click(funcmakerInsertLinkInToBox(box, model_name))
#    box.children += (b,)
#    return box


##############################################################################


def button_callback(function, *args):
    def callback(button):
        function(*args)

    return callback


class ModelListGridBox(widgets.GridspecLayout):
    def __init__(
        self,
        inspection_box,
        explicit_exclude_models: Set[str] = frozenset()
    ):
        self.inspection_box = inspection_box
        #self.names = list_target_models(frozenset((CompartmentalMatrix,)))
        self.names = list_models(explicit_exclude_models)
        super().__init__(len(self.names), 10)
        self.populate()

    def inspect_mvs(self, name):
        self.inspection_box.update(name)

    def populate(self):
        for i, name in enumerate(self.names):
            button_inspect_mvs = widgets.Button(description=name,)
            button_inspect_mvs.on_click(
                button_callback(self.inspect_mvs, name)
            )
            self[i, 0] = button_inspect_mvs
            # fixme mm:
            # Just temporarily disabled the output of any real model parameters.
            # the overview takes way too lang to appear in the notebook
            # The design is also flawed (ticket) by the fact that 
            # the rows are all equally tall, which creates a lot of whitespace
            mvs = CMTVS_from_model_name(name)
            target_class = CompartmentalMatrix
            if target_class in mvs.computable_mvar_types():
                res = mvs._get_single_value(target_class)
            else:
                res = f"can not compute {target_class}" 
            out = widgets.Output()
            with out:
                display(res)
            self[i, 1:9] = out

class GeneralMvarSetListGridBox(widgets.GridspecLayout):

    def __init__(
        self,
        inspection_box: 'MvarSetInspectionBox',
        target_classes: Tuple = (CompartmentalMatrix,StateVariableTuple),
        explicit_exclude_models: Set[str] = frozenset()
    ):
        self.inspection_box = inspection_box
        self.target_classes = target_classes
        self.names = list_target_models(frozenset(target_classes))
        super().__init__(len(self.names), 10)
        self.populate()

    def inspect_mvs(self, name):
        self.inspection_box.update(name)

    def populate(self):
        for i, name in enumerate(self.names):
            button_inspect_mvs = widgets.Button(description=name,)
            button_inspect_mvs.on_click(
                button_callback(self.inspect_mvs, name)
            )
            self[i, 0] = button_inspect_mvs
            mvs = CMTVS_from_model_name(name)
            results = [ 
                mvs._get_single_mvar_value(target_class)
                for target_class in self.target_classes
            ]
            out = widgets.Output()
            with out:
                for res in results:
                    display(res)
            self[i, 1:9] = out



class MvarSetInspectionBox(widgets.VBox):

    nb_link_box = None

    def create_notebook(self, model_name):
        file_name, nb_path = createSingleModelNbFile(model_name)
        self.nb_link_box.children = (
            widgets.HTML(
                value="""
                <a href="{path}" target="_blank">{text}</a>
                """.format(
                    path=nb_path.as_posix(), text=file_name,
                )
            ),
        )

    def update(self, model_name):
        mvs = CMTVS_from_model_name(model_name)

        self.children = (
            widgets.HTML(
                value="""
                <h1>{name}</h1>
                Overview
                """.format(
                    name=model_name
                )
            ),
            widgets.HTML(
                "computable_mvars( @Thomas perhaps as links to the docs or some graph ui ...)"
                + "<ol>\n"
                + "\n".join("<li>{}</li>".format(var) for var in mvs.computable_mvar_types())
                + "</ol>\n"
            ),
        )

        graph = compartmental_graph(mvs)
        if graph:
            graph_out = widgets.Output()
            with graph_out:
                display(graph)
            self.children += (graph_out,)

        rendered_vars = widgets.Output()
        with rendered_vars:
            for var in mvs.computable_mvar_types():
                latex_render(var,mvs._get_single_value(var))
        self.children += (rendered_vars,)

        b = widgets.Button(
            layout=widgets.Layout(width="auto", height="auto"),
            description="Create notebook from template",
        )
        b.on_click(button_callback(self.create_notebook, model_name))
        self.children += (b,)

        self.nb_link_box = widgets.VBox()
        self.children += (self.nb_link_box,)

def CMTVS_from_model_name(model_id):
        """
        convenience method to get the instance from a submodule of bgc.models
        by just giving the name of the submodule
        """
        sep = "."
        #model_mod_name = sep.join(__name__.split(sep)[:-1])
        models_mod_name = 'bgc_md2.models'
        # in case the module has been created or changed
        # after the current session started
        importlib.invalidate_caches()

        mod = importlib.import_module(
            sep + str(model_id) + sep + "source", package=models_mod_name
        )
        cls = CMTVS
        retVal= mod.mvs
        if type(retVal)==cls:
            return  retVal
        else:
            raise Exception(
                "The variable mvs in the target module is not of type {}.".format(cls.__name__)              )

def bgc_md2_computers():
    sep = "."
    pkg_name = __name__.split(sep)[0]
    cmod = importlib.import_module(".resolve.computers", package=pkg_name)
    def pred(a):
        return inspect.isfunction(a) or isinstance(a,_lru_cache_wrapper)
    return frozenset(
        [
            getattr(cmod, c)
            for c in cmod.__dir__()
            if pred(getattr(cmod, c)) 
        ]
    )

def compartmental_graph(mvs):
    # fixme mm 11-04 2021 
    # Could we have a computer for the graph
    # and keep only the display part here
    target_var = SmoothReservoirModel
    if target_var not in mvs.computable_mvar_types():
        return

    srm = mvs._get_single_value(target_var)
    fig = plt.figure()
    rect = (0, 0, 0.8, 1.2)  # l, b, w, h
    ax = fig.add_axes(rect)
    ax.clear()
    srm.plot_pools_and_fluxes(ax)
    plt.close(fig)
    return ax.figure


# fixme mm 11/04/2021
def latex_render(var,res,capture=True):
    out = widgets.Output()
    with out:
        display(
            "\\text{var.__name__}" + "=" +Math(latex(res))
        )
        # The latex could be filtered to display subscripts better
        # display(res)
    if capture:
        return out
    else:
        display(out)


