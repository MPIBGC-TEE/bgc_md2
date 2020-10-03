from typing import Dict, List, Set, TypeVar
from frozendict import frozendict
import networkx as nx
import inspect
import importlib

# from ..models.helpers import provided_mvar_values
from .helpers import bgc_md2_computers
from . import non_graph_helpers as ngh
from .graph_helpers import (
    sparse_powerset_graph,
    minimal_target_subgraph_for_single_var,
    minimal_startnodes_for_single_var,
    node_2_string,
    nodes_2_string,
)
# the next imports should not be necessary after the model specific part is factored out
import matplotlib.pyplot as plt
from IPython.display import Math
from IPython.display import display
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from sympy import latex

class MVarSet:
    # This class is explicitly instantiated (even in the model source.py) from a set.
    # The models source.py becomes A (as opposed to THE) way to create a MVarSet instance
    # which can also be created arbitrarily (and programmatically) somewhere else.
    # Not all MVarSets of interest are "Models" some are "ModelRuns" (=Simulations) 
    # The models submodule of bgc_md2 is just a convienience to assemble MVarSets that 
    # express models.
    # Actually in the case of Williams the MvarSet expresses a ModelRun with 
    # units (which motivates a possible second submodule modelruns 
    # that will contain modelruns which will refer to a particular model and
    # just add parameterizations and start values)
    
    # fixme: mm 9/15/2020 
    # 
    # Along with the graph and non_graph helpers the bigger part of the class
    # could be factored out into a more abstract class or meta class (possibly
    # living in a different python package)  that does not depend on the 
    # specific mvars and computers in this package.  
    # This more abstract version would also not
    # include references to specific variables like SmoothReservoirModel in the
    # graph method (and would not implement the graph_method)

    def __init__(self, s):
        self.provided_mvar_values=s




    @property
    def provided_mvar_types(self) -> Set[type]:
        return frozenset(type(v) for v in self.provided_mvar_values)

    #@property
    def computable_mvar_types(self) -> Set[type]:
        return ngh.computable_mvars(
            allComputers=bgc_md2_computers(),
            available_mvars=self.provided_mvar_types
        )

    @property
    def computable_mvar_names(self):
        return [var.__name__ for var in self.computable_mvar_types()]


    def __dir__(self):
        return super().__dir__() + [
            "get_{}".format(name) for name in self.computable_mvar_names
        ]

    def __getattr__(self, name):
        if name.startswith("get_"):
            var_name = name[4:]
            #for var in self.mvars:
            for var in self.computable_mvar_types():
                if var.__name__ == var_name:
                    return lambda: self._get_single_mvar_value(var)
        return super().__getattr__(name)

    def path_dict_to_single_mvar(
            self,
            mvar: type
        ) -> Dict[type, List[Set[type]]]:
        # fixme mm 09-15-2020:
        # should be deprecated since the class MVarSet implements a similar method now
        node = frozenset({mvar})
        spsg = sparse_powerset_graph(bgc_md2_computers())
        graph_min_nodes = minimal_startnodes_for_single_var(spsg, mvar)
        pmvs = self.provided_mvar_types
    
        model_min_nodes = list(filter(lambda n: n.issubset(pmvs), graph_min_nodes))
        if len(model_min_nodes) < 1:
            raise (
                Exception(
                    "The desired mvar can not be computed from the provided mvars:"
                    + node_2_string(pmvs)
                    + "Minimal sets to compute it are"
                    + nodes_2_string(graph_min_nodes)
                )
            )
    
        path_dict = frozendict(
            {
                n: list(nx.all_shortest_paths(spsg, source=n, target=node))
                for n in model_min_nodes
            }
        )
        return path_dict

    def _get_single_mvar_value(
            self,
            mvar: type,
            path: List[Set[type]] = []
        ):  # ->mvar:
        # fixme mm 03-07-2020:
        # This is interesting: The function actually returns
        # an instance of class mvar, I do not know yet how to express that with
        # the static type hint system.
        # (Obviously the return type  is a function of the input types)
    
        pvs = self.provided_mvar_values
        pv_dict = {type(v): v for v in pvs}
        if mvar in [type(v) for v in pvs]:
            return pv_dict[mvar]
    
        path_dict = self.path_dict_to_single_mvar(mvar)
        start_nodes = path_dict.keys()
    
        if path == []:
            default_start_node = sorted(start_nodes, key=lambda node: len(node))[0]
            path = path_dict[default_start_node][0]
        else:
            # check if the given path is among the possible paths
            start_node = path[0]
            if start_node not in path_dict.keys():
                raise (Exception("There are no path to the target with this startnode"))
            starting_here = path_dict[start_node]
            if not path in starting_here:
                raise (Exception("the given path is not possible"))
    
        # create results step by step along the graph
        spsg = sparse_powerset_graph(bgc_md2_computers())
        rg = spsg.subgraph(path).copy()
        rg.nodes[path[0]]["values"] = pvs
        for i in range(1, len(path)):
            computers = rg.get_edge_data(path[i - 1], path[i])[0][
                "computers"
            ]  # if we have more
    
            def apply(comp):
                arg_classes = [p.annotation for p in inspect.signature(comp).parameters.values()]
                arg_values = [pv_dict[cl] for cl in arg_classes]
                res = comp(*arg_values)
                return res
    
            pv_dict.update({ngh.output_mvar(c): apply(c) for c in computers})
    
        return pv_dict[mvar]


    # fixme mm 9/15/2020
    # This method is bgc or even Model specific and would not be part of the more general package
    def graph(self):
        target_var = SmoothReservoirModel
        if target_var not in self.computable_mvar_types():
            return

        srm = self._get_single_mvar_value(target_var)
        fig = plt.figure()
        rect = (0, 0, 0.8, 1.2)  # l, b, w, h
        ax = fig.add_axes(rect)
        ax.clear()
        srm.plot_pools_and_fluxes(ax)
        plt.close(fig)
        return ax.figure
        out = widgets.Output()
        with out:
            display(var.__name__ + "=")
            display(Math(latex(res)))
            # The latex could be filtered to display subscripts better
            # display(res)
        if capture:
            return out
        else:
            display(out)


    # fixme mm 9/15/2020
    # This method is bgc or even Model specific and would not be part of the more general package
    def render(self, var):
        res = self._get_single_mvar_value(var)
        display(Math("\\text{" + var.__name__ + "} =" + latex(res)))
        # The latex could be filtered to display subscripts better
        # display(res)
    

    # fixme mm 9/15/2020
    # This method is bgc or even Model specific and would not be part of the more general package
    @classmethod
    def from_model_name(cls,model_id):
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
        retVal= mod.mvs
        if type(retVal)==cls:
            return  retVal
        else:
            raise Exception(
                "The variable mvs in the target module is not of type {}.".format(cls.__name__)              )
