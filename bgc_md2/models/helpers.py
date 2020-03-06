import importlib,inspect
import bgc_md2.resolve.non_graph_helpers as ngh
from bgc_md2.resolve.graph_helpers import ( 
    sparse_powerset_graph
)

def available_mvars(model_id):
    sep="."
    parent_mod_name=sep.join(__name__.split(sep)[:-1])
    #parent_mod_name='bgc_md2.models'
    # in case the module has been created or changed
    # after the current session started
    importlib.invalidate_caches()
    
    mod=importlib.import_module(sep+str(model_id)+sep+"source",package=parent_mod_name)
    return frozenset(type(v) for v in mod.specialVars)

def computable_mvars(model_id):
    amvs=available_mvars(model_id)
    cmod=importlib.import_module('bgc_md2.resolve.computers')
    computers=frozenset(
        [
            getattr(cmod,c) for c in cmod.__dir__() 
            if inspect.isfunction(getattr(cmod,c))
        ]
    )
    return ngh.computable_mvars(
        allComputers=computers,
        available_mvars=available_mvars(model_id)
    )

def get_mvar(model_id):
    spgs=sparse_powerset_graph(computers)
    # create a resultgraph and replace the nodes with 
    # results for the instances
    rg=deepcopy(spgs)
    #return frozenset

