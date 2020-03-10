import importlib,inspect
from typing import Dict,List,Set,TypeVar
from frozendict import frozendict
from copy import deepcopy
from inspect import signature
import bgc_md2.resolve.non_graph_helpers as ngh
import networkx as nx
from bgc_md2.resolve.graph_helpers import ( 
    sparse_powerset_graph
    ,minimal_target_subgraph_for_single_var
    ,minimal_startnodes_for_single_var
    ,node_2_string
    ,nodes_2_string
)
def bgc_md2_computers():
    sep='.'
    pkg_name=__name__.split(sep)[0]
    cmod=importlib.import_module('.resolve.computers',package=pkg_name)
    return frozenset(
        [
            getattr(cmod,c) for c in cmod.__dir__() 
            if inspect.isfunction(getattr(cmod,c))
        ]
    )

def provided_mvar_values(model_id):
    sep="."
    parent_mod_name=sep.join(__name__.split(sep)[:-1])
    #parent_mod_name='bgc_md2.models'
    # in case the module has been created or changed
    # after the current session started
    importlib.invalidate_caches()
    
    mod=importlib.import_module(sep+str(model_id)+sep+"source",package=parent_mod_name)
    return mod.specialVars
    
def provided_mvars(model_id:str)->Set[type]:
    return frozenset(type(v) for v in provided_mvar_values(model_id))

def computable_mvars(model_id:str)->Set[type]:
    return ngh.computable_mvars(
        allComputers=bgc_md2_computers(),
        available_mvars=provided_mvars(model_id)
    )

# fixme: mm 03-07-2020
# - generalizing from single mvar to node should be easy
def path_dict_to_single_mvar(
        mvar:type
        ,model_id:str
    )->Dict[type,List[Set[type]]]:
    node=frozenset({mvar}) 
    spsg=sparse_powerset_graph(bgc_md2_computers())
    graph_min_nodes=minimal_startnodes_for_single_var(spsg,mvar) 
    model_min_nodes=list(
        filter(
            lambda n:n.issubset(provided_mvars(model_id))
            ,graph_min_nodes
        )
    )
    if len(model_min_nodes)<1:
        raise(
            Exception("The desired mvar can not be computed from the provided mvars:"+node_2_string(pms)+"Minimal sets to compute it are"+nodes_2_string(min_nodes)))

    path_dict=frozendict(
        {
            n : list(nx.all_shortest_paths(spsg,source=n,target=node)) 
            for n in model_min_nodes
        }
    )
    return path_dict

def get_single_mvar_value(
        mvar        :type
        ,model_id   :str
        ,path    :List[Set[type]]=[]
    ): # ->mvar: 
    # fixme mm 03-07-2020:
    # This is interesting: The function actually returns
    # an instance of class mvar, I do not know yet how to express that with
    # the static type hint system. 
    # (Obviously the return type  is a function of the input types)

    pvs=provided_mvar_values(model_id)
    pv_dict=frozendict({type(v):v for v in pvs})
    if mvar in [type(v) for v in pvs]:
        return pv_dict[mvar]
    
    path_dict= path_dict_to_single_mvar(mvar,model_id)
    start_nodes=path_dict.keys()

    if path==[]:
        default_start_node=sorted(start_nodes,key=lambda node:len(node))[0]
        path=path_dict[default_start_node][0]
    else:
        #check if the given path is among the possible paths
        start_node=path[0]
        if start_node not in path_dict.keys():
            raise(Exception("There are no path to the target with this startnode"))
        starting_here=path_dict[start_node]
        if not path in starting_here:
            raise(Exception("the given path is not possible"))

    # create a resultgraph from the path replace the nodes with 
    spsg=sparse_powerset_graph(bgc_md2_computers())
    rg=spsg.subgraph(path).copy()
    rg.nodes[path[0]]['values']=pvs
    for i in range(1,len(path)):
        computers= rg.get_edge_data(path[i-1],path[i])[0]['computers']
        for c in computers:
            arg_classes=[p.annotation for p in signature(c).parameters.values()]
            arg_values=[pv_dict[cl] for cl in arg_classes]
            print(arg_classes)
            print(arg_values)
            res=c(*arg_values)
            print(res)
    # Attach the values of the pvs as node data to the startnode
    return rg
    # now we go through the edges of the path and apply the stored computers 
    # to the compute instances
    

