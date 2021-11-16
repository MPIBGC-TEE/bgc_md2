from functools import lru_cache, reduce
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS, BASE_COLORS, TABLEAU_COLORS
from pygraphviz.agraph import AGraph
from typing import List, Set, Tuple, Callable
from copy import deepcopy
from frozendict import frozendict
from testinfrastructure.helpers import pp, pe

from .non_graph_helpers import (
    # computable_mvars
    # ,directly_computable_mvars
    input_mvars
    # ,output_mvar
    ,
    arg_set
    # ,arg_set_set
    ,
    all_mvars
    # ,applicable_computers
    ,
    all_computers_for_mvar,
    pretty_name,
)


def compset_2_string(compset):
    return "{" + ",".join([pretty_name(c) for c in compset]) + "}"


def node_2_string(node):
    return "{" + ",".join([pretty_name(v) for v in node]) + "}"


def nodes_2_string(node):
    return "[ " + ",".join([node_2_string(n) for n in node]) + " ]"


def edge_2_string(e):
    return "(" + node_2_string(e[0]) + "," + node_2_string(e[1]) + ")"


def immutable_edge(edge):
    s, d, dat = edge
    return (s, d, frozendict(dat))


def equivalent_singlegraphs(g1_single: nx.DiGraph, g2_single: nx.DiGraph) -> bool:
    return all(
        [
            g1_single.get_edge_data(*e) == g2_single.get_edge_data(*e)
            for e in g1_single.edges()
        ]
        + [
            g1_single.get_edge_data(*e) == g2_single.get_edge_data(*e)
            for e in g2_single.edges()
        ]
    ) & (g1_single.nodes() == g2_single.nodes())


def equivalent_multigraphs(
    g1_multi: nx.MultiDiGraph, g2_multi: nx.MultiDiGraph
) -> bool:
    # since we deal with a multigraph
    # it is possible to get several edges between two nodes.
    # The method get_edge_data returns a dictionary
    # with the numbers of these edges as keys.
    # But we want to consider two graphs equivalent if the resulting
    # SET is equal, in other words:
    # If a graph has two edges EACH from : A->B
    # we do not care which of the edges has which computerset
    # Therefore we compare the set of computersets belonging

    g1_single = toDiGraph(g1_multi)
    g2_single = toDiGraph(g2_multi)
    return all(
        [
            g1_single.get_edge_data(*e) == g2_single.get_edge_data(*e)
            for e in g1_single.edges()
        ]
        + [
            g1_single.get_edge_data(*e) == g2_single.get_edge_data(*e)
            for e in g2_single.edges()
        ]
    ) & (g1_single.nodes() == g2_single.nodes())


@lru_cache(maxsize=None)
def arg_set_graph(mvar: type, allComputers: Set[Callable]) -> nx.MultiDiGraph:
    # return the subgraph of arg_name_sets for all computers that
    # return this mvar
    # For compatibility we return a multigraph
    # although we do not have more than one edge between a pair
    # of nodes
    target = frozenset({mvar})
    g = nx.MultiDiGraph()
    g.add_node(target)
    for c in all_computers_for_mvar(mvar, allComputers):
        g.add_edge(arg_set(c), target, computers=frozenset({c}))
    return g


def product_graph(*graphs: Tuple[nx.MultiDiGraph]) -> nx.MultiDiGraph:
    return reduce(lambda u, v: product_graph_2(u, v), graphs)


def product_graph_2(
            g1: nx.MultiDiGraph,
            g2: nx.MultiDiGraph
        ) -> nx.MultiDiGraph:
    cp = nx.cartesian_product(g1, g2)
    prod = nx.MultiDiGraph()

def product_graph_2(g1: nx.MultiDiGraph, g2: nx.MultiDiGraph) -> nx.MultiDiGraph:
    cp = nx.cartesian_product(g1, g2)
    prod = nx.MultiDiGraph()

    for edge in cp.edges(data=True):
        s_tup, d_tup, data = edge
        prod.add_edge(
            reduce(lambda acc, n: acc.union(n), s_tup),
            reduce(lambda acc, n: acc.union(n), d_tup),
            computers=data["computers"],
        )

    # the cartesian product can also contain nodes that
    # are not connected.
    for cart_node in cp.nodes:
        prod_node = reduce(lambda acc, n: acc.union(n), cart_node)

        prod.add_node(prod_node)

    return prod


def initial_sparse_powerset_graph(computers: Set[Callable]) -> nx.MultiDiGraph:
    spsg = nx.MultiDiGraph()
    allMvars = all_mvars(computers)
    for v in allMvars:
        spsg.add_edges_from(arg_set_graph(v, computers).edges(data=True))
    return spsg


def update_step(spsg: nx.MultiDiGraph, computers: Set[Callable]) -> nx.MultiDiGraph:
    new = deepcopy(spsg)
    start_nodes = [n for n in new.nodes() if len(new.in_edges(n)) == 0]
    # fixme: the startnode choice is too exclusive
    for node in start_nodes:
        # print(tuple(v for v in node))
        pg = product_graph(*[arg_set_graph(v, computers) for v in node])
        # new.add_edges_from(pg.edges(data=True))
        new = nx.compose(new, pg)

    return new

@lru_cache
def sparse_powerset_graph(computers: Set[Callable]) -> nx.MultiDiGraph:
    old = initial_sparse_powerset_graph(computers)
    new = update_step(old, computers)
    while not (equivalent_multigraphs(old, new)):
        old = deepcopy(new)
        new = update_step(old, computers)
        # print(equivalent_multigraphs(old, new))
    return new


def update_generator(computers: Set[Callable], max_it: int) -> List[nx.MultiDiGraph]:

    if max_it < 0:
        raise IndexError("update sequence indices have to be larger than 0")

    val = initial_sparse_powerset_graph(computers)
    yield val
    old = deepcopy(val)
    val = update_step(val, computers)

    # print(equivalent_multigraphs(old, val))

    counter = 1
    while max_it > counter and not (equivalent_multigraphs(old, val)):
        yield val
        old = deepcopy(val)
        val = update_step(val, computers)
        counter += 1
        # print("counter", counter, "equivalent?", equivalent_multigraphs(old, val))


def toDiGraph(g_multi: nx.MultiDiGraph) -> nx.DiGraph:
    def edgeDict_to_set(ed):
        target = "computers"
        comp_set_set = frozenset([
            v[target]
            for v in ed.values() if target in v.keys()
        ])
        return comp_set_set

    g_single = nx.DiGraph()
    for e in g_multi.edges():
        s, t = e
        edgeDict = g_multi.get_edge_data(s, t)
        comp_set_set = edgeDict_to_set(edgeDict)
        if g_single.has_edge(s, t):
            comp_set_set = comp_set_set.union(
                g_single.get_edge_data(s, t)["computers"]
            )
        g_single.add_edge(s, t, computers=comp_set_set)
    return g_single


def minimal_startnodes_for_single_var(spg: nx.Graph, targetVar: type):
    """ spg is a sparse powerset Graph, which means that it only contains all one element sets as targets."""
    # We first create a graph with the direction of the edges reversed
    rev_spg = spg.reverse(copy=True)
    targetSet = frozenset({targetVar})
    res = nx.single_source_shortest_path(rev_spg, source=targetSet)
    # res=nx.shortest_path(spg,target=targetSet)
    possible_startnodes = frozenset(res.keys())
    minimal_startnodes = [
        n for n in filter(lambda n: not (n.issuperset(targetSet)), possible_startnodes)
    ]
    return frozenset(minimal_startnodes)


# def target_subgraph(
#         spg:nx.MultiDiGraph
#        ,targetNode:Set[type]
#    )->nx.MultiDiGraph:
#    def filter_func(n):
#        # remove every set that contains one of the variables we are looking for ...
#        return not(any([ (v in n) for v in targetNode]))
#
#    if spg.has_node(targetNode):
#        # The targetNode is already part of the spg.
#        # (because it has been added in one of the update
#        # steps as an argument set of a computer)
#        # in which case we can simply return the startnodes
#        # of the paths leading to it.
#        spg=spg.copy()
#    else:
#        # Although we do not find the node itself
#        # we can find the nodes for the single element sets
#        # of the mvars in the node, since we have built the graph
#        # starting wiht them. E.g if node {A,B,C} is not part of
#        # the graph we know at least that the nodes {A}, {B} amd {C}
#        # are in the graph.
#        # For each of them we can compute the subgraph of
#        # spg that leads to it. If we compute the product of these
#        # subgraphs it will contain the desired node.
#
#        # fixme: mm 02-26-2020
#        # We could make this more efficient by looking for all the
#        # disjoint unions of the targetNode Mvars and compute
#        # the product of the graphs leading to the subsets.
#        # If the subsets are not one element sets we need fewer
#        # multiplications.
#        #prod_g=product_graph(*[target_subgraph(spg,frozenset({v})) for v in targetNode])
#        spg=product_graph(*[minimal_target_subgraph_for_single_var(spg,v) for v in targetNode])
#
#    path_dict=nx.shortest_path(spg,target=targetNode)
#    possible_startnodes=frozenset(path_dict.keys())
#    minimal_startnodes=frozenset([n for n in filter(filter_func,possible_startnodes)])
#    path_list_list=[
#        nx.shortest_path(spg,target=targetNode, source=sn)
#        for sn in minimal_startnodes
#    ]
#    connected_nodes=reduce(lambda acc, pl : acc.union(pl),path_list_list,frozenset({}))
#    return spg.subgraph(connected_nodes).copy()


def minimal_target_subgraph_for_single_var(spg: nx.Graph, targetVar: type):
    # all minimal starting points
    targetNode = frozenset({targetVar})
    start_nodes = minimal_startnodes_for_single_var(spg, targetVar)
    path_list_list = [
        nx.shortest_path(spg, target=targetNode, source=sn) for sn in start_nodes
    ]
    connected_nodes = reduce(
        lambda acc, pl: acc.union(pl), path_list_list, frozenset({})
    )
    return spg.subgraph(connected_nodes).copy()


def minimal_startnodes_for_node(spg: nx.Graph, targetNode: Set[type]) -> Set[Set]:
    if spg.has_node(targetNode):
        # The targetNode is already part of the spg.
        # (because it has been added in one of the update
        # steps as an argument set of a computer)
        # in which case we can simply return the startnodes
        # of the paths leading to it.
        path_dict = nx.shortest_path(spg, target=targetNode)
        possible_startnodes = frozenset(path_dict.keys())
    else:
        # Although we do not find the node itself
        # we can find the nodes for the single element sets
        # of the mvars in the node, since we have built the graph
        # starting wiht them. E.g if node {A,B,C} is not part of
        # the graph we know at least that the nodes {A}, {B} amd {C}
        # are in the graph.
        # For each of them we can compute the subgraph of
        # spg that leads to it. If we compute the product of these
        # subgraphs it will contain the desired node.

        # fixme: mm 02-26-2020
        # We could make this more efficient by looking for all the
        # disjoint unions of the targetNode Mvars and compute
        # the product of the graphs leading to the subsets.
        # If the subsets are not one element sets we need fewer
        # multiplications.
        # prod_g=product_graph(*[target_subgraph(spg,frozenset({v})) for v in targetNode])
        prod_g = product_graph(
            *[minimal_target_subgraph_for_single_var(spg, v) for v in targetNode]
        )
        prod_path_dict = nx.shortest_path(prod_g, target=targetNode)
        possible_startnodes = frozenset(prod_path_dict.keys())

    def filter_func(n):
        # remove every set that contains one of the variables we are looking for ...
        return not (any([(v in n) for v in targetNode]))

    minimal_startnodes = frozenset(
        [n for n in filter(filter_func, possible_startnodes)]
    )
    return minimal_startnodes
