from pygraphviz.agraph import AGraph
from typing import Callable
from functools import reduce
from frozendict import frozendict
import networkx as nx
import numpy as np

from bgc_md2.resolve.graph_helpers import (
    update_generator,
    minimal_startnodes_for_node,
)
from .non_graph_helpers import (
    pretty_name,
)


def compset_2_string(compset, aliases=frozendict({})):
    return "{" + ",".join([pretty_name(c, aliases) for c in compset]) + "}"


def node_2_string(node, aliases=frozendict({})):
    return "{" + ",".join([pretty_name(v, aliases) for v in node]) + "}"


def nodes_2_string(node, aliases=frozendict({})):
    return "[ " + ",".join([node_2_string(n, aliases) for n in node]) + " ]"


def edge_2_string(e):
    return "(" + node_2_string(e[0]) + "," + node_2_string(e[1]) + ")"


def draw_sequence(fig, tups):
    axs = fig.subplots(len(tups), 1, sharex=True, sharey=True)
    for i, tup in enumerate(tups):
        g, title = tup
        axs[i].set(frame_on=True)
        # nx.draw(g,ax=axs[i],frame_on=True)
        draw_ComputerSetMultiDiGraph_matplotlib(axs[i], g)
        axs[i].set(title=title)


def draw_update_sequence(
    computers,
    max_it,
    fig,
    mvar_aliases=frozendict({}),
    computer_aliases=frozendict({})
):
    lg = [g for g in update_generator(computers, max_it=max_it)]
    nr = len(lg)
    fig.set_size_inches(20, 20 * nr)
    pos = nx.spring_layout(lg[-1])
    # layout alternatives
    # pos = nx.spring_layout(lg[-1], iterations=20)
    # pos = nx.circular_layout(lg[-1] )
    # pos = nx.kamada_kawai_layout (lg[-1])
    # pos = nx.planar_layout (lg[-1])
    # pos = nx.random_layout (lg[-1])
    # pos = nx.shell_layout (lg[-1])
    # pos = nx.spectral_layout (lg[-1])
    # pos = nx.spiral_layout (lg[-1])
    axs = fig.subplots(nr, 1, sharex=True, sharey=True)
    for i in range(nr):
        draw_ComputerSetMultiDiGraph_matplotlib(
            axs[i], lg[i], mvar_aliases, computer_aliases, pos=pos
        )


def draw_ComputerSetDiGraph_matplotlib(
    spsg: nx.DiGraph,
    ax,
    pos=None,
    **kwargs
):
    if pos is None:
        pos = nx.spring_layout(spsg)
        # pos = nx.circular_layout(spsg)

    nx.draw(
        spsg,
        labels={n: node_2_string(n) for n in spsg.nodes()},
        ax=ax,
        node_size=2000,
        node_shape="s",
        pos=pos,
        **kwargs
    )
    for e in spsg.edges():
        print(spsg.get_edge_data(*e))

    edge_labels = {
        e: compset_2_string(spsg.get_edge_data(*e)["computers"]) for e in spsg.edges()
    }
    nx.draw_networkx_edge_labels(spsg, ax=ax, edge_labels=edge_labels, pos=pos)


def draw_ComputerSetMultiDiGraph_matplotlib(
    ax,
    spsg,
    mvar_aliases=frozendict({}),
    computer_aliases=frozendict({}),
    targetNode=None,
    pos=None,
    **kwargs
):
    if pos is None:
        # layout alternatives
        # pos = nx.spring_layout(spsg)
        # pos = nx.circular_layout(spsg)
        # pos = nx.spring_layout(spsg, iterations=20)
        # pos = nx.circular_layout(spsg )
        pos = nx.kamada_kawai_layout(spsg)
        # pos = nx.planar_layout(spsg)
        # pos = nx.random_layout(spsg)
        # pos = nx.shell_layout(spsg)
        # pos = nx.spectral_layout(spsg)
        # pos = nx.spiral_layout(spsg)
    nx.draw(
        spsg,
        labels={n: node_2_string(n, mvar_aliases) for n in spsg.nodes()},
        ax=ax,
        node_size=1000,
        # node_shape='s',
        pos=pos,
        **kwargs
    )
    if targetNode is not None:
        res = minimal_startnodes_for_node(spsg, targetNode)
        nx.draw_networkx_nodes(
            spsg,
            pos,
            nodelist=[targetNode],
            node_color='r',
            alpha=0.8
        )
        nx.draw_networkx_nodes(
            spsg,
            pos,
            nodelist=list(res),
            node_color='r',
            alpha=0.4
        )

    ax.axis("On")
    # at the moment it is not possible to draw
    # more than one edge (egde_lables) between nodes
    # directly (no edgelabels for MultiDiGraphs)
    # therefore we draw only one line for all computersets
    # and assemble the label from the different edges
    def edgeDict_to_string(ed):
        target = "computers"
        comp_sets = [v[target] for v in ed.values() if target in v.keys()]
        comp_set_strings = [compset_2_string(cs, computer_aliases) for cs in comp_sets]
        res = "\n".join(comp_set_strings)
        # print(res)
        return res

    edge_labels = {e: edgeDict_to_string(spsg.get_edge_data(*e)) for e in spsg.edges()}

    nx.draw_networkx_edge_labels(spsg, ax=ax, edge_labels=edge_labels, pos=pos)
    mvar_aliases_inv = {val: key for key, val in mvar_aliases.items()}
    for i, k in enumerate(sorted(mvar_aliases_inv.keys())):
        ax.text(0, 0 - i / len(mvar_aliases), k + ": " + mvar_aliases_inv[k])


def AGraphComputerSetMultiDiGraph(spsg: nx.MultiDiGraph, cf: Callable) -> AGraph:
    A = nx.nx_agraph.to_agraph(spsg)
    A = AGraph(directed=True)
    A.node_attr["style"] = "filled"
    A.node_attr["shape"] = "rectangle"
    A.node_attr["fixedsize"] = "false"
    A.node_attr["fontcolor"] = "black"

    for node in spsg.nodes:
        A.add_node(node_2_string(node))
    edges = spsg.edges(data=True)
    for edge in edges:
        s, t, data_dict = edge
        computer_set = data_dict["computers"]
        ss, st = tuple(map(node_2_string, (s, t)))
        A.add_edge(ss, st)
        Ae = A.get_edge(ss, st)
        Ae.attr["label"] = "\n".join([c.__name__ for c in computer_set])
    return A


def AGraphComputerMultiDiGraph(spsg: nx.MultiDiGraph, cf: Callable) -> AGraph:
    A = nx.nx_agraph.to_agraph(spsg)
    A = AGraph(directed=True)
    A.node_attr["style"] = "filled"
    A.node_attr["shape"] = "rectangle"
    A.node_attr["fixedsize"] = "false"
    A.node_attr["fontcolor"] = "black"

    for node in spsg.nodes:
        A.add_node(node_2_string(node))
    edges = spsg.edges(data=True)
    for edge in edges:
        s, t, data_dict = edge
        computer_set = data_dict["computers"]
        for c in computer_set:
            ss, st = tuple(map(node_2_string, (s, t)))
            A.add_edge(ss, st)
            Ae = A.get_edge(ss, st)
            Ae.attr["color"] = cf(c)
            Ae.attr["fontcolor"] = cf(c)
            Ae.attr["label"] = c.__name__

    return A
