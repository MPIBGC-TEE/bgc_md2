"""
The purpose is of this module is to collect different vizualisations that at least work.
Some are very unlikely to be used later on but are interesting as prototypes.
The tests meryly ensure that these prototypes work.
"""
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS, BASE_COLORS, TABLEAU_COLORS
import networkx as nx
from unittest import skip
from string import ascii_lowercase, ascii_uppercase
from testinfrastructure.InDirTest import InDirTest
from computers_and_mvars_for_testing import (
    # classes
    A,
    B,
    C,
    D,
    E,
    F,
    G,
    H,
    I
    # computers
    ,
    a_from_i,
    b_from_c_d,
    b_from_e_f,
    c_from_b,
    d_from_b,
    d_from_g_h,
    e_from_b,
    f_from_b,
)
from bgc_md2.resolve.graph_plotting import (
    AGraphComputerSetMultiDiGraph,
    AGraphComputerMultiDiGraph,
    draw_update_sequence,
    draw_ComputerSetMultiDiGraph_matplotlib,
    # ,draw_Graph_with_computers_svg
)
from bgc_md2.resolve.graph_helpers import sparse_powerset_graph
from bgc_md2.models.helpers import (
    # provided_mvars,
    # computable_mvars,
    # path_dict_to_single_mvar,
    # get_single_mvar_value,
    bgc_md2_computers,
    bgc_md2_computer_aliases,
    bgc_md2_mvar_aliases,
)

from bgc_md2.resolve.non_graph_helpers import (
    all_mvars,
    all_computers_for_mvar,
    input_mvars,
    pretty_name,
)

# for easier debugging in ipython
computers = frozenset(
    {
        a_from_i,
        b_from_c_d,
        b_from_e_f,
        c_from_b,
        d_from_b,
        d_from_g_h,
        e_from_b,
        f_from_b,
    }
)


def computers_between_mvars(allComputers):
    # note that this is not a graph we can query for connectivity
    allMvars = all_mvars(allComputers)
    G = nx.DiGraph()
    for v in allMvars:
        for c in all_computers_for_mvar(v, allComputers):
            ans = input_mvars(c)
            for an in ans:
                G.add_edge(pretty_name(an), pretty_name(v), computer=c)
    return G


def computer_color_func(allComputers):
    colordict = TABLEAU_COLORS
    color_names = [n for n in colordict.keys()]
    color_dict = {c.__name__: color_names[i] for i, c in enumerate(allComputers)}

    def cf(c):
        col = colordict[color_dict[c.__name__]]
        return col

    return cf


class TestGraphVizualization(InDirTest):
    def setUp(self):
        # fixme mm 02-03-2020:
        # assemble this set from a file of class definitions
        self.mvars = {
            A,
            B,
            C,
            D,
            E,
            F,
            G,
            H,
            I,
        }
        # fixme mm 02-03-2020:
        # assemble this set from a file of annotated functions in
        # special file
        self.computers = computers
        self.spsg = sparse_powerset_graph(self.computers)
        self.cf = computer_color_func(self.computers)

    def test_AGraphComputerSetMultiDiGraph(self):
        A = AGraphComputerSetMultiDiGraph(self.spsg, self.cf)
        A.draw("MultiDiGraph.svg", prog="circo")  # draw using circo

    def test_AGraphComputerMultiDiGraph(self):
        A = AGraphComputerMultiDiGraph(self.spsg, self.cf)
        A.draw("MultiDiGraph.svg", prog="circo")  # draw using circo

    def test_ComputerSetMultiGraph_matplotlib(self):
        # The direct visualization of networkx (using matplotlib)
        # is very rudimentary.
        fig = plt.figure(figsize=(20, 20))
        ax = fig.add_subplot(1, 1, 1)
        draw_ComputerSetMultiDiGraph_matplotlib(ax, self.spsg)
        fig.savefig("SetMultiGraph.pdf")

    def test_update_generator(self):
        fig = plt.figure()
        draw_update_sequence(self.computers, max_it=8, fig=fig)
        fig.savefig("c1.pdf")

    def test_ComputerSetMultiGraph_matplotlib_bgc_md2(self):
        # The direct visualization of networkx (using matplotlib)
        # is very rudimentary.
        spsg = sparse_powerset_graph(bgc_md2_computers())
        fig = plt.figure(figsize=(20, 20))
        ax = fig.add_subplot(1, 1, 1)
        draw_ComputerSetMultiDiGraph_matplotlib(
            ax,
            spsg,
            bgc_md2_mvar_aliases(), 
            bgc_md2_computer_aliases(),
        )
        fig.savefig("SetMultiGraph.pdf")

    def test_update_generator_bgc_md2(self):
        fig = plt.figure()
        draw_update_sequence(
            bgc_md2_computers(),
            8,
            fig,
            bgc_md2_mvar_aliases(), 
            bgc_md2_computer_aliases(),
        )
        fig.savefig("c1.pdf")

    @skip(
        """ very immature and nearly manual, but maybe one of the possibilities
        to make the connections clickable?"""
    )
    def test_draw_multigraph_plotly(self):
        draw_multigraph_plotly(self.mvars, self.computers)
