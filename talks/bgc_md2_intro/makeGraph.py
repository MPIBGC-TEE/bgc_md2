import matplotlib.pyplot as plt
from copy import  deepcopy
from ComputabilityGraphs.fast_graph_helpers import (
        add_combi_arg_set_graph,
        add_combis_arg_set_graphs_to_decomp,
        add_all_arg_set_graphs_to_decomp,
)
from ComputabilityGraphs.graph_helpers import (
        minimal_startnodes_for_single_var,
)
import ComputabilityGraphs.fast_graph_helpers as fgh
from ComputabilityGraphs.graph_plotting import (
    draw_ComputerSetMultiDiGraph_matplotlib
)
from ComputabilityGraphs.FastGraph import FastGraph
import ComputabilityGraphs.helpers as h
from ComputabilityGraphs.Node import Node
from ComputabilityGraphs.Decomposition import Decomposition
from ComputabilityGraphs.ComputerSet import ComputerSet
from ComputabilityGraphs.ComputerSetSet import ComputerSetSet

from testComputers import (
    A, A2, A3, B, B1, B0, C, D, E, F, G, H, I,
    a_from_i,
    b_from_c_d,
    b_from_e_f,
    c_from_b,
    d_from_b,
    d_from_g_h,
    a3_from_a2,
    b1_from_b0,
    a3_from_b0,
    b1_from_a2,
    e_from_b,
    f_from_b,
)
computers = frozenset(
    [
        a_from_i,
        b_from_c_d,
        b_from_e_f,
        c_from_b,
        d_from_b,
        d_from_g_h,
        e_from_b,
        f_from_b,
    ]
)

### target B
fspsg_B = fgh.project_to_multiDiGraph(
    fgh.fast_graph(
        root_type=B,
        cs=computers,
        given=frozenset()
    )
)
res_B = minimal_startnodes_for_single_var(fspsg_B, B)

# plot
fig = plt.figure(figsize=(8, 8))
axs = fig.subplots(1, 1)
draw_ComputerSetMultiDiGraph_matplotlib(
    axs,
    fspsg_B,
    targetNode=Node({B})
)
fig.savefig("StartNodes.pdf")
