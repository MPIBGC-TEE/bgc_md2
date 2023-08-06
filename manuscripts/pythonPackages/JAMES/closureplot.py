import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import bgc_md2
from importlib import import_module

def _format_axes(ax):
    """Visualization options for the 3D axes."""
    # Turn gridlines off
    ax.grid(False)
    # Suppress tick labels
    for dim in (ax.xaxis, ax.yaxis):#, ax.zaxis):
        dim.set_ticks([])
    # Set axes labels
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    #ax.set_zlabel("z")
#
#
#plt.show()
mf="jon_yib"
model_mod=f'bgc_md2.models.{mf}'
##
mvs = import_module(f"{model_mod}.source").mvs
seq = mvs.closure_graph_sequence()
final_g = seq[-1]
given_nodes, computer_nodes, computable_nodes = final_g._nodes_tup()

pos = nx.bipartite_layout(final_g,given_nodes.union(computable_nodes))
n_p = 2 
fig = plt.figure(figsize=10*np.array((n_p,1)))
axs = fig.subplots(1, n_p)
# ax = fig.add_subplot(111)#, projection="3d")
for i, g in enumerate(seq[1:n_p]):
    g.draw_matplotlib(axs[i],pos=pos)
    #_format_axes(ax)
final_g.draw_matplotlib(axs[n_p-1],pos=pos)
fig.tight_layout()

from IPython import embed; embed()
# plt.show()
fig.savefig("closure.pdf")
