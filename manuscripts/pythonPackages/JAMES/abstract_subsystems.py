import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

import matplotlib.pyplot as plt

ifbs = ["a","c","d"]
ofbs = ["a","b"]
swoss = {"a","b","c"}
ss = {"d","e","f"}
s = swoss.union(ss)
# first add the inner node
G=nx.DiGraph()
G.add_nodes_from(ss)
f=plt.figure(figsize=(10,10))
ax=f.add_subplot(1,1,1)
#pos_ss=nx.circular_layout(G,center=np.array([0,0]))
pos_ss={
 'd': np.array([-1.00000000e+00, 1.98682151e-08]),
 'e': np.array([0.50000007,  0.86602542]),
 'f': np.array([0.49999993, -0.86602544]),
}
def project_left(arr):
    print(arr.shape)
    x,y=arr
    return np.array([-5,y])

G.add_nodes_from(swoss)
pos_swoss = {
    "a":2*pos_ss["d"],
    "b":2*pos_ss["e"],
    "c":2*pos_ss["f"],
}        
pos_s={**pos_ss,**pos_swoss}
pos_ext= {
    "to_a":project_left(pos_s["a"]),
    "to_b":project_left(pos_s["b"]),
    "to_e":project_left(pos_ss["e"])
    #"a_to":2*pos_ss["g"],
}        
#pos_ext = {k:3*v for k,v in pos_ss.items()}
in_ext=[
    ("to_a","a"),
    ("to_b","b"),
    ("to_e","e")
]
out_ext=[
    ("a_to","a"),
    ("b_to","b"),
]
pos={**pos_s,**pos_ext}
G.add_edges_from(in_ext)

nx.draw_networkx_nodes(
        G=G,
        pos=pos_ss,
        ax=ax,
        nodelist=ss,
        node_color="black"
)
nx.draw_networkx_labels(
        G=G,
        pos=pos_ss,
        ax=ax,
        labels={k:k for k in ss},
        font_color="white"
)
nx.draw_networkx_nodes(
        G=G,
        pos=pos_swoss,
        ax=ax,
        nodelist=swoss,
        node_color="black",
        alpha=0.2
)
nx.draw_networkx_labels(
        G=G,
        pos=pos_swoss,
        ax=ax,
        labels={k:k for k in swoss},
        font_color="white"
)
nx.draw_networkx_edges(
        G=G,
        pos=pos,
        ax=ax,
        edgelist=in_ext,
        edge_color="blue",
        alpha=0.3
)
circ_s=plt.Circle((0,0),2.1,color='red',alpha=0.1)
circ_ss=plt.Circle((0,0),1.1,color='red',alpha=0.1)
ax.add_patch(circ_s)
ax.add_patch(circ_ss)
f.savefig("abstract_subsystem.pdf")

