import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

import matplotlib.pyplot as plt

G=nx.DiGraph()
ifbs=["a","c","d"]
ofbs=["a","b"]
swss={"a","b","c"}
ss={"d","e","f"}
s=sss.union(ss)
# first add the inner node
G.add_nodes_from(ss)
#G.add_edges_from([("s_a","a")])
f=plt.figure(figsize=(10,10))
ax=f.add_subplot(1,1,1)
pos=nx.circular_layout(G,center=np.array([0,0]))
G.add_nodes_from(s)
#pos={
#    'a': np.array([1.00000000e+00, 2.45045699e-08]),
#    'b': np.array([0.49999998, 0.86602546]),
#    'c': np.array([-0.50000004,  0.8660254 ]),
#    'd': np.array([-9.99999970e-01, -6.29182054e-08]),
#    'e': np.array([-0.49999989, -0.86602541]),
#    'f': np.array([ 0.49999992, -0.86602541])
#}
pos_ext={
nx.draw_networkx_nodes(
        G=G,
        pos=pos,
        ax=ax,
        nodelist=nns,
        node_color="red"
)
#nx.draw_networkx_nodes(
#        ax=ax,G=G,
#        nodelist=[enns],
#        size=1,
#        node_color="blue"
#)
f.savefig("abstract_subsystem.pdf")

