import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from functools import reduce
from typing import List

ifbs = ["a","c","e"]
ofbs = ["f","b"]
intfbs=[
    ("a","d"),
    ("a","e"),
    ("d","f"),
    ("e","b"),
    ("e","f"),
    ("c","f"),
]
swoss = {"a","b","c"}
ss = {"d","e","f"}
s = swoss.union(ss)
# first add the inner node
G=nx.DiGraph()
G.add_nodes_from(ss)
f=plt.figure(figsize=(10,10))
axs=f.subplots(2,1)
ax=axs[0]
#pos_ss=nx.circular_layout(G,center=np.array([0,0]))
leftb = -3
rightb = 3
pos_ss={
 'd': np.array([-1.00000000e+00, 1.98682151e-08]),
 'e': np.array([0.50000007,  0.86602542]),
 'f': np.array([0.49999993, -0.86602544]),
}
def project_in(arr):
    x,y=arr
    return np.array([leftb, y])

def project_out(arr):
    x,y=arr
    return np.array([rightb, y])

def ext_in_name(n:str):
    return f"to_{n}"
    
def ext_out_name(n:str):
    return f"{n}_to"

def flatten(ll: List[List])->List:
    return reduce(
        lambda acc, el: acc + el,
        ll
    )

def plot_nodes_with_labels(
        G,
        ax,
        pos,
        nodelist,
        labels,
        alpha=1
    ):
    nx.draw_networkx_nodes(
            G=G,
            #pos=pos_ss,
            pos=pos,
            ax=ax,
            nodelist=nodelist,
            node_color="black",
            alpha=alpha
    )
    nx.draw_networkx_labels(
            G=G,
            #pos=pos_ss,
            pos=pos,
            ax=ax,
            labels=labels,
            font_color="white",
            alpha=1
            #alpha=alpha
    )

def plot(
        G,
	    ax,
	    pos,
        s,
	    ss,
        circs,
        in_flux_color="blue",
        out_flux_color="red",
        internal_flux_color="black",
    ):
    ax.set_aspect(1)
    ax.set_xlim(-3,3)
    ax.set_ylim(-3,3)
    # draw ss nodes
    plot_nodes_with_labels(G,ax,pos,ss,labels={k:k for k in ss})
    # draw swoss nodes
    swoss=s.difference(set(ss))
    plot_nodes_with_labels(
        G,
        ax,
        pos,
        swoss,
        labels={
            k:k for k in swoss
        },
        alpha=0.2
    )
    # draw ss influxes 
    ss_in_es = [
        (src,t) for (src,t) in flatten(
            [
                list(G.in_edges(n)) 
                for n in ss 
            ]
        ) 
        if src  not in ss
    ]
    nx.draw_networkx_edges(
            G=G,
            pos=pos,
            ax=ax,
            edgelist=ss_in_es,
            edge_color=in_flux_color,
            alpha=1
    )
    # draw ss outfluxes 
    ss_out_es = [
        (src,t) for (src,t) in flatten(
            [
                list(G.out_edges(n)) 
                for n in ss 
            ]
        ) 
        if t not in ss
    ]
    nx.draw_networkx_edges(
            G=G,
            pos=pos,
            ax=ax,
            edgelist=ss_out_es,
            edge_color=out_flux_color,
            alpha=1
    )
    # draw internal ss edges
    ss_internal_es=[
        (src,t) for (src,t) in G.edges() 
        if ((src in ss) and (t in ss))
    ]
    nx.draw_networkx_edges(
            G=G,
            pos=pos,
            ax=ax,
            edgelist=ss_internal_es,
            edge_color=internal_flux_color,
            alpha=1
    )
    # draw remaining edges
    rest_es = set(G.edges()).difference( ss_in_es + ss_out_es + ss_internal_es)
    nx.draw_networkx_edges(
            G=G,
            pos=pos,
            ax=ax,
            edgelist=rest_es,
            edge_color=internal_flux_color,
            alpha=0.3
    )
    for circ in circs:
        ax.add_patch(circ)

G.add_nodes_from(swoss)
pos_swoss = {
    "a":2*pos_ss["d"],
    "b":2*pos_ss["e"],
    "c":2*pos_ss["f"],
}        
pos_s={**pos_ss,**pos_swoss}
s={k for (k,v) in pos_s.items()}

pos_ext_in = {ext_in_name(k):project_in(pos_s[k]) for k in ifbs}
pos_ext_out = {ext_out_name(k):project_out(pos_s[k]) for k in ofbs}
pos={**pos_s, **pos_ext_in, **pos_ext_out}

G.add_edges_from([(ext_in_name(n),n) for n in ifbs])
G.add_edges_from([(n,ext_out_name(n)) for n in ofbs])
G.add_edges_from(intfbs)


black_alpha = 0.05
red_alpha = 0.1
r_ss = 1.3
r_s = 2.3
circ_s_b = plt.Circle((0,0),r_s,color='black',alpha=black_alpha)
circ_s = plt.Circle((0,0),r_s,color='red',alpha=red_alpha)
plot(
    G,
	axs[0],
	pos,
    s,
	s,
	[circ_s_b, circ_s]
)

circ_s_b=plt.Circle((0,0),r_s,color='black',alpha=black_alpha)
circ_ss=plt.Circle((0,0), r_ss, color='red', alpha=red_alpha)
plot(
    G,
	axs[1],
	pos,
    s,
	ss,
	[circ_s_b, circ_ss]
)
f.savefig("abstract_subsystems.pdf")

