from sympy import symbols,var, exp, srepr
import networkx as nx
var("x,y,z")
expr = exp(x**2) + x*y
par_dict={
    x: 1000,
    y: 2000
}
print(srepr(expr))

# build a subgraph a
sg=nx.DiGraph()
sub_expr=expr.args[0]
sf=sub_expr.func
sg.add_node(sf)
sg.add_edges_from(
    [(sf,a) for a in sub_expr.atoms()]
)

def subg(subex):
    if isinstance(subex,Symbol):
        return subex
    else:
        subg=nx.DiGraph()
        f=subex.func
        subg.add_edges_from(
            [(f,subg(a)) for a in subg.atoms()]
        )
        return subg
