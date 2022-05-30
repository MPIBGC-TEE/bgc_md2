from sympy import Symbol, diff, lambdify
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Callable

RT,u=map(Symbol,["RT","u"])
expr=RT*u

def map_cc(expr,syms,vals):
    RT,u=syms
    RTs,us=vals
    return np.array(
            [
                [
                    lambdify((RT,u),expr)(RTv,uv) for RTv in RTs
                ] for uv in us
            ]
    )

RTs=np.linspace(0,1,11)
us=np.linspace(0,2,21)
X,Y=np.meshgrid(RTs,us)
ax=plt.axes(projection='3d')
#ax.plot_surface(
ax.plot_wireframe(
        X,
        Y,
        map_cc(expr,(RT,u),(RTs,us)),
        color="blue"
)
prefix_fun=lambda sym :Symbol("Delta_{}".format(str(sym)))
Delta_RT,Delta_u=map(prefix_fun,(RT,u))

def taylor_delta(
        expr,
        syms: Tuple[Symbol],
        prefix: Callable
    ):
    return (
        expr
        +
        sum([
             diff(expr,v1)*prefix(v1)
            for v1 in syms
        ])
        +
        1/2*sum([
            diff(expr,v1,v2)*prefix(v1)*prefix(v2) 
            for v1 in (RT,u) for v2 in syms
        ])
    )

def delta(arr):
    return arr-arr[0]

Delta_RTs, Delta_us= map(delta,(RTs,us))

X,Y=np.meshgrid(Delta_RTs,Delta_us)
#ax.plot_wireframe(
ax.contour3D(
    X,
    Y,
    map_cc(
        taylor_delta(expr,(RT,u),prefix_fun).subs(
            {
                RT:RTs[0],
                u:us[0]
            }
        ),
        (Delta_RT,Delta_u),
        (Delta_RTs,Delta_us)),
    color="red"
)
plt.show()
