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

RTs=np.linspace(1,3,11)
us=np.linspace(2,4,21)
X,Y=np.meshgrid(RTs,us)
ax=plt.axes(projection='3d')
#ax.plot_surface(
ax.plot_wireframe(
        X,
        Y,
        map_cc(expr,(RT,u),(RTs,us)),
        color="blue"
)
#prefix_fun=lambda sym :Symbol("Delta_{}".format(str(sym)))
suffix_fun=lambda sym :Symbol("{}_0".format(str(sym)))
#Delta_RT,Delta_u=map(prefix_fun,(RT,u))
RT_0,u_0=map(suffix_fun,(RT,u))

# +
def first_order(
        expr,
        syms: Tuple[Symbol],
        suffix: Callable
    ):
    d = {v:suffix(v) for v in syms}
    return (
        expr.subs(d)
        +
        sum([
             diff(expr,v).subs(d)*(v - suffix(v))
            for v in syms
        ])
    )

def taylor(
        expr,
        syms: Tuple[Symbol],
        suffix: Callable
    ):
    d = {v:suffix(v) for v in syms}
    print(d)
    return (
        first_order( expr, syms, suffix)+
        1/2*sum([
            diff(expr,v1,v2).subs(d)*(v1-suffix(v1))*(v2-suffix(v2)) 
            for v1 in (RT,u) for v2 in syms
        ])
    )


# -

#expr
taylor(expr,(RT,u),suffix=suffix_fun)
#first_order(expr,(RT,u),suffix=suffix_fun)


# def delta_l(arr):
#    return arr-arr[0]
#
# Delta_RTs_l, Delta_us_l= map(delta_l,(RTs,us))


# +
def center(arr):
    return (arr[-1]+arr[0])/2

#def delta_c(arr):
#    return arr-center(arr)
#Delta_RTs_c, Delta_us_c= map(delta_c,(RTs,us))


# -

taylor_l= taylor(expr,(RT,u),suffix_fun).subs(
    {
        RT_0:RTs[0],
        u_0:us[0]
    }
)
first_order_l= first_order(expr,(RT,u),suffix_fun).subs(
    {
        RT_0:RTs[0],
        u_0:us[0]
    }
)
taylor_c= (taylor(expr,(RT,u),suffix_fun)).subs(
    {
        RT_0:center(RTs),
        u_0:center(us)
    }
)
first_order_c= (first_order(expr,(RT,u),suffix_fun)).subs(
    {
        RT_0:center(RTs),
        u_0:center(us)
    }
)
first_order_l,first_order_c

# +
#ax.plot_wireframe(
X,Y=np.meshgrid(RTs,us)
fig=plt.figure()
ax=plt.axes(projection="3d")
ax.set_xlabel("RT")
ax.set_ylabel("u")

ax.plot_surface(
    X,
    Y,
    map_cc(
        expr,
        (RT,u),
        (RTs,us)
    ),
    color="red"
)
#ax.plot_wireframe(
#    X,
#    Y,
#    map_cc(
#        taylor_l,
#        (RT,u),
#        (RTs,us)
#    ),
#    color="black"
#)
#ax.plot_wireframe(
#    X,
#    Y,
#    map_cc(
#        taylor_c,
#        (RT,u),
#        (RTs,us)
#    ),
#    color="blue"
#)
#ax.plot_wireframe(
#    X,
#    Y,
#    map_cc(
#        first_order_l,
#        (RT,u),
#        (RTs,us)
#    ),
#    color="green"
#)
ax.plot_wireframe(
    X,
    Y,
    map_cc(
        first_order_c,
        (RT,u),
        (RTs,us)
    ),
    color="blue"
)
#plt.show()
from collections import namedtuple
p2d=namedtuple("p2d",["RT","u"])
corners=[
    p2d(RTs[0],us[0]),
    p2d(RTs[0],us[-1]),
    p2d(RTs[-1],us[-1]),
    p2d(RTs[-1],us[0]),
    p2d(center(RTs),center(us)),
]    
#ax=plt.axes(projection="3d")
#ax.set_xlabel("RT")
#ax.set_ylabel("u")
for c in corners:
    ax.scatter(
        xs=[c.RT],
        ys=[c.u],
        zs=[expr.subs({RT:c.RT,u: c.u})],
        s=200,
        color="blue"
    )
# -

comparisons=[(c1,c2) for c1 in corners for c2 in corners if c1 != c2]
comparisons


# Now we can use our linear (green) approximation to attribute differences

# +
def unaccounted(c1,c2):
    c1_dict,c2_dict=({RT: c.RT, u: c.u} for c in [c1,c2])
    
    z1,z2=(expr.subs(d) for d in [c1_dict, c2_dict])
    Delta_z=z1-z2
    
    z1_fo,z2_fo = (first_order_c.subs(d) for d in [c1_dict, c2_dict])
    Delta_z_fo = z1_fo-z2_fo
    
    #z1_tc,z2_tc=(taylor_c.subs(d) for d in [c1_dict, c2_dict])
    #Delta_z_tc=z1_tc-z2_tc
    return Delta_z-Delta_z_fo

unaccounted(corners[0],corners[1])

# comparisons with zero approximation error
aes=[ c for c in comparisons if unaccounted(*c) == 0 ]
aes

# +
#c1,z1,c2,z2,

# +
from sympy import symbols, simplify,factor
RT_1, u_1, RT_2, u_2=symbols("RT_1 u_1 RT_2 u_2")
d1={RT: RT_1, u: u_1}
d2={RT: RT_2, u: u_2}

taylor_c.subs(d1)-taylor_c.subs(d2) - first_order_c.subs(d1)-first_order_c.subs(d2)
# -


