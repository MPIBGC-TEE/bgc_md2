
from sympy import Symbol, diff, lambdify
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Callable
from sklearn.linear_model import LinearRegression

RT,u=map(Symbol,["RT","u"])
#expr=RT*u
#expr=5+3*RT + 4*u
expr=5+0.3*RT + .4*u + RT*u

#def map_cc(expr,syms,vals):
#    RT,u=syms
#    f = lambdify(syms,expr)
#    RTs,us=vals
#    return np.array(
#            [
#                [
#                    f(RTv,uv) for RTv in RTs
#                ] for uv in us
#            ]
#    )

def map_cols(expr,syms,vals):
    f = lambdify(syms,expr)
    return np.array(
        [
            f(vals[i,0],vals[i,1]) for i in range(vals.shape[0])
        ]
    )

def tensor_prod_1d(arr1,arr2):
    # transform to a data.frame like structure
    # by enumerating the tupels in a way that the last index 
    # changes fastest
    n1, n2 = len(arr1), len(arr2)
    fc = np.concatenate(
        [
            np.ones((n2,))*v for v in arr1
        ]
    )
    sc = np.concatenate(
        [
            arr2 for i in range(n1)
        ]
    )
    return fc,sc

def double(arr1,arr2):    
    fc,sc=tensor_prod_1d(arr1,arr2)
    return np.stack([fc,sc],axis=1)

def triple(arr1,arr2):
    fc,sc=tensor_prod_1d(arr1,arr2)
    tc=fc*sc
    return np.stack([fc,sc,tc],axis=1)

RTs=np.linspace(1,30,31)
us=np.linspace(2,40,21)
#Y=map_cc(expr,(RT,u),(RTs,us))
XS2D=double(RTs,us)
Y=map_cols(expr,(RT,u),XS2D)
lin_model = LinearRegression().fit(XS2D, Y.flatten())
print(lin_model.intercept_)
lin_app_1d = lin_model.predict(XS2D)

XS3D=triple(RTs,us)
bi_lin_model = LinearRegression().fit(XS3D, Y.flatten())
bi_lin_app_1d = bi_lin_model.predict(XS3D)


X1, X2 = np.meshgrid(RTs,us)
s=X1.shape
ax=plt.axes(projection='3d')
#ax.plot_surface(
ax.scatter(
        XS2D[:,0],
        XS2D[:,1],
        Y, 
        color="blue"
)
ax.scatter(
        XS2D[:,0],
        XS2D[:,1],
        lin_app_1d,
        color="red"
)
ax.scatter(
        XS2D[:,0],
        XS2D[:,1],
        bi_lin_app_1d,
        color="green"
)
plt.show()
