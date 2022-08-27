import sys
sys.path.insert(0,'..')
from pathlib import Path
import json 
import netCDF4 as nc
import numpy as np
from sympy import sin,integrate,var,lambdify


with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

dataPath = Path(conf_dict['dataPath'])
nds=nc.Dataset(dataPath.joinpath('YIBs_S0_Annual_cSoil.nc'),decode_times=False)
longs=nds.variables['longitude'].__array__()
lats=nds.variables['latitude'].__array__()
# We assume that lat,lon actually specifies the center of the gridcell
# and that it extends from lat-0.5 to lat+0.5 and long-0.5
for v in ('theta','phi','delta_theta','delta_phi','r'):
    var(v)

# we compute the are of a delta_phi * delta_theta patch 
# on the unit ball (which depends also on theta but not
# on phi)
A_sym=integrate(
    integrate(
        sin(theta),
        (theta,theta-delta_theta/2,theta+delta_theta/2)
     ),
     (phi,phi-delta_phi/2,phi+delta_phi/2)
)

A_fun=lambdify((theta,delta_theta,delta_phi),A_sym)

def grad2rad(alpha):
    return np.pi/180.0 * alpha

def area(lat):
    delta_phi=1
    delta_theta=1
    r= 6378.1370 #km 
    # translate to normal spherical coordinates
    theta=(90.0-lat)
    theta, delta_theta, delta_phi = map(
        grad2rad, 
        (
            theta,
            delta_theta,
            delta_phi,
        )
    )
    return A_fun(theta,delta_theta,delta_phi)*r
    
    
weights=np.array(list(map(area,lats)))    
cSoil = nds.variables['cSoil'].__array__()
sum_over_all_lons=cSoil.sum(axis=2)

#res = sum([sum_over_all_lons[:,i]*weights[i] for i in range(weights.shape[0])])
# equivalent but faster (since with direct array operation without list)
res = (sum_over_all_lons*weights).sum(axis=1)
# the result has one value for every time..
