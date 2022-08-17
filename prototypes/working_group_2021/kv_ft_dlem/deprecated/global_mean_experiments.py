# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import json
from pathlib import Path
import netCDF4 as nc
import sys
import numpy as np
sys.path.insert(0,'..')
import general_helpers as gh
with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 
# %load_ext autoreload    
# %autoreload 2

# +
dataPath=Path(conf_dict["dataPath"])
import model_specific_helpers_2 as msh
def get_var(vn):
    fn="ISAM_S2_{}.nc".format(vn)
    path = dataPath.joinpath(fn)
    ds = nc.Dataset(str(path))
    #scale fluxes vs pools
    return ds.variables[vn]

    

names=msh.Observables._fields + msh.Drivers._fields
#names=["rh"]
    #prin
# -



svs,dvs=msh.get_globalmean_vars(dataPath=dataPath)

# +
#invalid_pixels=[
#    gh.get_nan_pixels(get_var(name))
#    for name in names
#]
# -

masks=[
    gh.get_nan_pixel_mask(get_var(name))
    for name in names
]

m=np.zeros((3,3))
#x=1 if m else 0
m=False
isinstance(m,bool)


from functools import reduce    
combined_mask= reduce(lambda acc,m: np.logical_or(acc,m),masks)
combined_mask


np.array(tuple(map(lambda x:x**2,np.linspace(0,10,11))))

fn="DLEM_S2_cSoil.nc"
path=Path(conf_dict["dataPath"])
ds=nc.Dataset(str(path.joinpath(fn)))
ds.variables['lat']
lats=ds.variables['lat'].__array__()
lons=ds.variables['lon'].__array__()

if get_var("cSoil")[:,0,0].mask

gh.global_mean_var(lats,lons,combined_mask,get_var("rh"))

# +
arr=var.__array__()
N_t,N_lat,N_lon = var.shape
cs=30
def f(I_lat,I_lon):
    n_lat = min(cs,N_lat-I_lat)
    n_lon = min(cs,N_lon-I_lon)
    chunk = var[:,I_lat:I_lat+n_lat,I_lon:I_lon+n_lon]
    return [
        (I_lat + i_lat, I_lon + i_lon) 
        for i_lat in range(n_lat) 
        for i_lon in range(n_lon)
        if np.isnan(
            chunk[:,i_lat,i_lon]
        ).any()
    ]
                
l=[
    f(I_lat,I_lon) 
    for I_lat in range(0,N_lat,cs)
    for I_lon in range(0,N_lon,cs)
]
# -

from functools import reduce
reduce(lambda x,y:x+y,l)

lats=ds.variables["lat"].__array__()
lons=ds.variables["lon"].__array__()
lats,lons

[i for i in range(0,11,2)]

gh.global_mean(lats,lons,arr)

for lat in range(len(lats)):
    for lon in range(len(lons)):
        if not arr.mask[0,lat,lon]:
            tl=arr[:,lat,lon]
            hn=np.isnan(tl).any()
            if hn: 
                print(lat,lon,hn) 
                #print(tl)


var[:,55,659]
#var[:,71,579]
#var[:,83,175]

# +
weight_mask=arr.mask[0,:,:] if  arr.mask.any() else False
delta_lat=(lats.max()- lats.min())/(len(lats)-1)
delta_lon=(lons.max() -lons.min())/(len(lons)-1)

pixel_area = gh.make_pixel_area_on_unit_spehre(delta_lat, delta_lon)

weight_mat= np.ma.array(
    np.array(
        [
                [   
                    pixel_area(lats[lat_ind]) 
                    for lon_ind in range(len(lons))    
                ]
            for lat_ind in range(len(lats))    
        ]
    ),
    mask = weight_mask 
)

# -

n_t=2
n_lats=3
n_lons=4
ref=np.zeros((n_t,n_lats,n_lons))
ref[0,2,3]=np.nan
ref[1,1,3]=np.nan
ref


# +
def get_nan_pixels_from_var(
    var 
    ):
    # we consider the  dimensions to be time,lat,lon
    n_t,n_lats,n_lons=var.shape
    return tuple(
        (
            (i_lat,i_lon) 
            for i_lat in range(n_lats) 
            for i_lon in range(n_lons)
            if np.isnan(var[:,i_lat,i_lon]).any()
        )
    )

#get_nan_pixels(ref)


# -

weight_mat.shape,arr.shape

w_arr=weight_mat*arr
w_arr.sum(-1)

weight_mat.sum()

masks=[arr.mask[i,:,:] for i in  range(arr.shape[0])] 

np.all(
    [
        (m==masks[0]).all()
        for m in masks
    ] 
)


lats.shape


lons.shape

arr.shape

(np.stack([weight_mat*arr[i,:,:] for i in range(arr.shape[0])])==weight_mat*arr).all()


for lat in range(len(lats)):
    for lon in range(len(lons)):
        if not arr.mask[0,lat,lon]:
            #print(lat,lon,arr[0,lat,lon],w_arr[lat,lon]) 
            print(lat,lon,w_arr[0,lat,lon]) 

weight_mat.min()

weight_mat.sum()

im=  np.array(
    tuple(
        (
            tuple((   
                pixel_area(lats[lat_ind]) 
                for lon_ind in range(len(lons))    
            ))
        for lat_ind in range(len(lats))
        )
    )
)
im2=np.zeros(arr.shape[1:])
weight_mask

weight_mat.sum()

t=np.array(((1,2,3),(4,5,6)))


# +



# -

gh.get_nan_pixels(get_var(name))  
for name in names

# import general_helpers as gh


get_var('rh').__array__().mask


