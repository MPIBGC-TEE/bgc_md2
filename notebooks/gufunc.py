# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from getpass import getuser
import xarray as xr
import numpy as np
import dask
from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
from sympy import Symbol,var

# +
# To be able to work independently we
# set up independent dask clusters for different users 
# with different dashboards on different ports 
# You can use the same cluster from different notebooks though ()
#  
# Note that only the dashboard port (in addition of your jupyterlab or jupyter port) 
# has to be forwarded to your local maschine since the scheduler port  
# will be used for communication between jupyter and the dask cluster who 
# both run on the same machine (matagorda or antakya)
# The portnumber has just to be different for every user so that we can all kill our clusters.
# and ports are a resource that has to be shared between users of a machine. 
# (In a real hpc scenarion we would most likely use ONE cluster but this would not run under user
# priveledges)

port_dict = {
    'mm':(8689,8789),       # first is the port of the actual scheduler, second the port for the dashboard
    'hmetzler':(8690,8790), # change at will to a port you forward via ssh to your local machine
    'cs':(8691,8791)        # change at will
}
my_user_name = getuser()
print(my_user_name)
addr = 'localhost:'+str(port_dict[my_user_name][0])
try:
    Client(addr)
except IOError:
    my_cluster = LocalCluster(
        scheduler_port= port_dict[my_user_name][0],   
        dashboard_address='localhost:'+str(port_dict[my_user_name][1])
    )
    Client(my_cluster)#same as Client(addr)

    
# -


