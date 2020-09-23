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
from dask.distributed import Client
from dask.distributed import LocalCluster

# +
# To be able to work independently.
# set up independent dask clusters for different users 
# with different dashboards on different ports 
#  
# The resulting port number has to be forwarded via ssh to look at the dashboard 
#(Alternativly we probably could use a common cluster though) 

port_dict = {
    'mm':8789,
    'hmetzler':8790, # change at will
    'cs':8791        # change at will
}
my_user_name = getuser()
print(my_user_name)

my_port = port_dict[my_user_name]
print(my_port)

my_cluster = LocalCluster(dashboard_address='localhost:'+str(my_port))

# -

Client(my_cluster)
