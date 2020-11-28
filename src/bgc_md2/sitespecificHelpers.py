# These functions help keeping the dask clusters managable 
# They have nothing to do  with the framework and are just included for convinience


# To be able to work independently we set up independent dask clusters for
# different users with different dashboards on different ports You can use the
# same cluster from different notebooks though ()
#  
# Note that only the dashboard port (in addition of your jupyterlab or jupyter
# port) has to be forwarded to your local maschine since the scheduler port
# will be used for communication between jupyter and the dask cluster who both
# run on the same machine (matagorda or antakya) The portnumber has just to be
# different for every user so that we can all kill our clusters.  and ports are
# a resource that has to be shared between users of a machine.  (In a real hpc
# scenarion we would most likely use ONE cluster but this would not run under
# user priveledges)
from dask.distributed import Client
from dask.distributed import LocalCluster
from getpass import getuser
from subprocess import run

port_dict = {
    'mm': (8689, 8789),       # first is the port of the actual scheduler,
                              # second the port for the dashboard
    'hmetzler': (8690, 8790), # change at will to a port you forward via ssh to your local machine
    'cs': (8691, 8791)        # change at will
}

my_user_name = getuser()
print(my_user_name)
user_specific_scheduler_port = port_dict[my_user_name][0]
user_specific_scheduler_addr = 'localhost:'+str(user_specific_scheduler_port)
user_specific_dashboard_addr = 'localhost:'+str(port_dict[my_user_name][1])


def get_client():
    try:
        client = Client(user_specific_scheduler_addr)
    except IOError:
        my_cluster = getCluster()
        client = Client(my_cluster)#same as Client(addr)
    print(client)
    return client


def getCluster():
    my_cluster = LocalCluster(
        scheduler_port=user_specific_scheduler_port,   
        dashboard_address=user_specific_dashboard_addr
        #,
        #n_workers=24,
        #threads_per_worker=4 
        #,
        #n_workers=48,
        #threads_per_worker=1 
    )
    return my_cluster



