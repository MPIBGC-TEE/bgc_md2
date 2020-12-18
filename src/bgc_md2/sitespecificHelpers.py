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
import dask
from dask.distributed import Client
from dask.distributed import LocalCluster
from getpass import getuser
from subprocess import run, check_output


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

hosts={
        #b'matagorda':(32,1,'1GB'),
        #b'matagorda':(16,2,'4GB'),
        b'matagorda':(8,4,'4GB'),
        b'antakya':(96,1,'2.5GB')
}
def get_client():
    try:
        client = Client(user_specific_scheduler_addr)
    except IOError:
        my_cluster = getCluster()
        client = Client(my_cluster)#same as Client(addr)
    print(client)
    return client


def getCluster():
    n_w,t_p_w,mem=hosts[check_output(['hostname']).strip()]

    # The next line allows the workers to start
    # subprocesses which is important to implement
    # timeouts for computations that take too long
    dask.config.set({'distributed.worker.daemon': False})

    my_cluster = LocalCluster(
        scheduler_port=user_specific_scheduler_port,
        dashboard_address=user_specific_dashboard_addr,
        #n_workers=24,
        #threads_per_worker=4 
        #,
        n_workers=n_w,
        threads_per_worker=t_p_w, 
        #worker_kwargs = {
        #    'memory_limit': mem,
        #    'memory_target_fraction': 0.6,
        #    'memory_spill_fraction': 0.7,
        #    'memory_pause_fraction': 0.8,
        #    #'memory_terminate_fraction': False, # leads to errors if commented in
        #}
    )
    return my_cluster


if __name__ == "__main__":
    pass
