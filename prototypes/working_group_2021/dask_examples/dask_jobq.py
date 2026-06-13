from os import environ
from datetime import datetime
from dask_mpi import initialize
from distributed import Client
from dask_jobqueue import SLURMCluster

def run_test(client):
    res=client.gather(
        client.map(
            lambda x:x**3,
            [i for i in range(100)]
        )
    )
    return(res)

def main():
    cluster = SLURMCluster(cores=2,
                     memory="10GB",
                     walltime='00:05:00')
    cluster.scale(5)  # Start 100 workers in 100 jobs that match the description above
    client = Client(cluster)    # Connect to that cluster    client = Client()
    print(client)
    start = datetime.now()
    results = run_test(client=client)
    end = datetime.now()

    print(f"Time taken: {end - start}")
    print(results)

if __name__ == '__main__':
    main()