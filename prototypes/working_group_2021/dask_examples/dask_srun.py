from os import environ
from datetime import datetime
from dask_mpi import initialize
from distributed import Client

def run_test(client):
    res=client.gather(
        client.map(
            lambda x:x**3,
            [i for i in range(100)]
        )
    )
    return(res)

def main():
    # Work out from the environment how many threads to allocate
    num_threads = int(environ.get(
        'SLURM_CPUS_PER_TASK',
        environ.get('OMP_NUM_THREADS', 1)
    ))

    print("num_threads: ", num_threads)

    # Create the Dask workers
    initialize(interface='ib0', nthreads=num_threads)
    # Create the Dask object that will manage the communications
    client = Client()

    start = datetime.now()
    results = run_test(client=client)
    end = datetime.now()

    print(f"Time taken: {end - start}")
    print(results)

if __name__ == '__main__':
    main()