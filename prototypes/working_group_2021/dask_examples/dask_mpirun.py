from os import environ
from datetime import datetime
from dask_mpi import initialize
initialize()
from distributed import Client
client=Client()

def run_test(client):
    res=client.gather(
        client.map(
            lambda x:x**3,
            [i for i in range(1000000)]
        )
    )
    return(res)

def main():
    print("start of main function")

    start = datetime.now()
    results = run_test(client=client)
    end = datetime.now()

    print(f"Time taken: {end - start}")
    print(results)

if __name__ == '__main__':
    main()