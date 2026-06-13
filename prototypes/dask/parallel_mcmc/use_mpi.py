from dask.distributed import Client
from os import environ
client = Client(scheduler_file='/home/'+environ.get('USER')+'/scheduler.json')

res=client.gather(
    client.map(
        lambda x:x**3,
        [i for i in range(100000)]
    )
)
print(res)
