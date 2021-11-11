#import dask client function and os to read environmental variables
from dask.distributed import Client
from os import environ
from pathlib import Path
import pandas as pd
import json

#read in user environment, read in dask scheduler file from user dir
user = '/scratch/'+ environ.get('USER') + '/scheduler.json'
print(user)

client=Client(scheduler_file=user)

with Path('config.json').open(mode='r') as f:
    conf_dict=json.load(f) 

dataPath = Path(conf_dict['dataPath'])

#function to distribute across dask cluster
res=client.gather(
        client.map(
            lambda x:x**3,
            [i for i in range(1000000)]
        )
)

map_cubes = dataPath.joinpath('map_cubes.csv')
pd.DataFrame(res).to_csv(map_cubes,sep=',')

print("Data saved")
#print output
#print(res)