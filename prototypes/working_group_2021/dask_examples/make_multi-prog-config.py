from os import environ
from sys import argv 
from pathlib import Path
_, user, ntasks, config_filename, pyscript = argv
config="""0	dask-scheduler --scheduler-file /scratch/{user}/scheduler.json
1-{w} 	dask-worker --scheduler-file /scratch/{user}/scheduler.json
{n}	python {s}""".format(
    user=user,
    w=int(ntasks)-2,
    n=int(ntasks)-1,
    s=pyscript
)
with Path(config_filename).open("w") as f:
    f.write(config)

