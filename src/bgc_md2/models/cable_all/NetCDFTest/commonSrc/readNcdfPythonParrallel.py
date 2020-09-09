# to run: mpirun -np 4 python mpi_example.py
import sys
from mpi4py import MPI
import numpy as np
from netCDF4 import Dataset

rank = MPI.COMM_WORLD.rank  # The process ID (integer 0-3 for 4-process run)
testfile = "simple_xy_par.nc"

nc = Dataset(testfile, parallel=True, comm=MPI.COMM_WORLD, info=MPI.Info())

if rank == 0:
    # obviously rank 0 can see the whole data set
    # whis is interesting but unnecessary and not scalable since
    # the data might not fit into the RAM of a single node
    for k, v in nc.variables.items():
        print("Variablename={0}".format(k))
        print("Value\n")
        print(v)
        print("Value as numpy.array \n")
        print(np.array(v))

    # The aim of the exercise is actually to read only the part relevant
    # for the computation the actual rank is responsible for

nc.close()
