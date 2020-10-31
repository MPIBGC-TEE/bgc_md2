import unittest
import dask.array
import numpy as np
from dask.distributed import Client
from dask.distributed import LocalCluster
class TestDaskBatch(unittest.TestCase):

    def setUp(self):
        self.cluster=LocalCluster()

    def testComputationResult(self):
        c=Client(self.cluster)
        def my_main():
            x=dask.array.stack(
                [   
                    i*np.ones((100,100))
                    for i in range(2000)
                ],
                2
            )
            return 2*x
        futures=c.submit(my_main)
        x=futures.result()

        print(x[1,1,2])
        print(x[1,1,2].compute())
