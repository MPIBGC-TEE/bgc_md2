import dask.array as da
import numpy as np
## output has the same size
nc=3
x = da.ones((16,nc), chunks=(1,nc))
#y = x.map_blocks(lambda x: x * 2).compute()
#print(y,y.shape)
# output is bigger in one of the dimensions the same size
y2 = x.map_blocks(lambda x: np.concatenate(
        [
            x *2,
            np.array([5]).reshape(1,1)
        ],
        axis=1
    ),
    meta=np.zeros((1,nc))
).compute()
print(y2,y2.shape)
