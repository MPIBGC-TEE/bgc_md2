
import dask.array as da
from dask.distributed import Client, LocalCluster
import numpy as np
from scipy.optimize import least_squares

def euklidian(x):
    return np.dot(x,x)

def make_single_func(compute_loss,bounds):
    def single_func(x0):
        res=least_squares(
            fun=compute_loss,
            x0=x0.reshape(-1),
            method="trf",
            bounds=bounds
        )
        if res.success:
            return  np.concatenate(
                [
                    res.x.reshape(1,-1), 
                    np.array(res.cost).reshape(1,1)
                ],
                axis=1
            ) 
        else: 
            return  np.concatenate(
                [
                    x0.reshape(1,-1), 
                    np.array(np.NAN).reshape(1,1)
                ],
                axis=1
            ) 
    return single_func

if __name__ == '__main__':
    cluster = LocalCluster(n_workers=8)
    client = cluster.get_client()
    guesses=np.random.uniform(-1,+1,(8,20))
    nc=guesses.shape[1]
    gs=da.from_array(guesses,chunks=(1,nc))
    single_func=make_single_func(
        compute_loss=euklidian,
        bounds=(
            -np.ones((nc,)),
            np.ones((nc,))
        )
    )
    res = gs.map_blocks(
        single_func,
        meta=np.zeros((1,nc))
    ).compute()
    print(res)
