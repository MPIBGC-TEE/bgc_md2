import numpy as np

def compute_start_values_14C(times, Bs_12C, us_12C, Fa_func, method):
    ## compute 14C start_values

    # prepare Bs_14C, us_14C in case they are needed by 'method'
    lamda = 0.0001209681 # daily decay_rate 
    len_Bs = Bs_12C.shape[0]
    nr_pools = Bs_12C.shape[1]

    Bs_14C = np.array(
        [
            np.matmul(
                Bs_12C[k],
                np.exp(-lamda*365.25/12)*np.identity(nr_pools)
            ) 
            for k in range(len_Bs)
        ]
    )
    with np.errstate(divide='ignore'):
        us_14C = us_12C * Fa_func(times[:-1]).reshape(-1,1)
    us_14C =  np.nan_to_num(us_14C, posinf=0)


    if 'method' == 'D1':
        # assume system at t0 in eq., use matrices and vectors
        B0_14C = Bs_14C[0]
        start_values_14C = np.linalg.solve(
            (np.eye(nr_pools)-B0_14C),
            us_14C[0]
        )

    if method == 'D2':
        # version from whiteboard with Paul in Berkeley, pool-wise
        # based on meas through time
        # problem for pools with no external input
        ks = np.diag(np.mean(Bs_12, axis=0))
        start_values_14C = np.mean(us_14C, axis=0)\
                                /(1-ks*np.exp(-lamda*365.25/12))

    if method == 'D3':
        # mean of 14C Bs and mean of 14C us, use as eq.
        B0_14C = np.mean(Bs_14C, axis=0)
        u0_14C = np.mean(us_14C, axis=0)
        start_values_14C = np.linalg.solve(
            (np.eye(nr_pools)-B0_14C),
            u0_14C
        )

    if method == 'C3':
        # mean of 14C Bs and mean of 14C us
        B0_14C = np.mean(Bs_14C, axis=0)
        u0_14C = np.mean(us_14C, axis=0)

        start_values_14C = np.linalg.solve(-B0_14C, u0_14C)

    return start_values_14C

