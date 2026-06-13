import matplotlib.pyplot as plt
import numpy as np
from collections import namedtuple

params = namedtuple('params',['A_1', 'omega_1', 'phi_1', 'A_2', 'omega_2','phi_2'])

def forward_solution(
        times: np.ndarray,
        A_1: float,
        omega_1: float,
        phi_1: float,
        A_2: float,
        omega_2: float,
        phi_2: float,
    ) -> np.ndarray:
    # The real forward model will be implemented by a difference equation
    # x_{i+1} = f(x_i,i) but to show how the parameter estimation works it is
    # enough to consider it to be a function of the parameters and returns a
    # solution on some time grid here we represent it by a 2 sin functions for
    # the contents of the two pools
    x1 = A_1 * np.sin(omega_1*times + phi_1)
    x2 = A_2 * np.sin(omega_2*times + phi_2)
    return np.stack([x1, x2], axis=1)


# run the forward model to make sure that it works
nt = 41
times = np.linspace(0,10,nt)
p_test = params(
    A_1=10,
    omega_1=2*np.pi/4.0,
    phi_1=0,
    A_2=20,
    omega_2=2*np.pi/4.1,
    phi_2=np.pi/4
 )._asdict()

sol_test = forward_solution(times, **p_test)

# now we use the forward solution to produce synthetic data  
# with a bit of random noise
err_1 = 4*np.random.rand(nt)
err_2 = 2*np.random.rand(nt)
err = np.stack([err_1,err_2],axis=1)
data = forward_solution(times, **p_test) + err

fig = plt.figure()
colors = ['blue','green']
ax = fig.add_subplot()
for i in (0,1):
    ax.scatter(times, data[:, i],color=colors[i])
    ax.plot( times, sol_test[:,i],color=colors[i])

fig.savefig('test.pdf')

