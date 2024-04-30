# This is a minimal example how to create a numeric rhs from a symbolic expression and solve it
# The central tool is the funcion lambdify of the sympy package, that translates sympy expressions into python functions.

# Note that this is not necessarry to solve ODEs we could write the rhs function ourselfes..
# It is just more versatile if we want to work on the symbolic side of things too..
# e.g. by checking if our guess for an analytic solution is true or to display the equation for the rhs which is easier to check
# than a piece of code...

from sympy import symbols, Function, Tuple, lambdify, exp, Abs, Matrix, diff
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

t,x_a,x_b,k,l=symbols("t,x_a,x_b,k,l")
state_var_tuple=Tuple(
        x_a,
        x_b
)    
rhs_sym=Tuple(
        k*x_a,
        l*x_b
)
start_value_dict={
        x_a:1,
        x_b:2
}
y0=state_var_tuple.subs(start_value_dict)
par_dict = {
        k: -5,
        l: -2
}        

# note that in the call to lamdify we specify the argument tuple including t since solve ivp wants a function 
# that has time as the first and the iterable y (in our case the tuple (x_a, x_b) as the second argument

sol2=solve_ivp(
    fun=lambdify((t,state_var_tuple),rhs_sym.subs(par_dict),"numpy"),
    t_span=(0,5),
    y0=y0,
    dense_output=True
).sol

plot_times=np.linspace(sol2.t_min,sol2.t_max,100)
fig=plt.figure()
ax=fig.add_subplot()
values=sol2(plot_times)
for i in range(values.shape[0]):
    ax.plot(
        plot_times,
        values[i,:],
        label=str(state_var_tuple[i])
    ) 
ax.legend()

