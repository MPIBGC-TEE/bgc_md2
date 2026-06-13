from sympy import latex
from scipy.integrate import solve_ivp
import numpy as np
import CompartmentalSystems.helpers_reservoir as hr
def fig_save(
    module,
    func_name,
    path,
    **kwargs
):
    func = getattr(module,func_name)
    fig=func(**kwargs)
    fig.savefig(path.joinpath(f"{func_name}.pdf"))
    return fig



def mass_balance_equation(mvs):
    try:
        ip = mvs.get_InputTuple()
        cm = mvs.get_CompartmentalMatrix()
        sv = mvs.get_StateVariableTuple()
    except Exception as e:
        print(
            """The provided argument mvs does not allow all necessary computations to show the mass 
            balance equation:
             mvs.get_InputTuple()
             mvs.get_CompartmentalMatrix()
             mvs.get_StateVariableTuple()
            """
        )

    eq_str = rf'\frac{{d}}{{dt}}{latex(sv)} = {latex(ip)} + {latex(cm)} {latex(sv)}'
    return eq_str


###############################################################################
def compute_trajectories(mvs, par_dict, func_dict, X_0, t_eval):
    # create the numerical system from the symbolic one
    # this part uses the symbolic framework I developed together with Holger
    # Althouhg in this case we could have written the matrices by hand...
    state_vector = mvs.get_StateVariableTuple()
    t = mvs.get_TimeSymbol()

    B_func = hr.numerical_array_func(
        state_vector=state_vector,
        time_symbol=t,
        expr=mvs.get_CompartmentalMatrix(),
        parameter_dict=par_dict,
        func_dict=func_dict,
    )
    u_func = hr.numerical_array_func(
        state_vector=state_vector,
        time_symbol=t,
        expr=mvs.get_InputTuple(),
        parameter_dict=par_dict,
        func_dict=func_dict,
    )
    n = len(mvs.get_StateVariableTuple())

    # uncomment to test that we got the function right
    B_0 = B_func(0, X_0)  # in our linear case B does not realy depend on X
    print(B_0)
    # I_0=u_func(0,X_0) # in our linear case u does not realy depend on X
    # np.linalg.inv(B_0)

    def X_dot(t, X):
        B = B_func(t, X)
        # code for checking Yiqi's claim that the inverse of the (negative) Compartmental
        # matrix has only positive entries is true in this instance
        tau = np.linalg.inv(-B)
        if not (
            all(
                [
                    np.sign(tau[i, j]) >= 0
                    for i in range(2)
                    for j in range(2)
                ]
            )
        ):
            print(f"non positive inverse: min(tau) = {np.min(tau)}")
        return u_func(t, X).reshape(n) + B @ X

    # uncomment to check our right hand side
    # X_dot(0,X_0)
    # Here you could also just write your own rhs in python without having to use
    # the symbolic translation and the packages that provide it.

    # actually solve the system
    sol = solve_ivp(
        X_dot,
        t_span=[t_eval.min(), t_eval.max()],
        t_eval=t_eval,
        y0=X_0,
        method="Radau",  # use an implicit solver for safety (overkill ...)
    )

    Xs = sol.y.transpose()
    times = sol.t

    # compute some diagnostic variables as arrays (time indexed)
    nt = len(times)
    rnt = range(nt)

    # the derivative vectors
    X_dots = np.array([X_dot(times[i], Xs[i, :]) for i in rnt])

    Bs = np.array([B_func(times[i], Xs[i]) for i in rnt])
    # in Yiqis nomenclature the system is written
    # as d/dt X=I(t)-B(t)X
    # whereas others (like me) would write
    # it as d/dt X = I(t) +B(t) X... hence with my B_func Tau=(-B)^-1
    Taus = np.array(list(map(np.linalg.inv, -Bs)))

    # input vectors
    Is = np.array(
        [
            # more general u_func for nonlinear pool dependent input (Xs)
            # which in the linear case are just ignored
            u_func(times[i], Xs[i])
            for i in rnt
        ]
    )

    X_cs = np.array([Taus[i] @ Is[i].reshape(n) for i in rnt])

    # in the paper the Carbon Storage Potential X_p
    # which in the case of a one pool system has the same sign as the
    # derivative \dot{X}  which is the argument on which calling X_c an attractor
    # for X rests.
    # The plots will show that this property does not persist in higher dimensions
    X_ps = X_cs - Xs

    return Xs, X_cs, X_ps, X_dots, Is, Bs


