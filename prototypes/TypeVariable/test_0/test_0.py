# assume a data assimilation procedure  like mcmc.
# One function it needs is param2res which receives a tuple of guessed parameters
# and produces the output by running the model which also depends on a set of 
# fixed parameters. This output can then be compared to observations.
# Now imagine that we have at least two different versions param2res_1 param2res_2
# which work for different tuples of parameters to be estimated.
# Imagine further that we have these two param2res versions for every model (here model_a and model_b)
# 
# To distinguish at compile time (or mypy check time) if we called the correct 
# version of param2res with with the right kind of estimated parameter we
# would implement different classes of estimated parameters.
# 
# We implement everything in a rudimentary way so that we can check the code with mypy# and run it but don't attemp to implement a real mcmc which would also need
# a costfuncion a filter and so on..
from typing import TypeVar, Callable 
import model_a.mod1
import model_a.mod2
import model_a.Obs
import model_b.mod1
import model_b.mod2
import model_b.Obs
# EP is a type VARIABLE that can represent either of the EP types for either models
EPT = TypeVar(
        'EPT',
        model_a.mod1.EP,
        model_a.mod2.EP,
        model_b.mod1.EP,
        model_b.mod2.EP
) 
ObsT = TypeVar(
        'ObsT',
        model_a.Obs.Obs,
        model_b.Obs.Obs,
) 


def mcmc(
    ep_0: EPT,
    obs: ObsT,
    param2res: Callable[[EPT], ObsT]
) -> EPT:
    # we don't do anything real
    # We run the model once (thereby making sure that it does run) 
    # but return the start parameters unimproved
    obs_mod = param2res(ep_0)
    return ep_0

print(mcmc(
    model_a.mod1.EP(a=3,b=5),
    model_a.mod1.param2res
))    
print(mcmc(
    model_a.mod2.EP(a=3,c=5,d=4),
    model_a.mod2.param2res
))    

# run also mypy on this file!
