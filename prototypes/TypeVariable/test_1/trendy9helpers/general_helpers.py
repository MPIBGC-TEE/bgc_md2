from typing import TypeVar, Callable
import numpy as np
'''
create mixin types to distinguish EstimatedParameter Tuples 
to define a TypeVariable for all its subtypes
'''


class EP():
    pass


EPT = TypeVar('EPT', bound=EP) 


# same for Observations
class Obs():
    pass


ObsT = TypeVar('ObsT', bound=Obs) 


def mcmc(
    ep_0: EPT,
    obs: ObsT,
    param2res: Callable[[EPT], ObsT]
) -> EPT:
    # we don't do anything real
    # We run the model once (thereby making sure that it does run) 
    # but return the start parameters unimproved
    obs_mod = param2res(ep_0)
    print(np.sum(np.array(obs)-np.array(obs_mod)))
    return ep_0
