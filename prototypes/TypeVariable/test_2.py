# Extension of test_0.py to an potentially unlimited number
from collections import namedtuple
from typing import TypeVar, Callable, NamedTuple, Generic

class NT_1(NamedTuple):
    a: float
    b: float

EP_1=EP[NT_1]

class NT_2(NamedTuple):
    a: float
    c: float
    d: float

EP_2=EP[NT_2]

def param2res_1(ep:EP_1) -> Obs:
    t=ep.t
    return Obs(x= t.a+t.b, y=t.a-t.b)

def param2res_2(ep:EP_2) -> Obs:
    t=ep.t
    return Obs(x= t.a+t.c+t.d, y=t.a-t.c-t.d)

def mcmc(
    ep_0:EP[T],
    param2res: Callable[[EP[T]],Obs]
    )->EP[T]:
    # we don't do anything real
    # We run the model once (thereby making sure that it does run) 
    # but return the start parameters unimproved
    obs=param2res(ep_0)
    return ep_0

res1=mcmc(
    EP(NT_1(a=3,b=5)),
    param2res_1
)    
print(res1)

res2=mcmc(
    EP(NT_2(a=3,c=5,d=3)),
    param2res_2
)    
print(res2)

res3=mcmc(
    #this should fail the typecheck (and also produce a runtime error)
    EP(NT_2(a=3,c=5,d=3)),
    param2res_1
)    
