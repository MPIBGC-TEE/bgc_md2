from itertools import islice, tee, takewhile
import sys 
from copy import copy, deepcopy
sys.setrecursionlimit(1000000)


def natural_numbers(start):
    """defines an infinite stream of numbers 
    >=start

    lt2=natural_numbers(2)
    next(lt2)
    """
    yield start
    #while True: # faster
    #    start+=1
    #    yield start
    yield from natural_numbers(start + 1) # recursive is slower

def stream_sieve(candidate_number_stream):
    d=next(candidate_number_stream)
    yield d # note the side effects on candidate_number_stream
    yield from stream_sieve((n for n in candidate_number_stream if n%d !=0))

def stream_sieve_f(candidate_number_stream):
    d=next(candidate_number_stream)
    yield d
    yield from stream_sieve_f(filter(lambda n:n%d !=0 , candidate_number_stream))

primes=stream_sieve(natural_numbers(2))

# show first 10 prime numbers
for i in range(10):
    print(next(primes))

# or more elegant
print(list(islice(stream_sieve(natural_numbers(2)),0,10)))
print(list(islice(stream_sieve_f(natural_numbers(2)),0,10)))

#find prime twins

def primes():
    return stream_sieve(natural_numbers(2))

def prime_twins():
    prs1,prs2=tee(primes(),2)
    next(prs2)
    tups=zip(prs1,prs2)
    return filter(lambda tup:tup[1]-tup[0] ==2,tups)

for t in takewhile(lambda t:t[0]<10,prime_twins()):
    print(t)

# summary:
# One a generator is created it is impossible to get an indipendent clone
# that does not have side effects on the original. We therefore don't pass 
# generators around.

# To be able to use side effect free functions with iterator arguments we can use our own iterator (in contrast to generator) classes

class Nats():
    def __init__(self,start):
        self.start=start
        self.current=start
    
    def __next__(self):
        val=self.current
        self.current+=1
        return val

    def __iter__(self):
        #return copy(self)
        return self.__class__(self.start) # a fresh one

class Primes():
    def __init__(self):
        self.primes=stream_sieve(Nats(2))
    
    def __next__(self):
        return next(self.primes)
    
    def __iter__(self):
        #return copy(self) # the present one
        return self.__class__() # a fresh one

# now completely without streams
from functools import reduce
from itertools import takewhile, dropwhile
import numpy as np

class PrimesLazyList():
    def __init__(self):
        self.cache=[2]
    
    def __next__(self):
        val=self.cache[-1]
        candidate=val+1
        #while reduce(lambda acc, p: acc or (candidate%p ==0),self.cache,False):
        #    candidate+=1
        #
        #self.cache.append(candidate)
        n_c=Nats(candidate)
        accepted = dropwhile(
            lambda candidate: reduce(lambda acc, p: acc or (candidate%p ==0),self.cache,False),
            n_c
        )
        
        self.cache.append(next(accepted))
        
        return val

    
    def __iter__(self):
        #return copy(self) # the present one
        return self.__class__() # a fresh one

prll=PrimesLazyList()
#for i in range(10):
#    print(prll.__next__())

print(list(islice(prll,0,10000)))
