from typing import Set,Callable
from inspect import signature 
from functools import lru_cache,reduce

@lru_cache(maxsize=None) 
def computable_mvars(
         allComputers:Set[Callable],
         available_mvars:Set[type]
    )->Set[type]:
    # fixme mm: 
    # if possible replace by the new graph based method
    # We can already compute the graph. We only have to do it once and can easyly inver the union
    # of all nodes reachable from the startset.
    #
    # this is the old bottom up approach: repeatedly compute all
    # directly (in the next step) reachable Mvars and use the enriched set for
    # the next iteration until the set stays constant 
    dcmvs=directly_computable_mvars(allComputers,available_mvars)
    
    if dcmvs.issubset(available_mvars):
        return available_mvars
    else:
        return computable_mvars(allComputers,available_mvars.union(dcmvs))

@lru_cache(maxsize=None) 
def directly_computable_mvars(
        allComputers:Set[Callable],
        available_mvars:Set[type]
    )->Set[type]:
    # find the computers that have a source_set contained in the available_set
    return frozenset([output_mvar(c) for c in applicable_computers(allComputers,available_mvars)])

@lru_cache(maxsize=None) 
def applicable_computers(
         allComputers:Set[Callable],
         available_mvars:Set[type]
        )->Set[Callable]:
    return frozenset(
        [c for c in allComputers if input_mvars(c).issubset(available_mvars)]
    )


@lru_cache(maxsize=None) 
def all_mvars(all_computers: Set[Callable])->Set[type]:
    # the set of all mvars is implicitly defined by the 
    # parameterlists and return values of the computers
    return reduce(
        lambda acc,c : acc.union(input_mvars(c),{output_mvar(c)}),
        all_computers,
        frozenset({})
    )


def input_mvars(computer:Callable)->Set[type]:
    params = signature(computer).parameters.values()
    return frozenset({param.annotation for param in params})


def output_mvar(computer:Callable)->type:
    return signature(computer).return_annotation
