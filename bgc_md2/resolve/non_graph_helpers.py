from typing import Set, Callable
from inspect import signature
from functools import lru_cache, reduce
from frozendict import frozendict


@lru_cache(maxsize=None)
def computable_mvars(
    allComputers: Set[Callable], available_mvars: Set[type]
) -> Set[type]:
    # fixme mm:
    # if possible replace by the new graph based method
    # We can already compute the graph. We only have to do it once and can easyly infer the union
    # of all nodes reachable from the startset.
    #
    # this is the old bottom up approach: repeatedly compute all
    # directly (in the next step) reachable Mvars and use the enriched set for
    # the next iteration until the set stays constant
    dcmvs = directly_computable_mvars(allComputers, available_mvars)

    if dcmvs.issubset(available_mvars):
        return available_mvars
    else:
        return computable_mvars(allComputers, available_mvars.union(dcmvs))


@lru_cache(maxsize=None)
def directly_computable_mvars(
    allComputers: Set[Callable], available_mvars: Set[type]
) -> Set[type]:
    # find the computers that have a source_set contained in the available_set
    return frozenset(
        [output_mvar(c) for c in applicable_computers(allComputers, available_mvars)]
    )


@lru_cache(maxsize=None)
def applicable_computers(
    allComputers: Set[Callable], available_mvars: Set[type]
) -> Set[Callable]:
    return frozenset(
        [c for c in allComputers if input_mvars(c).issubset(available_mvars)]
    )


@lru_cache(maxsize=None)
def all_computers_for_mvar(mvar: type, allComputers: Set[Callable]) -> Set[Callable]:
    return frozenset([c for c in allComputers if output_mvar(c) == mvar])


def arg_set(computer: Callable) -> Set[type]:
    params = signature(computer).parameters.values()
    return frozenset({param.annotation for param in params})


@lru_cache(maxsize=None)
def arg_set_set(mvar: type, allComputers: Set[Callable]) -> Set[Set[type]]:
    # return the set of arg_name_sets for all computers that
    # return this mvar
    return frozenset([arg_set(c) for c in all_computers_for_mvar(mvar, allComputers)])


@lru_cache(maxsize=None)
def all_mvars(all_computers: Set[Callable]) -> Set[type]:
    # the set of all mvars is implicitly defined by the
    # parameterlists and return values of the computers
    return reduce(
        lambda acc, c: acc.union(input_mvars(c), {output_mvar(c)}),
        all_computers,
        frozenset({}),
    )


def pretty_name(mvar: type, aliases: frozendict = frozendict({})) -> str:
    if len(aliases) == 0:
        s = mvar.__name__
        # return ((s.split('<')[1]).split('>')[0]).split('.')[-1]
    else:
        s = aliases[mvar.__name__]
    return s


# synonym for arg_set
def input_mvars(computer: Callable) -> Set[type]:
    params = signature(computer).parameters.values()
    return frozenset({param.annotation for param in params})


def output_mvar(computer: Callable) -> type:
    return signature(computer).return_annotation
