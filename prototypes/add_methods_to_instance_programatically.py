def f(self):
    return sum(self)


class ModelVarSet(frozenset):
    sum = f


inst = ModelVarSet(frozenset({1, 2, 3}))

inst.sum()
