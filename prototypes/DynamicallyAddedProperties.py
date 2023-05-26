class Dyn():
    def __init__(self,fd):
        self.fd=fd
        for k,v in fd.items():
            self.__setattr__(k,v)

d=Dyn({"a":5})
print(d.fd)
print(d.__dir__())
print(d.a)
