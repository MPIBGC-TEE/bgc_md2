class SequenceOfUpdates():
    def __iter__(self):
        self.a=1
        return self
    def __next__(self):
        x=self.a
        self.a+=1
        if x<4:
            return x  
        else:
            raise StopIteration

mySeq=SequenceOfUpdates()
for x in iter(mySeq):
    print(x)
