class SequenceOfUpdates_1:
    def __iter__(self):
        self.a = 1
        return self

    def __next__(self):
        x = self.a
        self.a += 1
        if x <= 4:
            return x
        else:
            raise StopIteration


mySeq_1 = SequenceOfUpdates_1()
for x in iter(mySeq_1):
    print(x)


######################## now with arguments
class SequenceOfUpdates_2:
    def __init__(self, startValue, maxValue):
        self.a = startValue
        self.max_a = maxValue

    def __iter__(self):
        return self

    def __next__(self):
        x = self.a
        self.a += 1
        if x <= self.max_a:
            return x
        else:
            raise StopIteration


mySeq_2 = SequenceOfUpdates_2(7, 8)
for x in iter(mySeq_2):
    print(x)


######################### generator version (using yield statements)


def update_generator(startValue, maxValue):
    val = startValue
    while val <= maxValue:
        yield val
        val += 1


for x in update_generator(10, 12):
    print(x)
