# example of numpy matrix-matrix multiplication with threads
from os import environ
environ['OMP_NUM_THREADS'] = '8'
from numpy.random import rand
# size of arrays
n = 8000
# create an array of random values
data1 = rand(n, n)
data2 = rand(n, n)
# matrix-matrix multiplication
result = data1.dot(data2)
