from os import environ
from time import sleep
from random import random
from datetime import datetime

from dask_mpi import initialize
from distributed import Client


def slow_is_prime(num):
    sleep(random())
    if num == 1:
        return False
    for test_factor in range(2, num // 2):
        if num % test_factor == 0:
            return False

    return True


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    #parser.add_argument('min_num', type=int)
    #parser.add_argument('max_num', type=int)
    args = parser.parse_args()

    num_threads = int(environ.get(
        'SLURM_CPUS_PER_TASK',
        environ.get('OMP_NUM_THREADS', 1)
    ))
    initialize(interface='ib0', nthreads=num_threads)
    client = Client()
    
    min_num = 10
    max_num = 100
    start_time = datetime.now()
    num_primes = sum(
        client.gather(
            client.map(slow_is_prime,
                       range(min_num, max_num + 1))
        )
    )
    end_time = datetime.now()

    print(f'{num_primes} primes between {min_num} and {max_num} '
          f'[{end_time - start_time}]')


if __name__ == '__main__':
    main()