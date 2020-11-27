#!/usr/bin/env python3
# vim:set ff=unix expandtab ts=4 sw=4:
# this is a pure python version
# run with pyhton3 run_tests.py in a venv

#from concurrencytest import ConcurrentTestSuite, fork_for_tests
import unittest
import sys
from dask.distributed import Client, LocalCluster

def main():
    print("\n###################### running tests ##########################\n")

    s = unittest.defaultTestLoader.discover("", pattern="localTestCableHelpers*")
    # p = unittest.defaultTestLoader.discover('', pattern='Pinned_Test*')
    # s.addTests(p)
    # concurrent_suite = s
    #concurrent_suite = ConcurrentTestSuite(s, fork_for_tests(64))
    #res = r.run(concurrent_suite)
    r = unittest.TextTestRunner()
    res = r.run(s)
    return res


if __name__ == "__main__":
    cl=LocalCluster()
    c=Client(cl)
    #fut = c.submit(main)
    #res=fut.result()
    #if len(res.errors) + len(res.failures) > 0:
    #    sys.exit(1)
    main()


