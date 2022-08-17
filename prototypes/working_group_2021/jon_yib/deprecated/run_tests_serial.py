#!/usr/bin/env python3
# vim:set ff=unix expandtab ts=4 sw=4:
# this is a pure python version 
# run with pyhton3 run_tests.py in a venv

import unittest
import sys
from pathlib import Path

def main():

    print("\n###################### running tests ##########################\n")

    s = unittest.defaultTestLoader.discover('', pattern='Test*')
    r = unittest.TextTestRunner()

    res = r.run(s)
    if len(res.errors) + len(res.failures) > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
