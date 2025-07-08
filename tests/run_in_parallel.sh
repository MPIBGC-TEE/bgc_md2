#!/bin/bash
# needs unittest-parallel installed
# -j 0 is the default and meanss all cores
# --level 'test' means parallelize on test level
unittest-parallel  --level 'test' --pattern 'Test*'

