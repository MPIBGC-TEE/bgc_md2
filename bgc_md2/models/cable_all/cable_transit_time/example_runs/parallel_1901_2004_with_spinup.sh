#!/bin/bash
set -e
source ./common_config

bN=$(basename $0) #this scripts name
runDirName="${bN%.sh}"
run_and_log_in_target_dir \
  ${cableDataDir} \
  "complete_run $runDirName 1901 1910 5 10 1901 2004" \
  "log_stdout_and_stderr_${runDirName}"
