#!/bin/bash
set -e
source ./common_config
source functions.sh

complete_run(){
  local timeRunId
  timeRunId=$1
  # see ../README about how to run cable in general for several years
  # Here we try to reproduce exactly the conditions that produced
  # the data of the paper.
  # First we check if the spinup has already been run or run it.
  spinStart=1901
  spinEnd=1910
  rep1=5
  rep2=10
  spinupDir=$(spinup_dir $spinStart $spinEnd $rep1 $rep2)
  echo $spinupDir

  if [ ! -d ${spinupDir} ]; then
    mkdir $spinupDir;
  fi
  if [ ! -f "../${spinupDir/$spinnupSuccessIndicatorFile}" ]; then
    # sometimes the spinup has to be increased.
    # Use Jianyang Xia 's method to compute.
    run_and_log_in_target_dir "${spinupDir}" "run_spinup $spinStart $spinEnd $rep1 $rep2" "log_stdout_and_stderr"
  fi
    run_and_log_in_target_dir "${timeRunId}" "after_spinup_run $spinupDir $spinStart 2004" "log_stdout_and_stderr"
}

bN=$(basename $0) #this scripts name
timeRunId="${bN%.sh}"
run_and_log_in_target_dir ${cableDataDir} complete_run "log_stdout_and_stderr_${timeRunId}"
