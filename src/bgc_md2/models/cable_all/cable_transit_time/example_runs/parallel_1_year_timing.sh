#!/bin/bash
set -e
source ./common_config
source functions.sh

complete_run(){
  spinStart=1901 
  spinEnd=1901 
  cp ${nml_dir}/rcp85_conc_kf.txt .
  mknml_files $spinStart $spinEnd
  echo copy pool file
  cp -p ${cable_run_test_lqy_dir}/spinup/poolcnp_in.csv ./poolcnp_in.csv

  echo copy restart file
  cp -p  ${cable_run_test_lqy_dir}/spinup/${filename_restart_in} ./${filename_restart_in}

  yr=$spinStart
  #prepare directory
  mkdir -p ${odir}
  cp -p cable_CN_spindump_${yr}.nml cable.nml #for the paper only CN
  timeLogFileName='time.log'

  # instead of running cable once we run it with different number of cpus
  # and capture the times in a log file
  if [ -f $timeLogFileName ]
  then
    mv $timeLogFileName ${timeLogFileName}.bak
  fi

  for np in 8 16 24 32 40 48 56 64;
   do
     odir="out/np_${np}"
     mkdir -p ${odir} #capture the output for later comparisons

     # Task:
     # 1. We want to capture the output of the time command in a single variable.
     # 2. We still want to see and log the output of cable.
     #
     # Solution:
     #    To isolate the output of time we redirect the output of the
     #    timed command to fd3 which we set before to the current value of stdout 
     #    of the outer process.
     exec 3>&1 4>&2
     #et=$( { time ls 1>&3  ;ls crap 2>&3; } 2>&1 )
     et=$( { time mpirun -n ${np} --oversubscribe cable-mpi 1>&3 2>&3; } 2>&1 )
     echo "################
  np: $np
  time: $et" | tee -a $timeLogFileName
     exec 3>&- 4>&-
  done
  }

bN=$(basename $0) #this scripts name
timeRunId="${bN%.sh}"
run_and_log_in_target_dir "${cableDataDir}/$timeRunId" complete_run "log_stdout_and_stderr_${timeRunId}"
