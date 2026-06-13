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
  mpirun -np ${np} --oversubscribe ${cable_exe}
}

bN=$(basename $0) #this scripts name
timeRunId="${bN%.sh}"
run_and_log_in_target_dir "${cableDataDir}/$timeRunId" complete_run "log_stdout_and_stderr_${timeRunId}"
