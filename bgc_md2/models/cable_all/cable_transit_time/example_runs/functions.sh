run_and_log_in_target_dir(){
  local targetDir comm logFileName
  targetDir=$1
  comm=$2
  logFileName=$3
  #echo "##################### running and logging in $targetDir"
  #echo "logfile = $targetDir/$logFileName"
  local oldDir=$(pwd)
  mkdir -p $targetDir
  cd $targetDir 
  cwd=$(pwd)
  echo "##################### running and logging in $(pwd)"
  echo "logfile = $cwd/$logFileName"
  eval ${comm}|tee $logFileName 2>&1
  cd $oldDir
}
