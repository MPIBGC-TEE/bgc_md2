source functions.sh
bN=$(basename $0)
targetDir=${bN%.sh}
echo $targetDir
runSomething(){
  touch "testfile"
  for i in $(seq 30)
  do 
    echo $i
  done
}
run_and_log_in_target_dir $targetDir runSomething "log"
