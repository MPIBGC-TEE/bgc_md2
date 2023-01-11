fail=0
in_dir_command(){
	cd $1
	echo "################################################################################"
	echo $(pwd)
	echo $2
	eval $2
	fail = $?
	cd -
}
  
in_dir_command tests "python -m unittest discover -t . -p 'Test*'"
in_dir_command prototypes/working_group_2021 "python -m unittest Test_general_helpers.py"
in_dir_command notebooks "pytest --nbmake './'"

# make sure that the script has non-zero exit status if one of the suites has failing test
test(fail)        
