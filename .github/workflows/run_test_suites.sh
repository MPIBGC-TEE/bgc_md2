#global variable for sum of return values
fail=0 
in_dir_command(){
	cd $1
	echo "################################################################################"
	echo $(pwd)
	echo ${2}
	eval $2
	fail=$((fail + $?))
	cd -
}

# to make shure that the python unittest module reports both failures and
# errors as exitstatus != 0 uncomment the next lines
# in_dir_command "." "python -m unittest PseudoTestError.py"
# in_dir_command "." "python -m unittest PseudoTestFail.py"

in_dir_command tests "python -m unittest discover -t . -p 'Test*'"
in_dir_command prototypes/working_group_2021 "python -m unittest Test_general_helpers.py"
in_dir_command notebooks "pytest --nbmake './'"
exit ${fail}
