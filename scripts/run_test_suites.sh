#global variable for sum of return values
fail=0 
#global variable for summary
log=""
com(){
	echo "################################################################################"
	echo $(pwd)
	echo ${1}
	eval ${1}
	es=$?
	log="${log}
	${1}: ${es}
	"
	fail=$((fail + ${es}))
}
in_dir_command(){
	cd "$1"
	com "${2}"
	cd -
}

# to make shure that the python unittest module reports both failures and
# errors as exitstatus != 0 uncomment the next lines
# in_dir_command "." "python -m unittest PseudoTestError.py"
# in_dir_command "." "python -m unittest PseudoTestFail.py"

com "python scripts/test_notebooks.py tests/notebooks/"
com "python scripts/test_notebooks.py binder_notebooks/"

in_dir_command tests "python -m unittest discover -t . -p 'Test*'"
in_dir_command prototypes/working_group_2021 "python -m unittest Test_general_helpers.py"
## in_dir_command notebooks "pytest --nbmake './'"
echo "#############
summary
${log}"
exit ${fail}
