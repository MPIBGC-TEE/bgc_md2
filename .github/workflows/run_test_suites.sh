in_dir_command(){
	cd $1
	echo "################################################################################"
	echo $(pwd)
	echo $2
	eval $2
	cd -
}
  
in_dir_command tests "python -m unittest discover -t . -p 'Test*'"
    
in_dir_command prototypes/working_group_2021 "python -m unittest Test_general_helpers.py"

in_dir_command notebooks "pytest --nbmake './'"

        
