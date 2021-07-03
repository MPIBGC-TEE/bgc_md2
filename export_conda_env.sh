suffix="_$(whoami)_$(date -Iminutes)_$(hostname)"
path=conda_environment_files/environment${suffix}.yml

conda env export > $path
echo "env saved in $path."
echo "To create an exact copy type:"
echo "conda env create -f $path"
