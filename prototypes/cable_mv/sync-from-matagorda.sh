src_dir="/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4/"
target_dir=${src_dir}
rsync -r --progress matagorda-from-home:${src_dir} ${target_dir} --exclude="zarr*"
