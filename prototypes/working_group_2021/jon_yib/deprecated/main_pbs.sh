#!/bin/bash
#PBS -A UOKL0017
#PBS -q economy 
#PBS -N yibs
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1
#PBS -o /glade/scratch/jonw/bgc_md2/prototypes/working_group_2021/jon_yib/dask.out
#PBS -e /glade/scratch/jonw/bgc_md2/prototypes/working_group_2021/jon_yib/dask.err

module load conda/latest 
module load openmpi
conda activate bgc_env

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

conda list
echo ""
echo "Conda location and version:"
which conda
conda --version
echo ""
echo "Subsequent python location and version:"
which python
python --version
echo ""
echo "Job information create by bash file:"
echo "Node number: ${PBS_NUM_NODES}"
echo "CPUs on node: ${PBS_NUM_PPN}"
echo "Task number: ${PBS_TASKNUM}"
echo "Initial node: ${PBS_O_HOST}"
echo "Allocated nodes: ${PBS_NODEFILE}"
echo "OMP thread count: ${OMP_NUM_THREADS}"
echo ""
echo "Output from main.py:"

# Get Slurm to run the Python program that uses Dask
# Get Slurm to run the Python program that uses Dask
#srun --overcommit \
#        --distribution=cyclic \
#        --nodes=${SLURM_NNODES} \
#	--ntasks=${SLURM_NTASKS} \
#	--cpus-per-task=${SLURM_CPUS_PER_TASK} \
#        python main.py
mpirun -n ${PBS_TASKNUM} python main_local.py

