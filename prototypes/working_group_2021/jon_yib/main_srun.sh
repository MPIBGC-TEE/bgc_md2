#!/bin/bash
#SBATCH --job-name=main_srun
#SBATCH --time=80:00:00
#SBATCH --mem=10GB
#SBATCH --nodes=4
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --chdir=/scratch/mm4967/bgc_md2/prototypes/working_group_2021/jon_yib
#SBATCH --output=/scratch/mm4967/bgc_md2/prototypes/working_group_2021/jon_yib/dask.out
#SBATCH --error=/scratch/mm4967/bgc_md2/prototypes/working_group_2021/jon_yib/dask.err

module purge
module load anaconda3/2021.05 #has to match the currently activated anaconda version in the shell
conda activate /scratch/mm4967/miniconda_envs/bgc_md2

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
echo "Node number: $SLURM_NNODES"
echo "CPUs on node: $SLURM_CPUS_ON_NODE"
echo "Task number: $SLURM_NTASKS"
echo "Initial node: $SLURMD_NODENAME"
echo "Allocated nodes: $SLURM_NODELIST"
echo "OMP thread count: $OMP_NUM_THREADS"
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
mpirun -n ${SLURM_NTASKS} python main.py

