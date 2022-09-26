#!/bin/bash
#SBATCH --job-name=dask_exp
#SBATCH --time=00:10:00
#SBATCH --mem=10GB
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=5
#SBATCH --chdir=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable
#SBATCH --output=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable/dask.out.%J
#SBATCH --error=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable/dask.err.%J

module purge
module load anaconda3/2020.11
conda activate /scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable/dask

conda list
which python
python --version
which conda
conda --version

# Get Slurm to run the Python program that uses Dask
srun --overcommit \
     --distribution=cyclic \
	 --nodes=${SLURM_NNODES} \
	 --ntasks=$[SLURM_NTASKS+2] \
	 --cpus-per-task=${SLURM_CPUS_PER_TASK} \
	 python dask_example.py

