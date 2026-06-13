#!/bin/bash
#SBATCH --job-name=dask_exp
#SBATCH --time=00:10:00
#SBATCH --mem=1GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --chdir=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/dask_examples
#SBATCH --output=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/dask_examples/dask.out.%J
#SBATCH --error=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/dask_examples/dask.err.%J

module purge
module load anaconda3/2021.05
conda activate /scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable/bgc_dask_2021.05

conda list
which python
python --version
which conda
conda --version

# Get Slurm to run the Python program that uses Dask
srun python dask_jobq.py

