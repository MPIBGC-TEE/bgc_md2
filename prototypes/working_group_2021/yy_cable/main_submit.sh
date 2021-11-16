#!/bin/bash
#SBATCH --job-name=bgc_yy_cable
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --chdir=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable
#SBATCH --output=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable/yy_cable_out.o

module purge
module load anaconda3/2020.11
conda activate /scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable/bgcEnv_2020.11

conda list
which python
python --version
which conda
conda --version
srun python3 main.py
