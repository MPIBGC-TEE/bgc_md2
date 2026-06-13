#!/bin/bash

#SBATCH --job-name=main_srun
#SBATCH --time=00:05:00
#SBATCH --mem=5GB
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1

module purge
module load anaconda3/2021.05
conda activate /scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable/bgc_dask_2021.05

rm /scratch/$USER/scheduler.json
configfile="silly.config"
python make_multi-prog-config.py $USER $SLURM_NTASKS $configfile silly_mpirun.py
 
srun -n ${SLURM_NTASKS} -l --multi-prog $configfile
