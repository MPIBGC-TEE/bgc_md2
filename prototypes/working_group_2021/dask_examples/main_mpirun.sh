#!/bin/bash
#SBATCH --job-name=main_mpirun
#SBATCH --time=00:35:00
#SBATCH --mem=3GB
#SBATCH --nodes=4
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --chdir=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/dask_examples
#SBATCH --output=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/dask_examples/dask.out.%J
#SBATCH --error=/scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/dask_examples/dask.err.%J

module purge
module load anaconda3/2021.05
conda activate /scratch/jw2636/CMIP6/bgc_md2/prototypes/working_group_2021/yy_cable/bgc_dask_2021.05

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
echo "Output from main_dask.py:"

# Get Slurm to run the Python program that uses Dask
mpirun --np ${SLURM_NTASKS} python main_dask.py

