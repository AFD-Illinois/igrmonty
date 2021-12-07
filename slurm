#!/bin/bash
#SBATCH --account=def-dhaggard
#SBATCH --job-name=S93MR40i50
#SBATCH --nodes=1         # Number of nodes
#SBATCH --ntasks-per-node=1  # Number of tasks per node
#SBATCH --cpus-per-task=40     # CPU cores per MPI process
#SBATCH --mem=120G        # memory per node
#SBATCH --time=0-12:00      # time (DD-HH:MM)

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load hdf5
module load gsl
#module load intel

cd /lustre04/scratch/astrodoo/iGRmonty_EHT/grmonty_SA93M/Rh40/inc50
make clean
make

srun -n1 -N1 -c40 --cpu-bind=threads ./run_grmonty.sh >& run_grmonty.err
