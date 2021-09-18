#!/bin/bash 
#SBATCH -A PHY129
#SBATCH -J grmonty
#SBATCH -N 1
#SBATCH -t 1:00:00

module load gsl
module load hdf5
module load intel

cd $MEMBERWORK/phy129/iGRmonty_kappa/
make clean
make 

srun -n1 -N1 -c32 --cpu-bind=threads ./grmonty 10000 /gpfs/alpine/proj-shared/phy129/NC3R256/RT/grmonty_150.hdf5 8.400000e-18 >& spectrum_150_sp1e5.err
mv spectrum.h5 spectrum_150_sp1e5.h5

