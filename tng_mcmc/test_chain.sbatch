#!/bin/bash

#SBATCH -n 16
#SBATCH -c 5
#SBATCH --mem-per-cpu=4G
#SBATCH -t 01:00:00

KAPPA_MIN=0.05
NWALKERS=$((2 * (SLURM_NTASKS-1)))

module load anaconda3/2021.11 gsl/2.6 hdf5/gcc/1.10.0 fftw/gcc/3.3.4 openmpi/gcc/3.1.5/64
conda activate galaxies

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -n $SLURM_NTASKS -W 0 \
  python -u chain.py $NWALKERS $KAPPA_MIN
