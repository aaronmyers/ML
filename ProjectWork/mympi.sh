#!/bin/bash
#SBATCH -J myMPI           # job name
#SBATCH -o myMPI.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 256              # total number of mpi tasks requested
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -t 04:00:00        # run time (hh:mm:ss) - 1.5 hours
# SBATCH --mail-user=aaron.myers@utexas.edu
# SBATCH --mail-type=begin  # email me when the job starts
# SBATCH --mail-type=end    # email me when the job finishes
ibrun -n 128 -o 0 mpigradientopen /scratch/01396/naga86/hw7/livejournal.dat .85 50 10 4 &      # run the MPI executable named a.out
ibrun -n 128 -o 128 mpigradientopen /scratch/01396/naga86/hw7/livejournal.dat .85 200 10 4          # run the MPI executable named a.out
