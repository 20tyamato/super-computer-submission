#!/bin/bash
#PBS -q lecture-c
#PBS -l select=1:mpiprocs=4
#PBS -l walltime=00:01:00
#PBS -W group_list=gt48

cd $PBS_O_WORKDIR

#export I_MPI_DEBUG=4

mpiexec.hydra ./sp02 > sp02.out
