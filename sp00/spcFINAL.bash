#!/bin/bash
#PBS -q lecture-c
#PBS -l select=2:mpiprocs=112
#PBS -l walltime=00:15:00
#PBS -W group_list=gt48

cd $PBS_O_WORKDIR

#export I_MPI_DEBUG=4

cp spcDEBUG1.h spc.h 
make clean
make
mpiexec.hydra ./spc >  spcFINAL1.out
cp spcDEBUG2.h spc.h 
make clean
make
mpiexec.hydra ./spc >  spcFINAL2.out
cp spcFINAL1.h spc.h 
make clean
make
mpiexec.hydra ./spc >  spcFINAL3.out
cp spcFINAL2.h spc.h 
make clean
make
mpiexec.hydra ./spc >  spcFINAL4.out

