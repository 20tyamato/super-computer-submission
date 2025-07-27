#!/bin/bash
#PBS -q lecture-c
#PBS -l select=2:mpiprocs=112
#PBS -l walltime=00:01:00
#PBS -W group_list=gt48
#PBS -o sp04.out
#PBS -j oe

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1

run_test() {
    local np=$1
    local grid_desc=$2
    
    local available_procs=$(cat $PBS_NODEFILE | wc -l)
    if [ $np -gt $available_procs ]; then
        echo "$np ($grid_desc): SKIP (not enough procs)"
        return
    fi
    
    if mpirun -np $np ./sp04 2>/dev/null; then
        echo "$np ($grid_desc): ✓"
    else
        echo "$np ($grid_desc): ✗"
    fi
}

echo "Testing Cannon's Algorithm:"
run_test 4 "2x2"
run_test 9 "3x3"
run_test 16 "4x4"
run_test 25 "5x5"
run_test 36 "6x6"
run_test 49 "7x7"
run_test 64 "8x8"
run_test 100 "10x10"
run_test 144 "12x12"
