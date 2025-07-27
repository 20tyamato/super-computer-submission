#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int local_value = rank + 1;
    int reduce_result = -999;
    int allreduce_result = -999;
    
    printf("=== BEFORE REDUCTION ===\n");
    printf("Rank %d: local_value = %d\n", rank, local_value);
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Reduce(&local_value, &reduce_result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    printf("\n=== AFTER MPI_Reduce ===\n");
    printf("Rank %d: reduce_result = %d\n", rank, reduce_result);
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Allreduce(&local_value, &allreduce_result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    printf("\n=== AFTER MPI_Allreduce ===\n");
    printf("Rank %d: allreduce_result = %d\n", rank, allreduce_result);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("\n=== SUMMARY ===\n");
        printf("MPI_Reduce: Only rank 0 gets the result\n");
        printf("MPI_Allreduce: ALL ranks get the same result\n");
    }
    
    MPI_Finalize();
    return 0;
}
