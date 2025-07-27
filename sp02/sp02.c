#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void matrix_transpose_scatter_gather(double *local_A, double *local_AT, 
                                   int N, int local_rows, int rank, int size) {
    double *global_A = NULL;
    double *global_AT = NULL;
    
    if (rank == 0) {
        global_A = (double*)calloc(N * N, sizeof(double));
        global_AT = (double*)calloc(N * N, sizeof(double));
        printf("Root: Allocated global matrices\n");
    }
    
    MPI_Gather(local_A, local_rows * N, MPI_DOUBLE,
               global_A, local_rows * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("Root: Gathered matrix:\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%6.1f ", global_A[i * N + j]);
            }
            printf("\n");
        }
        printf("\n");
        
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                global_AT[j * N + i] = global_A[i * N + j];
            }
        }
        
        printf("Root: Transposed matrix:\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%6.1f ", global_AT[i * N + j]);
            }
            printf("\n");
        }
        printf("\n");
    }
    
    MPI_Scatter(global_AT, local_rows * N, MPI_DOUBLE,
                local_AT, local_rows * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        free(global_A);
        free(global_AT);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) {
        printf("Starting matrix transpose with %d processes\n", size);
    }
    
    const int N = 8;
    const int local_rows = N / size;
    
    printf("Rank %d: local_rows = %d\n", rank, local_rows);
    
    double *local_A = (double*)malloc(local_rows * N * sizeof(double));
    double *local_AT = (double*)malloc(local_rows * N * sizeof(double));
    
    for (int i = 0; i < local_rows; i++) {
        for (int j = 0; j < N; j++) {
            int global_row = rank * local_rows + i;
            local_A[i * N + j] = global_row * 10 + j + 1;
        }
    }
    
    for (int p = 0; p < size; p++) {
        if (rank == p) {
            printf("Rank %d - Initial local matrix:\n", rank);
            for (int i = 0; i < local_rows; i++) {
                for (int j = 0; j < N; j++) {
                    printf("%6.1f ", local_A[i * N + j]);
                }
                printf("\n");
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    matrix_transpose_scatter_gather(local_A, local_AT, N, local_rows, rank, size);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int p = 0; p < size; p++) {
        if (rank == p) {
            printf("Rank %d - Transposed local matrix:\n", rank);
            for (int i = 0; i < local_rows; i++) {
                for (int j = 0; j < N; j++) {
                    printf("%6.1f ", local_AT[i * N + j]);
                }
                printf("\n");
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    free(local_A);
    free(local_AT);
    MPI_Finalize();
    return 0;
}
