#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#define BASE_N 224

int calculate_optimal_N(int grid_size, int base_N) {
    return ((base_N + grid_size - 1) / grid_size) * grid_size;
}

void matrix_multiply_local(double *A, double *B, double *C, int size) {
    int i, j, k;
    for (i = 0; i < size; i++) {
        for (k = 0; k < size; k++) {
            double a_ik = A[i * size + k];
            for (j = 0; j < size; j++) {
                C[i * size + j] += a_ik * B[k * size + j];
            }
        }
    }
}

void init_matrix(double *matrix, int size, double value) {
    int i;
    for (i = 0; i < size * size; i++) {
        matrix[i] = value;
    }
}

int verify_result(double *C, int local_size, int rank, int nprocs, int N) {
    int i, error_count = 0;
    double expected = (double)N;
    
    for (i = 0; i < local_size * local_size; i++) {
        if (fabs(C[i] - expected) > 1e-6) {
            error_count++;
        }
    }
    
    int total_errors;
    MPI_Allreduce(&error_count, &total_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    return total_errors;
}

int main(int argc, char **argv) {
    int rank, nprocs;
    int grid_size, local_size, N;
    int my_row, my_col;
    int left_rank, right_rank, up_rank, down_rank;
    double *A_local, *B_local, *C_local;
    double start_time, end_time;
    MPI_Status status;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    grid_size = (int)sqrt(nprocs);
    if (grid_size * grid_size != nprocs) {
        if (rank == 0) {
            printf("Error: Number of processes must be a perfect square\n");
        }
        MPI_Finalize();
        return 1;
    }
    
    N = calculate_optimal_N(grid_size, BASE_N);
    local_size = N / grid_size;
    
    my_row = rank / grid_size;
    my_col = rank % grid_size;
    
    left_rank = my_row * grid_size + ((my_col - 1 + grid_size) % grid_size);
    right_rank = my_row * grid_size + ((my_col + 1) % grid_size);
    up_rank = ((my_row - 1 + grid_size) % grid_size) * grid_size + my_col;
    down_rank = ((my_row + 1) % grid_size) * grid_size + my_col;
    
    if (rank == 0) {
        printf("Cannon's Algorithm: N=%d, Grid=%dx%d, Local=%dx%d, Procs=%d\n", 
               N, grid_size, grid_size, local_size, local_size, nprocs);
    }
    
    A_local = (double*)calloc(local_size * local_size, sizeof(double));
    B_local = (double*)calloc(local_size * local_size, sizeof(double));
    C_local = (double*)calloc(local_size * local_size, sizeof(double));
    
    if (!A_local || !B_local || !C_local) {
        printf("Rank %d: Memory allocation failed\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    init_matrix(A_local, local_size, 1.0);
    init_matrix(B_local, local_size, 1.0);
    
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    
    for (int shift = 0; shift < my_row; shift++) {
        MPI_Sendrecv_replace(A_local, local_size * local_size, MPI_DOUBLE,
                           left_rank, 100, right_rank, 100, MPI_COMM_WORLD, &status);
    }
    
    for (int shift = 0; shift < my_col; shift++) {
        MPI_Sendrecv_replace(B_local, local_size * local_size, MPI_DOUBLE,
                           up_rank, 200, down_rank, 200, MPI_COMM_WORLD, &status);
    }
    
    for (int step = 0; step < grid_size; step++) {
        matrix_multiply_local(A_local, B_local, C_local, local_size);
        
        if (step < grid_size - 1) {
            MPI_Sendrecv_replace(A_local, local_size * local_size, MPI_DOUBLE,
                               left_rank, 300 + step, right_rank, 300 + step, MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv_replace(B_local, local_size * local_size, MPI_DOUBLE,
                               up_rank, 400 + step, down_rank, 400 + step, MPI_COMM_WORLD, &status);
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    
    if (rank == 0) {
        double total_time = end_time - start_time;
        double gflops = (2.0 * N * N * N) / (total_time * 1e9);
        printf("Time: %.6f sec, Performance: %.3f GFLOPS\n", total_time, gflops);
    }
    
    int total_errors = verify_result(C_local, local_size, rank, nprocs, N);
    if (rank == 0) {
        printf("Verification: %s\n", (total_errors == 0) ? "PASSED" : "FAILED");
    }
    
    free(A_local);
    free(B_local);
    free(C_local);
    
    MPI_Finalize();
    return 0;
}
