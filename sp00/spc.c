#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "spc.h"

#define  NPROCS   224

#define  EPS    2.220446e-16

double  A[N][N];
double  b[M][N];
double  x[M][N];
double  c[N];


int     myid, numprocs;


void spc(double [N][N], double [M][N], double [M][N], int, int); 

void main(int argc, char* argv[]) {

     double  t0, t1, t2, t_w;
     double  dc_inv, d_mflops, dtemp, dtemp2, dtemp_t;

     int     ierr;
     int     i, j;
     int     ii;      
     int     ib;

     ierr = MPI_Init(&argc, &argv);
     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
     ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);



      /* matrix generation --------------------------*/
      for(j=0; j<N; j++) {
        ii = 0;
        for(i=j; i<N; i++) {
          A[j][i] = (N-j) - ii;
          A[i][j] = A[j][i];
          ii++;
        }
      }
      /* end of matrix generation -------------------------- */

     /* set vector b  -------------------------- */
      for (i=0; i<N; i++) {
        b[0][i] = 0.0;
        for (j=0; j<N; j++) {
          b[0][i] += A[i][j];
        }
      }
      for (i=0; i<M; i++) {
        for (j=0; j<N; j++) {
          b[i][j] = b[0][j];
        }
      }
     /* ----------------------------------------------------- */


     /* Start of spc routine ----------------------------*/
     ierr = MPI_Barrier(MPI_COMM_WORLD);
     t1 = MPI_Wtime();

     spc(A, b, x, N, M);

     //ierr = MPI_Barrier(MPI_COMM_WORLD);
     t2 = MPI_Wtime();
     t0 =  t2 - t1; 
     ierr = MPI_Reduce(&t0, &t_w, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
     /* End of spc routine --------------------------- */

     if (myid == 0) {

       printf("--------------------------- \n");
       printf("N = %d , M = %d \n",N,M);
       printf("LU solve time  = %lf [sec.] \n",t_w);

       d_mflops = 2.0/3.0*(double)N*(double)N*(double)N;
       d_mflops += 7.0/2.0*(double)N*(double)N;
       d_mflops += 4.0/3.0*(double)N*(double)M;
       d_mflops = d_mflops/t_w;
       d_mflops = d_mflops * 1.0e-6;
       printf(" %lf [MFLOPS] \n", d_mflops);

     }

     /* Verification routine ----------------- */
     ib = N / NPROCS;
     dtemp_t = 0.0;
     for(i=0; i<M; i++) {
       for(j=myid*ib; j<(myid+1)*ib; j++) {
         dtemp2 = x[i][j] - 1.0;
         dtemp_t += dtemp2*dtemp2;
       }
     }
     dtemp_t = sqrt(dtemp_t);
     /* -------------------------------------- */

     MPI_Reduce(&dtemp_t, &dtemp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

     /* do not modify follows. -------- */ 
     if (myid == 0) {
       dtemp2 = (double)N*(double)N;      
       dtemp_t = EPS*(double)N*dtemp2*(double)M;
       dtemp_t = sqrt(dtemp_t);
       printf("Pass value: %e \n", dtemp_t);
       printf("Calculated value: %e \n", dtemp);
       if (dtemp > dtemp_t) {
          printf("Error! Test is falled. \n");
          exit(1);
       } 
       printf(" OK! Test is passed. \n");
       printf("--------------------------- \n");
     }
     /* ----------------------------------------- */


     ierr = MPI_Finalize();

     exit(0);
}


void spc(double A[N][N], double b[M][N], double x[M][N], int n, int m) 
{
     int i, j, k, ne;
     double dtemp;
     register double temp1, temp2, temp3, temp4, temp5, temp6;
     
     for (k = 0; k < n; k++) {
         dtemp = 1.0 / A[k][k];
         
         #pragma vector aligned
         #pragma ivdep
         for (i = k + 1; i < n; i++) {
             A[i][k] *= dtemp;   
         }
         
         for (j = k + 1; j < n; j++) {
             temp1 = A[j][k];
             i = k + 1;
             for (; i + 5 < n; i += 6) {
                 temp2 = A[k][i];
                 temp3 = A[k][i+1]; 
                 temp4 = A[k][i+2];
                 temp5 = A[k][i+3];
                 temp6 = A[k][i+4];
                 dtemp = A[k][i+5];
                 
                 A[j][i]   -= temp2 * temp1;
                 A[j][i+1] -= temp3 * temp1;
                 A[j][i+2] -= temp4 * temp1;
                 A[j][i+3] -= temp5 * temp1;
                 A[j][i+4] -= temp6 * temp1;
                 A[j][i+5] -= dtemp * temp1;
             }
             for (; i < n; i++) {
                 A[j][i] -= A[k][i] * temp1; 
             }
         }
     }

     for (ne = 0; ne < m; ne++) {
         
         for (k = 0; k < n; k++) {
             temp1 = b[ne][k];
             j = 0;
             for (; j + 5 < k; j += 6) {
                 temp1 -= A[k][j]   * c[j] +
                          A[k][j+1] * c[j+1] +
                          A[k][j+2] * c[j+2] +
                          A[k][j+3] * c[j+3] +
                          A[k][j+4] * c[j+4] +
                          A[k][j+5] * c[j+5];
             }
             for (; j < k; j++) {
                 temp1 -= A[k][j] * c[j];
             }
             c[k] = temp1;
         }
         x[ne][n-1] = c[n-1] / A[n-1][n-1];
         
         for (k = n-2; k >= 0; k--) {
             temp1 = c[k];
             j = k + 1;
             for (; j + 5 < n; j += 6) {
                 temp1 -= A[k][j]   * x[ne][j] +
                          A[k][j+1] * x[ne][j+1] +
                          A[k][j+2] * x[ne][j+2] +
                          A[k][j+3] * x[ne][j+3] +
                          A[k][j+4] * x[ne][j+4] +
                          A[k][j+5] * x[ne][j+5];
             }
             for (; j < n; j++) {
                 temp1 -= A[k][j] * x[ne][j];
             }
             
             x[ne][k] = temp1 / A[k][k];
         }
     }
}
