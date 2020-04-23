#include <stdio.h>
#include <stdlib.h>
#ifdef USE_MPI
#include <mpi.h>                /* MPI header file */
#endif
#ifdef _OPENMP
#include "omp.h"            /* OpenMP header file */
#endif
#define NUM_STEPS 10000
#define MAX_THREADS 64

int main(int argc, char *argv[])
{
  int nprocs, myid;
  int tid, nthreads, nbin;
  double start_time, end_time;
  double pi = 0.0, Psum=0.0, sum[MAX_THREADS]={0.0};
  double step = 1.0/(double) NUM_STEPS;/* initialize for MPI */
  #ifdef USE_MPI
  MPI_Init(&argc, &argv);/* starts MPI *//* get number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);/* get this process's number (ranges from 0 to nprocs - 1) */
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  #else
  myid = 0;
  nprocs = 1;
  #endif

  nbin= NUM_STEPS/nprocs;
  // MPI_OpenMP version
#pragma omp parallel private(tid)
   {
    int i;
    double x;
#ifdef _OPENMP
    nthreads=omp_get_num_threads();
    tid=omp_get_thread_num();
#else
    nthreads=1;
    tid=0;
#endif
    for (i=nbin*myid+tid; i < nbin*(myid+1); i+= nthreads)
      {
	/* changed*/ x = (i+0.5)*step;
	sum[tid] += 4.0/(1.0+x*x);
      }
   }
  for(tid=0; tid<nthreads; tid++)   /*sum by each mpi process*/
    Psum += sum[tid]*step;
#ifdef USE_MPI
  MPI_Reduce(&Psum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);/* added */
#endif
  if (myid == 0)
    {
      printf("parallel program results with %d processes and %d threads:\n", nprocs, nthreads);
      #ifdef USE_MPI
      printf("pi = %g  (%17.15f)\n",pi, pi);
      #else
      printf("pi = %g  (%17.15f)\n",Psum, Psum);
      #endif
    }
  #ifdef USE_MPI
  MPI_Finalize();
  #endif
  return 0;
}
