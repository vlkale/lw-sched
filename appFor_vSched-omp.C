/* -- Libraries for C++  -- */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* -- library for parallelization of code -- */
// #include <pthread.h>
#include <omp.h>

/* variables to support parallelization of code */
#define NUMTHREADS 16 // default constant and set value for the number of threads
int numThreads;


/* -- Debugging -- */
#define VERBOSE 1

/* --  Performance Measurement -- */
double totalTime = 0.0;
FILE* myfile;// output file for experimental data

/* --  Hardware Profiling -- */
// #define USE_PAPI 
// #ifdef USE_PAPI
// #include <papi.h>  // Leave PAPI out for now to make fully portable. Some platforms don't have the L2 and L3 cache misses available. TODO: need a way to check the counters in config or programmatically. 
// #endif

/* --Library for scheduling strategy and variables and macros associated with the library -- */
#include "vSched.h"
double constraint;
double fs;
// in the below macros, strat is how we specify the library
#define FORALL_BEGIN(strat, s,e, start, end, tid, numThds )  loop_start_ ## strat (s,e ,&start, &end, tid, numThds);  do {
#define FORALL_END(strat, start, end, tid)  } while( loop_next_ ## strat (&start, &end, tid));
double static_fraction = 1.0; /* runtime param */ // name this differently than the actual strategy
int chunk_size  = 10;

/* -- Application specific #defines and variables -- */
#define MAX_ITER 1000 
#define PROBSIZE 16384 // Default values based on architecture. For a proper test, the problem size ought to be such that data goes out of cache.
int probSize;
int numIters;
double sum = 0.0;
double globalSum = 0.0;
int iter = 0;
float* a;
float* b;

// function testing the scheduler vSched
void dotProdFunc(void* arg)
{
/* Test correctness, and test threads */
 double mySum = 0.0;
 int threadNum;
 int numThreads;
 int i = 0;
/* initialization */
int startInd;
int endInd;

while(iter < numIters) // timestep, or outer iteration, loop
  {
    mySum = 0.0;
   sum = 0.0;
   setCDY(static_fraction, constraint, chunk_size); // set constraint parameter of scheduling strategy
#pragma omp parallel 
   {
   threadNum = omp_get_thread_num();
   numThreads = omp_get_num_threads();
   // The first parameter is the loop scheduling strategy.
   FORALL_BEGIN(statdynstaggered, 0, probSize, startInd, endInd, threadNum, numThreads)
   if(VERBOSE==1) printf("[%d] : iter = %d \t startInd = %d \t  endInd = %d \t\n", threadNum,iter, startInd, endInd);
    for (i = startInd; i < endInd; i++)
      {
          mySum += a[i]*b[i];
      }
    FORALL_END(statdynstaggered, startInd, endInd, threadNum)
    if(VERBOSE == 1) printf("[%d] out of iter\n", threadNum);
    #pragma omp critical
      sum += mySum; // the vSched library could support reductions too. However, better would be to just use user-defined schedules being proposed for OpenMP and then use OpenMP's reductions. Code is kept like this to make an illustrative point about this software design tradeoff.
    } // end omp paralllel
     iter++;
  } // end timestep loop 
} // end dotProdFunc

int main(int argc, char* argv[])
{
  long i;
  double totalTime = 0.0;
  int checkSum;
  numThreads = NUMTHREADS;
  probSize = PROBSIZE;
  numIters = MAX_ITER;
  if(argc < 3)
    printf("usage: appName [probSize][numThreads] (static_fraction) (constraint) (numIters) \n");
  else if(argc > 2)
  {
    probSize = atoi(argv[1]);
    numThreads = atoi(argv[2]);
  }
  if(argc > 3) static_fraction = atof(argv[3]);
  if(argc > 4) constraint = atof(argv[4]);
  if(argc > 5) numIters = atoi(argv[5]);
  if(argc > 6) chunk_size = atoi(argv[6]);

  printf("starting OpenMP application using vSched. threads = %d \t probSize = %d \t numIters = %d \n", numThreads, probSize, numIters);
  vSched_init(numThreads);
  a = (float*)malloc(sizeof(float)*probSize);
  b = (float*)malloc(sizeof(float)*probSize);
	
  // initialize input vectors, use standard worksharing here. 
  #pragma omp parallel for
  for (int i = 0 ; i < probSize ; i++)
  {
    a[i] = i*1.0;
    b[i] = 1.0;
  } // The input vectors are initialized in this way to simplify checking the correctness of the output: the sum of n numbers from 1..n is (n*(n+1))/2 
 
  totalTime = -nont_vSched_get_wtime(); 
  dotProdFunc(NULL);
  totalTime += nont_vSched_get_wtime(); 
  
  printf("totalTime: %f \n", totalTime);

  myfile = fopen("outFilePerf.dat","a+");
  fprintf(myfile, "\t%d\t%d\t%f\t%f\n", numThreads, probSize, static_fraction, totalTime);
  fclose(myfile);

  printf("Completed the program dot Prod for testing vSched. The solution of the program is: %f \n", sum);

  vSched_finalize(numThreads);
}
