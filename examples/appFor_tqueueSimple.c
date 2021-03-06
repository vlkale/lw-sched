/* -- Libraries for C++  -- */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* -- library for parallelization of code -- */
#include <pthread.h>
#define NUMTHREADS 16 // default constant and set value for the number of threads
int numThreads;
pthread_t callThread[NUMTHREADS];
pthread_barrier_t myBarrier;
pthread_mutex_t timerLock;
pthread_mutex_t myLock;
pthread_attr_t attr;

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

void* dotProdFunc(void* arg)
{
/* Test correctness, and test threads */
 double mySum = 0.0;
 long threadNum = (long) arg;
 int i = 0;
/* initialization */
int startInd = (probSize*threadNum)/numThreads;
int endInd = (probSize*(threadNum+1))/numThreads;
 while(iter < numIters) // timestep loop
  {
    mySum = 0.0; 
    if(threadNum == 0) sum = 0.0;
    if(threadNum == 0) setCDY(static_fraction, constraint, chunk_size); // set constraint parameter of scheduling strategy
    pthread_barrier_wait(&myBarrier);
    FORALL_BEGIN(statdynstaggered, 0, probSize, startInd, endInd, threadNum, numThreads)
    if(VERBOSE==1) printf("[%d] : iter = %d \t startInd = %d \t  endInd = %d \t\n", threadNum,iter, startInd, endInd);
    for (i = startInd; i < endInd; i++)
      {
	mySum += a[i]*b[i];
        //mySum += (sqrt(a[i])*sqrt(b[i])) / 4.0; // Uncomment this line and comment the line above this one if you'd like to increase the floating point operations per outer iteration done by the program.
      }
    FORALL_END(statdynstaggered, startInd, endInd, threadNum)
    if(VERBOSE == 1) printf("[%d] out of iter\n", threadNum);
    pthread_mutex_lock(&myLock);
    sum += mySum;
    pthread_mutex_unlock(&myLock);
    pthread_barrier_wait(&myBarrier);
    if(threadNum == 0) iter++;
    pthread_barrier_wait(&myBarrier);
  }
}

int main(int argc, char* argv[])
{
  int rcThread;
  long i;
  void *status;
  pthread_attr_t attr;
  double totalTime = 0.0;
  int checkSum;
  numThreads = NUMTHREADS;
  probSize = PROBSIZE;
  numIters = MAX_ITER;
  if(argc < 3)
    printf("usage: appName [probSize] [numThreads] (static_fraction) (constraint) (numIters) \n");
  else if(argc > 2)
  {
    probSize = atoi(argv[1]);
    numThreads = atoi(argv[2]);
  }
  if(argc > 3) static_fraction = atof(argv[3]);
  if(argc > 4) constraint = atof(argv[4]);
  if(argc > 5) numIters = atoi(argv[5]);
  if(argc > 6) chunk_size = atoi(argv[6]);

  printf("starting pthreads application.  threads = %d \t probSize = %d \t numIters = %d \n", numThreads, probSize, numIters);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_barrier_init(&myBarrier, NULL, numThreads);
  pthread_mutex_init(&myLock, NULL);
  vSched_init(numThreads);
  a = (float*)malloc(sizeof(float)*probSize);
  b = (float*)malloc(sizeof(float)*probSize);
  for (int i = 0 ; i < probSize ; i++)
  {
    a[i] = i*1.0;
    b[i] = 1.0;
  } // The input vectors are initialized in this way to simplify checking the correctness of the output: the sum of n numbers from 1..n is (n*(n+1))/2
  totalTime = -nont_vSched_get_wtime(); // set this to 0 because we are not in a threaded computation region
  for(i=0;i<numThreads;i++)
  {
    rcThread = pthread_create(&callThread[i], &attr, dotProdFunc, (void*)i);
    if(rcThread) printf("ERROR: return code from pthread_create() is %d \n", rcThread);
  }
  pthread_attr_destroy(&attr);
  for(i=0;i<numThreads;i++) pthread_join(callThread[i], &status);
  totalTime += nont_vSched_get_wtime(); // set this to 0 because we are not in a threaded computation region
  
  printf("totalTime: %f \n", totalTime);

  myfile = fopen("outFileVecSum.dat","a+");
  fprintf(myfile, "\t%d\t%d\t%f\t%f\n", numThreads, probSize, static_fraction, totalTime);
  fclose(myfile);
	
  printf("Completed the program vecSum. The solution is: %f \n", sum);

  pthread_barrier_destroy(&myBarrier);
  pthread_mutex_destroy(&myLock);
  vSched_finalize(numThreads);
}
