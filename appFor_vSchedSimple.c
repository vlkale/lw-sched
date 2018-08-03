// -- libraries for C++  --
#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#define VERBOSE 1
// --  Performance Measurement --
//#include <papi.h>  // leave PAPI out for now
double totalTime = 0.0;
//double get_wtime();
FILE* myfile;// output file for experimental data

// -- library for parallelization of code --
#include <pthread.h>
// default values for threads
#define NUMTHREADS 16
int numThreads;
pthread_t callThread[NUMTHREADS];
pthread_barrier_t myBarrier;
pthread_mutex_t timerLock;
pthread_mutex_t myLock;
pthread_attr_t attr;


// --sched. strategy library and variables --
#include "vSched.h"
double constraint;
double fs;
// in the below macros, strat is how we specify the library
#define FORALL_BEGIN(strat, s,e, start, end, tid, numThds )  loop_start_ ## strat (s,e ,&start, &end, tid, numThds);  do {
#define FORALL_END(strat, start, end, tid)  } while( loop_next_ ## strat (&start, &end, tid));
double static_fraction = 1.0; /* runtime param */ // name this differently than the actual strategy
int chunk_size  = 10;

// -- application specific defines and variables --
// default values for application
#define MAX_ITER 1000
#define PROBSIZE 16384
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
/* initialization  */
int startInd =  (probSize*threadNum)/numThreads;
int endInd = (probSize*(threadNum+1))/numThreads;
 while(iter < numIters) // timestep loop
  {
    mySum = 0.0; //reset sum to zero at the beginning of the product
    /* local computation */
    if(threadNum == 0) sum = 0.0;
    if(threadNum == 0) setCDY(static_fraction, constraint, chunk_size);
    pthread_barrier_wait(&myBarrier);
    // FORALL_BEGIN(cdy, 0, probSize, startInd, endInd, threadNum, numThreads)

#pragma omp parallel
     FORALL_BEGIN(statdynstaggered, 0, probSize, startInd, endInd, threadNum, numThreads)
     if(VERBOSE) printf("[%d] : iter = %d \t startInd = %d \t  endInd = %d \t\n", threadNum,iter, startInd, endInd);
/* get thread num and numThreads from functions */
     for (i = startInd ; i < endInd; i++)
      {
        //printf("[threadNum: %d] : iter = %d \t i = %d \t\n", threadNum, iter, i);
	mySum += a[i]*b[i];
        //mySum += (sqrt(a[i])*sqrt(b[i])) / 4.0;
      }
     //   FORALL_END(cdy, startInd, endInd, threadNum )
     FORALL_END(statdynstaggered, startInd, endInd, threadNum)
      /* thread reduction */
        //printf("[%d] out of iter\n", threadNum);
       pthread_mutex_lock(&myLock);
    sum += mySum;
    pthread_mutex_unlock(&myLock);
    pthread_barrier_wait(&myBarrier);
    if(threadNum == 0)
      iter++;
    pthread_barrier_wait(&myBarrier);
  } // end timestep loop
}

// TODO : barrier

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

  if(argc > 3)
    static_fraction = atof(argv[3]);
  if(argc > 4)
    constraint = atof(argv[4]);
  if(argc > 5)
    numIters = atoi(argv[5]);
  if(argc > 6)
    chunk_size = atoi(argv[6]);

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
  } // check this by looking at the sum of n numbers
  totalTime = -nont_vSched_get_wtime(); // set this to 0 because we are not in a threaded computation region
  for(i=0;i<numThreads;i++)
  {
    rcThread = pthread_create(&callThread[i], &attr, dotProdFunc, (void*)i);
    if(rcThread)
     printf("ERROR: return code from pthread_create() is %d \n", rcThread);
  }
  pthread_attr_destroy(&attr);
  for(i=0;i<numThreads;i++) {
   pthread_join(callThread[i], &status);
  }
  totalTime += nont_vSched_get_wtime(); // set this to 0 because we are not in a threaded computation region
  printf("totalTime: %f \n", totalTime);

  myfile = fopen("outFileVecSum.dat","a+");
  fprintf(myfile, "\t%d\t%d\t%f\t%f\n", numThreads, probSize, static_fraction, totalTime);
  fclose(myfile);
  printf("completed vecSum. Solution is : %f \n", sum);

  pthread_barrier_destroy(&myBarrier);
  pthread_mutex_destroy(&myLock);
  vSched_finalize(numThreads);
}
