/*  code for testing strategies with several different MPI+OpenMP apps
 */

/* ToDO:  
1.  make notes on how to add a new app
2.  make notes on how to add a new strategy
*/

//include files

// parallelization libraries 
#include <omp.h>
#include <pthread.h>
//#include "uSchedLib.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

FILE* perftuningData;

// Libraries for Performance Measurement
// for timers
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>

// for lower-level Perf. Profiling

//#define PAPI_PROFILING

#ifdef PAPI_PROFILING
#include <papi.h>
#endif

double get_wtime(void);

#define MAXTHRDS 16

//#define _GNU_SOURCE

// parallelism data structures
pthread_barrier_t myBarrier;
pthread_mutex_t myMutex;
pthread_attr_t attr;
pthread_t callThd[MAXTHRDS];

//strategies vars
double blockSize; //block size

int strat = 0 ;

// implementation  strategies
int numThreads;
int coresUsed;

int nodes;
int conf;

int schedLib = 1; // implementation of scheduler (as opposed to type of sched strategy)

char* algStrategy ;

//app problem size and other app considerations
int arrSize;

// app data structures for dot product
double* a;
double* b;
double* c;

// for scalar product
double* d;

// app data structures for spMV
double* x;
double* x_next;
double** solverMatrix;

// app data structures for stencil
double** grid1;
double** grid2;

// app data structures for indirect accesses
double* indAccArr1;
double* indAccArr2;


double globalProd = 0.0;
double scalarProdCheckSum = 0.0;

//hardware-specific
//TODO: need to find a way to obtain this from system.  Make sure this works for any number of specified threads.
int numCores = 16;

//strategy specification
// app specification
// TO DO: avoid repeated code in functions

void* iterativeSolve(void* tid)
{
  // pthread binding management - system-specific
  int s, cNum;
  cpu_set_t cpuset;
  pthread_t thread;
  thread = pthread_self();

  double startTime, endTime, totalTime;
  long prof1, prof2;

  long t = *((int*) tid);
  int myTid = (int) t;

  CPU_ZERO(&cpuset);
  CPU_SET(myTid, &cpuset);
  s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
  
  if (s != 0)
    printf("error: with setting pthread affinity. Error code = %d\n", s);

  double prod = 1.0;

  //warm up iters
#ifdef TRANSPOSE
  for(int i = myTid ; i < arrSize; i+=numThreads)
    {
      x_next[i] = b[i];
      for(int k = 0; k< arrSize ; k++)
        x_next[i] = x_next[i] + solverMatrix[i][k]*x[k];
    }
  pthread_barrier_wait(&myBarrier);
  for(int i = myTid ; i< arrSize; i+= numThreads)
    x[i] = x_next[i];
#endif

#ifdef BLOCK
  for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
    {
      x_next[i] = b[i];
      for(int k = 0; k< arrSize ; k++)
        x_next[i] = x_next[i] + solverMatrix[i][k]*x[k];
    }
  pthread_barrier_wait(&myBarrier);
  for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
    x[i] = x_next[i];
#endif
  // prod = 1.0;
  pthread_barrier_wait(&myBarrier);
  startTime = get_wtime();
  for(int compIter = 0; compIter < 10; compIter++ )
    {
      // transpose algorithm strategy
#ifdef TRANSPOSE
      for(int i = myTid ; i < arrSize; i+=numThreads)
	{
	  x_next[i] = b[i];
	  for(int k = 0; k< arrSize ; k++)
	    x_next[i] = x_next[i] + solverMatrix[i][k]*x[k];
	}
      pthread_barrier_wait(&myBarrier);
      for(int i = myTid ; i< arrSize; i+= numThreads)
	x[i] = x_next[i];
#endif

#ifdef BLOCK
      for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
	{
	  x_next[i] = b[i];
	  for(int k = 0; k< arrSize ; k++)
	    x_next[i] = x_next[i] + solverMatrix[i][k]*x[k];
	}
      pthread_barrier_wait(&myBarrier);
      for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
	x[i] = x_next[i];
#endif
    }
  // may need to take this out
  pthread_mutex_lock(&myMutex);
  globalProd += prod;
  pthread_mutex_unlock(&myMutex);

  pthread_barrier_wait(&myBarrier);
  endTime = get_wtime();
  totalTime = endTime - startTime;
  if(myTid == 0)
    {
      printf("Time(secs) to run spMV, reported from threadID %d:\t%f\n", myTid, totalTime);
      if((perftuningData = fopen("perfTuningData.dat", "a+")) != NULL)
	{
	  fprintf(perftuningData,  "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%f\t%ld\t%ld\n",
		  "spMV", arrSize, conf,
		  nodes, coresUsed, algStrategy, schedLib,
		  totalTime, prof1, prof2);
	  fclose(perftuningData);
	}
    }
}

void* daxpy(void* tid)
{
  // pthread binding management - system-specific
  int s, cNum;
  cpu_set_t cpuset;
  pthread_t thread;
  thread = pthread_self();

  int block;
  double startTime, endTime, totalTime;

  long prof1, prof2;

  long t = *((int*) tid);
  int myTid = (int) t;

  CPU_ZERO(&cpuset);
  CPU_SET(myTid, &cpuset);
  s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
  if (s != 0)
    printf("error: with setting pthread affinity. Error code = %d\n", s);

  double prod = 1.0;

  //warm up iters
#ifdef TRANSPOSE
  for(int ind = myTid ; ind < arrSize; ind+=numThreads)
      c[ind] +=
        a[ind]*b[ind];
#endif

#ifdef BLOCKCYCLIC
  for (int block = myTid; block*blockSize < arrSize; block+=numThreads)
    for (int i=block*blockSize; i < (block+1)*blockSize; i++)
      c[i] += a[i]*b[i];
#endif

#ifdef BLOCK
  for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
    c[i] += a[i]*b[i];
#endif

  // prod = 1.0;
  pthread_barrier_wait(&myBarrier);
  startTime = get_wtime();
  for(int compIter = 0; compIter < 10000; compIter++ )
    {
      // transpose algorithm strategy
#ifdef TRANSPOSE
      for(int ind = myTid ; ind < arrSize; ind+=numThreads)
	c[ind] += a[ind]*b[ind];
#endif

#ifdef BLOCKCYCLIC
      for (int block = myTid; block*blockSize < arrSize; block+=numThreads)
	for (int i=block*blockSize; i < (block+1)*blockSize; i++)
	  c[i] += a[i]*b[i];
#endif

#ifdef BLOCK
      for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
	c[i] += a[i]*b[i];
#endif
    }
  // may need to take this out
  pthread_mutex_lock(&myMutex);
  globalProd += prod;
  pthread_mutex_unlock(&myMutex);

  pthread_barrier_wait(&myBarrier);
  endTime = get_wtime();
  totalTime = endTime - startTime;

  if(myTid == 0)
    {
      printf("Time(secs) to run daxpy, reported from threadID %d:\t%f\n", myTid, totalTime);
      if((perftuningData = fopen("perfTuningData.dat", "a+")) != NULL)
	{
	  fprintf(perftuningData, "%s\t%d\t%d\t%d\t%d\t%s\t%d\t%f\t%ld\t%ld\n", "daxpy", arrSize, conf, nodes,
		  coresUsed, algStrategy, schedLib,
		  totalTime, prof1, prof2);

	  fclose(perftuningData);
	}
    }
}

void* stencil(void* tid)
{
  // pthread binding management - system-specific
  int s, cNum;
  cpu_set_t cpuset;
  pthread_t thread;
  thread = pthread_self();
  double startTime, endTime, totalTime;

  //profile vars
  long prof1;
  long prof2;

  long t = *((int*) tid);
  int myTid = (int) t;

  CPU_ZERO(&cpuset);
  CPU_SET(myTid, &cpuset);
  s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
  if (s != 0)
    printf("error: with setting pthread affinity. Error code = %d\n", s);

  double prod = 1.0;

#ifdef PAPI_PROFILING
  if(myTid==0)
    PAPI_library_init(PAPI_VER_CURRENT); // should be in the main program
  if(myTid == 0)
    PAPI_thread_init(pthread_self); // I don't understand this
#endif
  //printf("initialized\n" );
#ifdef PAPI_PROFILING
  // PAPI perf counters  (using PAPI high-level interface
  unsigned long int papi_tid;
  // int num_hwcntrs = 3;
  int Events[3]  = {PAPI_L2_DCM, PAPI_SR_INS, PAPI_LD_INS};
  int retval;
  int EventSet = PAPI_NULL;
  long_long values[3];
  if ((papi_tid = PAPI_thread_id()) == (unsigned long int)-1)
    {
      printf("error with assigning the papi thread ID\n");
      exit(1);
    }
  printf("Initial PAPI thread id is:\t%lu\n", papi_tid);
  //if (PAPI_start_counters(Events, num_hwcntrs) != PAPI_OK)
  if(PAPI_create_eventset(&EventSet) != PAPI_OK)
    printf("Error with creating event set \n");
  if(PAPI_add_event(EventSet, Events[0]) != PAPI_OK)
    printf("Error with creating event set 0 \n");
  if(PAPI_add_event(EventSet, Events[1]) != PAPI_OK)
    printf("Error with creating event set 1 \n");
  if(PAPI_add_event(EventSet, Events[2]) != PAPI_OK)
    printf("Error with creating event set 2 \n");
#endif

  //  printf(" starting user code \n ") ;
  //  return 1;
  //warm up iters
#ifdef TRANSPOSE
  for(int i = myTid+1; i < arrSize - 1; i+=numThreads)
    for(int j = 1; j < arrSize-1; j++)
      grid2[i][j] = (grid1[i][j] + grid1[i-1][j] + grid1[i+1][j]
                     + grid1[i][j+1]+grid1[i][j-1]) /5.0;

  pthread_barrier_wait(&myBarrier);
  for(int i = myTid+1; i < arrSize - 1; i+=numThreads)
    for(int j = 1; j < arrSize-1; j++)
      grid1[i][j] = (grid2[i][j] + grid2[i-1][j] + grid2[i+1][j]
		     + grid2[i][j+1]+grid2[i][j-1]) /5.0;
  pthread_barrier_wait(&myBarrier) ;
#endif

#ifdef BLOCK
  for(int i = ((arrSize-2)*myTid)/numThreads + 1; i < ((arrSize-2)*(myTid+1))/numThreads +  1; i++)
    for(int j = 1; j < arrSize - 1 ; j++)
      grid2[i][j] = (grid1[i][j] + grid1[i-1][j] + grid1[i+1][j]
                     + grid1[i][j+1]+grid1[i][j-1]) /5.0;
  pthread_barrier_wait(&myBarrier);
  for(int i = ((arrSize-2)*myTid)/numThreads + 1; i < ((arrSize-2)*(myTid+1))/numThreads +  1; i++)
    for(int j = 1; j < arrSize - 1 ; j++)
      grid1[i][j] = (grid2[i][j] + grid2[i-1][j] + grid2[i+1][j]
                     + grid2[i][j+1]+grid2[i][j-1]) /5.0;
  pthread_barrier_wait(&myBarrier);
#endif
  pthread_barrier_wait(&myBarrier);
  startTime = get_wtime();

  for(int compIter = 0; compIter < 10; compIter++ )
    {
#ifdef TRANSPOSE
      for(int i = myTid+1; i < arrSize - 1; i+=numThreads)
	for(int j = 1; j < arrSize-1; j++)
	  grid2[i][j] = (grid1[i][j] + grid1[i-1][j] + grid1[i+1][j]
			 + grid1[i][j+1]+grid1[i][j-1])/5.0;

      pthread_barrier_wait(&myBarrier);

      for(int i = myTid+1; i < arrSize - 1; i+=numThreads)
	for(int j = 1; j < arrSize-1; j++)
	  grid1[i][j] = (grid2[i][j] + grid2[i-1][j] + grid2[i+1][j]
			 + grid2[i][j+1]+grid2[i][j-1]) /5.0;
      pthread_barrier_wait(&myBarrier);

#endif

#ifdef BLOCK
      for(int i = ((arrSize-2)*myTid)/numThreads + 1; i < ((arrSize-2)*(myTid+1))/numThreads +  1; i++)
	for(int j = 1; j < arrSize - 1 ; j++)
	  grid2[i][j] = (grid1[i][j] + grid1[i-1][j] + grid1[i+1][j]
			 + grid1[i][j+1]+grid1[i][j-1]) /5.0;

      pthread_barrier_wait(&myBarrier);

      for(int i = ((arrSize-2)*myTid)/numThreads + 1; i < ((arrSize-2)*(myTid+1))/numThreads +  1; i++)
	for(int j = 1; j < arrSize - 1 ; j++)
	  grid1[i][j] = (grid2[i][j] + grid2[i-1][j] + grid2[i+1][j]
			 + grid2[i][j+1]+grid2[i][j-1]) /5.0;
      pthread_barrier_wait(&myBarrier);
#endif

    }

  pthread_barrier_wait(&myBarrier);
  endTime = get_wtime();
  totalTime = endTime - startTime;

#ifdef PAPI_PROFILING
  // end PAPI perf counters
  //  if (PAPI_read_counters(values, num_hwcntrs) != PAPI_OK)
  //    printf("Error with reading counters!\n");
  values = {0};
  if (PAPI_read(EventSet, values) != PAPI_OK)
    printf("Error with reading values from EventSet\n");
  // if (PAPI_stop_counters(values, num_hwcntrs) != PAPI_OK)
  //  printf("Error with stopping counters!\n");
  if (PAPI_stop(EventSet, values) != PAPI_OK)
    printf("Error with stopping counters!\n");
  //if(papi_tid == 0)
  printf("Total Prof: L1 data cache misses, reported from papi tid %lld:\t%lld\n", papi_tid, values[0]);
  printf("Total Prof: data TLB misses, reported from papi tid %lld:\t%lld\n", papi_tid, values[1]);
  printf("Total Prof: L1 data cache hit rate, reported from papi tid %lld:\t%f\n", papi_tid, 1.0 - 1.0*(values[0])/(1.0*(values[1] + values[2])));
#endif

  if(myTid == 0)
    {
      printf("Time (secs) to run stencil, reported from threadID %d:\t%f\n", myTid, totalTime);
      if((perftuningData = fopen("perfTuningData.dat", "a+")) != NULL)
	{
	  fprintf(perftuningData,  "%s\t%d\t%d\t%d\t%d\t%s\t%d\t%f\t%ld\t%ld\n",
		  "sten", arrSize, conf,
		  nodes, coresUsed, algStrategy, schedLib,
		  totalTime, prof1, prof2);
	  fclose(perftuningData);
	}
    }
}

void* scalarProd(void* tid)
{
  // pthread binding management - system-specific
  int s, cNum;
  cpu_set_t cpuset;
  pthread_t thread;
  thread = pthread_self();

  double startTime, endTime, totalTime;
  long prof1, prof2;

  long t = *((int*) tid);
  int myTid = (int) t;
  /* Set affinity mask to include CPUs 0 to numCores - TODO: done in app or in runtime ? */
  CPU_ZERO(&cpuset);
  CPU_SET(myTid, &cpuset);
  s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
  if (s != 0)
    printf("error: with setting pthread affinity. Error code = %d\n", s);
  double prod = 1.0;

  //warm up iters
#ifdef TRANSPOSE
  for(int i = myTid ; i < arrSize; i+=numThreads)
    {
      //  printf("thread %d \t  start %d \t end %d \n ", myTid,  myTid, arrSize - myTid );
      x[i] = x[i]*1.00001;
    }
#endif

#ifdef BLOCK
  for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
    {
      // printf("thread %d \t  start %d \t end %d \n ", myTid,  (arrSize*myTid)/numThreads, (arrSize*(myTid+1))/numThreads );
      x[i] = x[i]*1.00001;
    }
#endif
  // prod = 1.0;
  pthread_barrier_wait(&myBarrier);
  startTime = get_wtime();
  for(int compIter = 0; compIter < 10000; compIter++ )
    {
      // transpose algorithm strategy
#ifdef TRANSPOSE
      for (int i = myTid ; i < arrSize; i+=numThreads)
	x[i] = x[i]*1.00001;
#endif

      //TODO: print each thread.
      //TODO: print product 
      //block algorithmic strategy
#ifdef BLOCK
      for (int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
	x[i] = x[i]*1.00001;
#endif
    }
  // may need to take this out
  pthread_mutex_lock(&myMutex);
  globalProd += prod;
  pthread_mutex_unlock(&myMutex);
  pthread_barrier_wait(&myBarrier);
  endTime = get_wtime();
  totalTime = endTime - startTime;
  if(myTid == 0)
    {
      printf("Time(secs) to run Scalar Prod, reported from threadID %d:\t%f\n", myTid, totalTime);
      if((perftuningData = fopen("perfTuningData.dat", "a+")) != NULL)
	{
	  fprintf(perftuningData,  "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%f\t%ld\t%ld\n",
		  "scpd", arrSize, conf,
		  nodes, coresUsed, algStrategy, schedLib,
		  totalTime, prof1, prof2);
	  fclose(perftuningData);
	}
    }
}


//dot Product application
void* dotProd(void* tid)
{
  // pthread binding management - system-specific
  int s, cNum;//TODO: need a better name for cNum
  cpu_set_t cpuset;
  pthread_t thread;
  thread = pthread_self();
  double startTime, endTime, totalTime;
  long prof1, prof2;
  long t = *((int*) tid);
  int myTid = (int) t;
  /* Set affinity mask to include CPUs 0 to numCores */
  CPU_ZERO(&cpuset);
  CPU_SET(myTid, &cpuset);
  s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
  if (s != 0)
    printf("error: with setting pthread affinity. Error code = %d\n", s);
  double prod = 0.0;
  //warm up iters
#ifdef TRANSPOSE
  for(int i = myTid ; i < arrSize; i+=numThreads)
    {
      prod += a[i]*b[i];
    }
#endif

#ifdef BLOCK
  for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
    {
      prod += a[i]*b[i];
    }
#endif
  prod = 0.0;
  pthread_barrier_wait(&myBarrier);
  startTime = get_wtime();
  // transpose algorithm strategy
  for(int compIter = 0; compIter < 10000; compIter++ )
    {
#ifdef TRANSPOSE
      for(int i = myTid ; i < arrSize; i+=numThreads)
	prod += a[i]*b[i];
#endif

#ifdef BLOCK
      for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
	prod += a[i]*b[i];
#endif
    }
  pthread_mutex_lock(&myMutex);
  globalProd += prod;
  pthread_mutex_unlock(&myMutex);
  pthread_barrier_wait(&myBarrier);
  endTime = get_wtime();
  totalTime = endTime - startTime;
  if(myTid == 0)
    {
      printf("Time(secs) to run dot Prod, reported from threadID %d:\t%f\n", myTid, totalTime);
      if((perftuningData = fopen("perfTuningData.dat", "a+")) != NULL)
	{
	  fprintf(perftuningData,  "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%f\t%ld\t%ld\n",
		  "sten", arrSize, conf,
		  nodes, coresUsed, strat, schedLib,
		  totalTime, prof1, prof2);
	  fclose(perftuningData);
	}
    }
}


void* prefixSum(void* tid)
{ 


}


void* largestSubset(void* tid)
{


} 

//indirect accesses application
void* indAcc(void* tid)
{
  // pthread binding management - system-specific
  int s, cNum;//TODO: need a better name for cNum
  cpu_set_t cpuset;
  pthread_t thread;
  thread = pthread_self();

  double startTime, endTime, totalTime;

  long prof1, prof2;

  long t = *((int*) tid);
  int myTid = (int) t;

  /* Set affinity mask to include CPUs 0 to numCores */
  CPU_ZERO(&cpuset);
  CPU_SET(myTid, &cpuset);
  s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
  if (s != 0)
    printf("error: with setting pthread affinity. Error code = %d\n", s);
  double prod = 0.0;

  //warm up iters
#ifdef TRANSPOSE
  for(int i = myTid ; i < arrSize; i+=numThreads)
    indAccArr1[i]  += indAccArr1[indAccArr2[i]];
#endif

#ifdef BLOCK
  for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
    indAccArr1[i]  += indAccArr1[indAccArr2[i]];
#endif
  prod = 0.0;
  pthread_barrier_wait(&myBarrier);
  startTime = get_wtime();
  // transpose algorithm strategy
  for(int compIter = 0; compIter < 10000; compIter++ )
    {
#ifdef TRANSPOSE
      for(int i = myTid ; i < arrSize; i+=numThreads)
	indAccArr1[i]  += indAccArr1[indAccArr2[i]];
#endif

#ifdef BLOCK
      for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
	indAccArr1[i]  += indAccArr1[indAccArr2[i]];
#endif
    }

  pthread_mutex_lock(&myMutex);
  globalProd += prod;
  pthread_mutex_unlock(&myMutex);
  pthread_barrier_wait(&myBarrier);
  endTime = get_wtime();
  totalTime = endTime - startTime;

  if(myTid == 0)
    {
      printf("Time(secs) to run dot Prod, reported from threadID %d:\t%f\n", myTid, totalTime);
      if((perftuningData = fopen("perfTuningData.dat", "a+")) != NULL)
	{
	  fprintf(perftuningData,  "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%f\t%ld\t%ld\n",
		  "iacc", arrSize, conf,
		  nodes, coresUsed, strat, schedLib,
		  totalTime, prof1, prof2);
	  fclose(perftuningData);
	}
    }
}


//indirect accesses application
void* nbody(void* tid)
{
  // pthread binding management - system-specific
  int s, cNum;//TODO: need a better name for cNum
  cpu_set_t cpuset;
  pthread_t thread;
  thread = pthread_self();
  double startTime, endTime, totalTime;
  long prof1, prof2;

  long t = *((int*) tid);
  int myTid = (int) t;

  /* Set affinity mask to include CPUs 0 to numCores */
  CPU_ZERO(&cpuset);
  CPU_SET(myTid, &cpuset);
  s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
  if (s != 0)
    printf("error: with setting pthread affinity. Error code = %d\n", s);
  double prod = 0.0;

  //warm up iters
#ifdef TRANSPOSE
  for(int i = myTid ; i < arrSize; i+=numThreads)
    // indAccArr1[i] += indAccArr1[indAccArr2[i]];  <-- figure out how to do this
#endif

#ifdef BLOCK
    for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
      indAccArr1[i] += indAccArr1[indAccArr2[i]];
#endif

  prod = 0.0;
  pthread_barrier_wait(&myBarrier);
  startTime = get_wtime();

  // transpose algorithm strategy
  for(int compIter = 0; compIter < 10000; compIter++ )
    {
#ifdef TRANSPOSE
      for(int i = myTid ; i < arrSize; i+=numThreads)
	indAccArr1[i] += indAccArr1[indAccArr2[i]];
#endif

#ifdef BLOCK
      for(int i = (arrSize*myTid)/numThreads; i < (arrSize*(myTid+1))/numThreads; i++)
	indAccArr1[i]  += indAccArr1[indAccArr2[i]];
#endif
    }

  pthread_mutex_lock(&myMutex);
  globalProd += prod;
  pthread_mutex_unlock(&myMutex);

  pthread_barrier_wait(&myBarrier);
  endTime = get_wtime();
  totalTime = endTime - startTime;

  if(myTid == 0)
    {
      printf("Time(secs) to run dot Prod, reported from threadID %d:\t%f\n", myTid, totalTime);
      if((perftuningData = fopen("perfTuningData.dat", "a+")) != NULL)
	{
	  fprintf(perftuningData,  "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%f\t%ld\t%ld\n",
		  "nbod", arrSize, conf,
		  nodes, coresUsed, strat, schedLib,
		  totalTime, prof1, prof2);
	  fclose(perftuningData);
	}
    }
}



int main(int argc, char* argv[] )
{
  if(argc >= 2)
    {
      arrSize = atoi(argv[1]);
      numThreads = atoi(argv[2]);
    }
  else
    {
      printf("Usage: a.out [problem size] [<number of threads>]\n");
      exit(1); // can be system specific
    }
  //initialize app data structure
  a = (double*) malloc(sizeof(double)*arrSize);
  b = (double*) malloc(sizeof(double)*arrSize);
  c = (double*) malloc(sizeof(double)*arrSize);

  // for scalar prod
  d = (double*) malloc(sizeof(double)*arrSize);

  x = (double*) malloc(sizeof(double)*arrSize);
  x_next = (double*) malloc(sizeof(double)*arrSize);

  indAccArr1 = (double*) malloc(sizeof(double)*arrSize);
  indAccArr2 = (double*) malloc(sizeof(double)*arrSize);

  //nned for spMV
  solverMatrix = (double**) malloc(sizeof(double*)*arrSize);
  for(int i = 0; i < arrSize ; i++)
    solverMatrix[i] = (double*) malloc (sizeof(double)*arrSize);
  // needed for stencil
  grid1 = (double**) malloc (sizeof(double*)*arrSize);
  grid2 = (double**) malloc (sizeof(double*)*arrSize);

  for(int i = 0; i < arrSize; i++)
    {
      grid1[i] = (double*) malloc (sizeof(double)*arrSize);
      grid2[i] = (double*) malloc (sizeof(double)*arrSize);
    }

  for(int i = 0; i < arrSize; i++)
    {
      a[i] = 2.0;
      b[i] = 3.0;
      c[i] = 3.0 ;
      d[i] = 4.0 ;

      x[i] = 4.0;
      x_next[i] = x[i];
      srand(time(NULL));

      indAccArr1[i] = (double)rand()/(double)RAND_MAX;
      indAccArr2[i] = (double)rand()/(double)RAND_MAX;

      for(int j = 0; j < arrSize; j++)
	{
	  solverMatrix[i][j] = drand48();
	  grid1[i][j]=(double)rand()/(double)RAND_MAX;
	  grid2[i][j] = 0.0;
	}
    }
  globalProd = 0.0;

  //initialize pthread
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  //initialize barrier and mutex
  pthread_barrier_init(&myBarrier, NULL, numThreads);
  pthread_mutex_init(&myMutex, NULL);

  void* thread_status;
  int* thdsArr = (int*) malloc(sizeof(int)*numThreads);

  //print parallelism characteristics
  printf("Number of Threads per Process:\t%d\n", numThreads);

  //print application characteristics or strategy
  printf("Problem size:\t%d\n", arrSize);

  //print the perf. tuning technique used by runtime
#ifdef TRANSPOSE
  printf("Alg. Strategy: Transpose\n");
  algStrategy = "cyc";
#endif

#ifdef BLOCK
  printf("Alg. Strategy: Block\n");
  algStrategy = "blk";
#endif

#ifdef BLOCKCYCLIC
  printf("Alg. Strategy: Block-Cyclic\n");
  algStrategy = "bcy";
#endif

#ifdef CLA
  printf("Alg. Strategy: Block\n");
  algStrategy = "cla";
#endif

#ifdef PAPI_PROFILE
  PAPI_library_init(PAPI_VER_CURRENT); // should be in the main program
  PAPI_thread_init(pthread_self);
#endif
  //TODO: check if done for each core

  // goto daxpycode;
  goto stencilcode;


 dotProdcode:
  //run the dot product  (begin threaded computation region for the dot Product application)
  for(int threadID = 0; threadID < numThreads; threadID++)
    {
      thdsArr[threadID] = threadID; // add this to ensure we don't get duplicated threadID
      pthread_create(&callThd[threadID], &attr, dotProd, &thdsArr[threadID]);
    }
  pthread_attr_destroy(&attr);


  // TODO: check if this is correct
  for(int i=0;i < numThreads; i++)
    pthread_join(callThd[i], &thread_status);

 scalarProdCode:
  // run the scalar product in aparallel for multiple cores.  (begin threaded computation region for the scalar product applicaiton)
  for(int threadID = 0; threadID < numThreads; threadID++)
    {
      thdsArr[threadID] = threadID; // add this to ensure we don't get duplicated threadID
      pthread_create(&callThd[threadID], &attr, scalarProd, &thdsArr[threadID]);
    }

  pthread_attr_destroy(&attr);

  // TODO: check if this is correct, may need to do loop as above, with thdsArr.
  for(int threadID=0;threadID<numThreads;threadID++)
    pthread_join(callThd[threadID], &thread_status); 
  //TODO: do a profile for each app.

 daxpycode:
  // run the daxpy in aparallel for multiple cores.  (begin threaded computation region for the scibalar product applicaiton)
  for(int threadID = 0; threadID < numThreads; threadID++)
    {
      thdsArr[threadID] = threadID; // add this to ensure we don't get duplicated threadID
      pthread_create(&callThd[threadID], &attr, daxpy, &thdsArr[threadID]);
    }
  pthread_attr_destroy(&attr); 
  // TODO: check if this is correct, may need to do loop as above, with thdsArr.
  for(int threadID=0;threadID<numThreads;threadID++)
    pthread_join(callThd[threadID], &thread_status);

  // run the scalar product in aparallel for multiple cores.  (begin threaded computation region for the scalar product applicaiton)
  for(int threadID = 0; threadID < numThreads; threadID++)
    {
      thdsArr[threadID] = threadID; // add this to ensure we don't get duplicated threadID
      pthread_create(&callThd[threadID], &attr, iterativeSolve, &thdsArr[threadID]);
    }
  pthread_attr_destroy(&attr);
  //  TODO: get access to shared bus
  // TODO: check if this is correct, may need to do loop as above, with thdsArr.
  for(int threadID=0;threadID<numThreads;threadID++)
    pthread_join(callThd[threadID], &thread_status); 

  // run the stencil in aparallel for multiple cores.  (begin threaded computation region for the scalar product applicaiton)
 stencilcode:
  for(int threadID = 0; threadID < numThreads; threadID++)
    {
      thdsArr[threadID] = threadID; // add this to ensure we don't get duplicated threadID
      pthread_create(&callThd[threadID], &attr, stencil, &thdsArr[threadID]);
    }
  pthread_attr_destroy(&attr);

  // TODO: check if this is correct, may need to do loop as above, with thdsArr.
  for(int threadID=0;threadID<numThreads;threadID++)
    pthread_join(callThd[threadID], &thread_status);


  // run the stencil in aparallel for multiple cores.  (begin threaded computation region for the scalar product applicaiton)
 indArrAcccode:
  for(int threadID = 0; threadID < numThreads; threadID++)
    {
      thdsArr[threadID] = threadID; // add this to ensure we don't get duplicated threadID
      pthread_create(&callThd[threadID], &attr, indAcc, &thdsArr[threadID]);
    }
  pthread_attr_destroy(&attr);

  // TODO: check if this is correct, may need to do loop as above, with thdsArr.
  for(int threadID=0;threadID<numThreads;threadID++)
    pthread_join(callThd[threadID], &thread_status);
  // correctness checks for application
  printf("correctness check for Dot Product application: global prod = %f \n", globalProd);

  for (int i = 0 ; i < arrSize ; i++)
    scalarProdCheckSum += x[i];
  printf("correctness check for Scalar Product application: checkSum = %f \n", scalarProdCheckSum);
  //TODO: check to see why we can't get the right result for cache miss rate below

}


// functions for obtaining results
double get_wtime(void)
{
  struct rusage ruse;
  //  getrusage(RUSAGE_SELF, &ruse);
  //see documentation , need to use rusage thread
  getrusage(RUSAGE_THREAD, &ruse);
  return( (double)(ruse.ru_utime.tv_sec+ruse.ru_utime.tv_usec / 1000000.0) );
}

double get_profile(void)
{


}
