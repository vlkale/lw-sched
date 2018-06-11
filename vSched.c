#include <pthread.h>

//flag used for debugging output
//#define VERBOSE

#define PROFILING

#include <stdio.h> // use this for testing
pthread_mutex_t sched_lock;

#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
// could declare below as static
float f_s;
float f_d;
float constraint_;
int chunkSize;
int nextChunk;
int loopEnd;
int isLoopStarted;

int count ;
int threadCount;

// functions internal to the vSched library
int get_constraint();
double vSched_get_wtime();


typedef struct PossibleWork  // coem up with a better name
{
  int nextChunk;
  int chunkSize;
  int limit; //this is the iteration that we should not assign.
  pthread_mutex_t qLock;
} PossibleWork ;

 PossibleWork** dynwork; 
int selectAnotherThread(int tid, int numThreads);

int vSched_thread_init()
{ 
  // allocate the data structure for each thread
}

void vSched_init(int numThreads)
{
  pthread_mutex_init(&sched_lock, NULL);
  threadCount = numThreads;
  dynwork = (PossibleWork**) malloc(sizeof(void*)*numThreads);
  for (int i = 0 ; i < numThreads; i++)
  {
    dynwork[i] =  ( PossibleWork* ) malloc(sizeof(PossibleWork));
    pthread_mutex_init(&(dynwork[i]->qLock), NULL);
  }
}
void vSched_finalize(int numThreads)
{
  pthread_mutex_destroy(&sched_lock);
}

//int loop_start_static(int loopBegin, int loopEnd, int *pstart, int *pend, int threadID, int numThreads)
//{
// *pstart = loopBegin + ((loopEnd - loopBegin)*threadID)/numThreads; /* figure out algebra  here , based on loopBegin */
// *pend = loopBegin + ((loopEnd - loopBegin)*(threadID+1))/numThreads; /* figure out algebra  here , based on loopBegin ;  check */
// return 1;
//}


//int loop_next_static(int *pstart, int *pend)
//{
// return 0;
//}

void setStaticFraction(float f, int _chunkSize)
{
  f_s = f;
  chunkSize = _chunkSize;
  isLoopStarted = 0;
}

//int loop_start_static_fraction(int loopBegin, int _loopEnd, int *pstart, int *pend, int threadID, int numThreads)
//{
// pthread_mutex_lock(&sched_lock);
// if(!isLoopStarted)
// {
//   loopEnd = _loopEnd;
//   nextChunk = loopBegin + (loopEnd - loopBegin)*f_s;
//   isLoopStarted = 1;
// }
// pthread_mutex_unlock(&sched_lock);
//  *pstart = loopBegin + (((loopEnd - loopBegin)*threadID)*f_s)/numThreads; /* figure out algebra  here , based on loopBegin */
// *pend = loopBegin + (((loopEnd - loopBegin)*(threadID+1))*f_s)/numThreads; /* figure out algebra  here , based on loopBegin, check */
//return 1;
// }

// int loop_next_static_fraction(int *pstart, int *pend)
// {
// if(nextChunk  >=  loopEnd )
// {
//  isLoopStarted = 0;  /* protect with lock, or make only thread 0 do it */
    /* might be good place to do tasklet locality here */
//    return 0;
// }
// *pstart = nextChunk;
// nextChunk = nextChunk + chunkSize;
// *pend  = nextChunk;
//  if(*pend  > loopEnd)
//   *pend = loopEnd;
// return 1;
//}

void setCDY(float f, double c, int _chunkSize)
{
  f_s = f; f_d = 1.0 - f; constraint_ = c; chunkSize = _chunkSize; isLoopStarted = 0;
}

/*
  this is the initialization function for the constrained dynamic scheduling
*/
int loop_start_cdy(int loopBegin, int _loopEnd, int *pstart, int *pend, int threadID, int numThreads )  // think about adding to parameter list here
{
  pthread_mutex_lock(&sched_lock);
// printf("loop_start_cdy(): thread %d \n", threadID);
  if(!isLoopStarted)
  {
    loopEnd = _loopEnd;
    nextChunk = loopBegin + (loopEnd - loopBegin)*f_s;
    isLoopStarted = 1;
    #ifdef VERBOSE
      printf("loop_start_cdy(): thread %d : isLoopStarted = %d \t nextChunk = %d \t loopEnd = %d \n", threadID, isLoopStarted, nextChunk, loopEnd);
      #endif
  }
  pthread_mutex_unlock(&sched_lock);
  *pstart = loopBegin + (((loopEnd - loopBegin)*threadID)*f_s)/numThreads; /* figure out algebra  here , based on loopBegin */
  *pend = loopBegin + (((loopEnd - loopBegin)*(threadID+1))*f_s)/numThreads; /* figure out algebra  here , based on loopBegin, check */
// print ostart pend
#ifdef VERBOSE
    printf("thread%d:\t pstart =  %d \t pend = %d \t nextChunk = %d \n ", threadID, *pstart, *pend, nextChunk);
#endif
 return 1;
}

/*
    return 0 means that there is no more work
 */
int loop_next_cdy(int *pstart, int *pend, int tid)
{
#ifdef VERBOSE
  printf("starting loop_next_cdy() \t pstart %d \t pend %d \n", *pstart, *pend);
#endif
  if(isLoopStarted == 0)
    return 0;

#ifdef PROFILING
  double time_loop_next = 0.0;
  time_loop_next = - vSched_get_wtime();
#endif
  pthread_mutex_lock(&sched_lock);
  if((nextChunk  >=  loopEnd)  && (isLoopStarted == 1)) // if next chunk greater than end bound, and no one else noticed
  {
    isLoopStarted = 0;  /* protect with lock, or make only thread 0 do it */
    /* might be good place to do tasklet locality here */
#ifdef VERBOSE
    printf("loop ended\n");
#endif

#ifdef PROFILING
    time_loop_next += vSched_get_wtime();
    printf("loop_next_sds(): loop ended: thread %d \t dequeue time = %f  \n" , tid, time_loop_next);
#endif
    pthread_mutex_unlock(&sched_lock);
    return 0;
  }
  if(get_constraint())
  {
    *pstart = nextChunk;
    nextChunk = nextChunk + chunkSize;
    *pend  = nextChunk;
  }
  else // the constraint is not satisfied, so we make the thread do a dummy piece of work
  {
#ifdef VERBOSE
    printf("condition not passed. Thread working on dummy tasklet, and then dequeueing again.\n");
#endif
    *pstart = 0;  // this should generate bus traffic, it's only hitting registers
    *pend = 0;
  }
  pthread_mutex_unlock(&sched_lock);
#ifdef VERBOSE
    printf(" Loop_next_cdy(): \t pstart =  %d \t pend = %d \t nextChunk = %d \n ", *pstart, *pend, nextChunk);
#endif
  if(*pend  > loopEnd)
    *pend = loopEnd;

#ifdef PROFILING
  time_loop_next += vSched_get_wtime();
  printf("loop_next_sds(): thread %d \t dequeue time = %f  \n" ,tid, time_loop_next);
#endif

  return 1;
}

int loop_start_statdynstaggered(int loopBegin, int _loopEnd, int *pstart, int *pend, int threadID, int numThreads )  // think about adding to parameter list here
{
  pthread_mutex_lock(&sched_lock); // check if the lock is needed around the whole conditional clause
  if(count == 0)
    count = numThreads;
  pthread_mutex_unlock(&sched_lock);
  pthread_mutex_lock(&dynwork[threadID]->qLock);
// printf("loop_start_cdy(): thread %d \n", threadID);
  loopEnd = _loopEnd;
  dynwork[threadID]->limit = loopBegin + ((loopEnd - loopBegin)*(threadID+1))/numThreads;
  *pstart = loopBegin + (((loopEnd - loopBegin)*threadID))/numThreads; /* figure out algebra  here , based on loopBegin */
  *pend = *pstart + ((loopEnd - loopBegin)*f_s)/numThreads; /* figure out algebra  here , based on loopBegin, check */
  dynwork[threadID]->nextChunk = *pend;
  dynwork[threadID]->chunkSize = chunkSize;
  pthread_mutex_unlock(&(dynwork[threadID]->qLock));
  if (dynwork[threadID]->nextChunk >= dynwork[threadID]->limit) // this thread is done with its own queue
  {
    pthread_mutex_lock(&sched_lock);
    count--;
    pthread_mutex_unlock(&sched_lock);
  }
// print ostart pend
#ifdef VERBOSE
    printf("thread%d:\t pstart =  %d \t pend = %d \t nextChunk = %d \n ", threadID, *pstart, *pend, dynwork[threadID]->nextChunk);
#endif
  return 1;
}

/*
- This function chooses another thread to steal from, based on the threadId it is given.
*/
int selectAnotherThread(int tid, int numThreads)
{
  return 0;
  int another_tid = -1;
  #ifdef VERBOSE
    printf("given threadID %d \t thd to steal from is %d \n", tid, another_tid);
#endif

//for  now we choose the next tid
//  another_tid = (tid + 1)%numThreads;
  another_tid = (tid + 1)%numThreads ; // hardcode to 16 for now
// Strategy 1: randomized stealing
// Strategy 2: steal from thread with most elements in queue
  return another_tid;
}

/*
This is the staggered method for mixed static/dynamic scheduling .
*/
int loop_next_statdynstaggered(int *pstart, int *pend, int tid)
{
#ifdef PROFILING
  double time_loop_next = 0.0;
  time_loop_next = - vSched_get_wtime();
#endif

  int t_x  = -1;
  if (count == 0) return 0;
  pthread_mutex_lock(&(dynwork[tid]->qLock));
  if(dynwork[tid]->nextChunk < dynwork[tid]->limit) // the thread still has work to be done in its own queue
  {
    #ifdef VERBOSE
      printf("loop_next_sds(): thread %d . There is work in the local queue. nextChunk = %d \t limit = %d \n", tid, nextChunk, dynwork[tid]->limit);
      #endif
    *pstart = dynwork[tid]->nextChunk;
    dynwork[tid]->nextChunk = dynwork[tid]->nextChunk + dynwork[tid]->chunkSize;
    if(dynwork[tid]->nextChunk > dynwork[tid]->limit) dynwork[tid]->nextChunk = dynwork[tid]->limit;
    *pend  = dynwork[tid]->nextChunk;
    pthread_mutex_unlock(&(dynwork[tid]->qLock));
    if(dynwork[tid]->nextChunk >= dynwork[tid]->limit) // this thread is done with its own queue
    {
      pthread_mutex_lock(&sched_lock);
      count--;
      pthread_mutex_unlock(&sched_lock);
    }

#ifdef PROFILING
      time_loop_next += vSched_get_wtime();
      printf("loop_next_sds(): thread %d \t time = %f  \n" , tid, time_loop_next);
#endif
    return 1;
  }
  else // we steal from another thread
  {
//    printf("[%d]: loop_next, steal: count=%d\n", tid,count);
    pthread_mutex_unlock(&(dynwork[tid]->qLock));
    if(count == 0) return 0;
    t_x = selectAnotherThread(tid, threadCount);
    if (t_x == -1) // we couldn't steal from another thread, as no other thread has work
    {
      return 0;
    }
    else // there is work from another thread to be stolen
    {
      if(count == 0) return 0;
      pthread_mutex_lock(&(dynwork[t_x]->qLock));
#ifdef VERBOSE
      printf("loop_next_sds(): thread %d \t There is work from another thread to be stolen\n", tid);
#endif
      *pstart = dynwork[t_x]->nextChunk;
      if (*pstart >= dynwork[t_x]->limit ) {      pthread_mutex_unlock(&(dynwork[t_x]->qLock)); return 0; }
      dynwork[t_x]->nextChunk = dynwork[t_x]->nextChunk + dynwork[t_x]->chunkSize;
      *pend = dynwork[t_x]->nextChunk;
      if(*pend > dynwork[t_x]->limit) *pend = dynwork[t_x]->limit;
      pthread_mutex_unlock(&(dynwork[t_x]->qLock));
      return 1;
    }
  } // end condition for stealing
}

/*
 get_constraint() :  used internally by vSched to determine the constraint of the scheduling
  - returns 1 if constraint is satisfied and returns 0 if constraint is not satisfied
  - returns -1 on error in calculation (this  should never occurs, unless of course the implementer of this function made a mistake in the calculation. The implementer must take care to test their functions before using it. Erroneous behavior may occur. )
  - other parameters can be added here, but for simplicity we leave it without params
  - the constraint_ value is set at the beginning of the threaded computation using the function setCDY
  - note that this should be a fairly fast calculation, as this will be invoked many times.
*/

int get_constraint()
{
  double rand_val = (double)rand()/(double)RAND_MAX;

  if(rand_val < 0.0 || rand_val > 1.0)
    return -1;
  #ifdef VERBOSE
    printf("constraint_ =  %f \t rand= %f \n", constraint_, rand_val);
#endif
// below is a simple condition  used for constraint. More complicated functions can be used below, to implementor's liking */
  if(rand_val < constraint_)
    return 1;
  else
    return 0;
}


// TODO :  figure this out for both threaded and non-threaded regions

double vSched_get_wtime(void)
{
  struct rusage ruse;
// getrusage(RUSAGE_SELF, &ruse);
// see documentation , need to use rusage thread for proper timing
//  if(inThreadedRegion)
#ifdef RUSAGE_THREAD
  getrusage(RUSAGE_THREAD, &ruse);
#else
  getrusage(RUSAGE_CHILDREN, &ruse);
#endif

#ifdef VERBOSE
    printf("vSched_get_wtime(): \n") ;
#endif

  return( (double)(ruse.ru_utime.tv_sec+ruse.ru_utime.tv_usec / 1000000.0) );
}


double nont_vSched_get_wtime(void)
{
  struct rusage ruse;
getrusage(RUSAGE_SELF, &ruse);
// see documentation , need to use rusage thread for proper timing

#ifdef VERBOSE
    printf("vSched_get_wtime(): \n") ;
#endif
  return( (double)(ruse.ru_utime.tv_sec+ruse.ru_utime.tv_usec / 1000000.0) );
}
