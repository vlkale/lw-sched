#include <pthread.h>
using namespace std;

//flag used for debugging output - uncomment the line below if you want verbose debugging.
//#define VERBOSE

#define PROFILING
#include <stdio.h> 
// use this for testing
pthread_mutex_t sched_lock;

#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
// NOTE: one could declare the variable in line below as static so that we maintain state over loop iterations. 
float f_s;
float f_d;
float constraint_;
int taskSize;
int nextTask;
int taskNumEnd;
int isTaskStarted;
int count;
int threadCount;

// functions internal to this library
int get_constraint();
double vSched_get_wtime();

typedef struct ThreadedQueue // come up with a better name
{
  int nextChunk;
  int chunkSize;
  int limit; //this is the iteration that we should not assign.
  pthread_mutex_t qLock;
} ThreadedQueue;

ThreadedQueue** my_t_queue;

int selectAnotherThread(int tid, int numThreads);

int vSched_thread_init()
{ 
  // allocate the data structure for each thread
}

void vSched_init(int numThreads)
{
  pthread_mutex_init(&sched_lock, NULL);
  threadCount = numThreads;
  my_t_queue = (PossibleWork**) malloc(sizeof(void*)*numThreads);
  for (int i = 0 ; i < numThreads; i++)
    {
      my_t_queue[i] = (PossibleWork*) malloc(sizeof(PossibleWork));
      pthread_mutex_init(&(my_t_queue[i]->qLock), NULL);
    }
}

void vSched_finalize(int numThreads)
{
  pthread_mutex_destroy(&sched_lock);
}

int enqueue(int taskNumBegin, int _taskNumEnd, int *pstart, int *pend, int threadID, int numThreads )  // think about adding to parameter list here
{
  pthread_mutex_lock(&sched_lock); // check if the lock is needed around the whole conditional clause
  if(count == 0)
    count = numThreads;
  pthread_mutex_unlock(&sched_lock);
  pthread_mutex_lock(&my_t_queue[threadID]->qLock);
  // printf("loop_start_cdy(): thread %d \n", threadID);
  // loopEnd = _loopEnd;
  taskNumLast = _taskNumEnd;
  // The following code distributes the endpoint of each thread's queue. 
  my_t_queue[threadID]->limit = taskNumBegin + ((numTasks)*(threadID+1))/numThreads; // the limit is the begin point of the threadID + 1.  
  *pstart = taskNumBegin + (numTasks*threadID)/numThreads; /* figure out algebra here, based on taskBegin */
  *pend = *pstart + numTasks/numThreads; /* figure out algebra here, based on taskBegin */
  my_t_queue[threadID]->nextTask = *pend;
  my_t_queue[threadID]->chunkSize = chunkSize;
  pthread_mutex_unlock(&(my_t_queue[threadID]->qLock));
  if (my_t_queue[threadID]->nextTask >= my_t_queue[threadID]->limit) // this thread is done with its own queue
  {
    pthread_mutex_lock(&sched_lock);
    count--;
    pthread_mutex_unlock(&sched_lock);
  }
// print pstart and pend for VERBOSE debugging.
#ifdef VERBOSE
    printf("thread%d:\t pstart =  %d \t pend = %d \t nextTask = %d \n ", threadID, *pstart, *pend, my_t_queue[threadID]->nextTask);
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
// For now we choose the next tid
//  another_tid = (tid + 1)%numThreads;
  another_tid = (tid + 1)%numThreads ; // hardcode to 16 for now
// Strategy 1: randomized stealing
// Strategy 2: steal from thread with most elements in queue
  return another_tid;
}


/*
This is the staggered method for mixed static/dynamic scheduling.
*/
int dequeue(int *pstart, int *pend, int tid)
{
#ifdef PROFILING
  double time_loop_next = 0.0;
  time_loop_next = - vSched_get_wtime();
#endif
  int t_x  = -1;
  if (count == 0) return 0;
  pthread_mutex_lock(&(my_t_queue[tid]->qLock));
  if(my_t_queue[tid]->nextTask < my_t_queue[tid]->limit) // the thread still has work to be done in its own queue
  {
    #ifdef VERBOSE
    printf("task_next_sds(): thread %d . There are tasks in the local queue. nextTask = %d \t limit = %d \n", tid, nextTask, my_t_queue[tid]->limit);
      #endif
    *pstart = my_t_queue[tid]->nextTask;
    my_t_queue[tid]->nextTask = my_t_queue[tid]->nextTask;
    if(my_t_queue[tid]->nextTask > my_t_queue[tid]->limit) my_t_queue[tid]->nextTask = my_t_queue[tid]->limit;
    *pend  = my_t_queue[tid]->nextTask;
    pthread_mutex_unlock(&(my_t_queue[tid]->qLock));
    if(my_t_queue[tid]->nextTask >= my_t_queue[tid]->limit) // this thread is done with its own queue
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
      pthread_mutex_unlock(&(my_t_queue[tid]->qLock));
      if(count == 0) return 0;
      t_x = selectAnotherThread(tid, threadCount);
      if (t_x == -1) // we couldn't steal from another thread, as no other thread has work
	return 0;
      else // there is work from another thread to be stolen
	{
	  if(count == 0) return 0;
	  pthread_mutex_lock(&(my_t_queue[t_x]->qLock));
#ifdef VERBOSE
	  printf("loop_next_sds(): thread %d \t There is work from another thread to be stolen\n", tid);
#endif
	  *pstart = my_t_queue[t_x]->nextTask;
	  if (*pstart >= my_t_queue[t_x]->limit ) { pthread_mutex_unlock(&(my_t_queue[t_x]->qLock)); return 0; }
	  //  my_t_queue[t_x]->nextTask = my_t_queue[t_x]->nextTask + my_t_queue[t_x]->chunkSize;
	  *pend = my_t_queue[t_x]->nextTask;
	  if(*pend > my_t_queue[t_x]->limit) *pend = my_t_queue[t_x]->limit;
	  pthread_mutex_unlock(&(my_t_queue[t_x]->qLock));
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

// TODO: make the function vSched_get_wtime work for both threaded computation regions and serial computation regions.

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
