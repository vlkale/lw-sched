extern void vSched_init(int);
extern void vSched_finalize(int);
extern double vSched_get_wtime();
extern double nont_vSched_get_wtime();

extern int enqueue(int taskNumStart, int _taskNumEnd, int *pstart, int *pend, int threadID, int numThreads);
extern int dequeue(int *pstart, int *pend, int tid);
