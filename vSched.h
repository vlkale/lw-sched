


extern void vSched_init(int);
extern void vSched_finalize(int);
extern double vSched_get_wtime();
extern double nont_vSched_get_wtime();


// TODO: figure out whether below should be extern'd
extern int loop_start_static(int loopBegin, int loopEnd, int *pstart, int *pend, int threadID, int numThreads);
extern int loop_next_static(int *pstart, int *pend);

extern void setStaticFraction(float f, int _chunkSize);
extern void setCDY(float f, double constraint, int _chunkSize);

extern int loop_start_static_fraction(int loopBegin, int _loopEnd, int *pstart, int *pend, int threadID, int numThreads);
extern int loop_next_static_fraction(int *pstart, int *pend);

extern int loop_start_cdy(int loopBegin, int _loopEnd, int *pstart, int *pend, int threadID, int numThreads);
 extern int loop_next_cdy(int *pstart, int *pend);
// extern int loop_next_cdy(int *pstart, int *pend, double constraint_val);

extern int loop_start_cdy(int loopBegin, int _loopEnd, int *pstart, int *pend, int threadID, int numThreads);
extern int loop_next_cdy(int *pstart, int *pend, int tid);

extern int loop_start_statdynstaggered(int loopBegin, int _loopEnd, int *pstart, int *pend, int threadID, int numThreads);
extern int loop_next_statdynstaggered(int *pstart, int *pend, int tid);
