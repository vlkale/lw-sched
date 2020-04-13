#OPTS  = -O3 -lpthread -lrt -lm -openmp
OPTS  = -O3 -lpthread -lrt -lm
MPICC = mpicc
MPICXX = mpicxx
CC  = gcc $(OPTS)
GCC = gcc $(OPTS)
GXX = c++ $(OPTS)
CLANGXX = clang++
DEBUG=-g

all: test_vSched

test_vSched: appFor_vSchedSimple.c vSched.h vSched.c
	$(CLANGXX) -fPIC vSched.c appFor_vSchedSimple.c -DCDY_ $(OPTS) -o test_vSched

test_tQueue: appFor_vSchedSimple.C
	$(GXX) -fPIC -g threadedQueue.h threadedQueue.C appFor_threadedQueueSimple.C -DCDY_ $(OPTS) -o test_vSched

test_vSchedforMac: appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c
	$(GXX) -fPIC -g vSched.h vSched.c appFor_vSchedSimple.c pthBarrierforOSX.c -DCDY_ $(OPTS) -o test_vSchedforMac

test_tqueue_forMac: threadedQueue.h threadedQueue.C pthBarrierforOSX.c
	$(GXX) -fPIC -g threadedQueue.h threadedQueue.C appFor_tqueueSimple.c pthBarrierforOSX.c $(OPTS) -o test_tqueue_forMac

test_vSchedOpenMP: appFor_vSchedSimpleOpenMP.c vSched.h vSched.c
	$(GXX) -fPIC -g -I. vSched.h vSched.c appFor_vSchedSimpleOpenMP.c $(OPTS) -o test_vSchedOpenMP

test_vSchedomp: appFor_vSched-omp.C vSched.h vSched.c
	$(MPICXX) -fPIC -g $(OPTS) -I. vSched.h vSched.c appFor_vSched-omp.C -o test_vSchedomp

test: test_vSched
	./test_vSched 65536 1000 16 64 0.5 0.1

tgz: appFor_vSchedSimpleOpenMP.c appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c 
	tar -cvzf appFor_vSchedSimpleOpenMP.c appFor_vSchedSimple.c vSched.h vSched.c pthBarrierforOSX.c README

clean:
	rm -rf *.o test_vSched test_tqueue_forMac test_vSchedomp

realclean:
	rm -rf *.o core *.gch test_vSched test_vSchedforMac test_vSchedOpenMP test_vSchedomp
